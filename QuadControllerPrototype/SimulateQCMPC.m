%% Quadcopter MPC Simulation — Thanh 2022
% x (12x1): [x,y,z, phi,theta,psi, xdot,ydot,zdot, phidot,thetadot,psidot]
% u  (4x1): [fz, tau_phi, tau_th, tau_ps]
% Augmented state (40x1): [delta_x; x; u; r]
%   delta_x = x_actual - x_nominal (offset-free disturbance rejection)

clear all;

%% Physical Parameters
m  = 1.12;        % mass [kg]
g  = 9.81;        % gravitational acceleration [m/s^2]
l  = 0.23;        % rotor arm length [m]
Jr = 8.5e-4;      % rotor moment of inertia [kg·m^2]
Ix = 0.0019;      % roll  moment of inertia [kg·m^2]
Iy = 0.0019;      % pitch moment of inertia [kg·m^2]
Iz = 0.0223;      % yaw   moment of inertia [kg·m^2]
b  = 7.73212e-6;  % thrust coefficient [N/(rad/s)^2]
d  = 1.27513e-7;  % drag  coefficient [N·m/(rad/s)^2]
params = [m; g; Jr; Ix; Iy; Iz];

%% Simulation Time Grid (outer loop)
tend = 20;          % total simulation time (s)
dt = 0.1;              % outer loop time step (s)
tspan = 0:dt:tend;

%% MPC Parameters
N_mpc = tend;
dt_mpc = dt/5; 

%% Initial Conditions
% Quadcopter
x_start_qc = 0; % initial x position (m)
y_start_qc = 0; 
z_start_qc = 10;  
deltax = zeros(12, 1);                    % initial linearization error (zero at start)
x      = [x_start_qc; y_start_qc; z_start_qc; zeros(9, 1)];

% Target Vehicle
% reference trajectory vector for target vehicle, of same form as x for quadcopter (12x1)
x_start_tgt = 10; % initial x position of target vehicle (m)
y_start_tgt = 0; 
z_start_tgt = 0; 
r = [x_start_tgt; y_start_tgt; z_start_tgt; zeros(9,1)]; 

% Initial input (thrust and torques)
u      = zeros(4, 1);                     % initial input (no thrust, no torques)

% Preallocate arrays to store simulation history for plotting
Xs     = zeros(12, numel(tspan));         % full 12-state history
Xrefs  = zeros(3, numel(tspan));          % reference position history [x;y;z] of target vehicle at each tspan(k)

%% Cost Matrix Construction
% Weights = 1 / (acceptable_error)^2
q_pos = 1 / (0.05^2);
q_vel = 1 / (0.1^2);
r_fz  = 1 / (2.0^2);
r_tau = 1 / (0.5^2);

% Q encodes ||x_pos - r_pos||^2 via cross-terms (r at augmented indices 29:31)
Q = zeros(40, 40);
Q(13:15, 13:15) =  diag([q_pos, q_pos, q_pos]);
Q(13:15, 29:31) = -diag([q_pos, q_pos, q_pos]);
Q(29:31, 29:31) =  diag([q_pos, q_pos, q_pos]);
Q(29:31, 13:15) = -diag([q_pos, q_pos, q_pos]);
Q(19:21, 19:21) = diag([q_vel, q_vel, q_vel]);
R = diag([r_fz, r_tau, r_tau, r_tau]);
Q(25:28, 25:28) = R;
P = 100 * Q;  % terminal cost

%% MPC Control Loop
for k = 1 : numel(tspan)
    % Linearize at current operating point
    Ac = Jacobian(@(x) dynamics(0, x, u, params), x);
    Bc = Jacobian(@(u) dynamics(0, x + deltax, u, params), u);

    % Discretize linearized system over MPC horizon
    [Ads, Bds] = Discretize_dt([dt_mpc], ones(N_mpc, 1), Ac, Bc);

    % Build augmented discrete-time matrices (40x40):
    %   [delta_x; x; u; r]_{k+1} = A_aug*[...] + B_aug*delta_u
    Ads_augmented = zeros(N_mpc, 40, 40);
    Bds_augmented = zeros(N_mpc, 40, 4);
    for i = 1:20
        Ad_i = squeeze(Ads(i, :, :));
        Bd_i = squeeze(Bds(i, :, :));
        Ads_augmented(i, :, :) = [
                                Ad_i,           zeros(12, 12),  zeros(12, 4),   zeros(12, 12); ...
                                eye(12),        eye(12),        zeros(12, 4),   zeros(12, 12); ...
                                zeros(4, 12),   zeros(4, 12),   eye(4),         zeros(4, 12); ...
                                zeros(12, 12),  zeros(12, 12),  zeros(12, 4),   eye(12)
                                ];
        Bds_augmented(i, :, :) = [Bd_i;         zeros(12, 4);   eye(4);         zeros(12, 4)];
    end

    % Build prediction matrices (X = M*x_aug + S*U)
    [S, M] = StackedMatrix(Ads_augmented, Bds_augmented);

    x_aug     = [deltax; x; u; r];
    t_horizon = tspan(k) : dt_mpc : tspan(k)+N_mpc*dt_mpc;

    % Solve QP
    [H, q] = QPFormat(Q, R, P, S, M, t_horizon', x_aug, zeros(4*N_mpc, 1));
    u_mpc = quadprog(H, q);

    % Apply first control increment (receding horizon)
    u = u + u_mpc(1:4);

    % Advance nominal state
    x = x + deltax;

    Xs(:, k)    = x;
    Xrefs(:, k) = r(1:3);
    dist_diff   = sqrt(sum((Xs(1:3, k) - Xrefs(:, k)).^2));
    if dist_diff < 0.1
        disp('*'*50)
        fprintf('------> IMPACT @ %i seconds!! <---\n', tspan(k))
        disp('*'*50)
    end

    % Integrate nonlinear dynamics one outer step
    [T, X] = ode45(@(t, x) dynamics(t, x, u, params), 0:0.0001:0.01, x);

    % Linearization residual for next step
    deltax = X(end, :)' - x;

    % Update the trajectory of the target
    r(1) = r(1) + 0.01 * dt;
end

%% Reference Trajectory & Visualization
PlotResults(tspan, Xs, Xrefs);