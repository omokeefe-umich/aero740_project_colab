%% Quadcopter UAV Dynamics — Linearized MPC Prototype
% Reference: Thanh 2022
%
% State vector x (12x1):
%   [x, y, z, phi, theta, psi, xdot, ydot, zdot, phidot, thetadot, psidot]
%   Positions (m), Euler angles (rad), and their time derivatives.
%
% Input vector u (4x1):
%   [fz, tau_phi, tau_th, tau_ps]
%   Total thrust (N) and roll/pitch/yaw torques (N·m).
%
% Augmented state (28x1) used for offset-free MPC:
%   [delta_x (12x1); x (12x1); u (4x1)]
%   delta_x = x_actual - x_nominal captures linearization error so the
%   controller corrects it at the next step (integral-like disturbance
%   rejection without explicit integrators).

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

%% Simulation Time Grid (outer loop, dt = 0.01 s)
tspan = 0:0.01:2;   % 201 time steps over 2 seconds

%% Initial Conditions
% Quadcopter
x_start_qc = -10; % initial x position (m)
y_start_qc = 0; % initial y position (m)
z_start_qc = 100;   % initial z position (m)
deltax = zeros(12, 1);                    % initial linearization error (zero at start)
x      = [x_start_qc; y_start_qc; z_start_qc; zeros(9, 1)];  % initial state (hover at specified position)

% Target Vehicle
x_start_tgt = 0; 

u      = zeros(4, 1);                     % initial input (no thrust, no torques)
Xs     = zeros(12, numel(tspan));         % full 12-state history: [x;y;z;phi;…] at each tspan(k)
Xrefs  = zeros(3, numel(tspan));         % reference position history for plotting

%% Cost Matrix Construction
% C selects the absolute position states [x, y, z] from the 28-D augmented
% state. The first 12 elements are delta_x; elements 13-15 are [x, y, z].
C = [diag([1; 1; 1])*0.0, zeros(3, 9), diag([1; 1; 1]), zeros(3, 13)];

% Stage cost Q (28x28):
%   C'*C       — penalizes absolute position error (hard tracking weight)
%   + 0.01*I on the remaining delta_x components — mild regularization on
%     velocity and attitude errors inside the linearization residual.
Q = C' * C + blkdiag(zeros(3), eye(9) * 0.0, zeros(16));

%% MPC Control Loop
for k = 1 : 201

    % --- Linearize about the current operating point ---
    % Ac = df/dx evaluated at the current state x and input u
    Ac = Jacobian(@(x) dynamics(0, x, u, params), x);
    % Bc = df/du evaluated at the disturbed state (x + deltax) and input u
    Bc = Jacobian(@(u) dynamics(0, x + deltax, u, params), u);

    % --- Discretize the linearized system (dt = 0.01 s, 20-step horizon) ---
    % Returns Ads and Bds as 3-D arrays [N_steps x N_x x N_x/N_u]
    [Ads, Bds] = Discretize_dt([0.01], ones(20, 1), Ac, Bc);

    % --- Build augmented discrete-time matrices (28x28 state, 4 inputs) ---
    %   [delta_x_{k+1}]   [Ad  0 ] [delta_x_k]   [Bd]
    %   [    x_{k+1}  ] = [I  I  ] [  x_k    ] + [0 ] * u_k
    %   [    u_{k+1}  ]   [0  I  ] [  u_k    ]   [I ]
    % This appends the nominal state and input as "memory" states so the
    % QP can track absolute trajectories rather than deviations only.
    Ads_augmented = zeros(20, 28, 28);
    Bds_augmented = zeros(20, 28, 4);
    for i = 1:20
        Ad_i = squeeze(Ads(i, :, :));
        Bd_i = squeeze(Bds(i, :, :));
        Ads_augmented(i, :, :) = [Ad_i,           zeros(12, 16); ...
                                   [eye(12); zeros(4, 12)], eye(16)];
        Bds_augmented(i, :, :) = [Bd_i; zeros(12, 4); eye(4)];
    end

    % --- Compute prediction matrices S and M ---
    % Predicted state sequence: X = M*x_aug + S*U
    [S, M] = StackedMatrix(Ads_augmented, Bds_augmented);

    % --- Assemble current augmented state ---
    x_aug = [deltax; x; u];  % 28x1

    % --- Define desired trajectory over the prediction horizon ---
    t_horizon = (tspan(k) + 0.01) : 0.01 : (tspan(k) + 0.2);  % 20 future times
    [x_traj, x_ref] = ReferenceTrajectory(t_horizon, 'ground_vehicle');  % 36x1 stacked [x_1; y_1; z_1; ...; x_20; y_20; z_20]

    % --- Formulate the QP (H, q) ---
    % Minimizes:  (1/2)*U'*H*U + q'*U
    %   Q          — stage position/state cost (28x28)
    %   R = 0.001*I — control effort weight (4x4)
    %   P = Q*0.01  — terminal cost (28x28, scaled down)
    [H, q] = QPFormat(Q, diag([0.001; 0.001; 0.001; 0.001]), Q * 0.01, ...
        S, M, (0:0.01:0.2)', x_aug, x_traj, zeros(80, 1));

    % --- Solve QP for optimal input sequence ---
    % quadprog requires H to be symmetric; (H+H')/2 enforces this numerically
    u_mpc = quadprog((H + H') / 2, q);

    % Apply only the first control increment (receding-horizon principle)
    u = u + u_mpc(1:4);

    % --- Propagate the true state forward by one outer timestep ---
    % Shift nominal state by the residual from the previous step
    x = x + deltax;

    % Record the full 12-state at this outer time step
    Xs(:, k) = x;
    Xrefs(:, k) = x_ref(1:3, 1);  % reference position at the first prediction step

    % Integrate nonlinear dynamics over [0, 0.01] s at fine resolution
    [T, X] = ode45(@(t, x) dynamics(t, x, u, params), 0:0.00001:0.01, x);

    % Compute linearization residual for next step's augmented state
    deltax = X(end, :)' - x;
end

%% Reference Trajectory & Visualization
PlotResults(tspan, Xs, Xrefs);