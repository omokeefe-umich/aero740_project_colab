%% Quadcopter UAV Dynamics — Linearized MPC Prototype with Integral Error
% Reference: Thanh 2022
%
% Physical plant state x_phys (12x1):
%   [x, y, z, phi, theta, psi, xdot, ydot, zdot, phidot, thetadot, psidot]
%
% Input vector u (4x1):
%   [fz, tau_phi, tau_th, tau_ps]
%
% Augmented prediction state X_aug (40x1):
%   [delta_x (12x1); eI (12x1); xerr (12x1); u (4x1)]
%
% where
%   delta_x = x_actual,next - x_nominal,next   (linearization residual)
%   xerr    = x - r                            (instantaneous tracking error)
%   eI      = integral of xerr over time
%
% Notes:
% - The plant still propagates using the absolute physical state x.
% - The MPC prediction model is formulated in terms of the tracking-error
%   state xerr so that the reference does not need to be carried explicitly
%   inside the augmented state.
% - This keeps the total augmented size at 40, matching QPFormat.m.

clear all;
clc;

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
tend  = 20;        % total simulation time [s]
dt    = 0.1;       % outer-loop sample time [s]
tspan = 0:dt:tend;

%% MPC Parameters
N_mpc  = tend * 20;   % horizon length in MPC steps
dt_mpc = dt / 2;      % MPC discretization time [s]

%% Initial Conditions
% --- Physical plant state ---
x_start_qc = 0;
y_start_qc = 0;
z_start_qc = -10;
x = [x_start_qc; y_start_qc; z_start_qc; zeros(9,1)];   % absolute plant state

% --- Reference state ---
x_start_tgt = 10;
y_start_tgt = 0;
z_start_tgt = 0;
r = [x_start_tgt; y_start_tgt; z_start_tgt; zeros(9,1)];

% --- Augmented-state components ---
deltax = zeros(12,1);         % linearization residual
xerr   = x - r;               % instantaneous tracking error
eI     = zeros(12,1);         % integral of tracking error

% --- Initial input ---
u = zeros(4,1);

%% Storage for visualization
Xs    = zeros(12, numel(tspan));   % absolute plant state history
Xrefs = zeros(3,  numel(tspan));   % reference position history

%% Cost Matrix Construction
% Penalize:
%   - integral tracking error eI heavily
%   - instantaneous tracking error xerr moderately
%   - input usage lightly
%
% Augmented ordering:
%   1:12   -> delta_x
%   13:24  -> eI
%   25:36  -> xerr
%   37:40  -> u

smallcost = 1e-3;

Q = zeros(40,40);

% Integral error penalty (stronger)
QeI_pos = 50 * smallcost * diag([1 1 1]);
Q(13:15,13:15) = QeI_pos;

% Instantaneous tracking error penalty
Qxerr_pos = 1 * smallcost * diag([1 1 1]);
Q(25:27,25:27) = Qxerr_pos;

% Optional mild penalties on attitude / velocity tracking error if desired
% Q(28:36,28:36) = 0.1 * smallcost * eye(9);

R = diag([smallcost; smallcost; smallcost; smallcost]);

% Terminal cost
P = Q / smallcost;

%% MPC Control Loop
for k = 1:numel(tspan)

    % Current instantaneous tracking error
    xerr = x - r;

    % -------------------------------------------------------------
    % Linearize about the current PHYSICAL operating point (x, u)
    % -------------------------------------------------------------
    Ac = Jacobian(@(xlin) dynamics(0, xlin, u, params), x);
    Bc = Jacobian(@(ulin) dynamics(0, x + deltax, ulin, params), u);

    % -------------------------------------------------------------
    % Discretize the physical linearized system
    % -------------------------------------------------------------
    [Ads, Bds] = Discretize_dt([dt_mpc], ones(N_mpc,1), Ac, Bc);

    % -------------------------------------------------------------
    % Build augmented prediction model for:
    %   X_aug = [delta_x; eI; xerr; u]
    %
    % Model:
    %   delta_x(k+1) = Ad*delta_x(k) + Bd*du(k)
    %   xerr(k+1)    = xerr(k) + delta_x(k+1)
    %   eI(k+1)      = eI(k) + dt_mpc * xerr(k+1)
    %   u(k+1)       = u(k) + du(k)
    %
    % This treats the xerr block as the predicted tracking-error state.
    % -------------------------------------------------------------
    Ads_augmented = zeros(N_mpc, 40, 40);
    Bds_augmented = zeros(N_mpc, 40, 4);

    for i = 1:N_mpc
        Ad_i = squeeze(Ads(i,:,:));
        Bd_i = squeeze(Bds(i,:,:));

        % Block rows:
        % 1) delta_x_{k+1}
        % 2) eI_{k+1}
        % 3) xerr_{k+1}
        % 4) u_{k+1}
        Ads_augmented(i,:,:) = [ ...
            Ad_i,               zeros(12,12), zeros(12,12), zeros(12,4);   % delta_x+
            dt_mpc * Ad_i,      eye(12),      dt_mpc*eye(12), zeros(12,4); % eI+
            Ad_i,               zeros(12,12), eye(12),      zeros(12,4);   % xerr+
            zeros(4,12),        zeros(4,12),  zeros(4,12),  eye(4)         % u+
        ];

        Bds_augmented(i,:,:) = [ ...
            Bd_i;                 % delta_x+
            dt_mpc * Bd_i;        % eI+
            Bd_i;                 % xerr+
            eye(4)                % u+
        ];
    end

    % -------------------------------------------------------------
    % Prediction matrices
    % -------------------------------------------------------------
    [S, M] = StackedMatrix(Ads_augmented, Bds_augmented);

    % -------------------------------------------------------------
    % Assemble current augmented state
    % -------------------------------------------------------------
    x_aug = [deltax; eI; xerr; u];

    % Time nodes for QPFormat
    t_horizon = tspan(k):dt_mpc:(tspan(k) + N_mpc*dt_mpc);

    % Zero input-reference trajectory
    u_ref = zeros(4*N_mpc,1);

    % -------------------------------------------------------------
    % Form QP
    % -------------------------------------------------------------
    [H, q] = QPFormat(Q, R, P, S, M, t_horizon', x_aug, u_ref);

    % Solve QP
    u_mpc = quadprog(H, q);

    % Receding-horizon control increment
    du = u_mpc(1:4);
    u  = u + du;

    % -------------------------------------------------------------
    % Propagate physical nonlinear plant over one outer-loop step
    % -------------------------------------------------------------
    x_nominal = x + deltax;

    % Record current state/reference before propagation
    Xs(:,k)    = x_nominal;
    Xrefs(:,k) = r(1:3);

    % Integrate nonlinear plant
    [T, X] = ode45(@(tloc, xloc) dynamics(tloc, xloc, u, params), ...
                   0:0.0001:dt, x_nominal);

    x_next = X(end,:)';

    % Linearization residual for next iteration
    deltax = x_next - x_nominal;

    % Update absolute plant state
    x = x_nominal;

    % Update integral error using measured tracking error at the new step
    % (rectangle rule at end of interval)
    xerr_next = x_next - r;
    eI = eI + dt * xerr_next;

    % Shift plant state to the true nonlinear state for next iteration
    x = x_next;

    % Example moving reference, if desired:
    % r(1) = r(1) + 0.5 * dt;
end

%% Visualization
PlotResults(tspan, Xs, Xrefs);