%% Quadcopter UAV Dynamics
% Reference: Thanh 2022
%
% Simulates a 12-state quadcopter model by exciting each control input
% individually and plotting the resulting state trajectories.
%
% State vector (12x1):
%   [x, y, z, phi, theta, psi, xdot, ydot, zdot, phidot, thetadot, psidot]
%   positions (m), Euler angles (rad), linear velocities (m/s), angular rates (rad/s)
%
% Input vector (4x1):
%   [fz, tau_phi, tau_theta, tau_psi]
%   total thrust (N) and body-frame torques (N.m)

%% Model Parameters
m_kg  = 1.12;        % mass [kg]
g_ms2  = 9.81;        % gravitational acceleration [m/s^2]
l_m  = 0.23;        % arm length [m]
Jr_kgm2 = 8.5e-4;      % rotor moment of inertia [kg.m^2]
Ix_kgm2 = 0.0019;      % roll moment of inertia [kg.m^2]
Iy_kgm2 = 0.0019;      % pitch moment of inertia [kg.m^2]
Iz_kgm2 = 0.0223;      % yaw moment of inertia [kg.m^2]
b_nrads2  = 7.73212e-6;  % thrust coefficient [N/(rad/s)^2]
d_nrads2  = 1.27513e-7;  % drag (torque) coefficient [N.m/(rad/s)^2]

% Pack inertial parameters into a vector for passing to the dynamics function
params = [m_kg; g_ms2; Jr_kgm2; Ix_kgm2; Iy_kgm2; Iz_kgm2];

%% Simulation Setup
tspan = [0 10];         % simulation time window [s]
s0    = zeros(12,1);    % initial state: at rest at origin

% Magnitude of each isolated input perturbation
input_mag = [12; 0.001; 0.001; 0.001];   % [fz (N); tau_phi; tau_theta; tau_psi (N.m)]

% Human-readable labels for titling each figure
control_inputs = {'fz', 'tau phi', 'tau theta', 'tau psi'};

%% Single-Input Sweep: isolate each control input and simulate
for i = 1:4
    % Build input vector with only the i-th input active
    % (Fixed: previously ~i=0 was used incorrectly; now u is zeroed then set)
    u    = zeros(4,1);
    u(i) = input_mag(i);

    % For torque inputs (i > 1), apply hover thrust so gravity is balanced
    % and rotational effects are isolated from free-fall
    if i ~= 1
        u(1) = m * g;   % hover thrust: fz = mg => zddot = 0
    end

    % Integrate the nonlinear equations of motion with ode45
    [t, s] = ode45(@(t,s) dynamics(t, s, u, params), tspan, s0);

    %% Plotting
    figure;

    % --- Position ---
    subplot(4,1,1);
    plot(t, s(:,1:3));
    xlabel('Time (s)'); ylabel('Position (m)');
    legend('x','y','z');
    title(sprintf('Control Input: %s = %.3f, fz = %.2f', control_inputs{i}, u(i), u(1)));
    grid on;

    % --- Euler Angles (converted to degrees and wrapped to [0, 360)) ---
    subplot(4,1,2);
    plot(t, mod(rad2deg(s(:,4:6)), 360));
    xlabel('Time (s)'); ylabel('Orientation (deg)');
    legend('phi','theta','psi');
    grid on;

    % --- Linear Velocities ---
    subplot(4,1,3);
    plot(t, s(:,7:9));
    xlabel('Time (s)'); ylabel('Velocity (m/s)');
    legend('xdot','ydot','zdot');
    grid on;

    % --- Angular Rates (converted to degrees and wrapped to [0, 360)) ---
    subplot(4,1,4);
    plot(t, mod(rad2deg(s(:,10:12)), 360));
    xlabel('Time (s)'); ylabel('Angular Velocity (deg/s)');
    legend('phidot','thetadot','psidot');
    grid on;

    % Link all subplot x-axes and enforce the simulation time window
    linkaxes(findall(gcf,'type','axes'), 'x');
    xlim(tspan);
end

%% Dynamics Function
% Implements the nonlinear equations of motion (Thanh 2022, Eq. 6).
%
% Inputs:
%   ~      - time (unused; system is time-invariant)
%   s      - current state vector (12x1)
%   u      - control input vector [fz; tau_phi; tau_th; tau_ps] (4x1)
%   params - physical parameters [m; g; Jr; Ix; Iy; Iz]
%
% Output:
%   sdot   - state derivative vector (12x1)
function sdot = dynamics(~, s, u, params)
    % Unpack parameters
    m  = params(1); g  = params(2); Jr = params(3);
    Ix = params(4); Iy = params(5); Iz = params(6);

    % Unpack state: linear velocities
    xdot   = s(7);  ydot  = s(8);  zdot  = s(9);

    % Unpack state: Euler angles
    phi    = s(4);  th    = s(5);  ps    = s(6);

    % Unpack state: angular rates
    phidot = s(10); thdot = s(11); psdot = s(12);

    % Unpack control inputs
    fz      = u(1);
    tau_phi = u(2);
    tau_th  = u(3);
    tau_ps  = u(4);

    % Net rotor speed (gyroscopic effect term); set to 0 when not modelled
    Om = 0;

    % Translational accelerations in inertial frame (Eq. 6)
    % Derived from rotation matrix R(phi, theta, psi) applied to body-frame thrust [0; 0; fz]
    xddot = (cos(phi)*sin(th)*cos(ps) + sin(phi)*sin(ps)) * fz/m;
    yddot = (cos(phi)*sin(th)*sin(ps) - sin(phi)*cos(ps)) * fz/m;
    zddot = -g + cos(phi)*cos(th) * fz/m;   % gravity subtracted from vertical thrust

    % Rotational accelerations (Euler's rigid-body equations, Eq. 6)
    % Cross-product (gyroscopic) terms couple the three axes
    phiddot = thdot*psdot*(Iy-Iz)/Ix - Jr/Ix * thdot*Om + tau_phi/Ix;
    thddot  = phidot*psdot*(Iz-Ix)/Iy + Jr/Iy * phidot*Om + tau_th/Iy;
    psddot  = phidot*thdot*(Ix-Iy)/Iz + tau_ps/Iz;

    % Assemble state derivative: [pos_dot; angle_dot; vel_dot; rate_dot]
    sdot = [xdot; ydot; zdot; phidot; thdot; psdot;
            xddot; yddot; zddot; phiddot; thddot; psddot];
end