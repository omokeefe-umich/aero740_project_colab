%% Quadcopter Nonlinear Dynamics
% Reference: Thanh 2022, Eq. (6)
%
% Implements the continuous-time equations of motion for a quadrotor UAV
% under body-frame forces and torques.
%
% State s (12x1):
%   s(1:3)  = [x, y, z]            — inertial positions [m]
%   s(4:6)  = [phi, theta, psi]    — roll, pitch, yaw (ZYX Euler) [rad]
%   s(7:9)  = [xdot, ydot, zdot]   — inertial velocity [m/s]
%   s(10:12)= [phidot, thetadot, psidot] — Euler angle rates [rad/s]
%
% Input u (4x1):
%   u(1) = fz       — collective thrust force [N] (positive = upward in body frame)
%   u(2) = tau_phi  — roll  torque [N·m]
%   u(3) = tau_th   — pitch torque [N·m]
%   u(4) = tau_ps   — yaw   torque [N·m]
%
% Output sdot (12x1): time derivative of the state vector.

function sdot = dynamics(~, s, u, params)

    % Unpack physical parameters
    m  = params(1);   % mass [kg]
    g  = params(2);   % gravitational acceleration [m/s^2]
    Jr = params(3);   % rotor moment of inertia [kg·m^2]
    Ix = params(4);   % roll  inertia [kg·m^2]
    Iy = params(5);   % pitch inertia [kg·m^2]
    Iz = params(6);   % yaw   inertia [kg·m^2]

    % Unpack state
    xdot    = s(7);   ydot   = s(8);   zdot   = s(9);   % inertial velocities
    phi     = s(4);   th     = s(5);   ps     = s(6);   % Euler angles
    phidot  = s(10);  thdot  = s(11);  psdot  = s(12);  % Euler angle rates

    % Unpack inputs
    fz      = u(1);   % collective thrust
    tau_phi = u(2);   % roll  torque
    tau_th  = u(3);   % pitch torque
    tau_ps  = u(4);   % yaw   torque

    % Net rotor speed (Om = sum of signed rotor speeds for gyroscopic effect)
    % Set to zero here, meaning gyroscopic/propeller-reaction terms are neglected.
    Om = 0;

    % --- Translational Accelerations (Newton's 2nd law in inertial frame) ---
    % The rotation matrix R(phi,theta,psi) maps body-frame thrust [0;0;fz] to
    % inertial frame. Using ZYX Euler angle convention:
    xddot = (cos(phi)*sin(th)*cos(ps) + sin(phi)*sin(ps)) * fz/m;
    yddot = (cos(phi)*sin(th)*sin(ps) - sin(phi)*cos(ps)) * fz/m;
    zddot = -g + cos(phi)*cos(th) * fz/m;   % gravity opposes z, thrust assists

    % --- Rotational Accelerations (Euler's rigid-body equations) ---
    % Cross-axis inertia coupling (e.g., (Iy-Iz)*thdot*psdot / Ix) arises
    % from the off-diagonal terms in the full Euler equation I*alpha = tau - w×Iw.
    % The Jr*Om terms represent gyroscopic moments from spinning rotors (zero here).
    phiddot = thdot*psdot*(Iy-Iz)/Ix  - Jr/Ix * thdot*Om  + tau_phi/Ix;
    thddot  = phidot*psdot*(Iz-Ix)/Iy + Jr/Iy * phidot*Om + tau_th/Iy;
    psddot  = phidot*thdot*(Ix-Iy)/Iz                      + tau_ps/Iz;

    % Assemble state derivative
    sdot = [xdot; ydot; zdot; phidot; thdot; psdot; ...
            xddot; yddot; zddot; phiddot; thddot; psddot];
end