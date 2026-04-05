%% Quadcopter Nonlinear Dynamics — Thanh 2022, Eq. (6)
% s (12x1): [x,y,z, phi,theta,psi, xdot,ydot,zdot, phidot,thetadot,psidot]
% u  (4x1): [fz, tau_phi, tau_th, tau_ps]
% Returns sdot (12x1).

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

    Om = 0;  % net rotor speed (gyroscopic terms neglected)

    % Translational accelerations (ZYX Euler, body thrust mapped to inertial frame)
    xddot = (cos(phi)*sin(th)*cos(ps) + sin(phi)*sin(ps)) * fz/m;
    yddot = (cos(phi)*sin(th)*sin(ps) - sin(phi)*cos(ps)) * fz/m;
    zddot = -g + cos(phi)*cos(th) * fz/m;

    % Rotational accelerations (Euler rigid-body equations)
    phiddot = thdot*psdot*(Iy-Iz)/Ix  - Jr/Ix * thdot*Om  + tau_phi/Ix;
    thddot  = phidot*psdot*(Iz-Ix)/Iy + Jr/Iy * phidot*Om + tau_th/Iy;
    psddot  = phidot*thdot*(Ix-Iy)/Iz                      + tau_ps/Iz;

    % Assemble state derivative
    sdot = [xdot; ydot; zdot; phidot; thdot; psdot; ...
            xddot; yddot; zddot; phiddot; thddot; psddot];
end