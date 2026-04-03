%% Dynamics Function
% State: [x y z phi theta psi xdot ydot zdot phidot thetadot psidot]
% Input: [fz tau_phi tau_th tau_ps]
function sdot = dynamics(~, s, u, params)
    m  = params(1); g  = params(2); Jr = params(3);
    Ix = params(4); Iy = params(5); Iz = params(6);
    xdot    = s(7);  ydot   = s(8);  zdot   = s(9);
    phi     = s(4);  th     = s(5);  ps     = s(6);
    phidot  = s(10); thdot  = s(11); psdot  = s(12);
    fz      = u(1);
    tau_phi = u(2);
    tau_th  = u(3);
    tau_ps  = u(4);
    Om = 0;
    % Translational accelerations (Eq 6)
    xddot = (cos(phi)*sin(th)*cos(ps) + sin(phi)*sin(ps)) * fz/m;
    yddot = (cos(phi)*sin(th)*sin(ps) - sin(phi)*cos(ps)) * fz/m;
    zddot = -g + cos(phi)*cos(th) * fz/m;
    % Rotational accelerations (Eq 6)
    phiddot = thdot*psdot*(Iy-Iz)/Ix - Jr/Ix * thdot*Om + tau_phi/Ix;
    thddot  = phidot*psdot*(Iz-Ix)/Iy + Jr/Iy * phidot*Om + tau_th/Iy;
    psddot  = phidot*thdot*(Ix-Iy)/Iz + tau_ps/Iz;
    sdot = [xdot; ydot; zdot; phidot; thdot; psdot;
            xddot; yddot; zddot; phiddot; thddot; psddot];
end