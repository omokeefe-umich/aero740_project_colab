%% Dynamics Function
% State: [x y z phi theta psi xdot ydot zdot phidot thetadot psidot]
% Input: [fz tau_phi tau_th tau_ps]
function sdot = dynamics(~, s, u, params, wind_data)
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

    % Unpack the x and z states 
    x = s(1); z = s(3);

    % Use Lookup Table to get the disturbance values (wind velocity at
    % the quadcopters specific x and z coordinate)
    vw_x = interp2(wind_data.X,wind_data.Z,wind_data.U,x,z,"linear",0);
    vw_z = interp2(wind_data.X,wind_data.Z,wind_data.W,x,z,"linear",0);

    % Calculate the acceleration of the wind drag force ( F = ma )
    % D is a lumped paramater drag constant (D = 0.38 seems to be a
    %  moderate disturbance when running the "FixedTarget.m" file ) 
    D = 0.4;
    ax_wind = (D * (vw_x - xdot)) / m;
    az_wind = (D * (vw_z - zdot)) / m;
    
    
    %Translational accelerations (Eq 6)
    % Add the wind acceleration in the x direction to the xddot eq.
    xddot = (cos(phi)*sin(th)*cos(ps) + sin(phi)*sin(ps)) * fz/m + ax_wind;
    yddot = (cos(phi)*sin(th)*sin(ps) - sin(phi)*cos(ps)) * fz/m;
    % Add the wind acceleration in the z direction to the zddot eq.
    zddot = -g + cos(phi)*cos(th) * fz/m + az_wind;

    %Rotational accelerations (Eq 6)
    phiddot = thdot*psdot*(Iy-Iz)/Ix - Jr/Ix * thdot*Om + tau_phi/Ix;
    thddot  = phidot*psdot*(Iz-Ix)/Iy + Jr/Iy * phidot*Om + tau_th/Iy;
    psddot  = phidot*thdot*(Ix-Iy)/Iz + tau_ps/Iz;
    sdot = [xdot; ydot; zdot; phidot; thdot; psdot;
            xddot; yddot; zddot; phiddot; thddot; psddot];
end