clear; clc; close all;
%% Quadcopter UAV Dynamics
% from Thanh 2022
% State: [x y z phi theta psi xdot ydot zdot phidot thetadot psidot]
% Input: [fz tau_phi tau_th tau_ps]
%% Model Parameter
m  = 1.12;        % mass [kg]
g  = 9.81;        % gravity [m/s^2]
l  = 0.23;        % arm length [m]
Jr = 8.5e-4;      % rotor inertia [kg.m^2]
Ix = 0.0019;
Iy = 0.0019;
Iz = 0.0223;
b  = 7.73212e-6;  % thrust coeff
d  = 1.27513e-7;  % drag coeff

%% Dynamics Function 
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

% Residual rotor speed (0 for simplification)
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


