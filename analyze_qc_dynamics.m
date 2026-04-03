%  Compute the linearized EOM for our drone using symbolic differentiation
clear; clc;

% Symbolic variables
syms x y z phi th ps xdot ydot zdot phidot thdot psdot real
syms fz tau_phi tau_th tau_ps real
syms m g Jr Ix Iy Iz Om real

% State and input vectors
s = [x; y; z; phi; th; ps; xdot; ydot; zdot; phidot; thdot; psdot];
u = [fz; tau_phi; tau_th; tau_ps];

% Nonlinear dynamics
xddot = (cos(phi)*sin(th)*cos(ps) + sin(phi)*sin(ps)) * fz/m;
yddot = (cos(phi)*sin(th)*sin(ps) - sin(phi)*cos(ps)) * fz/m;
zddot = -g + cos(phi)*cos(th) * fz/m;

phiddot = thdot*psdot*(Iy-Iz)/Ix - Jr/Ix * thdot*Om + tau_phi/Ix;
thddot  = phidot*psdot*(Iz-Ix)/Iy + Jr/Iy * phidot*Om + tau_th/Iy;
psddot  = phidot*thdot*(Ix-Iy)/Iz + tau_ps/Iz;

sdot = [xdot; ydot; zdot; phidot; thdot; psdot; ...
        xddot; yddot; zddot; phiddot; thddot; psddot];

% Compute their Jacobians
A_sym = jacobian(sdot, s);
B_sym = jacobian(sdot, u);

disp('A = ');
pretty(A_sym)

disp('B = ');
pretty(B_sym)


% Assess poles and zeros present in the transfer function associated with
% our system
s0 = zeros(12,1);

u0 = [m*g; 0; 0; 0];   % symbolic still, because m and g are symbolic here

% Substitute hover point
subs_vars = [x y z phi th ps xdot ydot zdot phidot thdot psdot ...
             fz tau_phi tau_th tau_ps Om];
subs_vals = [0 0 0 0 0 0 0 0 0 0 0 0 ...
             m*g 0 0 0 0];

A_hover_sym = simplify(subs(A_sym, subs_vars, subs_vals));
B_hover_sym = simplify(subs(B_sym, subs_vars, subs_vals));

%% Numerical parameters
m_val  = 1.5;
g_val  = 9.81;
Jr_val = 0.02;
Ix_val = 0.02;
Iy_val = 0.025;
Iz_val = 0.03;

param_vars = [m g Jr Ix Iy Iz];
param_vals = [m_val g_val Jr_val Ix_val Iy_val Iz_val];

A_num = double(subs(A_hover_sym, param_vars, param_vals));
B_num = double(subs(B_hover_sym, param_vars, param_vals));

disp(A_num)
disp(B_num)