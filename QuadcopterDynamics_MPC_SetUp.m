clear; clc; close all;
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

params = [m; g; Jr; Ix; Iy; Iz];

%% Simulation
tspan = [0 10];
s0    = zeros(12,1);
input_mag = [12; 0.001; 0.001; 0.001];

control_inputs = {'fz', 'tau phi', 'tau theta', 'tau psi'};

for i = 1:4

%fixed ~i = 0 by replacing it with u = zeros(4,1) and setting u(i) = 12 to correctly isolate one input at a time. 
    u      = zeros(4,1);
    u(i) = input_mag(i); 

    % free-fall (fz=0 => zddot=-g)
    if i ~= 1
        u(1) = m * g;
    end

    % Simulate
    [t, s] = ode45(@(t,s) dynamics(t, s, u, params), tspan, s0);

    % Plot
    figure;

    subplot(4,1,1);
    plot(t, s(:,1:3));
    xlabel('Time (s)'); ylabel('Position (m)');
    legend('x','y','z');
    title(sprintf('Control Input: %s = 12, fz = %.2f', control_inputs{i}, u(1)));
    grid on;

    subplot(4,1,2);
    plot(t, rad2deg(s(:,4:6)));        % FIX 4: show in degrees
    xlabel('Time (s)'); ylabel('Orientation (deg)');
    legend('phi','theta','psi');
    grid on;

    subplot(4,1,3);
    plot(t, s(:,7:9));
    xlabel('Time (s)'); ylabel('Velocity (m/s)');
    legend('xdot','ydot','zdot');
    grid on;

    subplot(4,1,4);
    plot(t, rad2deg(s(:,10:12)));      % FIX 4: show in degrees
    xlabel('Time (s)'); ylabel('Angular Velocity (deg/s)');
    legend('phidot','thetadot','psidot');
    grid on;

    linkaxes(findall(gcf,'type','axes'), 'x');
    xlim(tspan);


end
%% MPC Setup 
% Design factors 
dt  = 0.1;              % discrete timestep [s]
T_H = 2.0;              % prediction horizon [s]
N   = round(T_H/dt);    % prediction steps 

%  state: [x z th xdot zdot thdot] 6 states
T_sim = 10;
nx = 6;  nu = 2; 

% Track
x0= zeros(nx,1);
xf = [50; 10; 0; 0; 0; 0]; % what we want to reach ? 

% Input Constraints 
u_min = [u1 ; u2 ];
u_max = [u1; u2];
z_bridge = 8.0;
rho = ;            % slack penalty

% Model mismatch (scenario parameter)
mismatch = 0.0;        % e.g., 0.1 = 10%
m_model  = m*(1-mismatch);

%% Linearize about hover + ZOH discretize
syms x1 x2 x3 x4 x5 x6 u1 u2 real
% State: [x z th xdot zdot thdot]
xs = [x1; x2; x3; x4; x5; x6];
us = [u1; u2];

fsys = [x4;
         x5;
         x6;
         -sin(x3)*u1/m_model;
         -g + cos(x3)*u1/m_model;
         u2/Iy];

% Jacobians
A_sym = jacobian(fsys, xs);
B_sym = jacobian(fsys, us);

% Evaluate at hover equilibrium: x=0, u=[mg; 0]
A_c = double(subs(A_sym, [xs; us], [zeros(nx,1); m_model*g; 0]));
B_c = double(subs(B_sym, [xs; us], [zeros(nx,1); m_model*g; 0]));

% Discretize
sys_d = c2d(ss(A_c, B_c, eye(nx), zeros(nx,nu)), dt, 'zoh');
A_d = sys_d.A;  B_d = sys_d.B;

%% Cost 
Q = ;
R = ;
[P,~,~] = dare(A_d, B_d, Q, R);



%% MPC Solver 



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
