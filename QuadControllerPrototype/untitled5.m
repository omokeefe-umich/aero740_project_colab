%% Quadcopter UAV Dynamics
% from Thanh 2022
% State: [x y z phi theta psi xdot ydot zdot phidot thetadot psidot]
% Input: [fz tau_phi tau_th tau_ps]

% Augmented State:
% [delta_x delta_y delta_z delta_phi delta_theta delta_psi...
%  delta_xdot delta_ydot delta_zdot delta_phidot delta_thetadot delta_psidot...
%  x y z phi theta psi xdot ydot zdot phidot thetadot psidot fz tau_phi tau_th tau_ps]
clear all;
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
tspan = 0:0.01:2;
deltax = zeros(12, 1);
x = zeros(12, 1);
u = zeros(4, 1);
Xs = [];

C = [diag([1; 1; 1]), zeros(3, 9), diag([1; 1; 1]), zeros(3, 13)];
Q = C' * C + blkdiag(zeros(3), eye(9) * 0.01, zeros(16));

for k = 1 : 201
    Ac = Jacobian(@(x) dynamics(0, x, u, params), x);
    Bc = Jacobian(@(u) dynamics(0, x + deltax, u, params), u);
    [Ads, Bds] = Discretize_dt([0.01], ones(20, 1), Ac, Bc);
    Ads_augmented = zeros(20, 28, 28);
    Bds_augmented = zeros(20, 28, 4);
    
    for i = 1:20
        Ads_augmented(i, :, :) = [squeeze(Ads(i, :, :)), zeros(12, 16); [eye(12); zeros(4, 12)], eye(16)];
        Bds_augmented(i, :, :) = [squeeze(Bds(i, :, :)); zeros(12, 4); eye(4)];
    end
    [S, M] = StackedMatrix(Ads_augmented, Bds_augmented);
    
    x_aug = [deltax; x; u];
    x_traj = [zeros(20, 12), cos((tspan(k) + 0.01):0.01:(tspan(k) + 0.2))' - 1, sin((tspan(k) + 0.01):0.01:(tspan(k) + 0.2))', zeros(20, 14)]';
    x_traj = x_traj(:);
    
    [H, q] = QPFormat(Q, diag([0.001; 0.001; 0.001; 0.001]), Q * 0.01,...
        S, M, (0:0.01:0.2)', x_aug, x_traj, zeros(80, 1));
    
    u_mpc = quadprog((H + H') / 2, q);
    u = u + u_mpc(1:4);
    x = x + deltax;
    [T, X] = ode45(@(t, x) dynamics(t, x, u, params), 0:0.00001:0.01, x);
    Xs = [Xs, X(:, 1:3)'];
    deltax = X(end, :)' - x;
end

plot(Xs(1, :), Xs(2, :));
hold on
plot(cos(0:0.01:2) - 1, sin(0:0.01:2));