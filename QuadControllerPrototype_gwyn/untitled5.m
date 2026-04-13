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
x_gate = 18;        % m

C = [diag([1; 1; 1]), zeros(3, 9), diag([1; 1; 1]), zeros(3, 13)];
Q = C' * C + blkdiag(zeros(3), eye(9) * 0.001, zeros(16));

% Target Vehicle trajectory
x_start_tgt = 10; % initial x position of target vehicle (m)
y_start_tgt = 0; 
z_start_tgt = 0; 
r = [x_start_tgt; y_start_tgt; z_start_tgt; zeros(9,1)]; 
e = x - r;  % initial tracking error

%% Bounds
% x (12x1): [x,y,z, phi,theta,psi, xdot,ydot,zdot, phidot,thetadot,psidot]
% u  (4x1): [fz, tau_phi, tau_th, tau_ps]
anglim      = deg2rad(45);   % rad
angratelim  = deg2rad(90);   % rad/s
torquelim   = 7;            % N-M
speedlim    = 45;            % m/s
Fzlim       = 75;            % N  > 5:1 T:W
accellim    = 5;             % m/s/s

xub = [inf*ones(3,1);    anglim*ones(3,1);   speedlim*ones(3,1);   angratelim*ones(3,1)];
xlb = [-inf*ones(3,1);  -anglim*ones(3,1);  -speedlim*ones(3,1);  -angratelim*ones(3,1)];
eub = inf*ones(12,1);
elb = -inf*ones(12,1);
uub = [ Fzlim;    torquelim;  torquelim;  (1/4)*torquelim];
ulb = [-Fzlim;   -torquelim; -torquelim; -(1/4)*torquelim];
dxub = [ speedlim*ones(3,1);    angratelim*ones(3,1);      accellim*ones(3,1);      (1/4)*angratelim*ones(3,1)];
dxlb = [-speedlim*ones(3,1);   -angratelim*ones(3,1);     -accellim*ones(3,1);     -(1/4)*angratelim*ones(3,1)];

Xub = [dxub; eub; xub; uub];
Xlb = [dxlb; elb; xlb; ulb];

duub = [(2/3)*Fzlim;    (1/2)*torquelim;    (1/2)*torquelim;  (1/2)*torquelim];
dulb = [-(2/3)*Fzlim;  -(1/2)*torquelim;   -(1/2)*torquelim; -(1/2)*torquelim];

% 5:1 weight ratio for performance drone -> 6kg thrust -> ~60N max thrust
% 4 motors -> 15 N force per motor
% 1m across drone -> 0.5 m moment arm -> 7.5N torque limit
% 
% https://wrekd.com/products/t-motor-v3115-3115-900kv-motor?_pos=6&_fid=1005fc1f0&_ss=c
% https://binsfeld.com/power-torque-speed-conversion-calculator/


for k = 1 : 201
    Ac = Jacobian(@(x) dynamics(0, x + deltax, u, params), x);
    Bc = Jacobian(@(u) dynamics(0, x + deltax, u, params), u);
    [Ads, Bds] = Discretize_dt([0.01], ones(20, 1), Ac, Bc);
    Ads_augmented = zeros(20, 28, 28);
    Bds_augmented = zeros(20, 28, 4);

    dt_impact = (x_gate - r(1)) / r(7); % delta time until ground vehicle impacts gate
    
    Axs_augmented = zeros(N_mpc, 40, 40);
    Bxs_augmented = zeros(N_mpc, 40, 1);
    for i = 1:20
        Ads_augmented(i, :, :) = [squeeze(Ads(i, :, :)), zeros(12, 16); [eye(12); zeros(4, 12)], eye(16)];
        Bds_augmented(i, :, :) = [squeeze(Bds(i, :, :)); zeros(12, 4); eye(4)];

        % Constraints
        % Constrain the altitude of the quadcopter to be above ground and the
        % velocities of the quadcopter to be sufficiently high as to intercept the
        % ground vehicle before it impacts the gate. 
        Ax_i = [1,  0,  0, 0, 0, 0, dt_impact,   0,          0,         0, 0, 0;
                0, -1,  0, 0, 0, 0, 0,           -dt_impact, 0,         0, 0, 0;
                0,  0, -1, 0, 0, 0, 0,           0,         -dt_impact, 0, 0, 0;
                zeros(9,12)];
        Axs_augmented(i, :, :) = [zeros(12,12),     zeros(12,12),  zeros(12,4);
                                  zeros(12,12),     Ax_i,          zeros(12,4);
                                  zeros(4,12),     zeros(4,12),    zeros(4,4);
                                ];
        Bxs_augmented(i, :, :) = [zeros(12,1); 
                                  x_gate;
                                  -(r(2) + r(8)*dt_impact);
                                  -(r(3) + r(9)*dt_impact);
                                  zeros(9,1);
                                  zeros(4,1)];
    end
    [S, M] = StackedMatrix(Ads_augmented, Bds_augmented);
    
    x_aug = [deltax; x; u];
    %x_traj = [zeros(20, 12), cos((tspan(k) + 0.01):0.01:(tspan(k) + 0.2))' - 1, 2 * sin((tspan(k) + 0.01):0.01:(tspan(k) + 0.2))', ((tspan(k) + 0.01):0.01:(tspan(k) + 0.2))', zeros(20, 13)]';
    x_traj = [zeros(20, 12), ones(20, 1), zeros(20, 15)]';
    x_traj = x_traj(:);
    P = Q * 0.01;
    [H, q] = QPFormat(Q, diag([0.001; 0.001; 0.001; 0.001]), P,...
        S, M, (0:0.01:0.2)', x_aug, x_traj, zeros(80, 1));
    
    u_mpc = quadprog((H + H') / 2, q);
    u = u + u_mpc(1:4);
    x = x + deltax;
    [T, X] = ode45(@(t, x) dynamics(t, x, u, params), 0:0.00001:0.01, x);
    Xs = [Xs, X(1:end - 1, 1:3)'];
    deltax = X(end, :)' - x;
    
    % r(1) = r(1) + 0.01 * dt;
end

plot(Xs(1, :));

%plot(Xs(1, :), Xs(2, :));
%hold on
%plot(cos(0:0.01:2.01) - 1, 2 * sin(0:0.01:2.01));