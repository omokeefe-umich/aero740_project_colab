clc; clear; close all;
%% Model Parameter
m  = 1.12; g  = 9.81; l  = 0.23;
Jr = 8.5e-4; Ix = 0.0019; Iy = 0.0019; Iz = 0.0223;
b  = 7.73212e-6; d  = 1.27513e-7;
params = [m; g; Jr; Ix; Iy; Iz];
tspan = 0:0.01:2;
deltax = zeros(12, 1);
x = zeros(12, 1);
u = zeros(4, 1);
Xs = [];
C = [ ...
    diag([1;1;1]), zeros(3,9), ...   % delta x (first 12)
    diag([1;1;1]), zeros(3,9), ...   % x_tilde (next 12)
    zeros(3,16) ...                      % x and u
    ];
Q = C' * C + blkdiag(zeros(3), eye(9)*0.001, zeros(28));
for k = 1 : 201
    Ac = Jacobian(@(xlin) dynamics(0, xlin + deltax, u, params), x);
    Bc = Jacobian(@(ulin) dynamics(0, x + deltax, ulin, params), u);
    [Ads, Bds] = Discretize_dt([0.01], ones(20,1), Ac, Bc);
    Ads_augmented = zeros(20, 40, 40);
    Bds_augmented = zeros(20, 40, 4);
    for i = 1:20
        Ad = squeeze(Ads(i,:,:));
        Bd = squeeze(Bds(i,:,:));
        Ads_augmented(i,:,:) = [
            Ad,            zeros(12,12), zeros(12,12), zeros(12,4);
            zeros(12,12),  Ad,           zeros(12,12), zeros(12,4);
            eye(12),       zeros(12,12), eye(12),      zeros(12,4);
            zeros(4,12),   zeros(4,12),  zeros(4,12),  eye(4)
            ];
        Bds_augmented(i,:,:) = [
            Bd;
            Bd;
            zeros(12,4);
            eye(4)
            ];
        end
    [S, M] = StackedMatrix(Ads_augmented, Bds_augmented);
    tk = tspan(k);
    x_ref = zeros(12,1);
    x_ref(1) = cos(tk) - 1;
    x_ref(2) = 2*sin(tk);
    x_ref(3) = tk;
    x_tilde = x - x_ref;
    x_aug = [deltax; x_tilde; x; u];
    tH = (tspan(k) + 0.01):0.01:(tspan(k) + 0.2);
    xRefH = cos(tH)' - 1;
    yRefH = 2*sin(tH)';
    zRefH = tH';
    x_traj_mat = [
        zeros(20,12), ...
        zeros(20,12), ...
        xRefH, yRefH, zRefH, zeros(20,9), ...
        zeros(20,4)
        ];
    x_traj = x_traj_mat';
    x_traj = x_traj(:);
    [H, q] = QPFormat(Q, diag([0.001;0.001;0.001;0.001]), Q*0.01,...
        S, M, (0:0.01:0.2)', x_aug, x_traj, zeros(80,1));
    u_mpc = quadprog((H + H') / 2, q);
    u = u + u_mpc(1:4);
    x = x + deltax;
    [T, X] = ode45(@(t, x) dynamics(t, x, u, params), 0:0.00001:0.01, x);
    Xs = [Xs, X(1:end - 1, 1:3)'];
    deltax = X(end, :)' - x;
end
plot(Xs(1, :));
figure;
plot(Xs(1, :), Xs(2, :));
hold on
plot(cos(0:0.01:2.01) - 1, 2 * sin(0:0.01:2.01));