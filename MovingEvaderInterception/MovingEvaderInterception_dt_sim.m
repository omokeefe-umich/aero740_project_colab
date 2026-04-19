%% Quadcopter UAV Dynamics with MPC pursuit of a moving evader
% from Thanh 2022
% State: [x y z phi theta psi xdot ydot zdot phidot thetadot psidot]
% Input: [fz tau_phi tau_th tau_ps]

% Augmented State:
% [x y z phi theta psi xdot ydot zdot phidot thetadot psidot 
%  xddot yddot zddot phiddot thetaddot psiddot fz tau_phi tau_th tau_ps]
%
% Lifted prediction state used in the QP:
%   z = [x; e; xddot; u], with e = x - x_ref
% This keeps physical dynamics in the compact 22-state model while adding
% explicit tracking-error coordinates to the optimization model.
clear; clc; close all;

%% Pursuit / evader setup
P0 = [0; 0];
E0 = [20; 0];
VE = 5;
VP = 6;
lambda = 50;
thetaE = pi / 2;
zTarget = 0.0;

dt = 0.04;
dt_sim = 0.04; 
horizonSteps = 40;
horizonTime = horizonSteps * dt;
simMaxTime = 60.0;
interceptRadius = 0.15;
leadTime = 0.0;
gateDistanceFromE0 = 16;

% The gate is a fixed point on the evader/ground-vehicle path. The evader
% moves from E0 toward this point, and the quadcopter must intercept the
% evader before the evader reaches the gate.
gateAxis = [cos(thetaE); sin(thetaE)];
gatePoint2D = E0 + gateDistanceFromE0 * gateAxis;
gateCoordinate = gateAxis' * gatePoint2D;
evaderInitialCoordinate = gateAxis' * E0;
gateTime = (gateCoordinate - evaderInitialCoordinate) / VE;
gatePoint = [gatePoint2D; zTarget];

%% Zero Wind Disturbance for Modelling
wind_data = DisturbanceModel(1); % 0 = no wind, 1 = wind on
no_wind = wind_data ;
no_wind.U = wind_data.U*0 ;
no_wind.W = wind_data.W*0 ;

%% Proportional pursuit baseline
R0 = norm(E0 - P0);
beta0 = atan2(E0(2) - P0(2), E0(1) - P0(1));
theta0PP = beta0;

tspanPP = 0:dt:100;
optsPP = odeset('Events', @(t, x) stopFcn(t, x, interceptRadius), ...
    'RelTol', 1e-6, 'AbsTol', 1e-6);
x0PP = [R0; beta0; theta0PP];

[tPP, xPP] = ode45(@(t, x) RB_PP_straight(t, x, VP, VE, lambda, thetaE), ...
    tspanPP, x0PP, optsPP);

RPP = xPP(:, 1);
betaPP = xPP(:, 2);
evaderPP = evaderPosition(tPP, E0, VE, thetaE, zTarget);
xEPP = evaderPP(:, 1);
yEPP = evaderPP(:, 2);
xPPP = xEPP - RPP .* cos(betaPP);
yPPP = yEPP - RPP .* sin(betaPP);

fprintf('PP lambda = %g -> tf = %.2f s, Rf = %.4f m\n', ...
    lambda, tPP(end), RPP(end));

%% Model parameters
m  = 1.12;        % mass [kg]
g  = 9.81;        % gravity [m/s^2]
Jr = 8.5e-4;      % rotor inertia [kg.m^2]
Ix = 0.0019;
Iy = 0.0019;
Iz = 0.0223;
params = [m; g; Jr; Ix; Iy; Iz];

nX = 12;
nA = 6;
nU = 4;
nAug = nX + nA + nU;      % [x; a; u]
nLift = nX + nX + nA + nU; % [x; e; a; u]

idxAug = struct('x', 1:nX, 'a', nX + (1:nA), 'u', nX + nA + (1:nU));
idxLift = struct( ...
    'x', 1:nX, ...
    'e', nX + (1:nX), ...
    'a', 2 * nX + (1:nA), ...
    'u', 2 * nX + nA + (1:nU));

x = zeros(12, 1);
x(1:2) = P0;
x(3) = 2;
u_hover = [m * g; 0; 0; 0];
u = u_hover;
% Assuming accurate acceleration measurement...
a = dynamics(0, x, u, params, wind_data);
a = a(7:12);

%% Bounds / constraints
% Gate reachability constraint:
%   gate < quad_position + quad_velocity * time_remaining
% Projected along the evader motion axis, this is a linear inequality in
% the predicted quad position and velocity.
%
% Build lifted mapping z_k = T_step * x_aug_k + d_k,
% where z_k = [x_k; e_k; a_k; u_k], e_k = x_k - x_ref_k.
T_step = zeros(nLift, nAug);
T_step(idxLift.x, idxAug.x) = eye(nX);
T_step(idxLift.e, idxAug.x) = eye(nX);
T_step(idxLift.a, idxAug.a) = eye(nA);
T_step(idxLift.u, idxAug.u) = eye(nU);
T_stacked = kron(eye(horizonSteps), T_step);

torquelim   = 0.4;
anglim      = deg2rad(70);

uub = [3.0 * m * g;  torquelim;  torquelim;  0.25 * torquelim];
ulb = [0.5 * m * g; -torquelim; -torquelim; -0.25 * torquelim];

duub = [4.0;  0.15;  0.15;  0.05] / dt;
dulb = [-4.0; -0.15; -0.15; -0.05] / dt;

duStackUB = repmat(duub, horizonSteps, 1);
duStackLB = repmat(dulb, horizonSteps, 1);

% Per-step hard constraints on the augmented predicted state:
% z >= 0, bounded attitude, and bounded physical input.
A_z = zeros(1, nLift);
A_z(3) = -1;
A_ang = zeros(3, nLift);
A_ang(:, idxLift.x(4:6)) = eye(3);
A_u = zeros(4, nLift);
A_u(:, idxLift.u) = eye(4);
A_step = [A_z; A_ang; -A_ang; A_u; -A_u];
b_step = [0; anglim * ones(3, 1); anglim * ones(3, 1); uub; -ulb];

stateLower = nan(12, 1);
stateUpper = nan(12, 1);
stateLower(3) = 0;
stateLower(4:6) = -anglim;
stateUpper(4:6) = anglim;
controlLower = ulb;
controlUpper = uub;
plotConstraints = struct( ...
    'stateLower', stateLower, ...
    'stateUpper', stateUpper, ...
    'controlLower', controlLower, ...
    'controlUpper', controlUpper);

qpOptions = optimoptions('quadprog', ...
         'Display', 'off', ...
         'Algorithm', 'interior-point-convex', ...
         'MaxIterations', 5000, ...
         'ConstraintTolerance', 1e-3, ...
         'OptimalityTolerance', 1e-3, ...
         'StepTolerance', 1e-8);

%% Penalty Matrices
Q = zeros(nLift);
Q(idxLift.e(1:3), idxLift.e(1:3)) = 0.1 * eye(3);
Q(idxLift.e(4:5), idxLift.e(4:5)) = 0.001 * eye(2);
Q(idxLift.e(7:9), idxLift.e(7:9)) = 0.01 * eye(3);
Q(idxLift.e(10:11), idxLift.e(10:11)) = 0.01 * eye(2);
Q(idxLift.e(12), idxLift.e(12)) = 0.001;
Q(idxLift.a(1:3), idxLift.a(1:3)) = 0.01 * eye(3);

P = zeros(nLift);
P(idxLift.e(1:3), idxLift.e(1:3)) = eye(3);
duub_norm = [4.0; 0.15; 0.15; 0.05] / dt;   % same as before, but now treat as 1.0 each
dulb_norm = -duub_norm;

% Scale R to reflect that all normalized rates are order-1:
R = diag([1e-4 / duub_norm(1)^2, ...
          1e-4 / duub_norm(2)^2, ...
          1e-4 / duub_norm(3)^2, ...
          1e-4 / duub_norm(4)^2]);

%% MPC pursuit of the moving evader
summedInput = u;
allXs = x;
times = 0;
eState = evaderPosition(0, E0, VE, thetaE, zTarget);
evaderState = eState(1,1:3);
inputHistory = u;

Xs = x(1:3);
Ts = 0;
mpcIntercepted = false;
maxSIMsteps = floor(simMaxTime / dt_sim);
wallTimes = zeros(maxSIMsteps,1);
for k = 1:maxSIMsteps
    currentTime = (k - 1) * dt_sim;

    J_xv = Jacobian(@(x) dynamics(0, x, u, params, no_wind), x);
    J_u  = Jacobian(@(u) dynamics(0, x, u, params, no_wind), u);
    Ac = [zeros(12, 6), eye(12), zeros(12, 4); zeros(6, 6), J_xv(7:12, :), zeros(6, 4); zeros(4, 22)];
    Bc = [zeros(12, 4); J_u(7:12, :); eye(4)];
    [Ads, Bds] = Discretize_dt(dt, ones(horizonSteps, 1), Ac, Bc);
    [S, M] = StackedMatrix(Ads, Bds);
    S_lift = T_stacked * S;
    M_lift = T_stacked * M;

    x_aug = [x; a; u];

    predictionTimes = (currentTime + leadTime + dt:dt:currentTime + leadTime + horizonTime)';
    targetXYZ = evaderPosition(predictionTimes, E0, VE, thetaE, zTarget);
    x_ref = [targetXYZ, zeros(horizonSteps, nX - 3)]';
    x_ref = x_ref(:);

    x_traj = zeros(horizonSteps * nLift, 1);
    d_lift = zeros(horizonSteps * nLift, 1);
    for i = 1:horizonSteps
        rowOffset = (i - 1) * nLift;
        refOffset = (i - 1) * nX;
        x_ref_i = x_ref(refOffset + (1:nX));

        x_traj(rowOffset + idxLift.x(1:3)) = x_ref_i(1:3);
        x_traj(rowOffset + idxLift.u) = u_hover;
        d_lift(rowOffset + idxLift.e) = -x_ref_i;
    end

    [H, q] = QPFormat(Q, R, P, S_lift, M_lift, (0:dt:horizonTime)', ...
        x_aug, x_traj, zeros(horizonSteps * 4, 1), d_lift);

    Hqp = (H + H') / 2;

    evaderNow = evaderPosition(currentTime, E0, VE, thetaE, zTarget);
    evaderCoordinate = gateAxis' * evaderNow(1, 1:2)';
    gateTimeRemainingNow = (gateCoordinate - evaderCoordinate) / VE;
    if gateTimeRemainingNow <= 0
        warning('Evader reached the gate before MPC intercept at t = %.2f s. Stopping MPC.', currentTime);
        break;
    end

    % Lifted state indices: [x(1:12); e(13:24); xddot(25:30); u(31:34)].
    Mx = M_lift * x_aug + d_lift;

    % Hard path constraints for each predicted augmented state.
    % This includes z >= 0, bounded attitude, and physical input bounds.
    A_state_stacked = kron(eye(horizonSteps), A_step);
    AineqBase = A_state_stacked * S_lift;
    bineqBase = repmat(b_step, horizonSteps, 1) - A_state_stacked * Mx;

    % Terminal gate reachability:
    % gate <= p_quad_N + v_quad_N * (time left after terminal prediction).
    terminalTime = currentTime + horizonTime;
    terminalTimeRemaining = gateTime - terminalTime;
    AineqGate = AineqBase;
    bineqGate = bineqBase;
    gateActive = false;
    if terminalTimeRemaining > 0
        rowIdxN = (horizonSteps - 1) * nLift + 1:horizonSteps * nLift;
        S_N = S_lift(rowIdxN, :);
        M_N = M_lift(rowIdxN, :);
        d_N = d_lift(rowIdxN);
        A_gate_N = zeros(1, nLift);
        A_gate_N(idxLift.x(1:2)) = -gateAxis';
        A_gate_N(idxLift.x(7:8)) = -terminalTimeRemaining * gateAxis';
        b_gate_N = -gateCoordinate;
        AineqGate = [A_gate_N * S_N; AineqBase];
        bineqGate = [b_gate_N - A_gate_N * (M_N * x_aug + d_N); bineqBase];
        gateActive = true;
    end

    tnow = tic; 
    if gateActive
        [u_mpc, ~, exitflag] = quadprog(Hqp, q, AineqGate, bineqGate, [], [], ...
            duStackLB, duStackUB, [], qpOptions);
    else
        [u_mpc, ~, exitflag] = quadprog(Hqp, q, AineqBase, bineqBase, [], [], ...
            duStackLB, duStackUB, [], qpOptions);
    end

    if exitflag <= 0 || isempty(u_mpc) || any(~isfinite(u_mpc))
        gateExitflag = exitflag;
        if gateActive
            [u_mpc, ~, exitflag] = quadprog(Hqp, q, AineqBase, bineqBase, [], [], ...
                duStackLB, duStackUB, [], qpOptions);
            if exitflag > 0 && ~isempty(u_mpc) && all(isfinite(u_mpc))
                fprintf(['Gate reachability rows infeasible at t = %.2f s ' ...
                    '(gate exitflag %d); using z-constrained tracking QP this step.\n'], ...
                    currentTime, gateExitflag);
                gateActive = false;
            end
        end
    end

    if exitflag <= 0 || isempty(u_mpc) || any(~isfinite(u_mpc))
        fprintf('QP fail debug: t = %.2f, z = %.3f, zdot = %.3f, u = [%.3f %.3f %.3f %.3f], gateActive = %d\n', ...
            currentTime, x(3), x(9), u(1), u(2), u(3), u(4), gateActive);
        warning('Constrained MPC quadprog failed at t = %.2f s with exitflag %d. Stopping MPC.', ...
            currentTime, exitflag);
        break;
    end
    wallTimes(k) = toc(tnow);
    u = u + u_mpc(1:4) * dt_sim;
    u = min(max(u, ulb), uub);
    [T, X] = ode45(@(t, state) dynamics(t, state, u, params, wind_data), 0:dt_sim/10:dt_sim, x);
    Xs = [Xs, X(2:end, 1:3)'];
    Ts = [Ts; currentTime + T(2:end)];
    a = (X(end, 7:12)' - x(7:12)) / dt_sim;
    x = X(end, :)';

    allXs = [allXs, X(end,:)'];
    inputHistory = [inputHistory, u];
    times = [times; times(end)+dt_sim];
    evaderState = [evaderState; evaderNow(1, 1:3)];

    evaderNow = evaderPosition(Ts(end), E0, VE, thetaE, zTarget);
    rangeNow = norm(Xs(1:3, end)' - evaderNow(1, 1:3));
    if rangeNow <= interceptRadius
        mpcIntercepted = true;
        numstepssaturated = sum(sum(inputHistory - ulb < 1e-2 | inputHistory - uub > -1e-2));
        time_at_saturation = dt_sim * numstepssaturated;
        fprintf('MPC intercepted at t = %.2f s. Range = %.2f m. Distance from gate = %.1f m. # of Steps = %i\n', Ts(end), rangeNow, (gateCoordinate - evaderCoordinate), k);
        fprintf('    input cost = %.2f Fz, %.5f tau_phi, %.5f tau_theta, %.5f tau_psi\n', u(1), u(2), u(3), u(4));
        fprintf('    time at saturation = %.2f s\n', time_at_saturation);
        fprintf('    average wall-clock-time = %.2f s\n', mean(wallTimes(wallTimes > 0)))
        break;
    end

    if mod(k, round(1 / dt_sim)) == 0 || k == maxSIMsteps
        fprintf('MPC progress: %.2f / %.2f s, range = %.2f m, gate active = %d\n', ...
            Ts(end), simMaxTime, rangeNow, gateActive);
    end
end

evaderMPC = evaderPosition(Ts, E0, VE, thetaE, zTarget);

rangeMPCFinal = norm(Xs(1:3, end)' - evaderMPC(end, 1:3));
if ~mpcIntercepted
    fprintf('MPC stopped at t = %.2f s without intercept, final range = %.4f m\n', ...
        Ts(end), rangeMPCFinal);
end

sameTimeEnd = min(Ts(end), tPP(end));
sameTimePP = tPP <= sameTimeEnd;
tEvaderPlot = (0:dt_sim:max(tPP(end), Ts(end)))';
evaderPlot = evaderPosition(tEvaderPlot, E0, VE, thetaE, zTarget);
evaderPlotCoarse = evaderPosition(times, E0, VE, thetaE, zTarget)';

%% Plot the results
comparisonData = struct( ...
    'evaderPlot', evaderPlot, ...
    'ppPath', [xPPP, yPPP], ...
    'ppIntercept', [xPPP(end), yPPP(end)], ...
    'gatePoint2D', gatePoint2D, ...
    'mpcIntercepted', mpcIntercepted);

PlotResults(times, allXs, evaderPlotCoarse, inputHistory, plotConstraints, wind_data, comparisonData)

%% Range comparison
figure;
hold on; grid on; box on;
rangeMPC = sqrt(sum((Xs(1:2, :)' - evaderMPC(:, 1:2)).^2, 2));
plot(tPP, RPP, 'b-', 'LineWidth', 2, 'DisplayName', 'Proportional pursuit');
plot(Ts, rangeMPC, 'r-', 'LineWidth', 2, 'DisplayName', 'MPC quadcopter');
xline(gateTime, 'k:', 'Gate time', 'DisplayName', 'Gate time');
xlabel('t [s]');
ylabel('planar range to evader [m]');
title('Range to Moving Evader');
legend('Location', 'best');

function targetXYZ = evaderPosition(t, E0, VE, thetaE, zTarget)
    t = t(:);
    targetXYZ = [E0(1) + VE * cos(thetaE) * t, ...
        E0(2) + VE * sin(thetaE) * t, ...
        zTarget * ones(size(t))];
end

function dx = RB_PP_straight(~, x, VP, VE, lambda, thetaE)
    R = x(1);
    beta = x(2);
    theta = x(3);

    dB = (-VE * sin(beta - thetaE) + VP * sin(beta - theta)) / R;
    dTheta = lambda * dB;
    dR = VE * cos(beta - thetaE) - VP * cos(beta - theta);

    dx = [dR; dB; dTheta];
end

function [value, isterminal, direction] = stopFcn(~, x, interceptRadius)
    value = x(1) - interceptRadius;
    isterminal = 1;
    direction = -1;
end
