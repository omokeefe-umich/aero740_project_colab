%% Quadcopter UAV Dynamics with MPC pursuit of a moving evader
% Nominal model vs mass-mismatch true plant
% from Thanh 2022
% State: [x y z phi theta psi xdot ydot zdot phidot thetadot psidot]
% Input: [fz tau_phi tau_th tau_ps]
%
% Augmented State (22):
% [x y z phi theta psi xdot ydot zdot phidot thetadot psidot
%  xddot yddot zddot phiddot thetaddot psiddot fz tau_phi tau_th tau_ps]
%
% MPC always linearizes with params_model (nominal).
% Case 1: true plant = nominal  (params_model)
% Case 2: true plant = 10% heavier (params_true_heavy)

clear; clc; close all;

%% Pursuit / evader setup
P0 = [0; 0];
E0 = [20; 0];
VE = 5;
VP = 6;
lambda = 15;
thetaE = pi / 2;
zTarget = 0.0;

dt = 0.04;
horizonSteps = 40;
horizonTime = horizonSteps * dt;
mpcMaxTime = 60.0;
interceptRadius = 0.15;
leadTime = 0.0;

gateAxis = [cos(thetaE); sin(thetaE)];
gateDistanceFromE0 = 50;
gatePoint2D = E0 + gateDistanceFromE0 * gateAxis;
gateCoordinate = gateAxis' * gatePoint2D;
evaderInitialCoordinate = gateAxis' * E0;
gateTime = (gateCoordinate - evaderInitialCoordinate) / VE;
gatePoint = [gatePoint2D; zTarget];

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

figure(1);
hold on; axis equal; grid on; box on;
plot(xEPP, yEPP, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Evader');
plot(xPPP, yPPP, 'b-', 'LineWidth', 2, 'DisplayName', 'Proportional pursuit');
plot(P0(1), P0(2), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'P_0');
plot(E0(1), E0(2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'E_0');
plot(gatePoint2D(1), gatePoint2D(2), 'md', ...
    'MarkerFaceColor', 'm', 'DisplayName', 'Gate');
plot(xPPP(end), yPPP(end), 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'PP intercept');
xlabel('x [m]'); ylabel('y [m]');
title('Moving Evader Pursuit: PP vs MPC');
legend('Location', 'best');

%% Model parameters
m  = 1.12;
g  = 9.81;
Jr = 8.5e-4;
Ix = 0.0019;
Iy = 0.0019;
Iz = 0.0223;

mass_err = 0.08;   % % heavier true plant

params_model      = [m; g; Jr; Ix; Iy; Iz];                % MPC model (both cases)
params_true_nom   = [m; g; Jr; Ix; Iy; Iz];                % Case 1: nominal
params_true_heavy = [m*(1+mass_err); g; Jr; Ix; Iy; Iz];   % Case 2: mismatch

u_hover = [m * g; 0; 0; 0];

%% Cost matrices (22x22)
Q = blkdiag(0.1 * eye(2), 10, 0.001 * eye(2), zeros(4), eye(3) * 0.01, ...
    eye(2) * 0.01, 0.0001, eye(3) * 0.01, zeros(4));
P = blkdiag(eye(3), zeros(19));
R = eye(4) * 1e-4;

%% Bounds / constraints (22-dim indices)
torquelim = 0.4;
anglim    = deg2rad(70);

uub = [3.0 * m * g;  torquelim;  torquelim;  0.25 * torquelim];
ulb = [0.5 * m * g; -torquelim; -torquelim; -0.25 * torquelim];

duub = [4.0;  0.15;  0.15;  0.05] / dt;
dulb = [-4.0; -0.15; -0.15; -0.05] / dt;

duStackUB = repmat(duub, horizonSteps, 1);
duStackLB = repmat(dulb, horizonSteps, 1);

A_z = zeros(1, 22);   A_z(3) = -1;           % z >= 0
A_ang = zeros(3, 22);  A_ang(:, 4:6) = eye(3); % attitude
A_u = zeros(4, 22);    A_u(:, 19:22) = eye(4);  % input bounds

A_step = [A_z; A_ang; -A_ang; A_u; -A_u];
b_step = [0; anglim*ones(3,1); anglim*ones(3,1); uub; -ulb];

qpOptions = optimoptions('quadprog', ...
    'Display', 'off', ...
    'Algorithm', 'interior-point-convex', ...
    'MaxIterations', 5000, ...
    'ConstraintTolerance', 1e-3, ...
    'OptimalityTolerance', 1e-3, ...
    'StepTolerance', 1e-8);

%% ========================= Run both cases =========================
% params_true_list{1} = nominal, params_true_list{2} = heavy
params_true_list = {params_true_nom, params_true_heavy};
caseNames = {'Nominal', sprintf('%.0f%% mass mismatch', mass_err*100)};
nCases = length(params_true_list);

% Storage
XsCell  = cell(1, nCases);
TsCell  = cell(1, nCases);
intCell = false(1, nCases);

for iCase = 1:nCases
    fprintf('\n--- Running %s case ---\n', caseNames{iCase});
    params_true = params_true_list{iCase};

    x = zeros(12, 1);
    x(1:2) = P0;
    x(3) = 0;
    u = u_hover;
    sdot0 = dynamics(0, x, u, params_model);
    a = sdot0(7:12);

    Xs = x(1:3);
    Ts = 0;
    mpcIntercepted = false;
    maxMPCSteps = floor(mpcMaxTime / dt);

    for k = 1:maxMPCSteps
        currentTime = (k - 1) * dt;

        % Linearize with NOMINAL model always
        J_xv = Jacobian(@(xin) dynamics(0, xin, u, params_model), x);
        J_u  = Jacobian(@(uin) dynamics(0, x, uin, params_model), u);

        Ac = [zeros(12, 6), eye(12), zeros(12, 4); ...
              zeros(6, 6),  J_xv(7:12, :), zeros(6, 4); ...
              zeros(4, 22)];

        Bc = [zeros(12, 4); ...
              J_u(7:12, :); ...
              eye(4)];

        [Ads, Bds] = Discretize_dt(dt, ones(horizonSteps, 1), Ac, Bc);
        [S, M] = StackedMatrix(Ads, Bds);

        x_aug = [x; a; u];

        predictionTimes = (currentTime + leadTime + dt:dt:currentTime + leadTime + horizonTime)';
        targetXYZ = evaderPosition(predictionTimes, E0, VE, thetaE, zTarget);
        x_traj = [targetXYZ, zeros(horizonSteps, 15), ...
            repmat(u_hover', horizonSteps, 1)]';
        x_traj = x_traj(:);

        [H, q] = QPFormat(Q, R, P, S, M, (0:dt:horizonTime)', ...
            x_aug, x_traj, zeros(horizonSteps * 4, 1));
        Hqp = (H + H') / 2;

        evaderNow = evaderPosition(currentTime, E0, VE, thetaE, zTarget);
        evaderCoordinate = gateAxis' * evaderNow(1, 1:2)';
        gateTimeRemainingNow = (gateCoordinate - evaderCoordinate) / VE;
        if gateTimeRemainingNow <= 0
            warning('%s: evader reached gate at t = %.2f s.', caseNames{iCase}, currentTime);
            break;
        end

        Mx = M * x_aug;

        A_state_stacked = kron(eye(horizonSteps), A_step);
        AineqBase = A_state_stacked * S;
        bineqBase = repmat(b_step, horizonSteps, 1) - A_state_stacked * Mx;

        terminalTime = currentTime + horizonTime;
        terminalTimeRemaining = gateTime - terminalTime;
        AineqGate = AineqBase;
        bineqGate = bineqBase;
        gateActive = false;

        if terminalTimeRemaining > 0
            S_N = S((horizonSteps-1)*22+1 : horizonSteps*22, :);
            M_N = M((horizonSteps-1)*22+1 : horizonSteps*22, :);
            A_gate_N = zeros(1, 22);
            A_gate_N(1:2) = -gateAxis';
            A_gate_N(7:8) = -terminalTimeRemaining * gateAxis';
            b_gate_N = -gateCoordinate;
            AineqGate = [A_gate_N * S_N; AineqBase];
            bineqGate = [b_gate_N - A_gate_N * M_N * x_aug; bineqBase];
            gateActive = true;
        end

        if gateActive
            [u_mpc, ~, exitflag] = quadprog(Hqp, q, AineqGate, bineqGate, [], [], ...
                duStackLB, duStackUB, [], qpOptions);
        else
            [u_mpc, ~, exitflag] = quadprog(Hqp, q, AineqBase, bineqBase, [], [], ...
                duStackLB, duStackUB, [], qpOptions);
        end

        % Fallback: drop gate row if infeasible
        if exitflag <= 0 || isempty(u_mpc) || any(~isfinite(u_mpc))
            gateExitflag = exitflag;
            if gateActive
                [u_mpc, ~, exitflag] = quadprog(Hqp, q, AineqBase, bineqBase, [], [], ...
                    duStackLB, duStackUB, [], qpOptions);
                if exitflag > 0 && ~isempty(u_mpc) && all(isfinite(u_mpc))
                    fprintf('%s: gate infeasible at t = %.2f s (exitflag %d); using base QP.\n', ...
                        caseNames{iCase}, currentTime, gateExitflag);
                    gateActive = false;
                end
            end
        end

        if exitflag <= 0 || isempty(u_mpc) || any(~isfinite(u_mpc))
            fprintf('%s QP fail: t=%.2f z=%.3f zdot=%.3f u=[%.2f %.3f %.3f %.3f] gate=%d\n', ...
                caseNames{iCase}, currentTime, x(3), x(9), u(1), u(2), u(3), u(4), gateActive);
            warning('%s: quadprog failed at t = %.2f s, exitflag %d.', ...
                caseNames{iCase}, currentTime, exitflag);
            break;
        end

        u = u + u_mpc(1:4) * dt;
        u = min(max(u, ulb), uub);

        % Propagate with TRUE plant (nominal or heavy)
        [T, X] = ode45(@(t, state) dynamics(t, state, u, params_true), 0:dt/10:dt, x);
        Xs = [Xs, X(2:end, 1:3)'];
        Ts = [Ts; currentTime + T(2:end)];
        a = (X(end, 7:12)' - x(7:12)) / dt;
        x = X(end, :)';

        evaderNow = evaderPosition(Ts(end), E0, VE, thetaE, zTarget);
        rangeNow = norm(Xs(1:2, end)' - evaderNow(1, 1:2));
        if rangeNow <= interceptRadius
            mpcIntercepted = true;
            fprintf('%s intercepted at t = %.2f s, range = %.4f m\n', ...
                caseNames{iCase}, Ts(end), rangeNow);
            break;
        end

        if mod(k, round(1/dt)) == 0
            fprintf('%s: %.2f s, range = %.2f m, gate = %d\n', ...
                caseNames{iCase}, Ts(end), rangeNow, gateActive);
        end
    end

    XsCell{iCase} = Xs;
    TsCell{iCase} = Ts;
    intCell(iCase) = mpcIntercepted;
end

%% ========================= Results =========================
fprintf('\n========== Summary ==========\n');
for iCase = 1:nCases
    evMPC = evaderPosition(TsCell{iCase}, E0, VE, thetaE, zTarget);
    rng = sqrt(sum((XsCell{iCase}(1:2,:)' - evMPC(:,1:2)).^2, 2));
    fprintf('%s: tf = %.2f s, final range = %.4f m, intercepted = %d\n', ...
        caseNames{iCase}, TsCell{iCase}(end), rng(end), intCell(iCase));
end

%%
colors = {'r', 'm'};
styles = {'-', '--'};

% 2D trajectory
figure(1);
tEvaderPlot = (0:dt:max(cellfun(@(t) t(end), TsCell)))';
evaderPlot = evaderPosition(tEvaderPlot, E0, VE, thetaE, zTarget);
plot(evaderPlot(:,1), evaderPlot(:,2), 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
for iCase = 1:nCases
    plot(XsCell{iCase}(1,:), XsCell{iCase}(2,:), [colors{iCase} styles{iCase}], ...
        'LineWidth', 2, 'DisplayName', ['MPC ' caseNames{iCase}]);
    if intCell(iCase)
        plot(XsCell{iCase}(1,end), XsCell{iCase}(2,end), [colors{iCase} 's'], ...
            'MarkerFaceColor', colors{iCase}, 'MarkerSize', 10, ...
            'DisplayName', [caseNames{iCase} ' intercept']);
    end
end
legend('Location', 'best');

% 3D trajectory
figure;
hold on; grid on; box on;
evMPC1 = evaderPosition(TsCell{1}, E0, VE, thetaE, zTarget);
plot3(evMPC1(:,1), evMPC1(:,2), evMPC1(:,3), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Evader');
for iCase = 1:nCases
    plot3(XsCell{iCase}(1,:), XsCell{iCase}(2,:), XsCell{iCase}(3,:), ...
        [colors{iCase} styles{iCase}], 'LineWidth', 2, 'DisplayName', ['MPC ' caseNames{iCase}]);
end
plot3(P0(1), P0(2), 0, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'P_0');
plot3(E0(1), E0(2), zTarget, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'E_0');
plot3(gatePoint(1), gatePoint(2), gatePoint(3), 'md', 'MarkerFaceColor', 'm', 'DisplayName', 'Gate');
zlim([0 5])
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('3D MPC: Nominal vs Mass Mismatch');
legend('Location', 'best'); view(3);

% Range vs time
figure;
hold on; grid on; box on;
plot(tPP, RPP, 'b-', 'LineWidth', 2, 'DisplayName', 'Proportional pursuit');
for iCase = 1:nCases
    evMPC = evaderPosition(TsCell{iCase}, E0, VE, thetaE, zTarget);
    rng = sqrt(sum((XsCell{iCase}(1:2,:)' - evMPC(:,1:2)).^2, 2));
    plot(TsCell{iCase}, rng, [colors{iCase} styles{iCase}], ...
        'LineWidth', 2, 'DisplayName', ['MPC ' caseNames{iCase}]);
end
xline(gateTime, 'k:', 'Gate time', 'DisplayName', 'Gate time');
xlabel('t [s]'); ylabel('planar range [m]');
title('Range to Evader');
legend('Location', 'best');

% Position vs time
figure;
dimLabels = {'x [m]', 'y [m]', 'z [m]'};

for iDim = 1:3
    subplot(3,1,iDim);
    hold on; grid on;
    for iCase = 1:nCases
        plot(TsCell{iCase}, XsCell{iCase}(iDim,:), [colors{iCase} styles{iCase}], ...
            'LineWidth', 1.5, 'DisplayName', caseNames{iCase});
    end

    evMPC1 = evaderPosition(TsCell{1}, E0, VE, thetaE, zTarget);
    % evMPC1 is assumed to be N×3 (time × dims). Plot column iDim.
    plot(TsCell{1}, evMPC1(:,iDim), 'k:', 'LineWidth', 1.2, 'DisplayName', 'Evader');

    ylabel(dimLabels{iDim});
    if iDim == 1
        legend('Location', 'best');
        title('Position vs time');
    end
end

xlabel('time [s]');

% --- Plotting: 3D on left, x/y components on right, 2D planar on bottom ---
figure;
tlay = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

% Precompute evader trajectories used in multiple plots
tMax = max(cellfun(@(t) t(end), TsCell));
tEvaderPlot = (0:dt:tMax)';
evaderPlot = evaderPosition(tEvaderPlot, E0, VE, thetaE, zTarget);
% Evader at first case times (for overlaying on component plots)
evaderCase1 = evaderPosition(TsCell{1}, E0, VE, thetaE, zTarget);

% Left: 3D trajectory spanning rows 1-2, col 1
ax3d = nexttile(1,[2 1]);
hold(ax3d,'on'); grid(ax3d,'on'); box(ax3d,'on');
plot3(ax3d, evaderPlot(:,1), evaderPlot(:,2), evaderPlot(:,3), 'k--', 'LineWidth', 1.5, 'DisplayName','Evader');
for iCase = 1:nCases
    plot3(ax3d, XsCell{iCase}(1,:), XsCell{iCase}(2,:), XsCell{iCase}(3,:), ...
        [colors{iCase} styles{iCase}], 'LineWidth', 2, 'DisplayName',['MPC ' caseNames{iCase}]);
    if intCell(iCase)
        plot3(ax3d, XsCell{iCase}(1,end), XsCell{iCase}(2,end), XsCell{iCase}(3,end), ...
            [colors{iCase} 's'], 'MarkerFaceColor', colors{iCase}, 'MarkerSize', 8, 'DisplayName',[caseNames{iCase} ' intercept']);
    end
end
plot3(ax3d, P0(1),P0(2),0,'bo','MarkerFaceColor','b','DisplayName','P_0');
plot3(ax3d, E0(1),E0(2),zTarget,'ro','MarkerFaceColor','r','DisplayName','E_0');
plot3(ax3d, gatePoint(1),gatePoint(2),gatePoint(3),'md','MarkerFaceColor','m','DisplayName','Gate');
zlim(ax3d,[0 5]);
xlabel(ax3d,'x [m]'); ylabel(ax3d,'y [m]'); zlabel(ax3d,'z [m]');
title(ax3d,'3D MPC: Nominal vs Mass Mismatch');
legend(ax3d,'Location','best');
view(ax3d,3);
hold(ax3d,'off');

% Right-top: x(t)
ax_x = nexttile(2);
hold(ax_x,'on'); grid(ax_x,'on');
for iCase = 1:nCases
    plot(ax_x, TsCell{iCase}, XsCell{iCase}(1,:), [colors{iCase} styles{iCase}], 'LineWidth', 1.5, 'DisplayName', caseNames{iCase});
end
plot(ax_x, TsCell{1}, evaderCase1(:,1), 'k:', 'LineWidth', 1.2, 'DisplayName', 'Evader');
ylabel(ax_x,'x [m]');
title(ax_x,'x(t)');
legend(ax_x,'Location','best');
hold(ax_x,'off');

% Right-middle: z(t)
ax_z = nexttile(4);
hold(ax_z,'on'); grid(ax_z,'on');
for iCase = 1:nCases
    plot(ax_z, TsCell{iCase}, XsCell{iCase}(3,:), [colors{iCase} styles{iCase}], ...
        'LineWidth', 1.5, 'DisplayName', caseNames{iCase});
end
plot(ax_z, TsCell{1}, evaderCase1(:,3), 'k:', 'LineWidth', 1.2, 'DisplayName', 'Evader');
ylabel(ax_z,'z [m]');
title(ax_z,'z(t)');
legend(ax_z,'Location','best');
hold(ax_z,'off');

% Bottom full-width: planar 2D trajectory spanning row 3 cols 1-2
ax2d = nexttile(5,[1 2]);
hold(ax2d,'on'); axis(ax2d,'equal'); grid(ax2d,'on'); box(ax2d,'on');
plot(ax2d, evaderPlot(:,1), evaderPlot(:,2), 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');
for iCase = 1:nCases
    plot(ax2d, XsCell{iCase}(1,:), XsCell{iCase}(2,:), [colors{iCase} styles{iCase}], ...
        'LineWidth', 2, 'DisplayName', ['MPC ' caseNames{iCase}]);
    if intCell(iCase)
        plot(ax2d, XsCell{iCase}(1,end), XsCell{iCase}(2,end), [colors{iCase} 's'], ...
            'MarkerFaceColor', colors{iCase}, 'MarkerSize', 10, 'DisplayName', [caseNames{iCase} ' intercept']);
    end
end
xlabel(ax2d,'x [m]'); ylabel(ax2d,'y [m]');
title(ax2d,'Planar Trajectories');
legend(ax2d,'Location','best');
hold(ax2d,'off');

%% ========================= Local functions =========================
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