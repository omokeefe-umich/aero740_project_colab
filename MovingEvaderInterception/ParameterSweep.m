%% ParameterSweep.m
% Sweep over quadcopter start position (P0), ground-vehicle speed (VE),
% and gate distance (gateDist).
%
% Metrics collected per case:
%   timeIntercept   – simulation time at intercept [s] (inf if not intercepted)
%   distToGate      – remaining evader distance to gate at intercept [m] (NaN if not intercepted)
%   fzControlUsage  – integral of absolute Fz deviation from hover over time [N·s]
%   torqueUsage     – integral of sum of torque magnitudes over time [N·m·s]
%   meanCompTime    – mean wall-clock time per MPC update [s]
%   minCtrlMargin   – smallest margin between any control channel and its bound [N or N·m]
%   minStateMargin  – smallest margin between any bounded state and its bound [m or rad]
%   timeAtSaturation– number of time steps where controls are at their limits 
%
% Results are printed as a table and saved to ParameterSweepResults.csv.

clear; clc; close all;

%% ── Fixed quadcopter physical parameters ───────────────────────────────────
m  = 1.12;
g  = 9.81;
Jr = 8.5e-4;
Ix = 0.0019;
Iy = 0.0019;
Iz = 0.0223;
params   = [m; g; Jr; Ix; Iy; Iz];
u_hover  = [m * g; 0; 0; 0];

%% ── Wind ────────────────────────────────────────────────────────────────────
wind_multiplier = 1;
wind_data = DisturbanceModel(1, wind_multiplier);
no_wind   = wind_data;
no_wind.U = wind_data.U * 0;
no_wind.W = wind_data.W * 0;

%% ── Control bounds (shared across cases) ───────────────────────────────────
torquelim = 0.4;
anglim    = deg2rad(70);
uub  = [ 3.0*m*g;  torquelim;  torquelim;  0.25*torquelim];
ulb  = [ 0.0*m*g; -torquelim; -torquelim; -0.25*torquelim];
duub = [ 4.0;  0.15;  0.15;  0.05];   % per dt_mpc, scaled in loop
dulb = [-4.0; -0.15; -0.15; -0.05];

%% ── MPC / simulation settings (shared) ─────────────────────────────────────
dt_sim       = 0.2;
dt_mpc       = 0.2;
horizonSteps = 40;
horizonTime  = horizonSteps * dt_mpc;
interceptRadius = 0.5;
leadTime     = 0.0;

%% ── Cost matrices ───────────────────────────────────────────────────────────
Q = blkdiag(0.1  * eye(3), ...
            0.001* eye(2), ...
            zeros(4), ...
            eye(3)   * 0.01, ...
            eye(2)   * 0.01, ...
            0.001, ...
            eye(3)   * 0.01, ...
            zeros(4));
P_term = blkdiag(eye(3), zeros(19));
R_mat  = [1e-6, 0,    0,    0;
          0,    1e-4, 0,    0;
          0,    0,    1e-4, 0;
          0,    0,    0,    1e-4];

%% ── Fixed evader parameters ─────────────────────────────────────────────────
E0      = [20; 0];
thetaE  = pi / 2;
zTarget = 0.0;

%% ── QP solver options ───────────────────────────────────────────────────────
qpOpts = optimoptions('quadprog', ...
    'Display',              'off', ...
    'Algorithm',            'interior-point-convex', ...
    'MaxIterations',        5000, ...
    'ConstraintTolerance',  1e-3, ...
    'OptimalityTolerance',  1e-3, ...
    'StepTolerance',        1e-8);

%% ── Per-step augmented-state constraint matrices (shared) ───────────────────
A_z    = zeros(1, 22);  A_z(3)    = -1;
A_ang  = zeros(3, 22);  A_ang(:, 4:6)   = eye(3);
A_u    = zeros(4, 22);  A_u(:, 19:22)   = eye(4);
A_step = [A_z; A_ang; -A_ang; A_u; -A_u];
b_step = [0; anglim*ones(3,1); anglim*ones(3,1); uub; -ulb];

%% ══════════════════════════════════════════════════════════════════════════
%% Parameter grid
%% ══════════════════════════════════════════════════════════════════════════
%  P0 cases  : three qualitatively different starting geometries.
%    [0;0]   – default (south-west of evader)
%    [30;0]  – far to the right (east of evader)
%    [0;10]  – ahead along the evader's path (north-west)
P0_cases    = {[0;0], [30;0], [15;0]};
VE_cases    = [3, 5, 7];          % evader speeds [m/s]
gDist_cases = [20, 25, 30];       % gate distances [m]

nP0    = numel(P0_cases);
nVE    = numel(VE_cases);
nGD    = numel(gDist_cases);
nTotal = nP0 * nVE * nGD;

%% ── Pre-allocate result storage ─────────────────────────────────────────────
P0x_res          = nan(nTotal, 1);
P0y_res          = nan(nTotal, 1);
VE_res           = nan(nTotal, 1);
gateDist_res     = nan(nTotal, 1);
timeIntercept_res = inf(nTotal, 1);
distToGate_res   = nan(nTotal, 1);
fzControlUsage_res = nan(nTotal, 1);
torqueUsage_res = nan(nTotal, 1);
meanCompTime_res  = nan(nTotal, 1);
minCtrlMargin_res = nan(nTotal, 1);
minStateMargin_res = nan(nTotal, 1);
timeAtSaturation_res = nan(nTotal, 1);
dt_mpc_res       = nan(nTotal, 1);
PosError_res    = nan(nTotal, 1);

%% ══════════════════════════════════════════════════════════════════════════
%% Main sweep loop
%% ══════════════════════════════════════════════════════════════════════════
idx = 0;
for iP = 1:nP0
  for iV = 1:nVE
    for iG = 1:nGD
        idx      = idx + 1;
        P0       = P0_cases{iP};
        VE       = VE_cases(iV);
        gateDist = gDist_cases(iG);
        posError = 0;

        fprintf('[%2d/%2d]  P0=[%5.1f,%5.1f]  VE=%g m/s  gateDist=%g m  ...', ...
            idx, nTotal, P0(1), P0(2), VE, gateDist);

        % ── Gate geometry ────────────────────────────────────────────────
        gateAxis     = [cos(thetaE); sin(thetaE)];
        gatePoint2D  = E0 + gateDist * gateAxis;
        gateCoord    = gateAxis' * gatePoint2D;               % scalar projection
        evaderInitCoord = gateAxis' * E0;
        simMaxTime   = ceil(gateDist / VE);
        maxSIMsteps  = floor(simMaxTime / dt_sim);
        
        dt_mpc_res(idx) = dt_mpc;
        duStackUB = repmat(duub / dt_mpc, horizonSteps, 1);
        duStackLB = repmat(dulb / dt_mpc, horizonSteps, 1);

        % ── Initial conditions ───────────────────────────────────────────
        x = zeros(12, 1);
        x(1:2) = P0;
        x(3)   = zTarget;
        u = u_hover;
        a = dynamics(0, x, u, params, wind_data);
        a = a(7:12);

        % ── History storage ──────────────────────────────────────────────
        allXs        = x;
        inputHistory = u;
        times        = 0;         %#ok<NASGU>
        Xs           = x(1:3);
        Ts           = 0;
        wallTimes    = zeros(maxSIMsteps, 1);
        mpcIntercepted = false;

        % ── Simulation loop ──────────────────────────────────────────────
        for k = 1:maxSIMsteps
            currentTime = (k - 1) * dt_sim;

            tnow = tic;
            J_xv = Jacobian(@(xj) dynamics(0, xj, u, params, no_wind), x);
            J_u  = Jacobian(@(uj) dynamics(0, x, uj, params, no_wind), u);
            Ac = [zeros(12,6), eye(12), zeros(12,4); ...
                  zeros(6,6),  J_xv(7:12,:), zeros(6,4); ...
                  zeros(4,22)];
            Bc = [zeros(12,4); J_u(7:12,:); eye(4)];
            [Ads, Bds] = Discretize_dt(dt_mpc, ones(horizonSteps,1), Ac, Bc);
            [S, Mpred] = StackedMatrix(Ads, Bds);

            x_aug = [x; a; u];

            predTimes = (currentTime + leadTime + dt_mpc : dt_mpc : ...
                         currentTime + leadTime + horizonTime)';
            targetXYZ = evaderPosition(predTimes, E0, VE, thetaE, zTarget);
            x_traj = [targetXYZ, zeros(horizonSteps,15), ...
                      repmat(u_hover', horizonSteps, 1)]';
            x_traj = x_traj(:);

            [H, q] = QPFormat(Q, R_mat, P_term, S, Mpred, ...
                              (0:dt_mpc:horizonTime)', x_aug, x_traj, ...
                              zeros(horizonSteps*4, 1));
            Hqp = (H + H') / 2;

            evaderNow       = evaderPosition(currentTime, E0, VE, thetaE, zTarget);
            evaderCoord     = gateAxis' * evaderNow(1,1:2)';
            gateTimeLeft    = (gateCoord - evaderCoord) / VE;
            if gateTimeLeft <= 0
                warning('Evader reached gate before intercept at t = %.2f s.', currentTime);
                break;
            end

            Mx = Mpred * x_aug;
            A_ss = kron(eye(horizonSteps), A_step);
            AineqBase = A_ss * S;
            bineqBase = repmat(b_step, horizonSteps, 1) - A_ss * Mx;

            terminalLeft = gateTimeLeft - (leadTime + horizonTime);
            AineqFinal   = AineqBase;
            bineqFinal   = bineqBase;
            gateActive   = false;
            if terminalLeft > 0
                S_N = S((horizonSteps-1)*22+1 : horizonSteps*22, :);
                M_N = Mpred((horizonSteps-1)*22+1 : horizonSteps*22, :);
                A_gate_N      = zeros(1, 22);
                A_gate_N(1:2) = -gateAxis';
                A_gate_N(7:8) = -terminalLeft * gateAxis';
                b_gate_N      = -gateCoord;
                AineqFinal = [A_gate_N * S_N; AineqBase];
                bineqFinal = [b_gate_N - A_gate_N * M_N * x_aug; bineqBase];
                gateActive = true;
            end

            % ── Solve QP ────────────────────────────────────────────────
            if gateActive
                [u_mpc, ~, exitflag] = quadprog(Hqp, q, AineqFinal, bineqFinal, ...
                    [], [], duStackLB, duStackUB, [], qpOpts);
            else
                [u_mpc, ~, exitflag] = quadprog(Hqp, q, AineqBase, bineqBase, ...
                    [], [], duStackLB, duStackUB, [], qpOpts);
            end

            % Fall back: drop gate constraint if infeasible
            if (exitflag <= 0 || isempty(u_mpc) || any(~isfinite(u_mpc))) && gateActive
                [u_mpc, ~, exitflag] = quadprog(Hqp, q, AineqBase, bineqBase, ...
                    [], [], duStackLB, duStackUB, [], qpOpts);
            end

            if exitflag <= 0 || isempty(u_mpc) || any(~isfinite(u_mpc))
                warning('QP failed at t = %.2f s (case %d/%d). Stopping.', currentTime, idx, nTotal);
                break;
            end

            wallTimes(k) = toc(tnow);

            u = u + u_mpc(1:4) * dt_mpc;
            u = min(max(u, ulb), uub);
            [T_ode, X_ode] = ode45(@(t,s) dynamics(t, s, u, params, wind_data), ...
                                   0:dt_sim/10:dt_sim, x);
            Xs   = [Xs,   X_ode(2:end, 1:3)'];   %#ok<AGROW>
            Ts   = [Ts;   currentTime + T_ode(2:end)]; %#ok<AGROW>
            a    = (X_ode(end, 7:12)' - x(7:12)) / dt_sim;
            x    = X_ode(end, :)';

            x_error = squeeze(Ads(1,:,:))*x_aug + squeeze(Bds(1,:,:))*u_mpc(1:4) - [x;a;u];
            posError = max(posError, norm(x_error(1:3)));
            % fprintf('Norm of augmented state error: %.3f\n', norm(x_error(1:3)));

            allXs        = [allXs,        x];          %#ok<AGROW>
            inputHistory = [inputHistory, u];           %#ok<AGROW>

            evaderAt = evaderPosition(Ts(end), E0, VE, thetaE, zTarget);
            rangeNow = norm(Xs(1:2, end)' - evaderAt(1, 1:2));
            if rangeNow <= interceptRadius
                mpcIntercepted = true;
                break;
            end
        end  % simulation loop

        % ── Metrics ─────────────────────────────────────────────────────
        P0x_res(idx)      = P0(1);
        P0y_res(idx)      = P0(2);
        VE_res(idx)       = VE;
        gateDist_res(idx) = gateDist;

        if mpcIntercepted
            timeIntercept_res(idx) = Ts(end);
            evaderAtInt = evaderPosition(Ts(end), E0, VE, thetaE, zTarget);
            evaderCoordInt = gateAxis' * evaderAtInt(1, 1:2)';
            distToGate_res(idx) = gateCoord - evaderCoordInt;
            PosError_res(idx) = posError;
        else
            timeIntercept_res(idx) = inf;
            distToGate_res(idx)    = 0;
        end

        % Control usage split into thrust and torque contributions.
        uDev = inputHistory - u_hover;
        fzControlUsage_res(idx) = sum(abs(uDev(1, :))) * dt_sim;
        torqueUsage_res(idx) = sum(sum(abs(inputHistory(2:4, :)), 1)) * dt_sim;

        % Mean computation time per MPC step
        wt = wallTimes(wallTimes > 0);
        meanCompTime_res(idx) = mean(wt);

        % Minimum margin to control bounds (positive = inside bounds)
        ctrlMarginUB = uub  - inputHistory;   % margin to upper bound
        ctrlMarginLB = inputHistory - ulb;    % margin to lower bound
        minCtrlMargin_res(idx) = min([ctrlMarginUB(:); ctrlMarginLB(:)]);

        % Minimum margin to state bounds (z >= 0, |angle| <= anglim)
        zMargin   = allXs(3, :);
        angMargin = anglim - abs(allXs(4:6, :));
        minStateMargin_res(idx) = min([zMargin(:); angMargin(:)]);
        numstepssaturated = sum(sum(inputHistory - ulb < 1e-2 | inputHistory - uub > -1e-2));
        timeAtSaturation_res(idx) = dt_sim * numstepssaturated;

        if mpcIntercepted
            fprintf('  intercepted  t=%.2f s  distToGate=%.2f m\n', ...
                timeIntercept_res(idx), distToGate_res(idx));
        else
            fprintf('  NOT intercepted\n');
        end
    end
  end
end

%% ══════════════════════════════════════════════════════════════════════════
%% Assemble and display results table
%% ══════════════════════════════════════════════════════════════════════════
ResultsTable = table( ...
    dt_mpc_res, ...
    P0x_res,           P0y_res,          VE_res,          gateDist_res, ...
    timeIntercept_res, distToGate_res,   fzControlUsage_res, torqueUsage_res, meanCompTime_res, ...
    minCtrlMargin_res, minStateMargin_res, timeAtSaturation_res, PosError_res, ...
    'VariableNames', { ...
        'dt_mpc_s', ...
        'P0_x', 'P0_y', 'VE_mps', 'gateDist_m', ...
        'timeIntercept_s', 'distToGate_m', 'fzControlUsage', 'torqueUsage', 'meanCompTime_s', ...
        'minCtrlMargin', 'minStateMargin', 'timeAtSaturation', 'MaxPositionError_m'});

fprintf('\n%s\n', repmat('=', 1, 110));
fprintf('PARAMETER SWEEP RESULTS\n');
fprintf('%s\n\n', repmat('=', 1, 110));
disp(ResultsTable)

writetable(ResultsTable, fullfile(pwd, 'ParameterSweepResults.csv'), "WriteMode", 'append');
fprintf('Results saved to ParameterSweepResults.csv\n');

%% ══════════════════════════════════════════════════════════════════════════
%% Plot parameter sweep results (single figure, tiled 2D scatter views)
%% ══════════════════════════════════════════════════════════════════════════
failMask = ~isfinite(timeIntercept_res);
boundMargin_res = min(minCtrlMargin_res, minStateMargin_res);

distToGate_plot = distToGate_res;
distToGate_plot(failMask) = nan;

metricSeries = {
    timeIntercept_res,      'Time To Intercept [s]';
    distToGate_plot,        'Distance To Gate At Intercept [m]';
    fzControlUsage_res,     'Fz Usage [N s]';
    torqueUsage_res,        'Torque Usage [N m s]';
    1e3 * meanCompTime_res, 'Computation Time [ms]';
    minCtrlMargin_res,      'Min Control Bound Margin';
    minStateMargin_res,     'Min State Bound Margin';
    boundMargin_res,        'Min Overall Bound Margin';
    timeAtSaturation_res,   'Time At Saturation [s]'
    };

plotFig = figure('Name', 'Parameter Sweep Response Maps', 'Color', 'w');
tl = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

markerList = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
legendAx = [];

% Plot/legend order: highest-to-lowest P0 starting range.
p0Range = zeros(nP0, 1);
for iP = 1:nP0
    p0Range(iP) = hypot(P0_cases{iP}(1), P0_cases{iP}(2));
end
[~, pOrder] = sort(p0Range, 'descend');

% Small visual offsets keep different P0 markers from fully overlapping.
if nP0 > 1
    veOffsets = linspace(-0.18, 0.18, nP0);
    gateOffsets = linspace(0.24, -0.24, nP0);
else
    veOffsets = 0;
    gateOffsets = 0;
end

for iMetric = 1:size(metricSeries, 1)
    metricVals = metricSeries{iMetric, 1};
    metricName = metricSeries{iMetric, 2};
    validMaskAll = isfinite(metricVals);

    ax = nexttile(tl);
    if isempty(legendAx)
        legendAx = ax;
    end
    hold(ax, 'on');

    finiteVals = metricVals(validMaskAll);
    hasFinite = ~isempty(finiteVals);
    if hasFinite
        cLim = [min(finiteVals), max(finiteVals)];
        if abs(cLim(2) - cLim(1)) < 1e-12
            cLim(2) = cLim(1) + 1;
        end
    end

    for kP = 1:nP0
        iP = pOrder(kP);
        idxCase = (P0x_res == P0_cases{iP}(1)) & (P0y_res == P0_cases{iP}(2));
        validMask = idxCase & validMaskAll;
        invalidMask = idxCase & ~validMaskAll;
        markerSym = markerList{mod(kP - 1, numel(markerList)) + 1};

        if any(validMask)
            vePlot = VE_res(validMask) + veOffsets(kP);
            gatePlot = gateDist_res(validMask) + gateOffsets(kP);
            h = scatter(ax, vePlot, gatePlot, 130, metricVals(validMask), ...
                markerSym, 'filled', 'MarkerEdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.6);
        end

        if any(invalidMask)
            vePlotBad = VE_res(invalidMask) + veOffsets(kP);
            gatePlotBad = gateDist_res(invalidMask) + gateOffsets(kP);
            scatter(ax, vePlotBad, gatePlotBad, 100, 'kx', 'LineWidth', 1.6);
        end
    end

    if hasFinite
        caxis(ax, cLim);
    end
    colormap(ax, bone);
    cb = colorbar(ax);
    cb.Label.String = metricName;
    cb.FontSize = 8;

    xlabel(ax, 'VE [m/s]');
    ylabel(ax, 'Gate Distance [m]');
    title(ax, metricName, 'FontSize', 10);
    xticks(ax, unique(VE_res));
    yticks(ax, unique(gateDist_res));
    grid(ax, 'on');
    box(ax, 'on');
end

% Build a marker-type legend that is always present and independent of data availability.
if ~isempty(legendAx)
    markerHandles = gobjects(nP0 + 1, 1);
    markerNames = cell(nP0 + 1, 1);
    for kP = 1:nP0
        iP = pOrder(kP);
        markerSym = markerList{mod(kP - 1, numel(markerList)) + 1};
        markerHandles(kP) = scatter(legendAx, nan, nan, 130, 0.55, markerSym, ...
            'filled', 'MarkerEdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.6);
        markerNames{kP} = sprintf('P0 = [%g, %g]', P0_cases{iP}(1), P0_cases{iP}(2));
    end
    markerHandles(nP0 + 1) = scatter(legendAx, nan, nan, 100, 'kx', 'LineWidth', 1.6);
    markerNames{nP0 + 1} = 'Undefined / no intercept';
    legend(legendAx, markerHandles, markerNames, 'Location', 'southoutside', 'Orientation', 'horizontal');
end

sgtitle(tl, 'Parameter Sweep Metrics (black X = undefined / no intercept)');
set(plotFig, 'Position', [80, 60, 1650, 980]);

plotImagePath = fullfile(pwd, 'ParameterSweepResultsFigure.png');
exportgraphics(plotFig, plotImagePath, 'Resolution', 300);
fprintf('Saved sweep figure to %s\n', plotImagePath);

%% ══════════════════════════════════════════════════════════════════════════
%% Helper: evader position as a function of time
%% ══════════════════════════════════════════════════════════════════════════
function targetXYZ = evaderPosition(t, E0, VE, thetaE, zTarget)
    t = t(:);
    targetXYZ = [E0(1) + VE * cos(thetaE) * t, ...
                 E0(2) + VE * sin(thetaE) * t, ...
                 zTarget * ones(size(t))];
end
