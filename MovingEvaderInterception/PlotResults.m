function PlotResults(t, Xs, x_ref, Us, constraints, windData, comparisonData)
%PLOTRESULTS  Visualize quadcopter MPC simulation results against a reference.
%
% Generates figures for:
%   1. XY trajectory vs reference path
%   2. 3-D trajectory vs reference path with optional wind overlay
%   3. Altitude vs time with bounds
%   4. All state histories with lower/upper constraints
%   5. Translational velocities and numerical accelerations vs time
%   6. Euler angles and angular rates vs time with bounds
%   7. Control inputs vs time with lower/upper constraints
%   8. Tracking error vs time
%
% Inputs:
%   t           - time vector (1×N or N×1) [s]
%   Xs          - 12×N state history, rows ordered as:
%                 [x; y; z; phi; theta; psi; xdot; ydot; zdot; ...
%                  phidot; thetadot; psidot]
%   x_ref       - 3×M or M×3 reference trajectory [x_ref, y_ref, z_ref]
%   Us          - 4×N or N×4 control input history (optional)
%   constraints - struct with fields:
%                 .stateLower   (12×1)
%                 .stateUpper   (12×1)
%                 .controlLower (4×1)
%                 .controlUpper (4×1)
%   windData    - optional wind-field struct from DisturbanceModel
%   comparisonData - optional struct for PP/MPC trajectory comparison

    if nargin < 4
        Us = [];
    end
    if nargin < 5 || isempty(constraints)
        constraints = struct();
    end
    if nargin < 6
        windData = [];
    end
    if nargin < 7 || isempty(comparisonData)
        comparisonData = struct();
    end

    t = t(:)';

    if size(Xs, 1) ~= 12
        error('Xs must be a 12xN state history.');
    end
    if size(Xs, 2) ~= numel(t)
        error('The number of columns in Xs must match the length of t.');
    end

    [stateLower, stateUpper, controlLower, controlUpper] = unpackConstraints(constraints);

    % ------------------------------------------------------------------ %
    %  Unpack state rows
    % ------------------------------------------------------------------ %
    x    = Xs(1,:);   y    = Xs(2,:);   z    = Xs(3,:);
    phi  = Xs(4,:);   th   = Xs(5,:);   ps   = Xs(6,:);
    xd   = Xs(7,:);   yd   = Xs(8,:);   zd   = Xs(9,:);
    phid = Xs(10,:);  thd  = Xs(11,:);  psd  = Xs(12,:);

    % Numerical accelerations via central finite differences.
    xdd = gradient(xd, t);
    ydd = gradient(yd, t);
    zdd = gradient(zd, t);

    % Reference trajectory components (accept either 3xM or Mx3).
    refXYZ = normalizeReferenceTrajectory(x_ref);
    xrPath = refXYZ(1,:);
    yrPath = refXYZ(2,:);
    zrPath = refXYZ(3,:);

    if size(refXYZ, 2) == numel(t)
        xr = xrPath;
        yr = yrPath;
        zr = zrPath;
    else
        tRef = linspace(t(1), t(end), size(refXYZ, 2));
        xr = interp1(tRef, xrPath, t, 'linear', 'extrap');
        yr = interp1(tRef, yrPath, t, 'linear', 'extrap');
        zr = interp1(tRef, zrPath, t, 'linear', 'extrap');
    end

    % Control history (accept either 4xN or Nx4).
    Uplot = [];
    if ~isempty(Us)
        if size(Us, 1) == 4
            Uplot = Us;
        elseif size(Us, 2) == 4
            Uplot = Us';
        else
            error('Us must be 4xN or Nx4.');
        end

        if size(Uplot, 2) ~= numel(t)
            tU = linspace(t(1), t(end), size(Uplot, 2));
            Uplot = interp1(tU, Uplot', t, 'linear', 'extrap')';
        end
    end

    clr = {'b', 'r', [0 0.55 0]};
    stateData = {x, y, z, phi, th, ps, xd, yd, zd, phid, thd, psd};
    stateLabels = {'x [m]', 'y [m]', 'z [m]', '$\phi$ [rad]', '$\theta$ [rad]', '$\psi$ [rad]', ...
                   '$\dot{x}$ [m/s]', '$\dot{y}$ [m/s]', '$\dot{z}$ [m/s]', ...
                   '$\dot{\phi}$ [rad/s]', '$\dot{\theta}$ [rad/s]', '$\dot{\psi}$ [rad/s]'};

    % ================================================================== %
    %  Figure 1 – Trajectory comparison: PP and MPC against the same moving evader
    % ================================================================== %
    figure(1); clf;
    hold on; axis equal; grid on; box on;

    if isfield(comparisonData, 'evaderPlot') && ~isempty(comparisonData.evaderPlot)
        plot(comparisonData.evaderPlot(:, 1), comparisonData.evaderPlot(:, 2), 'k--', ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
    else
        plot(xrPath, yrPath, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end

    if isfield(comparisonData, 'ppPath') && ~isempty(comparisonData.ppPath)
        plot(comparisonData.ppPath(:, 1), comparisonData.ppPath(:, 2), 'b-', ...
            'LineWidth', 2, 'DisplayName', 'Proportional pursuit');
    end

    plot(x, y, 'r-', 'LineWidth', 2, 'DisplayName', 'MPC quadcopter');

    if isfield(comparisonData, 'gatePoint2D') && numel(comparisonData.gatePoint2D) >= 2
        plot(comparisonData.gatePoint2D(1), comparisonData.gatePoint2D(2), 'md', ...
            'MarkerFaceColor', 'm', 'HandleVisibility', 'off');
    end

    if isfield(comparisonData, 'ppIntercept') && ~isempty(comparisonData.ppIntercept)
        plot(comparisonData.ppIntercept(1), comparisonData.ppIntercept(2), 'ko', ...
            'MarkerFaceColor', 'k', 'DisplayName', 'PP intercept');
    end

    if isfield(comparisonData, 'mpcIntercepted') && comparisonData.mpcIntercepted
        plot(x(end), y(end), 'ms', 'MarkerFaceColor', 'm', 'DisplayName', 'MPC intercept');
    end

    xlabel('x [m]', 'Interpreter', 'latex');
    ylabel('y [m]', 'Interpreter', 'latex');
    title('Moving Evader Pursuit: PP vs MPC', 'Interpreter', 'latex');
    legend('Location', 'best');

    % ================================================================== %
    %  Figure 2 – 3-D Trajectory
    % ================================================================== %
    figure('Name', '3-D Trajectory', 'NumberTitle', 'off');
    hold on; grid on;
    yWindSlices = linspace(min([y, yrPath]), max([y, yrPath]), 3);
    xyzLimits = [min([x, xrPath]), max([x, xrPath]); ...
                 min([y, yrPath]), max([y, yrPath]); ...
                 min([z, zrPath]), max([z, zrPath])];
    addWindOverlay3D(windData, yWindSlices, xyzLimits);
    plot3(xrPath, yrPath, zrPath, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Moving evader');
    plot3(x, y, z, 'k-', 'LineWidth', 1.5, 'DisplayName', 'MPC quadcopter');
    scatter3(x(1), y(1), z(1), 110, 'g', 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'Start');
    scatter3(x(end), y(end), z(end), 110, 'r', 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'End');
    xlabel('x [m]', 'Interpreter', 'latex');
    ylabel('y [m]', 'Interpreter', 'latex');
    zlabel('z [m]', 'Interpreter', 'latex');
    title('MPC Quadcopter Trajectory Toward Moving Evader', 'Interpreter', 'latex');
    legend('Location', 'best');
    view(45, 30);

    % ================================================================== %
    %  Figure 3 – Altitude vs Time
    % ================================================================== %
    figure('Name', 'z vs Time', 'NumberTitle', 'off');
    hold on; grid on;
    plot(t, zr, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
    plot(t, z, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Simulated');
    plotBounds(t, stateLower(3), stateUpper(3));
    scatter(t(1), z(1), 110, 'g', 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'Start');
    scatter(t(end), z(end), 110, 'r', 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'End');
    xlabel('t [s]', 'Interpreter', 'latex');
    ylabel('z [m]', 'Interpreter', 'latex');
    title('Altitude vs Time', 'Interpreter', 'latex');
    legend('Location', 'best');

    % ================================================================== %
    %  Figure 4 – State Histories with Constraints
    % ================================================================== %
    figure('Name', 'State Histories with Constraints', 'NumberTitle', 'off');
    for i = 1:12
        subplot(4, 3, i);
        hold on; grid on;
        plot(t, stateData{i}, 'k-', 'LineWidth', 1.2);
        plotBounds(t, stateLower(i), stateUpper(i));
        xlabel('t [s]', 'Interpreter', 'latex');
        ylabel(stateLabels{i}, 'Interpreter', 'latex');
        title(stateLabels{i}, 'Interpreter', 'latex');
        xlim([t(1), t(end)]);
    end
    sgtitle('State Histories and Active Constraints', 'Interpreter', 'latex');

    % ================================================================== %
    %  Figure 5 – Translational Velocities & Accelerations
    % ================================================================== %
    figure('Name', 'Translational Kinematics', 'NumberTitle', 'off');
    velData   = {xd, yd, zd};
    accelData = {xdd, ydd, zdd};
    velLabels = {'$\dot{x}$ [m/s]', '$\dot{y}$ [m/s]', '$\dot{z}$ [m/s]'};
    accLabels = {'$\ddot{x}$ [m/s^2]', '$\ddot{y}$ [m/s^2]', '$\ddot{z}$ [m/s^2]'};

    for i = 1:3
        subplot(2, 3, i);
        hold on; grid on;
        plot(t, velData{i}, 'Color', clr{i}, 'LineWidth', 1.2);
        plotBounds(t, stateLower(6 + i), stateUpper(6 + i));
        xlabel('t [s]', 'Interpreter', 'latex');
        ylabel(velLabels{i}, 'Interpreter', 'latex');
        title(velLabels{i}, 'Interpreter', 'latex');
        xlim([t(1), t(end)]);

        subplot(2, 3, i + 3);
        hold on; grid on;
        plot(t, accelData{i}, 'Color', clr{i}, 'LineWidth', 1.2);
        xlabel('t [s]', 'Interpreter', 'latex');
        ylabel(accLabels{i}, 'Interpreter', 'latex');
        title(accLabels{i}, 'Interpreter', 'latex');
        xlim([t(1), t(end)]);
    end
    sgtitle('Translational Velocities (top row) and Accelerations (bottom row)', 'Interpreter', 'latex');

    % ================================================================== %
    %  Figure 6 – Euler Angles & Angular Rates
    % ================================================================== %
    figure('Name', 'Rotational Kinematics', 'NumberTitle', 'off');
    angleData = {phi, th, ps};
    rateData  = {phid, thd, psd};
    angleLabels = {'$\phi$ [rad]', '$\theta$ [rad]', '$\psi$ [rad]'};
    rateLabels  = {'$\dot{\phi}$ [rad/s]', '$\dot{\theta}$ [rad/s]', '$\dot{\psi}$ [rad/s]'};

    for i = 1:3
        subplot(2, 3, i);
        hold on; grid on;
        plot(t, angleData{i}, 'Color', clr{i}, 'LineWidth', 1.2);
        plotBounds(t, stateLower(3 + i), stateUpper(3 + i));
        xlabel('t [s]', 'Interpreter', 'latex');
        ylabel(angleLabels{i}, 'Interpreter', 'latex');
        title(angleLabels{i}, 'Interpreter', 'latex');
        xlim([t(1), t(end)]);

        subplot(2, 3, i + 3);
        hold on; grid on;
        plot(t, rateData{i}, 'Color', clr{i}, 'LineWidth', 1.2);
        plotBounds(t, stateLower(9 + i), stateUpper(9 + i));
        xlabel('t [s]', 'Interpreter', 'latex');
        ylabel(rateLabels{i}, 'Interpreter', 'latex');
        title(rateLabels{i}, 'Interpreter', 'latex');
        xlim([t(1), t(end)]);
    end
    sgtitle('Euler Angles (top row) and Angular Rates (bottom row)', 'Interpreter', 'latex');

    % ================================================================== %
    %  Figure 7 – Control Inputs with Constraints
    % ================================================================== %
    if ~isempty(Uplot)
        figure('Name', 'Control Inputs', 'NumberTitle', 'off');
        inputLabels = {'$f_z$ [N]', '$\tau_{\phi}$ [N m]', '$\tau_{\theta}$ [N m]', '$\tau_{\psi}$ [N m]'};
        for i = 1:4
            subplot(2, 2, i);
            hold on; grid on;
            plot(t, Uplot(i, :), 'k-', 'LineWidth', 1.2);
            plotBounds(t, controlLower(i), controlUpper(i));
            xlabel('t [s]', 'Interpreter', 'latex');
            ylabel(inputLabels{i}, 'Interpreter', 'latex');
            title(inputLabels{i}, 'Interpreter', 'latex');
            xlim([t(1), t(end)]);
        end
        sgtitle('Control Inputs and Active Constraints', 'Interpreter', 'latex');
    end

    % ================================================================== %
    %  Figure 8 – Tracking Error vs Time
    % ================================================================== %
    figure('Name', 'Tracking Error', 'NumberTitle', 'off');
    subplot(2, 1, 1);
    trackingError = sqrt((x - xr).^2 + (y - yr).^2 + (z - zr).^2);
    plot(t, trackingError, 'm-', 'LineWidth', 1.5);
    grid on;
    xlabel('t [s]', 'Interpreter', 'latex');
    ylabel('Tracking Error [m]', 'Interpreter', 'latex');
    title('Tracking Error vs Time', 'Interpreter', 'latex');

    subplot(2, 1, 2);
    hold on; grid on;
    plot(t, x, 'b-', 'LineWidth', 1.2);
    plot(t, y, 'r-', 'LineWidth', 1.2);
    plot(t, z, 'Color', [0 0.55 0], 'LineWidth', 1.2);
    plotBounds(t, stateLower(1), stateUpper(1));
    plotBounds(t, stateLower(2), stateUpper(2));
    plotBounds(t, stateLower(3), stateUpper(3));
    xlabel('t [s]', 'Interpreter', 'latex');
    ylabel('Position [m]', 'Interpreter', 'latex');
    title('Position States', 'Interpreter', 'latex');
    legend('x', 'y', 'z', 'Location', 'best');

end

function refXYZ = normalizeReferenceTrajectory(x_ref)
    if size(x_ref, 1) == 3
        refXYZ = x_ref;
    elseif size(x_ref, 2) == 3
        refXYZ = x_ref';
    else
        error('x_ref must be 3xM or Mx3.');
    end
end

function [stateLower, stateUpper, controlLower, controlUpper] = unpackConstraints(constraints)
    stateLower = nan(12, 1);
    stateUpper = nan(12, 1);
    controlLower = nan(4, 1);
    controlUpper = nan(4, 1);

    if isfield(constraints, 'stateLower')
        stateLower = constraints.stateLower(:);
    end
    if isfield(constraints, 'stateUpper')
        stateUpper = constraints.stateUpper(:);
    end
    if isfield(constraints, 'controlLower')
        controlLower = constraints.controlLower(:);
    end
    if isfield(constraints, 'controlUpper')
        controlUpper = constraints.controlUpper(:);
    end
end

function plotBounds(t, lowerBound, upperBound)
    if ~isempty(lowerBound) && isfinite(lowerBound)
        plot([t(1), t(end)], [lowerBound, lowerBound], 'r--', 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
    end
    if ~isempty(upperBound) && isfinite(upperBound)
        plot([t(1), t(end)], [upperBound, upperBound], 'r--', 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
    end
end

function addWindOverlay3D(windData, ySlices, xyzLimits)
    if isempty(windData) || ~isstruct(windData)
        return;
    end

    requiredFields = {'X', 'Z', 'U', 'W'};
    for iField = 1:numel(requiredFields)
        if ~isfield(windData, requiredFields{iField})
            return;
        end
    end

    if ~any(abs(windData.U(:)) > 0 | abs(windData.W(:)) > 0)
        return;
    end

    nRows = size(windData.X, 1);
    nCols = size(windData.X, 2);
    rowIdx = unique(round(linspace(1, nRows, min(6, nRows))));
    colIdx = unique(round(linspace(1, nCols, min(9, nCols))));

    Xq = windData.X(rowIdx, colIdx);
    Zq = windData.Z(rowIdx, colIdx);
    Uq = windData.U(rowIdx, colIdx);
    Wq = windData.W(rowIdx, colIdx);

    if isfield(windData, 'Umag') && ~isempty(windData.Umag)
        Umag = windData.Umag;
    else
        Umag = sqrt(windData.U.^2 + windData.W.^2);
    end

    ySlices = unique(ySlices(:)');
    if isempty(ySlices)
        ySlices = 0;
    end

    yMid = ySlices(ceil(numel(ySlices) / 2));
    surf(windData.X, yMid * ones(size(windData.X)), windData.Z, Umag, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.18, 'DisplayName', 'Wind speed slice');
    colormap(turbo);
    cb = colorbar;
    cb.Label.String = 'Wind speed [m/s]';

    spanX = max(1, xyzLimits(1, 2) - xyzLimits(1, 1));
    spanY = max(1, xyzLimits(2, 2) - xyzLimits(2, 1));
    spanZ = max(1, xyzLimits(3, 2) - xyzLimits(3, 1));
    arrowScale = 0.08 * max([spanX, spanY, spanZ]);
    maxMag = max(sqrt(Uq(:).^2 + Wq(:).^2));
    if maxMag <= eps
        return;
    end

    Uplot = arrowScale * Uq / maxMag;
    Wplot = arrowScale * Wq / maxMag;

    for iSlice = 1:numel(ySlices)
        Yq = ySlices(iSlice) * ones(size(Xq));
        qh = quiver3(Xq, Yq, Zq, Uplot, zeros(size(Uplot)), Wplot, 0, ...
            'Color', [0.10 0.45 0.85], 'LineWidth', 1.0, 'MaxHeadSize', 0.8);
        if iSlice == 1
            qh.DisplayName = 'Wind vectors';
        else
            qh.HandleVisibility = 'off';
        end
    end
end


