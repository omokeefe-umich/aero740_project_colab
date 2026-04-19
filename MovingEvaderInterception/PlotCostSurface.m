%% Plot Cost Surface
% Plot six 2D contour slices of the MPC quadratic cost over selected pairs
% of the first-step control increments du(1:4), while all other decision
% variables are held at the current optimizer solution u_star.
%
% Optional movie capture fields in plotState:
%   captureMovie (logical)
%   moviePath (char/string)
%   movieFrameRate (scalar)
function plotState = PlotCostSurface(H, q, u_traj, u_star, plotState, k, currentTime)
    if nargin < 4 || isempty(u_star)
        error('u_star is required and must be the stacked decision vector.');
    end
    if nargin < 5 || isempty(plotState)
        plotState = struct();
    end
    if nargin < 6
        k = NaN;
    end
    if nargin < 7
        currentTime = NaN;
    end

    if numel(u_traj) ~= 4
        error('u_traj must be 4x1 (current physical control state).');
    end
    if numel(u_star) < 4
        error('u_star must contain at least the first 4 entries.');
    end

    pairList = [1 2; 2 3; 3 4; 1 3; 1 4; 2 4];
    pairNames = {'F_z vs \tau_{\phi}', '\tau_{\phi} vs \tau_{\theta}', '\tau_{\theta} vs \tau_{\psi}', ...
                 'F_z vs \tau_{\theta}', 'F_z vs \tau_{\psi}', '\tau_{\phi} vs \tau_{\psi}'};
    uNames = {'F_z', '\tau_{\phi}', '\tau_{\theta}', '\tau_{\psi}'};

    nGrid = 35;
    halfWidth = 5;

    if ~isfield(plotState, 'fig') || ~isgraphics(plotState.fig)
        plotState.fig = figure('Name', 'MPC Cost Contour Slices', 'NumberTitle', 'off');
        plotState.tiled = tiledlayout(plotState.fig, 2, 3, ...
            'TileSpacing', 'compact', 'Padding', 'compact');
        plotState.ax = gobjects(size(pairList, 1), 1);
        for p = 1:size(pairList, 1)
            plotState.ax(p) = nexttile(plotState.tiled, p);
        end
    end

    for p = 1:size(pairList, 1)
        idxA = pairList(p, 1);
        idxB = pairList(p, 2);

        centerA = u_star(idxA);
        centerB = u_star(idxB);
        rangeA = linspace(centerA - halfWidth, centerA + halfWidth, nGrid);
        rangeB = linspace(centerB - halfWidth, centerB + halfWidth, nGrid);
        [UA, UB] = meshgrid(rangeA, rangeB);

        Cost = zeros(size(UA));
        for i = 1:size(UA, 1)
            for j = 1:size(UA, 2)
                du_test = u_star;
                du_test(idxA) = UA(i, j);
                du_test(idxB) = UB(i, j);
                Cost(i, j) = 0.5 * du_test' * H * du_test + q' * du_test;
            end
        end

        ax = plotState.ax(p);
        cla(ax);
        contourf(ax, UA, UB, Cost, 20, 'LineStyle', 'none');
        hold(ax, 'on');
        plot(ax, centerA, centerB, 'rx', 'LineWidth', 2, 'MarkerSize', 10);
        hold(ax, 'off');
        colormap(ax, 'parula');
        colorbar(ax);
        grid(ax, 'on');
        xlabel(ax, uNames{idxA});
        ylabel(ax, uNames{idxB});
        title(ax, pairNames{p});
    end

    J_star = 0.5 * u_star' * H * u_star + q' * u_star;
    title(plotState.tiled, sprintf(['MPC cost contour slices at step %d, t = %.2f s, ' ...
        'J(U*_MPC) = %.4g, u = [%.3f %.3f %.3f %.3f]'], ...
        k, currentTime, J_star, u_traj(1), u_traj(2), u_traj(3), u_traj(4)));

    drawnow limitrate nocallbacks;

    if isfield(plotState, 'captureMovie') && plotState.captureMovie
        if ~isfield(plotState, 'moviePath') || isempty(plotState.moviePath)
            plotState.moviePath = fullfile(pwd, 'cost_surface_movie.mp4');
        end
        if ~isfield(plotState, 'movieFrameRate') || isempty(plotState.movieFrameRate)
            plotState.movieFrameRate = 10;
        end

        if ~isfield(plotState, 'writer') || isempty(plotState.writer)
            writerObj = VideoWriter(plotState.moviePath, 'MPEG-4');
            writerObj.FrameRate = plotState.movieFrameRate;
            open(writerObj);
            plotState.writer = writerObj;
        end

        frame = getframe(plotState.fig);
        writeVideo(plotState.writer, frame);
    end
end
