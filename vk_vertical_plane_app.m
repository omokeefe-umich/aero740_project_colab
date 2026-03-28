function vk_vertical_plane_app
% Interactive vertical-plane von Karman wind visualization
% Shows a 2D slice in the x-z plane with horizontal/along-plane wind u
% and vertical wind w.
%
% This uses a spectral-synthesis style approximation:
%   u(x,z) = sum A_u(k) cos(2*pi*(kx*x + kz*z) + phi)
%   w(x,z) = sum A_w(k) cos(2*pi*(kx*x + kz*z) + psi)
%
% with amplitudes shaped by von Karman-like spatial PSDs as a function of
% radial spatial frequency k = sqrt(kx^2 + kz^2).
%
% Intended for intuition / controls / simulation input generation,
% not as a strict high-fidelity atmospheric turbulence model.

    rng(2);

    %% Default parameters
    p.sigma_u = 1.0;     % m/s
    p.sigma_w = 0.5;     % m/s
    p.Lu      = 60;      % m
    p.Lw      = 20;      % m

    % spatial frequencies in cycles/m
    p.kmin    = 1/200;   % 0.005 cyc/m  -> wavelength 200 m
    p.kmax    = 1/8;     % 0.125 cyc/m  -> wavelength 8 m

    p.nModes  = 120;
    p.Lx      = 120;     % m
    p.Lz      = 60;      % m
    p.Nx      = 41;
    p.Nz      = 21;
    p.scale   = 1.0;     % quiver scale factor

    %% Build UI
    fig = uifigure('Name','Von Karman Vertical Plane Wind Viewer',...
        'Position',[80 80 1250 720]);

    gl = uigridlayout(fig,[1 2]);
    gl.ColumnWidth = {320,'1x'};

    left = uipanel(gl,'Title','Controls');
    right = uipanel(gl,'Title','Visualization');

    leftGrid = uigridlayout(left,[14 2]);
    leftGrid.RowHeight = repmat({50},1,14);
    leftGrid.ColumnWidth = {100,'1x'};

    ax = uiaxes(right);
    ax.Position = [30 40 860 620];
    title(ax,'Vertical Plane Wind Field (u,w)');
    xlabel(ax,'x [m]');
    ylabel(ax,'z [m]');
    grid(ax,'on');

    % --- Controls ---
    row = 1;

    uilabel(leftGrid,'Text','\sigma_u [m/s]','HorizontalAlignment','right');
    s_sigma_u = uislider(leftGrid,...
        'Limits',[0.05 5], 'Value',p.sigma_u,...
        'MajorTicks',[0.05 1 2 3 4 5], 'MinorTicks',[]);
    s_sigma_u.Layout.Row = row; s_sigma_u.Layout.Column = 2; row = row + 1;

    uilabel(leftGrid,'Text','\sigma_w [m/s]','HorizontalAlignment','right');
    s_sigma_w = uislider(leftGrid,...
        'Limits',[0.05 5], 'Value',p.sigma_w,...
        'MajorTicks',[0.05 1 2 3 4 5], 'MinorTicks',[]);
    s_sigma_w.Layout.Row = row; s_sigma_w.Layout.Column = 2; row = row + 1;

    uilabel(leftGrid,'Text','L_u [m]','HorizontalAlignment','right');
    s_Lu = uislider(leftGrid,...
        'Limits',[2 300], 'Value',p.Lu,...
        'MajorTicks',[2 25 50 100 200 300], 'MinorTicks',[]);
    s_Lu.Layout.Row = row; s_Lu.Layout.Column = 2; row = row + 1;

    uilabel(leftGrid,'Text','L_w [m]','HorizontalAlignment','right');
    s_Lw = uislider(leftGrid,...
        'Limits',[2 300], 'Value',p.Lw,...
        'MajorTicks',[2 25 50 100 200 300], 'MinorTicks',[]);
    s_Lw.Layout.Row = row; s_Lw.Layout.Column = 2; row = row + 1;

    uilabel(leftGrid,'Text','k_{min} [cyc/m]','HorizontalAlignment','right');
    s_kmin = uislider(leftGrid,...
        'Limits',[1/1000 1/5], 'Value',p.kmin,...
        'MajorTicks',[1/1000 1/200 1/100 1/50 1/20 1/10 1/5],...
        'MajorTickLabels',{'0.001','0.005','0.01','0.02','0.05','0.1','0.2'},...
        'MinorTicks',[]);
    s_kmin.Layout.Row = row; s_kmin.Layout.Column = 2; row = row + 1;

    uilabel(leftGrid,'Text','k_{max} [cyc/m]','HorizontalAlignment','right');
    s_kmax = uislider(leftGrid,...
        'Limits',[1/1000 1], 'Value',p.kmax,...
        'MajorTicks',[1/1000 1/200 1/100 1/50 1/20 1/10 1/5 1/2 1],...
        'MajorTickLabels',{'0.001','0.005','0.01','0.02','0.05','0.1','0.2','0.5','1.0'},...
        'MinorTicks',[]);
    s_kmax.Layout.Row = row; s_kmax.Layout.Column = 2; row = row + 1;

    uilabel(leftGrid,'Text','Number of modes','HorizontalAlignment','right');
    s_nModes = uislider(leftGrid,...
        'Limits',[10 300], 'Value',p.nModes,...
        'MajorTicks',[10 50 100 150 200 250 300], 'MinorTicks',[]);
    s_nModes.Layout.Row = row; s_nModes.Layout.Column = 2; row = row + 1;

    uilabel(leftGrid,'Text','Domain L_x [m]','HorizontalAlignment','right');
    s_Lx = uislider(leftGrid,...
        'Limits',[20 400], 'Value',p.Lx,...
        'MajorTicks',[20 50 100 200 300 400], 'MinorTicks',[]);
    s_Lx.Layout.Row = row; s_Lx.Layout.Column = 2; row = row + 1;

    uilabel(leftGrid,'Text','Domain L_z [m]','HorizontalAlignment','right');
    s_Lz = uislider(leftGrid,...
        'Limits',[10 200], 'Value',p.Lz,...
        'MajorTicks',[10 25 50 100 150 200], 'MinorTicks',[]);
    s_Lz.Layout.Row = row; s_Lz.Layout.Column = 2; row = row + 1;

    uilabel(leftGrid,'Text','Quiver scale','HorizontalAlignment','right');
    s_scale = uislider(leftGrid,...
        'Limits',[0.2 4], 'Value',p.scale,...
        'MajorTicks',[0.2 0.5 1 2 3 4], 'MinorTicks',[]);
    s_scale.Layout.Row = row; s_scale.Layout.Column = 2; row = row + 1;

    btnReseed = uibutton(leftGrid,'push','Text','Reseed Random Phases');
    btnReseed.Layout.Row = row; btnReseed.Layout.Column = [1 2]; row = row + 1;

    lblInfo = uitextarea(leftGrid,...
        'Value',{'Move sliders to update the field.',...
                 'u = along-plane wind', ...
                 'w = vertical wind', ...
                 'Magnitude contour is sqrt(u^2+w^2).', ...
                 'This is a 2D spectral visualization, not a full 3D atmosphere model.'},...
        'Editable','off');
    lblInfo.Layout.Row = [row 14];
    lblInfo.Layout.Column = [1 2];

    % Store random phases/directions
    phaseState = struct();
    refreshRandomState();

    % Initial plot
    updatePlot();

    % Callbacks
    s_sigma_u.ValueChangingFcn = @(src,event) updatePlot(event.Value, []);
    s_sigma_w.ValueChangingFcn = @(src,event) updatePlot([], event.Value);
    s_Lu.ValueChangingFcn      = @(~,~) updatePlot();
    s_Lw.ValueChangingFcn      = @(~,~) updatePlot();
    s_kmin.ValueChangingFcn    = @(~,~) updatePlot();
    s_kmax.ValueChangingFcn    = @(~,~) updatePlot();
    s_nModes.ValueChangingFcn  = @(~,~) updatePlot();
    s_Lx.ValueChangingFcn      = @(~,~) updatePlot();
    s_Lz.ValueChangingFcn      = @(~,~) updatePlot();
    s_scale.ValueChangingFcn   = @(~,~) updatePlot();

    s_sigma_u.ValueChangedFcn = @(~,~) updatePlot();
    s_sigma_w.ValueChangedFcn = @(~,~) updatePlot();
    s_Lu.ValueChangedFcn      = @(~,~) updatePlot();
    s_Lw.ValueChangedFcn      = @(~,~) updatePlot();
    s_kmin.ValueChangedFcn    = @(~,~) updatePlot();
    s_kmax.ValueChangedFcn    = @(~,~) updatePlot();
    s_nModes.ValueChangedFcn  = @(~,~) refreshAndPlot();
    s_Lx.ValueChangedFcn      = @(~,~) updatePlot();
    s_Lz.ValueChangedFcn      = @(~,~) updatePlot();
    s_scale.ValueChangedFcn   = @(~,~) updatePlot();

    btnReseed.ButtonPushedFcn = @(~,~) refreshAndPlot();

    %% nested functions
    function refreshAndPlot()
        refreshRandomState();
        updatePlot();
    end

    function refreshRandomState()
        nModes = round(s_nModes.Value);

        phaseState.theta = 2*pi*rand(nModes,1);  % direction angle in x-z plane
        phaseState.phi_u = 2*pi*rand(nModes,1);  % phase for u
        phaseState.phi_w = 2*pi*rand(nModes,1);  % phase for w
        phaseState.signs = sign(randn(nModes,2));
        phaseState.signs(phaseState.signs==0) = 1;
    end

    function updatePlot(varargin)
        if ~isempty(varargin)
            if numel(varargin) >= 1 && ~isempty(varargin{1})
                s_sigma_u.Value = varargin{1};
            end
            if numel(varargin) >= 2 && ~isempty(varargin{2})
                s_sigma_w.Value = varargin{2};
            end
        end

        pnow.sigma_u = s_sigma_u.Value;
        pnow.sigma_w = s_sigma_w.Value;
        pnow.Lu      = s_Lu.Value;
        pnow.Lw      = s_Lw.Value;
        pnow.kmin    = s_kmin.Value;
        pnow.kmax    = s_kmax.Value;
        pnow.nModes  = round(s_nModes.Value);
        pnow.Lx      = s_Lx.Value;
        pnow.Lz      = s_Lz.Value;
        pnow.scale   = s_scale.Value;
        pnow.Nx      = p.Nx;
        pnow.Nz      = p.Nz;

        if pnow.kmax <= pnow.kmin
            pnow.kmax = pnow.kmin * 1.05;
        end

        [X,Z,U,W,Umag] = synthesizeVKField2D(pnow, phaseState);

        cla(ax);

        % magnitude background
        contourf(ax, X, Z, Umag, 18, 'LineColor', 'none');
        hold(ax,'on');
        colormap(ax, turbo);
        caxis(ax, [0 3]);   % <-- set constant magnitude range in m/s
        cb = colorbar(ax);
        cb.Label.String = '|V_{wind}| [m/s]';

        % quiver
        quiver(ax, X, Z, U, W, pnow.scale, 'k', 'LineWidth', 1.0);

        xlabel(ax,'x [m]');
        ylabel(ax,'z [m]');
        title(ax, sprintf(['Vertical Plane VK Wind Field | ', ...
            '\\sigma_u = %.2f m/s, \\sigma_w = %.2f m/s, ', ...
            'L_u = %.1f m, L_w = %.1f m'], ...
            pnow.sigma_u, pnow.sigma_w, pnow.Lu, pnow.Lw));

        axis(ax,'equal');
        xlim(ax,[0 pnow.Lx]);
        ylim(ax,[0 pnow.Lz]);
        grid(ax,'on');
        hold(ax,'off');
    end
end

%% ------------------------------------------------------------------------
function [X,Z,U,W,Umag] = synthesizeVKField2D(p, phaseState)
% Build a vertical-plane (x-z) wind field from VK-shaped spatial spectra.

    x = linspace(0, p.Lx, p.Nx);
    z = linspace(0, p.Lz, p.Nz);
    [X,Z] = meshgrid(x,z);

    U = zeros(size(X));
    W = zeros(size(X));

    nModes = p.nModes;

    % Sample radial spatial frequencies (cycles/m)
    k = linspace(p.kmin, p.kmax, nModes).';

    % Convert to angular spatial frequency (rad/m)
    kap = 2*pi*k;

    % Random direction in x-z plane
    theta = phaseState.theta(1:nModes);
    kx = k .* cos(theta);
    kz = k .* sin(theta);

    % Approximate bin width
    if numel(k) > 1
        dk = mean(diff(k));
    else
        dk = max(1e-6, 0.05 * max(k,1e-3));
    end

    % VK spatial PSDs as function of angular spatial frequency kap
    % Longitudinal-like component
    Phi_u = p.sigma_u^2 .* (2*p.Lu/pi) .* (1 + (1.339*p.Lu*kap).^2).^(-5/6);

    % Vertical-like component
    Phi_w = p.sigma_w^2 .* (p.Lw/pi) .* ...
        (1 + (8/3)*(1.339*p.Lw*kap).^2) ./ (1 + (1.339*p.Lw*kap).^2).^(11/6);

    % Convert PSD to discrete mode amplitude
    % Practical visualization scaling:
    % amplitude ~ sqrt(2 * PSD * d(kappa))
    dkap = 2*pi*dk;
    Au = sqrt(2 * max(Phi_u,0) * dkap);
    Aw = sqrt(2 * max(Phi_w,0) * dkap);

    phi_u = phaseState.phi_u(1:nModes);
    phi_w = phaseState.phi_w(1:nModes);

    % Build field
    for i = 1:nModes
        phase_i = 2*pi*(kx(i)*X + kz(i)*Z);

        U = U + Au(i) * cos(phase_i + phi_u(i));
        W = W + Aw(i) * cos(phase_i + phi_w(i));
    end

    % Remove means
    U = U - mean(U(:));
    W = W - mean(W(:));

    % Normalize to approximately match target std devs
    stdU = std(U(:));
    stdW = std(W(:));

    if stdU > 1e-10
        U = U * (p.sigma_u / stdU);
    end
    if stdW > 1e-10
        W = W * (p.sigma_w / stdW);
    end

    Umag = sqrt(U.^2 + W.^2);
end