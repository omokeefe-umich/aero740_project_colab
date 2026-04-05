function PlotResults(t, Xs, x_ref)
%PLOTRESULTS  Visualize quadcopter MPC simulation results against a reference.
%
% Generates five figures:
%   1. Top-down XY trajectory vs reference path
%   2. 3-D trajectory vs reference path
%   3. Altitude z vs ground-plane range  r = sqrt(x^2 + y^2)
%   4. Translational velocities and numerical accelerations vs time
%   5. Euler angles and angular rates vs time
%
% A filled green circle marks the starting state and a filled red circle
% marks the ending state on all spatial trajectory figures (1-3).
%
% Inputs:
%   t     - time vector (1×N or N×1)  [s]
%   Xs    - 12×N state history, rows ordered as:
%             [x; y; z; phi; theta; psi; xdot; ydot; zdot; phidot; thetadot; psidot]
%   x_ref - 3×M reference trajectory  [x_ref; y_ref; z_ref]
%             M may differ from N; the reference is plotted as a continuous
%             spatial path independent of its time parameterisation.

    t = t(:)';    % enforce row vector

    % ------------------------------------------------------------------ %
    %  Unpack state rows
    % ------------------------------------------------------------------ %
    x    = Xs(1,:);   y    = Xs(2,:);   z    = Xs(3,:);
    phi  = Xs(4,:);   th   = Xs(5,:);   ps   = Xs(6,:);
    xd   = Xs(7,:);   yd   = Xs(8,:);   zd   = Xs(9,:);
    phid = Xs(10,:);  thd  = Xs(11,:);  psd  = Xs(12,:);

    % Ground-plane range (horizontal distance from origin)
    r = sqrt(x.^2 + y.^2);

    % Numerical accelerations via central finite differences.
    % gradient() uses central differences at interior points and one-sided
    % differences at the endpoints, respecting non-uniform spacing in t.
    xdd = gradient(xd, t);
    ydd = gradient(yd, t);
    zdd = gradient(zd, t);

    % Reference trajectory components
    xr = x_ref(1,:);
    yr = x_ref(2,:);
    zr = x_ref(3,:);
    rr = sqrt(xr.^2 + yr.^2);

    % Colour scheme for per-axis time-series (x=blue, y=red, z=green)
    clr = {'b', 'r', [0 0.55 0]};

    % ================================================================== %
    %  Figure 1 – XY Trajectory (Top View)
    % ================================================================== %
    figure('Name', 'XY Trajectory', 'NumberTitle', 'off');
    hold on;  grid on;  axis equal;

    plot(xr, yr, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
    plot(x,  y,  'k-',  'LineWidth', 1.5, 'DisplayName', 'Simulated');

    % Start: filled green circle
    scatter(x(1), y(1), 110, 'g', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'Start');

    % End: filled red circle
    scatter(x(end), y(end), 110, 'r', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'End');

    xlabel('x  [m]');  ylabel('y  [m]');
    title('Top-Down Trajectory (XY Plane)');
    legend('Location', 'best');

    % ================================================================== %
    %  Figure 2 – 3-D Trajectory
    % ================================================================== %
    figure('Name', '3-D Trajectory', 'NumberTitle', 'off');
    hold on;  grid on;

    plot3(xr, yr, zr, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
    plot3(x,  y,  z,  'k-',  'LineWidth', 1.5, 'DisplayName', 'Simulated');

    % Start: filled green circle
    scatter3(x(1), y(1), z(1), 110, 'g', 'filled', ...
             'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'Start');

    % End: filled red circle
    scatter3(x(end), y(end), z(end), 110, 'r', 'filled', ...
             'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'End');

    xlabel('x  [m]');  ylabel('y  [m]');  zlabel('z  [m]');
    title('3-D Trajectory');
    legend('Location', 'best');
    view(45, 30);

    % ================================================================== %
    %  Figure 3 – Altitude vs Ground-Plane Range
    % ================================================================== %
    figure('Name', 'z vs Range', 'NumberTitle', 'off');
    hold on;  grid on;

    plot(rr, zr, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
    plot(r,  z,  'k-',  'LineWidth', 1.5, 'DisplayName', 'Simulated');

    scatter(r(1), z(1), 110, 'g', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'Start');
    scatter(r(end), z(end), 110, 'r', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'End');

    xlabel('r = \surd(x^2 + y^2)  [m]');  ylabel('z  [m]');
    title('Altitude vs Ground-Plane Range');
    legend('Location', 'best');

    % ================================================================== %
    %  Figure 4 – Translational Velocities & Accelerations
    % ================================================================== %
    figure('Name', 'Translational Kinematics', 'NumberTitle', 'off');

    vel_data   = {xd;  yd;  zd};
    accel_data = {xdd; ydd; zdd};
    vel_labels = {'\dot{x}  [m/s]',    '\dot{y}  [m/s]',    '\dot{z}  [m/s]'};
    acc_labels = {'\ddot{x}  [m/s^2]', '\ddot{y}  [m/s^2]', '\ddot{z}  [m/s^2]'};

    for i = 1:3
        subplot(2, 3, i);
        plot(t, vel_data{i}, 'Color', clr{i}, 'LineWidth', 1.2);
        grid on;
        xlabel('t  [s]');  ylabel(vel_labels{i});
        title(vel_labels{i});
        xlim([t(1) t(end)]);

        subplot(2, 3, i+3);
        plot(t, accel_data{i}, 'Color', clr{i}, 'LineWidth', 1.2);
        grid on;
        xlabel('t  [s]');  ylabel(acc_labels{i});
        title(acc_labels{i});
        xlim([t(1) t(end)]);
    end

    sgtitle('Translational Velocities (top row) and Accelerations (bottom row)');

    % ================================================================== %
    %  Figure 5 – Euler Angles & Angular Rates
    % ================================================================== %
    figure('Name', 'Rotational Kinematics', 'NumberTitle', 'off');

    angle_data = {phi;  th;  ps};
    rate_data  = {phid; thd; psd};
    angle_labels = {'\phi  [rad]',         '\theta  [rad]',         '\psi  [rad]'};
    rate_labels  = {'\dot{\phi}  [rad/s]', '\dot{\theta}  [rad/s]', '\dot{\psi}  [rad/s]'};

    for i = 1:3
        subplot(2, 3, i);
        plot(t, angle_data{i}, 'Color', clr{i}, 'LineWidth', 1.2);
        grid on;
        xlabel('t  [s]');  ylabel(angle_labels{i});
        title(angle_labels{i});
        xlim([t(1) t(end)]);

        subplot(2, 3, i+3);
        plot(t, rate_data{i}, 'Color', clr{i}, 'LineWidth', 1.2);
        grid on;
        xlabel('t  [s]');  ylabel(rate_labels{i});
        title(rate_labels{i});
        xlim([t(1) t(end)]);
    end

    sgtitle('Euler Angles (top row) and Angular Rates (bottom row)');

end


