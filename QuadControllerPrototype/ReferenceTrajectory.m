% Generates a reference trajectory for the MPC to follow.
% The trajectory is defined as a 3×N matrix of desired positions [x_ref; y_ref; z_ref]
% at each of the N prediction steps. The MPC will track these absolute positions
% by penalizing the corresponding states in the cost function.
% Inputs:
%   t_horizon (1×N) — time points corresponding to the reference trajectory steps
% Outputs:
%   x_traj (28*N x 1) — stacked reference trajectory for all states (only x,y,z are nonzero)
%   x_ref  (3×N)     — reference positions for plotting and analysis
function [x_traj, x_ref] = ReferenceTrajectory(t_horizon, traj_type)
    if nargin < 2
        traj_type = 'circle';  % default trajectory type
    end
    switch traj_type
        case 'circle'
            x_traj = [zeros(20, 12), ...
                        cos(t_horizon)' - 1, ...    % x reference
                        sin(t_horizon)',     ...    % y reference
                        zeros(20, 14)]';            % all other states = 0
            x_traj = x_traj(:);  % stack into column vector (28*20 x 1)

            % The circular reference used by the MPC (unit circle in XY, z = 0).
            % Replace with any 3×N matrix [x_ref; y_ref; z_ref] for a different path.
            x_ref = [cos(t_horizon) - 1; sin(t_horizon); zeros(1, numel(t_horizon))];
        case 'ground_vehicle'
            % Example: straight line along x-axis at constant speed, z = 0
            x_start = 0; % initial x position (m)
            y_start = 0; % initial y position (m)
            z_start = 0; % initial z position (m)
            vx_start = 2; % starting speed along x-axis (m/s)
            speed = 5; % m/s
            x_traj = [zeros(20, 12), ...
                        (t_horizon - t_horizon(1))' * speed, ...  % x reference (speed m/s)
                        zeros(20, 1), ...                      % y reference
                        zeros(20, 14)]';                       % all other states = 0
            x_traj = x_traj(:);  % stack into column vector (28*20 x 1)

            x_ref = [(t_horizon - t_horizon(1)) * speed; zeros(1, numel(t_horizon)); zeros(1, numel(t_horizon))];    
        otherwise
            error('Unknown trajectory type: %s', traj_type);
    end