% Discretize a continuous-time LTI system at multiple time-step sizes.
%
% This function supports multi-rate MPC horizons where different prediction
% steps may use different sample intervals (e.g., fine grid near present,
% coarser grid further ahead).  Each unique interval is discretized once
% via zero-order hold (ZOH) using MATLAB's c2d, then copied to the output
% arrays according to the per-step level assignments.
%
% Continuous model:  x_dot = Ac*x + Bc*u
% Discrete model:    x_{k+1} = Ad_k*x_k + Bd_k*u_k
%
% Inputs:
%   dt_by_level    (N_levels x 1) — list of unique sample intervals [s]
%   timestep_levels (N_steps x 1) — index into dt_by_level for each prediction step
%   Ac              (N_x x N_x)   — continuous-time A matrix
%   Bc              (N_x x N_u)   — continuous-time B matrix
%
% Outputs:
%   Ads (N_steps x N_x x N_x) — discrete A matrices for each prediction step
%   Bds (N_steps x N_x x N_u) — discrete B matrices for each prediction step

function [Ads, Bds] = Discretize_dt(dt_by_level, timestep_levels, Ac, Bc)

    N_levels = size(dt_by_level, 1);   % number of distinct sample intervals

    % Validate level indices
    if(any(timestep_levels <= 0) || any(timestep_levels > N_levels))
        error("Invalid Timestep Level");
    end

    N_x = size(Ac, 1);   % state dimension
    N_u = size(Bc, 2);   % input dimension

    % Check that Ac is square and Bc rows match Ac
    if(size(Ac, 2) ~= N_x || size(Bc, 1) ~= N_x)
        error("Dimension Error");
    end

    N_steps = size(timestep_levels, 1);   % total prediction horizon length

    % Pre-allocate output arrays
    Ads = zeros(N_steps, N_x, N_x);
    Bds = zeros(N_steps, N_x, N_u);

    % Pre-allocate per-level discretized matrices (compute each level once)
    Ads_by_level = zeros(N_levels, N_x, N_x);
    Bds_by_level = zeros(N_levels, N_x, N_u);

    % Build a state-space object (dummy C/D since only A,B are needed)
    sys = ss(Ac, Bc, zeros(1, N_x), zeros(1, N_u));

    % --- Discretize at each unique sample interval using ZOH ---
    for i = 1 : N_levels
        sysd = c2d(sys, dt_by_level(i));   % zero-order hold discretization
        Ads_by_level(i, :, :) = sysd.A;
        Bds_by_level(i, :, :) = sysd.B;
    end

    % --- Assign discretized matrices to each prediction step by level ---
    for i = 1 : N_steps
        Ads(i, :, :) = Ads_by_level(timestep_levels(i), :, :);
        Bds(i, :, :) = Bds_by_level(timestep_levels(i), :, :);
    end
end