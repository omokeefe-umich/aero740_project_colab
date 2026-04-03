% LinTransitionMatrix - Builds MPC prediction matrices M and S for a multi-rate
% linear system xdot = Ac*x + Bc*u.
%
% The prediction model is: X = M*x0 + S*U
% where X is the stacked vector of predicted states and U is the stacked
% vector of control inputs over the horizon.
%
% Inputs:
%   dt_by_level    - (N_levels x 1) timestep size for each subdivision level
%   timestep_levels - (N_steps x 1) subdivision level index for each prediction step
%   Ac             - continuous-time state matrix (N_x x N_x)
%   Bc             - continuous-time input matrix (N_x x N_u)
%
% Outputs:
%   M - state transition matrix mapping x0 to predicted states (N_steps*N_x x N_x)
%   S - input-to-state matrix (N_steps*N_x x N_steps*N_u), lower block-triangular
function [S, M] = LinTransitionMatrix(dt_by_level, timestep_levels, Ac, Bc)
    % Number of distinct timestep sizes available
    N_levels = size(dt_by_level, 1);

    % Validate that all requested timestep levels are within the defined range
    if(any(timestep_levels <= 0) || any(timestep_levels > N_levels))
        error("Invalid Timestep Level")
    end

    % Extract state and input dimensions from system matrices
    N_x = size(Ac, 1);
    N_u = size(Bc, 2);

    % Ensure Ac is square and Bc is compatible with Ac
    if(size(Ac, 2) ~= N_x || size(Bc, 1) ~= N_x)
        error("Dimension Error");
    end

    % Total number of prediction steps in the horizon
    N_steps = size(timestep_levels, 1);

    % Preallocate M and S as zero matrices
    % M: maps initial state x0 to each predicted state over the horizon
    % S: maps stacked input sequence U to predicted states (lower block-triangular)
    M = zeros(N_steps * N_x, N_x);
    S = zeros(N_steps * N_x, N_steps * N_u);

    % Preallocate storage for the discretised A and B matrices at each timestep level
    Ads = zeros(N_levels, N_x, N_x);
    Bds = zeros(N_levels, N_x, N_u);

    % Build a continuous-time state-space model (dummy C/D since only A,B are needed)
    sys = ss(Ac, Bc, zeros(1, N_x), zeros(1, N_u));

    % Discretise the system at each timestep level using zero-order hold (c2d default)
    for i = 1:N_levels
        sysd = c2d(sys, dt_by_level(i));
        Ads(i, :, :) = sysd.A;
        Bds(i, :, :) = sysd.B;
    end

    % Initialise the first block row of M and S for step 1
    % x(1) = Ad(1)*x0 + Bd(1)*u(0)
    M(1:N_x, :) = Ads(timestep_levels(1), :, :);
    S(1:N_x, 1:N_u) = Bds(timestep_levels(1), :, :);
    
    % Recursively build M and S for each subsequent prediction step
    % x(i+1) = Ad(i+1)*x(i) = Ad(i+1)*M(i)*x0
    for i = 1 : N_levels - 1
        % M block row i+1: propagate previous M block row through Ad(i+1)
        M(N_x * i + 1:N_x * (i + 1), :) = Ads(timestep_levels(i + 1), :, :) * M(N_x * (i - 1) + 1:N_x * i, :);

        % S block row i+1: propagate previous M block row through Ad(i+1) (off-diagonal blocks)
        S(N_x * i + 1:N_x * (i + 1), :) = Ads(timestep_levels(i + 1), :, :) * M(N_x * (i - 1) + 1:N_x * i, :);

        % S diagonal block i+1: direct feedthrough of Bd(i+1) for input u(i)
        S(N_x * i + 1:N_x * (i + 1), N_u * i + 1:N_u * (i + 1)) = Bds(timestep_levels(i + 1), :, :);
    end
end