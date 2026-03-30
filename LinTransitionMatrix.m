% dt_by_level: defines the timestep size for each subdivision level.
% timestep_levels: subdivision level of each timestep.
% Discretises xdot = Ac x + Bc u for various timesteps, and returns S and M
% matrices
function [S, M] = LinTransitionMatrix(dt_by_level, timestep_levels, Ac, Bc)
    N_levels = size(dt_by_level, 1);
    if(any(timestep_levels <= 0) || any(timestep_levels > N_levels))
        error("Invalid Timestep Level")
    end
    N_x = size(Ac, 1);
    N_u = size(Bc, 2);
    if(size(Ac, 2) ~= N_x || size(Bc, 1) ~= N_x)
        error("Dimension Error");
    end
    N_steps = size(timestep_levels, 1);
    M = zeros(N_steps * N_x, N_x);
    S = zeros(N_steps * N_x, N_steps * N_u);

    Ads = zeros(N_levels, N_x, N_x);
    Bds = zeros(N_levels, N_x, N_u);

    sys = ss(Ac, Bc, zeros(1, N_x), zeros(1, N_u));
    for i = 1:N_levels
        sysd = c2d(sys, dt_by_level(i));
        Ads(i, :, :) = sysd.A;
        Bds(i, :, :) = sysd.B;
    end

    M(1:N_x, :) = Ads(timestep_levels(1), :, :);
    S(1:N_x, 1:N_u) = Bds(timestep_levels(1), :, :);
    
    for i = 1 : N_levels - 1
        M(N_x * i + 1:N_x * (i + 1), :) = Ads(timestep_levels(i + 1), :, :) * M(N_x * (i - 1) + 1:N_x * i, :);
        S(N_x * i + 1:N_x * (i + 1), :) = Ads(timestep_levels(i + 1), :, :) * M(N_x * (i - 1) + 1:N_x * i, :);
        S(N_x * i + 1:N_x * (i + 1), N_u * i + 1:N_u * (i + 1)) = Bds(timestep_levels(i + 1), :, :);
    end
end