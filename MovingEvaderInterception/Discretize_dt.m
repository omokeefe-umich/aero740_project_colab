function [Ads, Bds] = Discretize_dt(dt_by_level, timestep_levels, Ac, Bc)
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
    Ads = zeros(N_steps, N_x, N_x);
    Bds = zeros(N_steps, N_x, N_u);
    Ads_by_level = zeros(N_levels, N_x, N_x);
    Bds_by_level = zeros(N_levels, N_x, N_u);

    sys = ss(Ac, Bc, zeros(1, N_x), zeros(1, N_u));
    for i = 1 : N_levels
        sysd = c2d(sys, dt_by_level(i));
        Ads_by_level(i, :, :) = sysd.A;
        Bds_by_level(i, :, :) = sysd.B;
    end

    for i = 1 : N_steps
        Ads(i, :, :) = Ads_by_level(timestep_levels(i), :, :);
        Bds(i, :, :) = Bds_by_level(timestep_levels(i), :, :);
    end
end