% Returns S and M from Ad and Bd matrices, stored in 3D array Ads and Bds.
function [S, M] = StackedMatrix(Ads, Bds)
    N_steps = size(Ads, 1);
    N_x = size(Ads, 2);
    N_u = size(Bds, 3);
    if(size(Bds, 1) ~= N_steps || size(Bds, 2) ~= N_x || size(Ads, 3) ~= N_x)
        error("Wrong dimension");
    end
    M = zeros(N_steps * N_x, N_x);
    S = zeros(N_steps * N_x, N_steps * N_u);

    M(1:N_x, :) = squeeze(Ads(1, :, :));
    S(1:N_x, 1:N_u) = squeeze(Bds(1, :, :));
    
    for i = 1 : N_steps - 1
        M(N_x * i + 1:N_x * (i + 1), :) = squeeze(Ads(i + 1, :, :)) * M(N_x * (i - 1) + 1:N_x * i, :);
        S(N_x * i + 1:N_x * (i + 1), :) = squeeze(Ads(i + 1, :, :)) * S(N_x * (i - 1) + 1:N_x * i, :);
        S(N_x * i + 1:N_x * (i + 1), N_u * i + 1:N_u * (i + 1)) = squeeze(Bds(1, :, :));
    end
end