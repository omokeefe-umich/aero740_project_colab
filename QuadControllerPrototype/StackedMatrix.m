% Build the MPC prediction matrices M and S for a time-varying linear system.
%
% Given the discrete-time prediction model
%   x_{k+1} = Ad_k * x_k + Bd_k * u_k
% with a different (Ad_k, Bd_k) at each step (e.g., linearization points or
% multi-rate discretizations), this function assembles the batch form:
%
%   X = M * x0 + S * U
%
% where X = [x_1; x_2; ...; x_N]  and  U = [u_0; u_1; ...; u_{N-1}].
%
% The M and S matrices are built recursively:
%   M_1 = Ad_1
%   M_k = Ad_k * M_{k-1}       (k = 2..N)
%
%   S_k = Ad_k * S_{k-1}  then  S_k(col k) += Bd_k   <- direct-feedthrough column
%
% Inputs:
%   Ads (N_steps x N_x x N_x) — stacked discrete A matrices
%   Bds (N_steps x N_x x N_u) — stacked discrete B matrices
%
% Outputs:
%   S (N_steps*N_x x N_steps*N_u) — input sensitivity matrix
%   M (N_steps*N_x x N_x)         — initial-state sensitivity matrix

function [S, M] = StackedMatrix(Ads, Bds)

    N_steps = size(Ads, 1);   % prediction horizon length
    N_x     = size(Ads, 2);   % state dimension
    N_u     = size(Bds, 3);   % input dimension

    % Dimension consistency check
    if(size(Bds, 1) ~= N_steps || size(Bds, 2) ~= N_x || size(Ads, 3) ~= N_x)
        error("Wrong dimension");
    end

    M = zeros(N_steps * N_x, N_x);            % initial-condition propagation
    S = zeros(N_steps * N_x, N_steps * N_u);  % input-to-state sensitivity

    % --- Initialise first prediction step ---
    M(1:N_x, :)      = squeeze(Ads(1, :, :));  % x_1 = Ad_1 * x_0
    S(1:N_x, 1:N_u)  = squeeze(Bds(1, :, :));  % x_1 += Bd_1 * u_0

    % --- Recursively fill remaining prediction steps ---
    for i = 1 : N_steps - 1
        Ad_i = squeeze(Ads(i + 1, :, :));

        % x_{i+1} = Ad_{i+1} * x_i  ->  M_{i+1} = Ad_{i+1} * M_i
        M(N_x*i + 1 : N_x*(i+1), :) = Ad_i * M(N_x*(i-1) + 1 : N_x*i, :);

        % Propagate all previous input columns: S_{i+1} = Ad_{i+1} * S_i
        S(N_x*i + 1 : N_x*(i+1), :) = Ad_i * S(N_x*(i-1) + 1 : N_x*i, :);

        % Add direct influence of u_i on x_{i+1}: S_{i+1, col i} = Bd_1
        % (Note: Bds(1,:,:) is used here; if Bd varies by step replace with Bds(i+1,:,:))
        S(N_x*i + 1 : N_x*(i+1), N_u*i + 1 : N_u*(i+1)) = squeeze(Bds(i+1, :, :));
    end
end