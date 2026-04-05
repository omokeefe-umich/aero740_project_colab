% Build QP cost matrices for MPC: min (1/2)*U'*H*U + q'*U
%
% J = sum_{k=1}^{N-1} x_k'*Q*x_k + x_N'*P*x_N + sum_{k=1}^{N} u_k'*R*u_k
%
% Inputs:
%   Q_density, R_density — stage state/input cost matrices
%   P    — terminal state cost
%   S, M — prediction matrices (X = M*x0 + S*U)
%   t    — (N+1)x1 time nodes
%   x0   — current augmented state [dx; x; u; r]
%   u_traj — stacked nominal input trajectory
%
% Outputs:
%   H — QP Hessian (N*n_u x N*n_u)
%   q — QP linear term (N*n_u x 1)

function [H, q] = QPFormat(Q_density, R_density, P, S, M, t, x0, u_traj)

    N   = size(t, 1) - 1;         % number of prediction steps
    n_x = size(Q_density, 1);     % state dimension
    n_u = size(R_density, 1);     % input dimension

    % --- Dimension checks ---
    if size(u_traj, 1) ~= N * n_u  
        error("Input Trajectory Size Mismatch. Expected %d, got %d", N * n_u, size(u_traj, 1));
    elseif size(Q_density, 2) ~= n_x 
        error("State Cost Matrix Size Mismatch. Expected %d, got %d", n_x, size(Q_density, 2));
    elseif size(R_density, 2) ~= n_u  
        error("Input Cost Matrix Size Mismatch. Expected %d, got %d", n_u, size(R_density, 2));
    elseif size(P, 1) ~= n_x 
        error("Terminal Cost Matrix Size Mismatch. Expected %d, got %d", n_x, size(P, 1));
    elseif size(P, 2) ~= n_x 
        error("Terminal Cost Matrix Size Mismatch. Expected %d, got %d", n_x, size(P, 2));
    elseif size(x0, 1) ~= n_x 
        error("Initial State Size Mismatch. Expected %d, got %d", n_x, size(x0, 1));
    elseif size(S, 1) ~= N * n_x 
        error("State Sensitivity Matrix Size Mismatch. Expected %d, got %d", N * n_x, size(S, 1));
    elseif size(S, 2) ~= N * n_u 
        error("State Sensitivity Matrix Size Mismatch. Expected %d, got %d", N * n_u, size(S, 2));
    elseif size(M, 1) ~= N * n_x 
        error("Initial-State-to-Prediction Matrix Size Mismatch. Expected %d, got %d", N * n_x, size(M, 1));
    elseif size(M, 2) ~= n_x
        error("Initial-State-to-Prediction Matrix Size Mismatch. Expected %d, got %d", n_x, size(M, 2));
    end

    delta_t = ones(length(t(2:end)), 1);

    Q_bar = blkdiag(kron(diag(delta_t(1:end-1)), Q_density), P);
    R_bar = kron(diag(delta_t), R_density);

    H = S' * Q_bar * S + R_bar;
    q = S' * Q_bar * M * x0 - R_bar * u_traj;

    H = (H + H') / 2;  % enforce symmetry

    % Regularize if not PSD
    [~, p] = chol(H);
    if p > 0
        fprintf("Warning: H not PSD, adding regularization.\n");
        H = H + 1e-6 * eye(size(H));
    end
end