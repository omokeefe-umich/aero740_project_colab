% x_traj: "Desired trajectory" from x_1 to x_N
% u_traj: "Desired input" from u_0 to u_N-1
% d_traj: Optional affine prediction offset from x_1 to x_N
%         (e.g., lifted mapping offsets such that X = S*du + M*x0 + d_traj)
% How do we weigh each timestep when it's discretised into different
% intervals?
function [H, q] = QPFormat(Q_density, R_density, P, S, M, t, x0, x_traj, u_traj, d_traj)
    N = size(t, 1) - 1;
    n_x = size(Q_density, 1);
    n_u = size(R_density, 1);
    
    if nargin < 10 || isempty(d_traj)
        d_traj = zeros(N * n_x, 1);
    end
    
    % --- Dimension checks ---
    if size(x_traj, 1) ~= N * n_x
        error("State Trajectory Size Mismatch. Expected %d, got %d", N * n_x, size(x_traj, 1));
    elseif size(u_traj, 1) ~= N * n_u
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
    elseif size(d_traj, 1) ~= N * n_x || size(d_traj, 2) ~= 1
        error("Affine Prediction Offset Size Mismatch. Expected %d, got %d", N * n_x, size(d_traj, 1));
    end
    delta_t = t(2:end) - t(1:end-1);
    Q_bar = blkdiag(kron(diag(delta_t(2:end)), Q_density), P);
    R_bar = kron(diag(delta_t), R_density);
    H = S' * Q_bar * S + R_bar;
    q = S' * Q_bar * (M * x0 + d_traj - x_traj) - R_bar * u_traj;
end