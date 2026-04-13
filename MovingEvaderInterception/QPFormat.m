% x_traj: "Desired trajectory" from x_1 to x_N
% u_traj: "Desired input" from u_0 to u_N-1
% How do we weigh each timestep when it's discretised into different
% intervals?
function [H, q] = QPFormat(Q_density, R_density, P, S, M, t, x0, x_traj, u_traj)
    N = size(t, 1) - 1;
    n_x = size(Q_density, 1);
    n_u = size(R_density, 1);
    if(size(x_traj, 1) ~= N * n_x || size(u_traj, 1) ~= N * n_u ...
            || size(Q_density, 2) ~= n_x || size(R_density, 2) ~= n_u ...
            || size(P, 1) ~= n_x || size(P, 2) ~= n_x || size(x0, 1) ~= n_x...
            || size(S, 1) ~= N * n_x || size(S, 2) ~= N * n_u...
            || size(M, 1) ~= N * n_x || size(M, 2) ~= n_x)
        error("Sample Number Mismatch");
    end
    delta_t = t(2:end) - t(1:end-1);
    Q_bar = blkdiag(kron(diag(delta_t(2:end)), Q_density), P);
    R_bar = kron(diag(delta_t), R_density);
    H = S' * Q_bar * S + R_bar;
    q = S' * Q_bar * (M * x0 - x_traj) - R_bar * u_traj;
end