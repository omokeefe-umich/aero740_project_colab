% Build QP cost matrices (H, q) for a Model Predictive Control problem.
%
% The MPC optimal control problem minimises over the input sequence U:
%
%   J = (1/2) * U' * H * U  +  q' * U
%
% where the continuous cost is approximated by a time-weighted sum:
%
%   J = sum_{k=1}^{N-1} delta_t_k * (x_k - x_ref_k)' * Q * (x_k - x_ref_k)
%     + (x_N - x_ref_N)' * P * (x_N - x_ref_N)          <- terminal cost
%     + sum_{k=0}^{N-1} delta_t_k * (u_k - u_ref_k)' * R * (u_k - u_ref_k)
%
% The prediction model is  X = M * x0 + S * U  (from StackedMatrix).
%
% Inputs:
%   Q_density (n_x x n_x) — stage state cost per unit time
%   R_density (n_u x n_u) — stage input cost per unit time
%   P         (n_x x n_x) — terminal state cost (not time-weighted)
%   S         (N*n_x x N*n_u) — input-to-state sensitivity matrix
%   M         (N*n_x x n_x)   — initial-state-to-prediction matrix
%   t         ((N+1) x 1) — time nodes; intervals are delta_t = diff(t)
%   x0        (n_x x 1)   — current augmented state
%   x_traj    (N*n_x x 1) — stacked desired state trajectory x_1…x_N
%   u_traj    (N*n_u x 1) — stacked desired input trajectory u_0…u_{N-1}
%
% Outputs:
%   H (N*n_u x N*n_u) — positive semi-definite QP cost matrix
%   q (N*n_u x 1)     — linear QP cost vector

function [H, q] = QPFormat(Q_density, R_density, P, S, M, t, x0, x_traj, u_traj)

    N   = size(t, 1) - 1;         % number of prediction steps
    n_x = size(Q_density, 1);     % state dimension
    n_u = size(R_density, 1);     % input dimension

    % --- Dimension checks ---
    if(size(x_traj, 1) ~= N * n_x || size(u_traj, 1) ~= N * n_u ...
            || size(Q_density, 2) ~= n_x || size(R_density, 2) ~= n_u ...
            || size(P, 1) ~= n_x || size(P, 2) ~= n_x || size(x0, 1) ~= n_x ...
            || size(S, 1) ~= N * n_x || size(S, 2) ~= N * n_u ...
            || size(M, 1) ~= N * n_x || size(M, 2) ~= n_x)
        error("Sample Number Mismatch");
    end

    % Inter-sample intervals (used as quadrature weights)
    delta_t = t(2:end) - t(1:end-1);  % N x 1

    % Block-diagonal state cost: time-weight stages 2…N, then terminal cost P
    % Stage 1 is implicitly handled via M*x0 below.
    Q_bar = blkdiag(kron(diag(delta_t(2:end)), Q_density), P);  % (N*n_x x N*n_x)

    % Block-diagonal input cost: time-weight all N stages
    R_bar = kron(diag(delta_t), R_density);                      % (N*n_u x N*n_u)

    % QP Hessian:  H = S'*Q_bar*S + R_bar
    H = S' * Q_bar * S + R_bar;

    % QP linear term: q arises from expanding (X - x_ref)'*Q*(X - x_ref)
    % = (M*x0 + S*U - x_traj)'*Q_bar*(…) and collecting the U-linear terms,
    % minus the input reference contribution from R_bar.
    q = S' * Q_bar * (M * x0 - x_traj) - R_bar * u_traj;
end