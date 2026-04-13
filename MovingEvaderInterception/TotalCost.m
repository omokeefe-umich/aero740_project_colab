% T: Sampling time ((N_steps + 1) x 1)
% X: States ((N_steps + 1) x d_x)
% U: Inputs (N_steps x d_u)
% Note that the state cost and input cost density functions are multiplied
% by the timestep for adaptive meshing.
% Not the most efficient.
function J = TotalCost(T, X, U, StateCostDensity, InputCostDensity, TerminalStateCost)
    if(size(T, 1) ~= size(X, 1) || size(T, 1) - 1 ~= size(U, 1))
        error("Dimension Mismatch")
    end
    N_steps = size(T, 1) - 1;
    J = 0;
    for i = 1:N_steps
        J = J + (StateCostDensity(X(i, :)) + InputCostDensity(U(i, :))) * (T(i) - T(i - 1));
    end
    J = J + TerminalStateCost(X(N_steps + 1, :));
end
