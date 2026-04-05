% Numerical Jacobian via central finite differences
%
% Computes the Jacobian matrix H = df/dx at the point x using two-sided
% (central) differencing.  Central differences achieve O(epsilon^2) accuracy
% vs. O(epsilon) for one-sided differences, at the cost of two function
% evaluations per column.
%
% Inputs:
%   f         — function handle f: R^n -> R^m
%   x         — point at which to evaluate the Jacobian (n x 1)
%   epsilon_x — optional perturbation size(s);
%               scalar applies the same epsilon to all directions;
%               vector (n x 1) allows per-direction tuning (default: 1e-6)
%
% Output:
%   H (m x n) — Jacobian matrix where H(i,j) = df_i/dx_j

function H = Jacobian(f, x, epsilon_x)
    % Default perturbation magnitude
    if(nargin < 3)
        epsilon_x = 1e-6;
    end

    n_x = size(x, 1);         % dimension of the input space
    n_f = size(f(x), 1);      % dimension of the output space
    H   = zeros(n_f, n_x);   % pre-allocate Jacobian

    % Expand scalar epsilon to a per-direction vector
    if(size(epsilon_x, 1) == 1)
        epsilon_x = ones(n_x, 1) * epsilon_x;
    else
        if(size(epsilon_x, 1) ~= size(x, 1))
            error("Dimension Mismatch");
        end
    end

    % Fill each column via central difference:
    %   dH(:,i) ≈ [f(x + e_i*eps) - f(x - e_i*eps)] / (2*eps)
    for i = 1:n_x
        dx    = zeros(n_x, 1);
        dx(i) = epsilon_x(i);        % unit perturbation along direction i
        H(:, i) = (f(x + dx) - f(x - dx)) / (2 * epsilon_x(i));
    end
end