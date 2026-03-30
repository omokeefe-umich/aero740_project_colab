function H = Jacobian(f, x, epsilon_x)
    if(nargin < 3)
        epsilon_x = 1e-6;
    end
    n_x = size(x, 1);
    n_f = size(f(x), 1);
    H = zeros(n_f, n_x);
    if(size(epsilon_x, 1) == 1)
        epsilon_x = ones(n_x, 1) * epsilon_x;
    else 
        if(size(epsilon_x, 1) ~= size(x, 1))
            error("Dimension Mismatch");
        end
    end
    dx = zeros(n_x, 1);
    for i = 1:n_x
        dx(i) = epsilon_x(i);
        H(:, i) = (f(x + dx) - f(x - dx)) / (2 * epsilon_x(i));
    end
end