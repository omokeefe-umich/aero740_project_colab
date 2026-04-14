function wind_data = DisturbanceModel(wind_on)
    % This function generates the wind field and returns it as a structure
    if nargin < 1
        wind_on = 1;
    end

    % Domain Parameters
    Lx = 120; Lz = 60; Nx = 41; Nz = 21;

    if ~wind_on
        x = linspace(0, Lx, Nx);
        z = linspace(0, Lz, Nz);
        [X, Z] = meshgrid(x, z);

        wind_data.X = X;
        wind_data.Z = Z;
        wind_data.U = zeros(size(X));
        wind_data.W = zeros(size(X));
        wind_data.Umag = zeros(size(X));
        return
    end

    % Chosen Disturbance Parameters - These can be changed if needed
    p.sigma_u = 1.50; p.sigma_w = 0.9;
    p.Lu = 60.0; p.Lw = 50.0;
    p.kmin = 0.005; p.kmax = 0.2;
    p.nModes = 120;

    % Fixed Random Seed - Ensures the disturbance profile remains constant
    % between runs
    rng(2);

    % Pre-calculate Random States
    phaseState.theta = 2*pi*rand(p.nModes,1);
    phaseState.phi_u = 2*pi*rand(p.nModes,1);
    phaseState.phi_w = 2*pi*rand(p.nModes,1);

    % Synthesize Field
    [X, Z, U, W, Umag] = WindField(p, phaseState, Lx, Lz, Nx, Nz);

    % Pack into a struct for the MPC script
    wind_data.X = X;
    wind_data.Z = Z;
    wind_data.U = U;
    wind_data.W = W;
    wind_data.Umag = Umag;
end

function [X, Z, U, W, Umag] = WindField(p, phaseState, Lx, Lz, Nx, Nz)
    x = linspace(0, Lx, Nx);
    z = linspace(0, Lz, Nz);
    [X, Z] = meshgrid(x,z);
    U = zeros(size(X)); W = zeros(size(X));
    k = linspace(p.kmin, p.kmax, p.nModes).';
    kap = 2*pi*k;
    kx = k .* cos(phaseState.theta);
    kz = k .* sin(phaseState.theta);
    dkap = 2*pi*mean(diff(k));
    
    Phi_u = p.sigma_u^2 * (2*p.Lu/pi) * (1 + (1.339*p.Lu*kap).^2).^(-5/6);
    Phi_w = p.sigma_w^2 * (p.Lw/pi) * (1 + (8/3)*(1.339*p.Lw*kap).^2) ./ (1 + (1.339*p.Lw*kap).^2).^(11/6);
    
    Au = sqrt(2 * max(Phi_u,0) * dkap);
    Aw = sqrt(2 * max(Phi_w,0) * dkap);
    
    for i = 1:p.nModes
        phase_i = 2*pi*(kx(i)*X + kz(i)*Z);
        U = U + Au(i) * cos(phase_i + phaseState.phi_u(i));
        W = W + Aw(i) * cos(phase_i + phaseState.phi_w(i));
    end
    
    U = (U - mean(U(:))) * (p.sigma_u / std(U(:)));
    W = (W - mean(W(:))) * (p.sigma_w / std(W(:)));
    Umag = sqrt(U.^2 + W.^2);
end