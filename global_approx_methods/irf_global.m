function [out1, out2, out3, out4] = irf_global(shock_scal, shock_sigma,T, Sims, first, params, theta_val)
%Impulse Response Function (Global Approximation)
%   Computes IRF for global approximation
% Create Exogenous variable matrices with S1 having first row set to shock
rng(12345)
S = shock_sigma*randn(T, Sims);
S1 = S;
S1(1,:) = shock_scal*shock_sigma;
S2 = S;

% Create endogenous variable matrices
Z_lin_1 = [repmat(first(1),1,Sims); zeros(T-1,Sims)];
K_lin_1 = [repmat(first(2),1,Sims); zeros(T-1,Sims)];
C_lin_1 = [repmat(first(3),1,Sims); zeros(T-1,Sims)];
H_lin_1 = [repmat(first(4),1,Sims); zeros(T-1,Sims)];

Z_lin_2 = [repmat(first(1),1,Sims); zeros(T-1,Sims)];
K_lin_2 = [repmat(first(2),1,Sims); zeros(T-1,Sims)];
C_lin_2 = [repmat(first(3),1,Sims); zeros(T-1,Sims)];
H_lin_2 = [repmat(first(4),1,Sims); zeros(T-1,Sims)];

% Simulate technology shocks
for i=1:Sims
    for j=2:T
        Z_lin_1(j,i) = params.rho*Z_lin_1(j-1,i) + S1(j-1,i);
        Z_lin_2(j,i) = params.rho*Z_lin_2(j-1,i) + S2(j-1,i);
    end
end

% Loops to simulate endogenous choice variables for S1
for i=1:Sims
    for j=2:T
        k_use = K_lin_1(j-1,i);
        z_use = Z_lin_1(j-1,i);
        % Project k and z onto [-1,1]
        basis_k = 2 * (k_use - params.k_min) / (params.k_max - params.k_min) -1;
        basis_z = 2 * (z_use - params.Z_min) / (params.Z_max - params.Z_min) -1;
        % Create basis functions for k and z and kronecker product
        psi_use_k = cos([0:(params.n-1)]'*acos(basis_k));
        psi_use_z = cos([0:(params.m-1)]'*acos(basis_z));
        psi_use = kron(psi_use_k,psi_use_z);
        
        % Compute endogenous choice variables
        H_lin_1(j,i) = psi_use' * theta_val;
        C_lin_1(j,i) = ((1-params.alpha)/params.alpha)*(1-params.theta)* exp(z_use)*k_use^params.theta*H_lin_1(j,i)^(-params.theta)*(1-H_lin_1(j,i));
        K_lin_1(j,i)= (exp(z_use) * (k_use .^ params.theta) * (H_lin_1(j,i).^(1-params.theta)) + ((1-params.delta).*k_use) - C_lin_1(j,i))/((1+params.gamma)*(1+params.eta));
        
        % Restrict k_prime to be within [k_min, k_max]
        if K_lin_1(j,i) < params.k_min
            K_lin_1(j,i) = 1.001*params.k_min;
        elseif K_lin_1(j,i) > params.k_max
            K_lin_1(j,i) = 0.999*params.k_max;
        end
    end
end

% Loops to simulate endogenous choice variables for S2
for i=1:Sims
    for j=2:T
        k_use = K_lin_2(j-1,i);
        z_use = Z_lin_2(j-1,i);
        % Project k and z onto [-1,1]
        basis_k = 2 * (k_use - params.k_min) / (params.k_max - params.k_min) -1;
        basis_z = 2 * (z_use - params.Z_min) / (params.Z_max - params.Z_min) -1;
        % Create basis functions for k and z and kronecker product
        psi_use_k = cos([0:(params.n-1)]'*acos(basis_k));
        psi_use_z = cos([0:(params.m-1)]'*acos(basis_z));
        psi_use = kron(psi_use_k,psi_use_z);
        
        % Compute endogenous choice variables
        H_lin_2(j,i) = psi_use' * theta_val;
        C_lin_2(j,i) = ((1-params.alpha)/params.alpha)*(1-params.theta)* exp(z_use)*k_use^params.theta*H_lin_2(j,i)^(-params.theta)*(1-H_lin_2(j,i));
        K_lin_2(j,i)= (exp(z_use) * (k_use .^ params.theta) * (H_lin_2(j,i).^(1-params.theta)) + ((1-params.delta).*k_use) - C_lin_2(j,i))/((1+params.gamma)*(1+params.eta));
        
        % Restrict k_prime to be within [k_min, k_max]
        if K_lin_2(j,i) < params.k_min
            K_lin_2(j,i) = 1.001*params.k_min;
        elseif K_lin_2(j,i) > params.k_max
            K_lin_2(j,i) = 0.999*params.k_max;
        end
        
    end
end

% Compute IRFs
t1 = mean(Z_lin_1') - mean(Z_lin_2');
t2 = mean(K_lin_1') - mean(K_lin_2');
t3 = mean(C_lin_1') - mean(C_lin_2');
t4 = mean(H_lin_1') - mean(H_lin_2');

out1 = t1(2:length(t1));
out2 = t2(2:length(t2));
out3 = t3(2:length(t3));
out4 = t4(2:length(t4));
end

