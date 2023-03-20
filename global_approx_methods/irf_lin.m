function [out1, out2, out3, out4] = irf_lin(shock_scal, shock_sigma,T, Sims,first, coeff, param)
%Impulse Response Function (Linear solution)
%   Computes IRF for linear solution
% Create SIGMA matrix for shocks - bifurcated into SIGMA1 which has shock
% in first row
S = shock_sigma*randn(T, Sims);
S1 = S;
S1(1,:) = shock_scal*shock_sigma;
S2 = S;

% Simulate endogenous variables using S1
Z_lin_1 = [repmat(first(1),1,Sims); zeros(T-1,Sims)];
K_lin_1 = [repmat(first(2),1,Sims); zeros(T-1,Sims)];
C_lin_1 = [repmat(first(3),1,Sims); zeros(T-1,Sims)];
H_lin_1 = [repmat(first(4),1,Sims); zeros(T-1,Sims)];

% Simulate endogenous variables using S2
Z_lin_2 = [repmat(first(1),1,Sims); zeros(T-1,Sims)];
K_lin_2 = [repmat(first(2),1,Sims); zeros(T-1,Sims)];
C_lin_2 = [repmat(first(3),1,Sims); zeros(T-1,Sims)];
H_lin_2 = [repmat(first(4),1,Sims); zeros(T-1,Sims)];

% Loop over Endogenous Matrices 1
for i=1:Sims
    for j=2:T
        Z_lin_1(j,i) = first(1) + param*Z_lin_1(j-1,i) + S1(j-1,i);
        K_lin_1(j,i) = first(2) + coeff(1)*(K_lin_1(j-1,i)-first(2)) + ...
            coeff(2)*(Z_lin_1(j-1,i)-first(1)) + coeff(3)*S1(j-1,i);
        C_lin_1(j,i) = first(3) + coeff(4)*(K_lin_1(j-1,i)-first(2)) + ...
            coeff(5)*(Z_lin_1(j-1,i)-first(1)) + coeff(6)*S1(j-1,i);
        H_lin_1(j,i) = first(4) + coeff(7)*(K_lin_1(j-1,i)-first(2)) + ...
            coeff(8)*(Z_lin_1(j-1,i)-first(1)) + coeff(9)*S1(j-1,i);
    end
end

% Loop over Endogenous Matrices 2
for i=1:Sims
    for j=2:T
        Z_lin_2(j,i) = first(1) + param*Z_lin_2(j-1,i) + S2(j-1,i);
        K_lin_2(j,i) = first(2) + coeff(1)*(K_lin_2(j-1,i)-first(2)) + ...
            coeff(2)*(Z_lin_2(j-1,i)-first(1)) + coeff(3)*S2(j-1,i);
        C_lin_2(j,i) = first(3) + coeff(4)*(K_lin_2(j-1,i)-first(2)) + ...
            coeff(5)*(Z_lin_2(j-1,i)-first(1)) + coeff(6)*S2(j-1,i);
        H_lin_2(j,i) = first(4) + coeff(7)*(K_lin_2(j-1,i)-first(2)) + ...
            coeff(8)*(Z_lin_2(j-1,i)-first(1)) + coeff(9)*S2(j-1,i);
    end
end

out1 = mean(Z_lin_1') - mean(Z_lin_2');
out2 = mean(K_lin_1') - mean(K_lin_2');
out3 = mean(C_lin_1') - mean(C_lin_2');
out4 = mean(H_lin_1') - mean(H_lin_2');

end

