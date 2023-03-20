%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECON630 PS3 Part B
% Global Approximation Method
% Stochastic Neoclassical Growth Model:
% - Consumption & Labor-Leisure
% - Population & Technology growth
% Author: Daniel Schwindt
% Date: 10/22/2022

%% 0. Housekeeping
clear all;
close all;
clc;
tic;

%% 1. Solve Linear Approximation w/Dynare
dynare cooley_prescott_linear

%% 2. Set Parameters
parameters; % Call parameters script

%% 3. Set initial theta guess using Dynare linear solution
%theta_0 = zeros(n*m,2);
theta_0 = zeros(n*m,1);
%c_0 = css + coef_ck*(K-kss) + coef_cz*(Z-zss); % initial c from linear decision rule
h_0 = hss + coef_hk*(K-kss) + coef_hz*(Z-zss); % initial h from linear decision rule
%dr = [c_0, h_0];
[theta_0] = fsolve(@(theta_0) psi'*theta_0-h_0, theta_0, options_use);

%% 4. Solve Residual System
[theta_use,fval,exitflag,output1] = fsolve (@(theta_0) resid_system(theta_0, par, psi, k, z, P),theta_0, options_use);
toc;

%% 5. Create Data for Decision Rule Plots
n_k = 100;
n_z = 25;

% Create productivity grid and find transition probabilities
[grid_z,transition_grid,dist] = tauchen(rho,sigma_eps,n_z);

grid_z = grid_z';

% Define capital grid
grid_k = [k_min:(k_max-k_min)/n_k:k_max]';

size_k = size(grid_k,1);

size_z = size(grid_z,1);

% Loop over capital and productivity to construct Chebyshev polynomials and
% to approximate capital, labor, and consumption solutions from decision
% rule
for i = 1:size_k  % k
    for j = 1:size_z  % a   
        
        k_use = grid_k(i);
        z_use = grid_z(j);
        
        basis_k = 2 * (k_use - k_min) / (k_max - k_min) -1; % These are grid_k projected into [-1,1]

        basis_z = 2 * (z_use - Z_min) / (Z_max - Z_min) -1; % These are grid_z projected into [-1,1]
        
        psi_use_k = cos([0:n-1]'*acos(basis_k));

        psi_use_z = cos([0:m-1]'*acos(basis_z));

        psi_use = kron(psi_use_k,psi_use_z);
        
        H(i,j) = psi_use' * theta_use;
        C(i,j) = ((1-alpha)/alpha)*(1-theta)* exp(z_use)*k_use^theta*H(i,j)^(-theta)*(1-H(i,j));
        
        K_P(i,j)= (exp(z_use) * (k_use .^ theta) * (H(i,j).^(1-theta)) + ((1-delta).*k_use) - C(i,j))/((1+gamma)*(1+eta));
        
        % Compute C_PRIME and the Euler Equation Error
        
        x_prime_k = 2 * (K_P(i,j) - k_min) / (k_max - k_min) -1; % This is k_prime projected into [-1,1]
        psi_k_prime = cos([0:n-1]'*acos(x_prime_k));  % These are n*1      
        
        for count = 1:size_z
    
            z_prime = grid_z(count);    
           
            % Construct psi_z_prime
            
            x_prime_z = 2 * (z_prime - Z_min) / (Z_max - Z_min) -1; % This is z_prime projected into [-1,1]
    
            psi_z_prime = cos([0:m-1]'*acos(x_prime_z)); % These are m*1

            psi_prime = kron(psi_k_prime,psi_z_prime);
            
            h_prime = psi_prime' * theta_use;  
           
            % Error Check        
            if h_prime < 0
                h = 10e-6;
                disp('Warning : At EE1, h_prime < 0 encountered!!!')
            elseif h_prime > 1
                h = 1-10e-6;
                disp('Warning: At EE1, h_prime > 1 encountered!!!')
            end   
            
            % Create consumption and leisure choice vars
            c_prime = ((1-alpha)/(alpha))*(1-theta)* ...
                exp(z_prime)*K_P(i,j)^theta*h_prime^(-theta)*(1-h_prime);
            l_prime = 1-h_prime;
            
            % Define u' and f' for EE
            u_prime = (1-alpha)*(c_prime^(1-alpha) * ...
                l_prime ^ alpha)^(-sigma)* ...
                c_prime^-alpha * l_prime^alpha;            
            f_prime = theta * exp(z_prime) * (K_P(i,j))^(theta-1)* ...
                (h_prime)^(1-theta);
            
            integrand(:,count) = u_prime * (f_prime + 1 - delta);               
                
        end

        sol_integ = integrand * transition_grid(j,:)';    

    end
end

% Create Values for Linear Solution
K_P_lin = zeros(size_k, size_z);
C_lin = zeros(size_k, size_z);
H_lin = zeros(size_k, size_z);
for i=1:size_k
    for j=1:size_z
        K_P_lin(i,j) = kss + coef_kk*(grid_k(i)-kss) + coef_kz*(grid_z(j)-zss);
        C_lin(i,j) = css + coef_ck*(grid_k(i)-kss) + coef_cz*(grid_z(j)-zss);
        H_lin(i,j) = hss + coef_hk*(grid_k(i)-kss) + coef_hz*(grid_z(j)-zss);
    end
end

%% 6. Create Data for IRF Plots
T = 60;
N = 1000;
%% (a) Start at Steady State
% Create vector of steady state values
ss = [zss, kss, css, hss];
% Create vector of linear decision rule coefficients
coeffs = [coef_kk, coef_kz, coef_ke, coef_ck, coef_cz, coef_ce, coef_hk, coef_hz, coef_he];
% Compute Linear IRFs using irf_lin function
[irf_z_a1, irf_k_a1, irf_c_a1, irf_h_a1] = irf_lin(1, sigma_eps, T, N, ss, coeffs, rho);
[irf_z_a2, irf_k_a2, irf_c_a2, irf_h_a2] = irf_lin(2, sigma_eps, T, N, ss, coeffs, rho);
[irf_z_a3, irf_k_a3, irf_c_a3, irf_h_a3] = irf_lin(3, sigma_eps, T, N, ss, coeffs, rho);

% Compute Global IRFs
[irf_zg_a1, irf_kg_a1, irf_cg_a1, irf_hg_a1] = irf_global(1, sigma_eps, T+1, N, ss, par, theta_use);
[irf_zg_a2, irf_kg_a2, irf_cg_a2, irf_hg_a2] = irf_global(2, sigma_eps, T+1, N, ss, par, theta_use);
[irf_zg_a3, irf_kg_a3, irf_cg_a3, irf_hg_a3] = irf_global(3, sigma_eps, T+1, N, ss, par, theta_use);
%% (b) Start at 1.1xKss and 2xsigma_eps
start = [2*sigma_eps, 1.1*kss, css, hss];
% Compute Linear IRFs using irf_lin function
[irf_z_b1, irf_k_b1, irf_c_b1, irf_h_b1] = irf_lin(1, sigma_eps, T, N, start, coeffs, rho);
[irf_z_b2, irf_k_b2, irf_c_b2, irf_h_b2] = irf_lin(2, sigma_eps, T, N, start, coeffs, rho);
[irf_z_b3, irf_k_b3, irf_c_b3, irf_h_b3] = irf_lin(3, sigma_eps, T, N, start, coeffs, rho);

% Compute Global IRFs
[irf_zg_b1, irf_kg_b1, irf_cg_b1, irf_hg_b1] = irf_global(1, sigma_eps, T+1, N, start, par, theta_use);
[irf_zg_b2, irf_kg_b2, irf_cg_b2, irf_hg_b2] = irf_global(2, sigma_eps, T+1, N, start, par, theta_use);
[irf_zg_b3, irf_kg_b3, irf_cg_b3, irf_hg_b3] = irf_global(3, sigma_eps, T+1, N, start, par, theta_use);

%% 7. Plot Decision Rule Figures
figure(1);
surf(grid_z,grid_k,K_P);
shading interp;
axis tight;
ylabel('Capital Today')
xlabel('Productivity Shock')
zlabel('Capital Tomorrow')
title('Decision Rule for Capital');
saveas(gcf, '../output/K_P_3d','epsc');

figure(2);
surf(grid_z,grid_k,C);
shading interp;
axis tight;
ylabel('Capital Today')
xlabel('Productivity Shock')
zlabel('Consumption')
title('Decision Rule for Consumption');
saveas(gcf, '../output/C_3d','epsc');

figure(3);
surf(grid_z,grid_k,H);
shading interp;
axis tight;
ylabel('Capital Today')
xlabel('Productivity Shock')
zlabel('Labor')
title('Decision Rule for Labor');
saveas(gcf, '../output/H_3d','epsc');

figure(4)
plot(grid_k, K_P(:,13), 'b');
hold on
plot(grid_k, K_P_lin(:,13), 'r');
axis tight;
xlabel('Capital Today')
ylabel('Capital Tomorrow')
title('Decision Rule for Capital, z=0');
legend('Global', 'Linear');
hold off
saveas(gcf, '../output/K_P_2d','epsc');

figure(5)
plot(grid_k, C(:,13), 'b');
hold on
plot(grid_k, C_lin(:,13), 'r');
axis tight;
xlabel('Capital Today')
ylabel('Consumption')
title('Decision Rule for Consumption, z=0');
legend('Global', 'Linear');
hold off
saveas(gcf, '../output/C_2d','epsc');

figure(6)
plot(grid_k, H(:,13), 'b');
hold on
plot(grid_k, H_lin(:,13), 'r');
axis tight;
xlabel('Capital Today')
ylabel('Labor')
title('Decision Rule for Labor, z=0');
legend('Global', 'Linear');
hold off
saveas(gcf, '../output/H_2d','epsc');

%% 8. Plot IRFs
x = [0:59];

% IRFs at the Steady State
% Capital IRF
figure (7)
hold on
subplot(3,1,1);
plot(x,irf_kg_a1, 'b', x,irf_k_a1, 'r');
xlabel('Time')
ylabel('Capital')
title('Shock = 1SD');
legend('Global', 'Linear', 'Location', 'southeast');
subplot(3,1,2);
plot(x, irf_kg_a2, 'b', x, irf_k_a2, 'r');
xlabel('Time')
ylabel('Capital')
title('Shock = 2SD');
subplot(3,1,3);
plot(x, irf_kg_a3, 'b', x, irf_k_a3, 'r');
xlabel('Time')
ylabel('Capital')
title('Shock = 3SD');
hold off
saveas(gcf, '../output/K_irf_ss','epsc');

% Consumption IRF
figure (8)
hold on
subplot(3,1,1);
plot(x,irf_cg_a1, 'b', x,irf_c_a1, 'r');
xlabel('Time')
ylabel('Consumption')
title('Shock = 1SD');
legend('Global', 'Linear', 'Location', 'southeast');
subplot(3,1,2);
plot(x, irf_cg_a2, 'b', x, irf_c_a2, 'r');
xlabel('Time')
ylabel('Consumption')
title('Shock = 2SD');
subplot(3,1,3);
plot(x, irf_cg_a3, 'b', x, irf_c_a3, 'r');
xlabel('Time')
ylabel('Consumption')
title('Shock = 3SD');
hold off
saveas(gcf, '../output/C_irf_ss','epsc');

% Labor IRF
figure (9)
hold on
subplot(3,1,1);
plot(x,irf_hg_a1, 'b', x,irf_h_a1, 'r');
xlabel('Time')
ylabel('Labor')
title('Shock = 1SD');
legend('Global', 'Linear', 'Location', 'southeast');
subplot(3,1,2);
plot(x, irf_hg_a2, 'b', x, irf_h_a2, 'r');
xlabel('Time')
ylabel('Labor')
title('Shock = 2SD');
subplot(3,1,3);
plot(x, irf_hg_a3, 'b', x, irf_h_a3, 'r');
xlabel('Time')
ylabel('Labor')
title('Shock = 3SD');
hold off
saveas(gcf, '../output/H_irf_ss','epsc');

% IRFs at 1.1KSS and 2*SD for Z
% Capital IRF
figure (10)
hold on
subplot(3,1,1);
plot(x,irf_kg_b1, 'b', x,irf_k_b1, 'r');
xlabel('Time')
ylabel('Capital')
title('Shock = 1SD');
legend('Global', 'Linear', 'Location', 'southeast');
subplot(3,1,2);
plot(x, irf_kg_b2, 'b', x, irf_k_b2, 'r');
xlabel('Time')
ylabel('Capital')
title('Shock = 2SD');
subplot(3,1,3);
plot(x, irf_kg_b3, 'b', x, irf_k_b3, 'r');
xlabel('Time')
ylabel('Capital')
title('Shock = 3SD');
hold off
saveas(gcf, '../output/K_irf_nonss','epsc');

% Consumption IRF
figure (11)
hold on
subplot(3,1,1);
plot(x,irf_cg_b1, 'b', x,irf_c_b1, 'r');
xlabel('Time')
ylabel('Consumption')
title('Shock = 1SD');
legend('Global', 'Linear', 'Location', 'southeast');
subplot(3,1,2);
plot(x, irf_cg_b2, 'b', x, irf_c_b2, 'r');
xlabel('Time')
ylabel('Consumption')
title('Shock = 2SD');
subplot(3,1,3);
plot(x, irf_cg_b3, 'b', x, irf_c_b3, 'r');
xlabel('Time')
ylabel('Consumption')
title('Shock = 3SD');
hold off
saveas(gcf, '../output/C_irf_nonss','epsc');

% Labor IRF
figure (12)
hold on
subplot(3,1,1);
plot(x,irf_hg_b1, 'b', x,irf_h_b1, 'r');
xlabel('Time')
ylabel('Labor')
title('Shock = 1SD');
legend('Global', 'Linear', 'Location', 'southeast');
subplot(3,1,2);
plot(x, irf_hg_b2, 'b', x, irf_h_b2, 'r');
xlabel('Time')
ylabel('Labor')
title('Shock = 2SD');
subplot(3,1,3);
plot(x, irf_hg_b3, 'b', x, irf_h_b3, 'r');
xlabel('Time')
ylabel('Labor')
title('Shock = 3SD');
hold off
saveas(gcf, '../output/H_irf_nonss','epsc');