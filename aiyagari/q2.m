%% ECON630 PS#6 Q2
% Transition Path w/Heterogeneous Agents a la Aiyagari (1994)
% Tom Boesche & Daniel Schwindt
% 12/21/2022
% Script for Q2 in PS#6 (ECON630 - Fall 2022)

%% 0. Housekeeping
clear all
close all
clc

%% 1. Set Parameters
% Consumer
params.bbeta    = 0.96;
params.ssigma   = 2;
params.n_eff    = 5;              % Number of efficiency units possible
params.rrho     = 0.95;
params.ssigma_u = 0.05;

% Production
params.aalpha   = 0.36;
params.ddelta   = 0.06;
Z               = [1,0.9];             % Set technology parameters

% Assets
params.n_asset      = 500;
params.n_asset_c    = 100;
params.a_upper_scale = 12;
params.xxi          = 0;               % borrowing limits

% Efficiency Units grid
params.tauchen_sd   = 3;

% Iteration
params.tol          = 1e-7;
params.tol_paths    = 1.1e-3;                 
params.max_iter_A    = 100;             % max iterations for SCRE
params.max_iter_mu     = 10000;           % max iterations for stationary dist
params.max_iter        = 10000;         % max iterations for transition paths
params.ppsi         = 0.005;           % Parameter to control weighted average interest rate iterations
params.gamma        = 0.01;
%% 2. Compute SRCE for Z=1 and Z=0.9
params.Z = Z(1);
[r1, w1, K1, A1, N1, mu1, v_asset1, v_eff1, m_pol_fun1, m_V_fun1] = srce(params);
params.Z = Z(2);
[r2, w2, K2, A2, N2, mu2, v_asset2, v_eff2, m_pol_fun2, m_V_fun2] = srce(params);

%% 3. Transition path from SS1 to SS2
[v_log_eff, m_trans, pi] = tauchen(params.rrho, params.ssigma_u, ...
                                    params.n_eff, params.tauchen_sd);
v_eff = exp(v_log_eff);
N = v_eff*pi;

T = 50;
% Create vectors of guesses for r, K, and w
r_guess = repmat(r2, 1,T);
K_guess = ((params.aalpha*Z(2))./(r_guess + params.ddelta)).^(1/(1-params.aalpha)).*N;
w_guess = (1-params.aalpha)*Z(2)*(K_guess./N).^params.aalpha;

% Create empty matrices to track value and policy functions over time
m_V_fun = nan(params.n_asset, params.n_eff, 50);
m_pol_fun = nan(params.n_asset, params.n_eff, 50);
m_pol_idx = nan(params.n_asset, params.n_eff, 50);

m_V_fun(:,:,50) = m_V_fun2;
m_pol_fun(:,:,50) = m_pol_fun2;
[~,m_pol_idx(:,:,50)] = ismember(m_pol_fun2, v_asset2);

%% While loop over guesses until convergence
diff = 100;
iter = 0;
while (diff > params.tol_paths) && (iter < params.max_iter) 
%% Backward induction step
for t=(T-1):-1:1
    m_V_fun_exp = m_V_fun(:,:,t+1)*m_trans';
    income = ((1+r_guess(t))*v_asset2 + w_guess(t)*v_eff);
    % Find optimal asset choice
    for i=1:params.n_asset
        for j=1:params.n_eff
            cons = max(income(i,j) - v_asset2,0); % don't allow consumption to go negative
            utils = cons.^(1-params.ssigma)./(1-params.ssigma);
            value = (1-params.bbeta)*utils + params.bbeta*m_V_fun_exp(:,j);
            [m_V_fun(i,j,t), m_pol_idx(i,j,t)] = max(value);
            m_pol_fun(i,j,t) = v_asset2(m_pol_idx(i,j,t));
        end
    end
end

%% Forward iteration of asset distribution
mu = zeros(params.n_asset, params.n_eff, T);
mu(:,:,1) = mu1;
for t=2:T
    for j=1:params.n_eff
        for i=1:params.n_asset
            for j_nxt = 1:params.n_eff
                mu(m_pol_idx(i,j,t),j_nxt,t) = mu(m_pol_idx(i,j,t),j_nxt,t) + ...
                    m_trans(j,j_nxt)*mu(i,j,t-1);
            end
        end
    end
end

%% Compute Asset supply for each t
A = nan(1,T);
for t=1:T
    A(t) = sum(sum(mu(:,:,t).*m_pol_fun(:,:,t)));
end

%% Compute implied market clearing interest rate
r_implied = params.aalpha*Z(2)*(N./A).^(1-params.aalpha) - params.ddelta;
diff = max(log(r_guess) - log(r_implied));

r_newguess = params.gamma*r_implied + (1-params.gamma)*r_guess;

iter = iter + 1;

if (mod(iter,10)==0 || iter ==1)
    fprintf('-------------------------------------------- \n')
    fprintf('Iterations: %d, \n diff_A: %e \n',iter, diff)
    fprintf('-------------------------------------------- \n')
end

%% Re-initialize guesses for next iteration
r_guess = r_newguess;
K_guess = ((params.aalpha*Z(2))./(r_guess + params.ddelta)).^(1/(1-params.aalpha)).*N;
w_guess = (1-params.aalpha)*Z(2)*(K_guess./N).^params.aalpha;

end

%% Output Charts

r_path = [r1, r_guess, r2];
K_path = [K1, K_guess, K2];
w_path = [w1, w_guess, w2];

figure(1)
subplot(3,1,1)
plot(r_path, 'k', 'LineWidth', 2)
title('Transition Path, r')
ylabel('Interest Rate')
xlabel('Time')

subplot(3,1,2)
plot(K_path, 'k', 'LineWidth', 2)
title('Transition Path, K')
ylabel('Capital Stock')
xlabel('Time')

subplot(3,1,3)
plot(w_path, 'k', 'LineWidth', 2)
title('Transition Path, w')
ylabel('Wage')
xlabel('Time')

saveas(gcf, '../output/transition_path', 'epsc');