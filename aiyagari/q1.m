%% ECON630 PS#6 Q1
% General Equilibrium model w/Heterogeneous Agents a la Aiyagari (1994)
% Tom Boesche & Daniel Schwindt
% 12/21/2022
% Script for Q1 in PS#6 (ECON630 - Fall 2022)

%% 0. Housekeeping
clear all
close all
clc

%% 1. Set Parameters
% Consumer
params.bbeta = 0.96;
params.ssigma = 2;
params.n_eff = 5;              % Number of efficiency units possible
params.rrho = 0.95;
params.ssigma_u = 0.05;

% Production
params.aalpha = 0.36;
params.ddelta = 0.06;
params.Z = 1;

% Assets
params.n_asset = 500;
params.n_asset_c = 100;
params.a_upper_scale = 12;
xxi = [0, 0.5];                 % borrowing limits

% Efficiency Units grid
params.tauchen_sd = 3;

% Iteration
params.tol = 1e-7;
params.max_iter_A = 100;                  % max iterations for scre
params.max_iter_mu = 10000;               % max iterations for stationary dist.
params.ppsi = 0.005;                   % Parameter to control weighted average interest rate iterations

%% 2. Compute SRCE for different values of xxi
params.xxi = xxi(1);                % Set borrowing limit to 0
[r1, w1, K1, A1, N1, mu1, v_asset1, v_eff1, m_pol_fun1, m_v_fun1] = srce(params);
params.xxi = xxi(2);                % Set borrowing limit to 0.5*natural limit
[r2, w2, K2, A2, N2, mu2, v_asset2, v_eff2, m_pol_fun2, m_v_fun2] = srce(params);

%% 4. Construct Lorenz Curve and Compute Gini Coefficient
mu = nan(params.n_asset, params.n_eff, 2);
mu(:,:,1) = mu1;
mu(:,:,2) = mu2;

v_asset = nan(params.n_asset,2);
v_asset(:,1) = v_asset1;
v_asset(:,2) = v_asset2;

for i=1:2
    H = sum(mu(:,:,i),2);
    F = [0;cumsum(H)];
    S = cumsum(v_asset(:,i).*H);
    S_prior = [0; S(1:end-1)];
    L = [0; S./S(end)];
    gini = 1-sum(H.*(S_prior+S))/S(end); 

    % Plot Lorenz Curve with box for Gini Coefficient
    figure(1)
    plot(F,L, 'b');
    hold on
    plot(F, F,'k');
    title('Lorenz Curve')
    xlabel('Cumul. Share of People, Low to High Assets')
    ylabel('Cumul. Share of Assets')
    legend('Lorenz Curve', 'Equality Line', 'Location', 'Northwest')
    txt = strcat('Gini Coeff:', num2str(round(gini,2)));
    text(0.02,0.8,txt)
    hold off
    saveas(gcf, strcat('../output/lorenz_curve',num2str((i))),'epsc')
end
%% Plot asset distributions by xxi
for i=1:2
    a_dist = sum(mu(:,:,i),2); % marginal dist of assets
    figure(2)
    plot(v_asset(:,i), a_dist, 'k', 'LineWidth', 2)
    title('Marginal Distribution of Assets')
    xlabel('Assets')
    ylabel('Probability')
    saveas(gcf, strcat('../output/a_marg_dist',num2str(i)),'epsc')
    
    a_cumdist = cumsum(a_dist);  % cumulative dist of assets
    figure(3)
    plot(v_asset(:,i), a_cumdist, 'k', 'LineWidth', 2);
    title('Cumul. Distribution of Assets')
    xlabel('Assets')
    ylabel('Probability')
    saveas(gcf, strcat('../output/a_cum_dist',num2str(i)),'epsc')

end

%% Plot decision rules by xxi
cash =nan(params.n_asset, params.n_eff, 2);
cash(:,:,1) = (1+r1)*v_asset(:,1) + w1*v_eff1;
cash(:,:,2) = (1+r2)*v_asset(:,2) + w2*v_eff2;
m_pol_fun = nan(params.n_asset, params.n_eff, 2);
m_pol_fun(:,:,1) = m_pol_fun1;
m_pol_fun(:,:,2) = m_pol_fun2;
titles = ['Borrowing Limit: 0', 'Borrowing Limit: Half of Natural Limit'];
figure(4)
for i=1:2
    subplot(2,1,i)
    plot(cash(:,:,i), m_pol_fun(:,:,i))
    xlabel('Resources')
    ylabel('Asset Choice')
    title('Borrowing Limit: 0')
end

saveas(gcf, '../output/decision_rules', 'epsc');