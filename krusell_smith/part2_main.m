%% ECON630 Final Project, Part II
% Heterogeneous-Agent Model with Aggregate Uncertainty, Krusell-Smith 1998
% Daniel Schwindt
% 12/25/2022
% Updated: 1/19/2023
% Script to run all programs for Part II

%% 0. Housekeeping
clear all;
close all;
clc;

%% 1. Set Parameters
rng(12345);
z       = [0.99 1.01];   % Aggregate productivity states
v_eff   = [0 1];                 % Employed & Unemployed

params = struct;
% Grids
params.n_k              = 100;
params.n_kbar           = 4;
params.n_agg            = length(z);
params.n_eff            = length(v_eff);
params.k_min            = 0;
params.k_max            = 125;
params.kbar_min         = 10;
params.kbar_max         = 50;
params.theta            = 7;
% Utility & Production
params.bbeta            = 0.99;
params.ssigma           = 1;
params.aalpha           = 0.36;
params.ddelta           = 0.025;
% Sampling
params.T                = 11000;
params.I                = 5000;
params.burn             = 1000;
% Transitions
params.u_g              = 0.04;
params.u_b              = 0.1;
params.dur_g            = 1.5;
params.dur_b            = 2.5;
params.d_g              = 8;
params.d_b              = 8;
lbar                    = [1-params.u_b 1-params.u_g]; % aggregate labor
% Convergence
params.tol              = 1e-8;
params.tol_b            = 1e-6;
params.max_iter         = 10000;
params.psi_k            = 0.7;
params.psi_b            = 0.3;

%% 2. Set Up
[v_k, v_kbar, m_k, m_kbar] = k_grid_create(params, 1); % Create grids
[P_z, P_ze, P_1, P_2, P_11, P_12, P_21, P_22] = trans_prob(params); % Assign transition probabilities
[shocks_z,shocks_z_idx] = agg_shocks_gen(params.T,z,P_z);
[shocks_e,~] = ind_shocks_gen(params.T,params.I,params.u_b,P_z, P_ze, shocks_z_idx);
[shocks_e, shocks_e_idx] = adjust_ind_shocks(shocks_e, shocks_z_idx, params.u_g, params.u_b);
eps_z = shocks_z - [nan; shocks_z(1:end-1)];
eps_l = mean(shocks_e,2) - [nan; mean(shocks_e(1:end-1,:),2)];
b0 = [0.09 0.09];       % Intercept: bad, good states respectively
b1 = [1 1];             % Slope: bad, good states respectively

%% 3. Solve Incomplete Markets Benchmark Model
[b0_bm, b1_bm, K_bm, A_bm, C_bm, m_k_nxt_bm, m_c_star_bm, iter_bm,...
    max_diff_bm, iter_k_bm, max_diff_k_bm] = solve_ks_im(params, z, v_eff, lbar, ...
    b0, b1, P_11, P_12, P_21, P_22,v_k, v_kbar, m_k, m_kbar,...
    shocks_z, shocks_z_idx, shocks_e, shocks_e_idx);

%% 4. Solve Incomplete Markets SIGMA=5
params.ssigma = 5;
params.psi_b = 0.05;
params.tol_b = 3e-3; % Have to increase tolerance threshold b/c of overshooting
[b0_sm, b1_sm, K_sm, A_sm, C_sm, m_k_nxt_sm, m_c_star_sm, iter_sm,...
    max_diff_sm, iter_k_sm, max_diff_k_sm] = solve_ks_im(params, z, v_eff, lbar, ...
    b0, b1, P_11, P_12, P_21, P_22,v_k, v_kbar, m_k, m_kbar,...
    shocks_z, shocks_z_idx, shocks_e, shocks_e_idx);

%% 5. Solve Complete Markets Benchmark Model
dynare ks_complete_s1
%simul_cm = oo_.endo_simul';
cm_endo_rules = oo_.dr.ghx;
cm_exo_rules = oo_.dr.ghu;
cm_cons_rules = oo_.dr.ys;
K_cm = nan(params.T, 1);
K_cm(1) = params.kbar_min;
for i=2:params.T
    K_cm(i) = (1-cm_endo_rules(5))*cm_cons_rules(3) + cm_endo_rules(5)*K_cm(i-1) + cm_exo_rules(5,1)*eps_z(i) + cm_exo_rules(5,2)*eps_l(i);
end
%% 6. Solve Complete Markets SIGMA=5
dynare ks_complete_s5
%simul_cs = oo_.endo_simul';
cs_endo_rules = oo_.dr.ghx;
cs_exo_rules = oo_.dr.ghu;
cs_cons_rules = oo_.dr.ys;
K_cs = nan(params.T, 1);
K_cs(1) = params.kbar_min;
for i=2:params.T
    K_cs(i) = (1-cs_endo_rules(5))*cs_cons_rules(3) + cs_endo_rules(5)*K_cs(i-1) + cs_exo_rules(5,1)*eps_z(i) + cs_exo_rules(5,2)*eps_l(i);
end
%% 7. Individual Decision Rules Charts (Benchmark Model)
% Decision Rule charts
figure(1)
subplot(2,2,1);
plot(v_k, m_k_nxt_bm(:,:,1,1));
title('Individual Decision Rule: Z=0.99, E=0');
xlabel('Assets Today');
ylabel('Assets Tomorrow');

subplot(2,2,2);
plot(v_k, m_k_nxt_bm(:,:,1,2));
title('Individual Decision Rule: Z=0.99, E=1');
xlabel('Assets Today');
ylabel('Assets Tomorrow');

subplot(2,2,3);
plot(v_k, m_k_nxt_bm(:,:,2,1));
title('Individual Decision Rule: Z=1.01, E=0');
xlabel('Assets Today');
ylabel('Assets Tomorrow');

subplot(2,2,4);
plot(v_k, m_k_nxt_bm(:,:,2,2));
title('Individual Decision Rule: Z=1.01, E=1');
xlabel('Assets Today');
ylabel('Assets Tomorrow');
saveas(gcf, '../../output/ind_dec_rules', 'epsc');
saveas(gcf, '../../output/ind_dec_rules', 'png');

%% 8. ALM Charts (Benchmark Model)
v_kbar_b = exp(b0_bm(1) + b1_bm(1)*log(v_kbar));
v_kbar_g = exp(b0_bm(2) + b1_bm(2)*log(v_kbar));
figure(2)
plot(v_kbar, v_kbar_b, 'r');
hold on
plot(v_kbar, v_kbar_g, 'b');
xlim([10,50]);
xlabel('Capital Today');
ylabel('Capital Tomorrow');
legend('Z=0.99', 'Z=1.01');
hold off
saveas(gcf, '../../output/agg_cap_rules', 'epsc');
saveas(gcf, '../../output/agg_cap_rules', 'png');

%% 9. Table 1 Replication: Distribution of Wealth
mean_A_bm = mean(A_bm(params.burn:end, :),1);
mean_A_bm = sort(mean_A_bm, 'descend');
pct_w = [0.01 0.05 0.1 0.2 0.3];
tops = pct_w*params.I;
tot_A_bm = sum(mean_A_bm);
top_w = nan(1,length(tops));
for i=1:length(tops)
    top_w(i) = 100*sum(mean_A_bm(1:tops(i)))/tot_A_bm;
end
writematrix(top_w, '../../output/part2_table1.xlsx');

%% 10. Figure 3 Replication: Lorenz Curve (w/ Gini Coeff)
mean_A_bm = sort(mean_A_bm, 'ascend');
frac_pop = linspace(0,1,11);
frac_pop_w = frac_pop*params.I;
frac_w = nan(1, length(frac_pop_w));
for i=1:length(frac_w)
    frac_w(i) = sum(mean_A_bm(1:frac_pop_w(i)))/tot_A_bm;
end

% Plot
figure(3)
plot(frac_pop, frac_pop, 'k')
hold on
plot(frac_pop, frac_w, 'bo-')
xlabel('Fraction of population')
ylabel('Fraction of wealth')
ylim([0,1])
xlim([0,1])
saveas(gcf, '../../output/part2_fig3', 'epsc');

% Compute Gini Coeff
S = cumsum(mean_A_bm)/tot_A_bm;
x = linspace(0,1,params.I);
dx = x - [0, x(1:end-1)];
gini = 1 - 2*sum(S.*dx);

%% 11. Table 2 Replication: Aggregate Time Series
%K_all = [K_bm simul_cm(:,3) K_sm simul_cs(:,3)];
K_all = [K_bm K_cm K_sm K_cs];
labor = lbar(shocks_z_idx)';
Y_all = shocks_z.*K_all.^params.aalpha.*labor.^(1-params.aalpha);
%Y_all(:,[2 4]) = [simul_cm(:,1) simul_cs(:,1)]; 
Y1_all = [nan nan nan nan; Y_all(1:end-1,:)];
Y4_all = [repmat(nan,4); Y_all(1:end-4,:)];
r_all = params.aalpha*shocks_z.*K_all.^(params.aalpha-1).*labor.^(1-params.aalpha);
%r_all(:,[2 4]) = [simul_cm(:,4) simul_cs(:,4)];
w_all = (1-params.aalpha)*shocks_z.*K_all.^params.aalpha.*labor.^-params.aalpha;
%w_all(:,[2 4]) = [simul_cm(:,5) simul_cs(:,5)];
C_all = Y_all + (1-params.ddelta)*[K_all(1,:); K_all(1:end-1,:)] - K_all;
%C_all(:, [2 4]) = [simul_cm(:,2) simul_cs(:,2)];

% Comput aggregate statistics
mean_K_all = mean(K_all(params.burn+1:end,:));
corr_C_Y = nan(1,4);
corr_Y_Y1 = nan(1,4);
corr_Y_Y4 = nan(1,4);
for i=1:4
    corr_C_Y(i) = corr(C_all(params.burn+1:end,i), Y_all(params.burn+1:end,i));
    corr_Y_Y1(i) = corr(Y_all(params.burn+1:end,i), Y1_all(params.burn+1:end,i));
    corr_Y_Y4(i) = corr(Y_all(params.burn+1:end,i), Y4_all(params.burn+1:end,i));
end
sd_r = 100*std(r_all(params.burn+1:end,:));

table2_out = [mean_K_all' corr_C_Y', sd_r', corr_Y_Y4'];
writematrix(table2_out, '../../output/part2_table2.xlsx');

%% 13. Extra Appendix Charts
% Decision Rule charts, SIGMA = 5
figure(1)
subplot(2,2,1);
plot(v_k, m_k_nxt_sm(:,:,1,1));
title('Individual Decision Rule: Z=0.99, E=0');
xlabel('Assets Today');
ylabel('Assets Tomorrow');

subplot(2,2,2);
plot(v_k, m_k_nxt_sm(:,:,1,2));
title('Individual Decision Rule: Z=0.99, E=1');
xlabel('Assets Today');
ylabel('Assets Tomorrow');

subplot(2,2,3);
plot(v_k, m_k_nxt_sm(:,:,2,1));
title('Individual Decision Rule: Z=1.01, E=0');
xlabel('Assets Today');
ylabel('Assets Tomorrow');

subplot(2,2,4);
plot(v_k, m_k_nxt_sm(:,:,2,2));
title('Individual Decision Rule: Z=1.01, E=1');
xlabel('Assets Today');
ylabel('Assets Tomorrow');
saveas(gcf, '../../output/ind_dec_rules_sig5', 'epsc');

% ALM Chart, SIGMA = 5
v_kbar_b = exp(b0_sm(1) + b1_sm(1)*log(v_kbar));
v_kbar_g = exp(b0_sm(2) + b1_sm(2)*log(v_kbar));
figure(2)
plot(v_kbar, v_kbar_b, 'r');
hold on
plot(v_kbar, v_kbar_g, 'b');
xlim([10,50]);
xlabel('Capital Today');
ylabel('Capital Tomorrow');
legend('Z=0.99', 'Z=1.01');
hold off
saveas(gcf, '../../output/agg_cap_rules_sig5', 'epsc');