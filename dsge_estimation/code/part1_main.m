%% ECON630 Final Project, Part I
% An & Schorfheide (2007)
% Run all sub-programs to produce final output for report
% Author: Daniel Schwindt
% Date: 12/23/2022

%% 0. Housekeeping
clear all; clc;
options  = optimset('Display', 'off');
min_options = optimset('Display', 'iter', 'TolFun', 1e-7, 'TolX', 1e-7, 'MaxFunEvals', 10000);
n = 100000;
%% 1. Set parameters for initial model solve
v_params = [2.00;       % tau
            0.15;       % kappa
            1.50;       % psi1
            1.00;       % psi2
            0.60;       % rho_R
            0.95;       % rho_G
            0.65;       % rho_Z
            0.40;       % r_A
            4.00;       % dp_A
            0.50;       % gamma_Q
            0.0020;     % sig_R
            0.0080;     % sig_G
            0.0045];    % sig_Z
params = struct;
params.tau =      v_params(1);
params.kappa =    v_params(2);
params.psi1 =     v_params(3);
params.psi2 =     v_params(4);
params.rho_R =    v_params(5);
params.rho_G =    v_params(6);
params.rho_Z =    v_params(7);
params.r_A =      v_params(8);
params.dp_A =     v_params(9);
params.gamma_Q =  v_params(10);
params.sig_R =    v_params(11);
params.sig_G =    v_params(12);
params.sig_Z =    v_params(13);
save('params.mat', 'params');
% Set distribution types and scalings
dist_type = {'gamma', 'gamma', 'gamma', 'gamma',...
             'beta', 'beta', 'beta', 'gamma', ...
             'gamma', 'normal', 'invgamma', 'invgamma', ...
             'invgamma'};
scal = [1 1 1 1 1 1 1 1 1 1 100 100 100];

%% 2. Solve log-linearized model
dynare an_schorfheide_model noclearall;

%% 3. Load Data
load('../../data/data_final_project.mat')
y = [dy_obs, infl_obs, R_obs];

%% 4. Convert An & Schorfheide (2007) hyperparameters into matlab format
% Enables use of built-in MATLAB pdf functions
hparam_map;

%% 5. Maximize posterior likelihood
fun_post = @(x) -1*(model_like(x, y) + prior(x,dist_type, scal, hyper_mat));
[opt_params, fval, exitflag, output, grad, hessian] = fminunc(fun_post, v_params, min_options);
vcv = hessian^-1;
% Hard-coded posterior to compare vs my output
params_post = [1.4161, 0.15149, 1.5182, 0.73615, 0.5398, 0.90801, 0.58091, ...
               0.26503, 4.0308, 0.52453, 0.0020484, 0.007467, 0.0048373]; 

%% Random Walk Metropolis-Hastings
rng(72256, 'twister');
% Set initial parameter vector
c0 = 1;
c  = 0.3;
theta0 = [3 opt_params(2:3)' 0.1 opt_params(5:13)';
          3 opt_params(2:3)' 2.0 opt_params(5:13)';
          0.4 opt_params(2:3)' 0.1 opt_params(5:13)';
          0.6 opt_params(2:3)' 2.0 opt_params(5:13)']; %mvnrnd(opt_params', (c0^2).*vcv, 4);
chains = nan(n+1, length(opt_params), 4);
rej_rates = nan(1,4);
times = nan(1,4);
% for i=1:1
%     t = tic;
%     [chains(:,:,i), rej_rates(i)] = rwmh(theta0(i,:), vcv, 0.3, n, y, dist_type, scal, hyper_mat, M_, options_, oo_);
%     times(i) = toc(t);
% end
% save '../../data/rwmh_dat.mat' chains rej_rates times

for j=1:4
    % Initialize
    accept = 0;
    reject = 0;
    chains(:,:,j) = [theta0(j,:); nan(n,length(theta0(j,:)))];
    t = tic;
    % Random draws
    for i=2:n+1
        theta = mvnrnd(chains(i-1,:,j),c^2.*vcv); % candidate theta
        % Check if negative parameters
        theta(theta<0) = chains(i-1,theta<0,j);
        % Model likelihoods
        try % to account for draws where there is no unique equil (i.e., psi1 <1)
            ml_c = model_like(theta', y);
        catch
            ml_c = -Inf;
        end
            ml_b = model_like(chains(i-1,:,j)', y);
        % Prior likelihoods
        pl_c = prior(theta', dist_type, scal, hyper_mat);
        pl_b = prior(chains(i-1,:,j)', dist_type, scal, hyper_mat);
        % Posterior likelihoods
        p_c = ml_c + pl_c;
        p_b = ml_b + pl_b;
        p = min(1, exp(p_c - p_b));

        if p == 1
            chains(i,:,j) = theta;
            accept = accept + 1;
        else
            u = unifrnd(0,1);
            if p < u
                chains(i,:,j) = chains(i-1,:,j);
                reject = reject + 1;
            else
                chains(i,:,j) = theta;
                accept = accept + 1;
            end
        end
        rej_rate = reject/(i-1);
        acc_rate = accept/(i-1);

        if (mod(i,10000)==0)
            fprintf('Iteration: %d \n Rejection Rate: %d \n', i, rej_rate);
        end
    end
    rej_rates(j) = rej_rate;
    times(j) = toc(t);
end

save '../../data/part1_data.mat'
%% Diagnostic Figures
% Simulations Figures
series = {'tau', 'kappa', 'psi1', 'psi2', 'rho_r', 'rho_g', 'rho_z',...
            'r_A', 'dp_A', 'gamma_Q', 'sig_R', 'sig_G', 'sig_Z'};
% Full Chain Figures
for i = 1:4
    figure (i)
    for j=1:13
        subplot(4,4,j)
        plot(chains(:,j,i))
        xline(50001)
        title(series(j))
    end
    saveas(gcf, strcat('../../output/rwmh_results_',num2str(i)), 'png');
end
close all

% Recursive Means Figures
rmeans = nan(n+1, length(opt_params), 4);
for i = 1:4
    for j=1:(n+1)
        rmeans(j,:,i) = mean(chains(1:j,:,i));
    end
end

% Plot figures
for j=1:10
    figure (j)
    hold on
    for i=1:4
        plot(rmeans(:,j,i))
    end
    xline(50001)
    title(series(j))
    hold off
    saveas(gcf, strcat('../../output/rmeans_',char(series(j))), 'png');
end
for j=11:13
    figure (j)
    hold on
    for i=1:4
        plot(100*rmeans(:,j,i))
    end
    xline(50001)
    ylim([0, 1])
    title(strcat('100', char(series(j))))
    hold off
    saveas(gcf, strcat('../../output/rmeans_',char(series(j))), 'png');
end
close all

%% Recursive Means figures grouped
figure(1)
for j = 1:4
    subplot(2,2,j)
    hold on
    for i=1:4
        plot(rmeans(:,j,i))
    end
    xline(50001)
    xlim([0 100000])
    title(series(j))
    hold off
end
saveas(gcf, '../../output/rmeans_diag1', 'epsc');

figure(2)
for j = 5:7
    subplot(3,1,j-4)
    hold on
    for i=1:4
        plot(rmeans(:,j,i))
    end
    xline(50001)
    xlim([0 100000])    
    title(series(j))
    hold off
end
saveas(gcf, '../../output/rmeans_diag2', 'epsc');

figure(3)
ylims = [0 1;
         3 4.5;
         0.4 0.9];
for j = 8:10
    subplot(3,1,j-7)
    hold on
    for i=1:4
        plot(rmeans(:,j,i))
    end
    xline(50001)
    title(series(j))
    xlim([0 100000])    
    ylim(ylims(j-7,:))
    hold off
end
saveas(gcf, '../../output/rmeans_diag3', 'epsc');

figure(4)
ylims = [0 0.5;
         0.5 1;
         0.25 0.75];
for j = 11:13
    subplot(3,1,j-10)
    hold on
    for i=1:4
        plot(100*rmeans(:,j,i))
    end
    xline(50001)
    title(strcat('100', char(series(j))))
    xlim([0 100000])
    ylim(ylims(j-10,:))
    hold off
end
saveas(gcf, '../../output/rmeans_diag4', 'epsc');
close all
%% Plot prior vs posterior distributions
% Create prior distribution draws
chains_pr = nan(50000,length(v_params));
for i = 1:length(v_params)
    chains_pr(:,i) = draw_prior(50000,dist_type(i),hyper_mat(i,:));
end

% Keep last 50000 draws from posteriors
chain_fp = chains(50002:100001,:,1);

figure(1)
%idx = [1 3 5 7 9 11];
idx = [1 3 5];
xlims = [0 5;
         0.5 3;
         0 1];
ylims = [0 1;
         0 2;
         0 1];
for i=1:length(idx)
    % Prior Plot
    subplot(3,2,idx(i))
    plot(chains_pr(:,idx(i)), chains_pr(:,idx(i)+1), 'o', 'MarkerSize',0.4);
    xline(opt_params(idx(i)));
    yline(opt_params(idx(i) + 1));
    xlabel(series(idx(i)));
    ylabel(series(idx(i) + 1));
    xlim(xlims(i,:));
    ylim(ylims(i,:));
    if i==1
        title('Prior');
    end
    
    % Posterior Plot
    subplot(3,2,idx(i)+1)
    plot(chain_fp(:,idx(i)), chain_fp(:,idx(i)+1),'o', 'MarkerSize',0.4);
    xline(opt_params(idx(i)));
    yline(opt_params(idx(i) + 1));
    xlabel(series(idx(i)));
    ylabel(series(idx(i) + 1));
    xlim(xlims(i,:));
    ylim(ylims(i,:));
    if i==1
        title('Posterior')
    end
end
saveas(gcf, '../../output/prior_v_post_dist1', 'epsc');

figure(2)
idx = [7 9 11];
xlims = [0 1;
         0 20;
         0 0.1;
         0 0.3];
ylims = [0 5;
         0 1;
         0 0.3;
         0 0.1];
for i=1:length(idx)
    % Prior Plot
    subplot(4,2,idx(i)-6)
    plot(chains_pr(:,idx(i)), chains_pr(:,idx(i)+1), 'o', 'MarkerSize',0.4);
    xline(opt_params(idx(i)));
    yline(opt_params(idx(i) + 1));
    xlabel(series(idx(i)));
    ylabel(series(idx(i) + 1));
    if i==1
        title('Prior');
        xlim(xlims(i,:));
        ylim(ylims(i,:));
    end
    % Posterior Plot
    subplot(4,2,idx(i)-5)
    plot(chain_fp(:,idx(i)), chain_fp(:,idx(i)+1),'o', 'MarkerSize',0.4);
    xline(opt_params(idx(i)));
    yline(opt_params(idx(i) + 1));
    xlabel(series(idx(i)));
    ylabel(series(idx(i) + 1));     
    if i==1
        title('Posterior')
        xlim(xlims(i,:));
        ylim(ylims(i,:)); 
    end
end
subplot(4,2,7)
plot(chains_pr(:,12), chains_pr(:,13),'o', 'MarkerSize',0.4);
xline(opt_params(12));
yline(opt_params(13));
xlabel(series(12));
ylabel(series(13));

subplot(4,2,8)
plot(chain_fp(:,12), chain_fp(:,13),'o', 'MarkerSize',0.4);
xline(opt_params(12));
yline(opt_params(13));
xlabel(series(12));
ylabel(series(13));
saveas(gcf, '../../output/prior_v_post_dist2', 'epsc');
close all
%% Table values for posterior means
% Use second half of chain
chains_f = chains(50002:100001,:,:);
mean_chains_f = mean(chains_f, 1);
p05_chains_f = prctile(chains_f, 5, 1);
p95_chains_f = prctile(chains_f, 95, 1);
% Reshape
means = reshape(mean_chains_f, 13, 4);
p05 = reshape(p05_chains_f, 13, 4);
p95 = reshape(p95_chains_f, 13, 4);
% Export to excel
writematrix(means, '../../output/part1_post_means.xlsx');
writematrix(p05, '../../output/part1_post_p05.xlsx');
writematrix(p95, '../../output/part1_post_p95.xlsx');