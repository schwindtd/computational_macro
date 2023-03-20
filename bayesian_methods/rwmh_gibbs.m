%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECON630 PS4 Part B
% Markov Chain Monte-Carlo
% Author: Daniel Schwindt
% Date: 10/29/2022

%% 0. Housekeeping
clear all;
close all;
clc;

%% 1. Set Parameters
mu = [0,0];
sigma = [1, 0.5; 0.5,1];
mu_eta = [0,0];
sigma_eta = eye(2);
c = 0.1;
lwr = 0.2;
upp = 0.4;
N = 100000;
burn_out = 5000;

%% 2. Set initial value
init = zeros(1,2);
tic;
%% 3. RW Metropolis-Hastings Algorithm
% Preliminary run using sigma_eta and c=0.1
% Initial value for x = [0,0]
[x_1,acc_rate_1] = rwmh(init,N, mu, sigma, mu_eta, sigma_eta, c);
% Estimate VCV of draws and use 
sigma_est = 1/(N+1)*(x_1'*x_1);
mu_est = 1/(N+1)*sum(x_1);
% Final run using sigma_est and c=3 to keep
% acceptance rate within [0.2, 0.4]
[x_f,acc_rate_f] = rwmh(init,N, mu, sigma, mu, sigma_est, 3);
time_rwmh = toc;
%% 4. Compute Estimated Moments (RWMH)
x_f_b = x_f((burn_out+1):N,:);
mx = sum(x_f_b)/(N-burn_out);
mx2 = sum(x_f_b.^2)/(N-burn_out);
mx3 = sum(x_f_b.^3)/(N-burn_out);
mx4 = sum(x_f_b.^4)/(N-burn_out);
mxy = sum(x_f_b(:,1).*x_f_b(:,2))/(N-burn_out);
%% 5. Gibbs Sampling
x_g = zeros(N+1, 2);
for i=2:N
    % Define mean and variance of x1 conditional on  x2
    mu_cond1 = mu(1) + sigma(1,2)*sigma(2,2)^-1*(x_g(i-1,2)-mu(2));
    sigma_cond1 = sqrt(sigma(1,1) - sigma(1,2)*sigma(2,2)^-1*sigma(2,1));
    % Draw x1_i from distribution conditional on x2_(i-1)
    x_g(i,1) = normrnd(mu_cond1, sigma_cond1);
    % Define mean and variance of x2 conditional on x1
    mu_cond2 = mu(2) + sigma(2,1)*sigma(1,1)^-1*(x_g(i,1) - mu(1));
    sigma_cond2 = sqrt(sigma(2,2) - sigma(2,1)*sigma(1,1)^-1*sigma(1,2));
    % Draw x2_i from distbution conditional on x1_i
    x_g(i,2) = normrnd(mu_cond2, sigma_cond2);
end
time_gibbs = toc - time_rwmh;
%% 6. Compute Estimated Moments (Gibbs)
x_g_b = x_g((burn_out+1):N,:);
mx_g = sum(x_g_b)/(N-burn_out);
mx2_g = sum(x_g_b.^2)/(N-burn_out);
mx3_g = sum(x_g_b.^3)/(N-burn_out);
mx4_g = sum(x_g_b.^4)/(N-burn_out);
mxy_g = sum(x_g_b(:,1).*x_g_b(:,2))/(N-burn_out);

%% 7. Analytical Results
syms m1 m2 t1 t2 s11 s12 s22
mgf = @(t1, t2, m1, m2, s11, s12, s22) ...
    exp(m1*t1 + m2*t2 + 0.5*(s11*t1^2+2*s12*t1*t2 + s22*t2^2));
a_mom = zeros(2,5);
for i=1:4
    df = matlabFunction(diff(mgf, t1, i));
    a_mom(1,i) = df(mu(1), mu(2), sigma(1,1), sigma(1,2), sigma(2,2),0, 0);
    df = matlabFunction(diff(mgf, t2, i));
    a_mom(2,i) = df(mu(1), mu(2), sigma(1,1), sigma(1,2), sigma(2,2),0, 0);
end

df = matlabFunction(diff(diff(mgf, t1),t2));
a_mom(1,5) = df(mu(1), mu(2), sigma(1,1), sigma(1,2), sigma(2,2),0, 0);
a_mom(2,5) = df(mu(1), mu(2), sigma(1,1), sigma(1,2), sigma(2,2),0, 0);

%% 7. Output data
header = {'Moment','Analytical','RWMH','Gibbs'};
labs = {'E[X]'; 'E[X^2]'; 'E[X^3]'; 'E[X^4]'; 'E[XY]'};
% Concatenate Data and Labels
x_dat = [a_mom(1,:);
         mx(1), mx2(1), mx3(1), mx4(1), mxy;
         mx_g(1), mx2_g(1), mx3_g(1), mx4_g(1), mxy_g]';
y_dat = [a_mom(2,:);
         mx(2), mx2(2), mx3(2), mx4(2), mxy;
         mx_g(2), mx2_g(2), mx3_g(2), mx4_g(2), mxy_g]';
x_out = [header; [labs,num2cell(x_dat)]];
y_out = [header; [labs,num2cell(y_dat)]];
% Write to CSV
writecell(x_out, '../output/x_moments.csv');
writecell(y_out, '../output/y_moments.csv');

%% 8. Create Histograms
% Metropolis-Hastings
figure(1);
hist3(x_f_b, 'NBins',[50,50], 'CDataMode', 'auto','FaceColor','interp');
title('Histogram: RW Metropolis-Hastings');
xlabel('X1');
ylabel('X2');
saveas(gcf, '../output/hist_rwmh','epsc');

% Gibbs Sampling
figure(2);
hist3(x_g_b, 'NBins',[50,50],'CDataMode', 'auto','FaceColor','interp');
title('Histogram: Gibbs Sampling');
xlabel('X1');
ylabel('X2');
saveas(gcf, '../output/hist_gibbs','epsc');

% Metropolis-Hastings
figure(3);
hist3(x_f_b, 'NBins',[50,50],'CDataMode', 'auto');
title('Histogram: RW Metropolis-Hastings');
xlabel('X1');
ylabel('X2');
colorbar;
view(2);
saveas(gcf, '../output/cont_rwmh','epsc');

% Gibbs Sampling
figure(4);
hist3(x_g_b, 'NBins',[50,50],'CDataMode', 'auto');
title('Histogram: Gibbs Sampling');
xlabel('X1');
ylabel('X2');
colorbar;
view(2);
saveas(gcf, '../output/cont_gibbs','epsc');
close all;