clear all; clc;

load('../../data/data_final_project.mat')
y = [dy_obs, infl_obs, R_obs];
Z = [1 0; 0.5 0.5; 0.5 1];
m=2;
H0 = eye(3);
%d = [0.9, 0.9];
T = eye(m); %[0.5 0.25; 0.25 0.5]; 
R = eye(m);
Q0 = eye(m);
a0 = zeros(m,1);
P0 = 1000*eye(m);
c = zeros(1,3); %[0.2 0.3 0.1];

% Run Kalman Filter
[v,F, K, L, a, P] = kalmanfilter(Z, H0, T, R, Q0, y, a0, m, c);

% compute predicted y's
y_p = nan(150,3);
for i=1:150
    y_p(i,:) = a(i,:)*Z';
end

figure(1)
for i=1:3
    subplot(3,1,i)
    plot(y(:,i))
    hold on
    plot(y_p(:,i))
    hold off
end

% Run Kalman Smoother
%[r, alpha_h, N, V] = kalmansmooth(Z,v,F, K, L, a, P);

% Output data
%writematrix(a, '../output/kalman_test_a.csv');

%% Check Stock & Watson
obs = xlsread('../../data/stock_watson_data.xlsx');
%obs = obs(1:(end-1),:);
coeffs = xlsread('../../data/stock_watson_coeffs.xlsx');
% State equations
m = 4;
T = [coeffs(1,5) 0 0 0;
    0 coeffs(1,4) 0 0;
    0 0 coeffs(2,4) 0;
    0 0 0 coeffs(3,4)];
R = eye(m);
Q0 = [1-coeffs(1,5)^2 0 0 0;
      0 exp(coeffs(1,3)) 0 0;
      0 0 exp(coeffs(2,3)) 0;
      0 0 0 exp(coeffs(3,3))];
% Observed equations
Z = [exp(coeffs(1,2)) 1 0 0;
     coeffs(2,2) 0 1 0;
     coeffs(3,2) 0 0 1];
 c = coeffs(:,1)'; %zeros(1,3); 
 H0 = zeros(3,3);

% Initial state vector
a0 = zeros(m,1);

% Run Kalman Filter (original)
[v,F, K, L, a, P] = kalmanfilter(Z, H0, T, R, Q0, obs, a0, m, c);
ll = kloglike(F,v);
% Run Kalman Filter (test)
[v_test,F_test, K_test, L_test, a_p, a_u, P_p, P_u] = kalmanfilter_test(Z, H0, T, R, Q0, obs, a0, m, c);
ll_test = kloglike(F_test, v_test);