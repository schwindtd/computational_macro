%% State Space Create
% Creates matrices and vectors defining the state space system
% Used as inputs to model_pm function

%% Define system of equations
% State Space Model:
% State Equations:
% x_t = [R_t, g_t, z_t, y_t, y_t-1, dp_t]'
% eps_t = [eps_R_t, eps_G_t, eps_Z_t]'
% x_t = T*x_t-1 + R*eps_t
% 
% Observed Equations:
% y_t = [YGR_t, INFL_t, INT_t]'
% D   = [1 0 0
%        0 1 0
%        4 1 1]
% con = [gamma_Q, dp_A, r_A]'
% y_t = Z*x_t + D*con
 
ssystem = struct;
ssystem.m = 6;
ssystem.Z = [0 0 100 100 -100 0;
     0 0 0 0 0 400;
     400 0 0 0 0 0];
% ssystem.Z = [0 0 100 100 0 0;
%       0 0 0 0 0 400;
%       400 0 0 0 0 0];

ssystem.T = [oo_.dr.ghx(2:5,:), zeros(4,3);
    zeros(1,3), 1, zeros(1,2);
    oo_.dr.ghx(6,:), zeros(1,3)];

ssystem.R = [oo_.dr.ghu(2:5,:);zeros(1,3);oo_.dr.ghu(6,:)];

d_q = [params.sig_R, params.sig_G, params.sig_Z];
ssystem.Q0 = diag((d_q).^2);
%d_h = [params.sig_ygr, params.sig_infl, params.sig_int]; % Use if
%measurement error is non-zero
ssystem.H0 = zeros(3);%diag(d_h);

ssystem.a0 = zeros(ssystem.m,1);
ssystem.D = [1 0 0; 0 1 0; 4 1 1];
cons_par = [params.gamma_Q, params.dp_A, params.r_A]';
ssystem.cons = (ssystem.D*cons_par)';