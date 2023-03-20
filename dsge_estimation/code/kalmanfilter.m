function [v,F, K, L, a, P] = kalmanfilter(Z, H0, T, R, Q0, y, a0, m, c)
%% Matlab Implementation of Kalman Filter
%   Based on notes from ECON630 and Durbin & Koopman
%%   General State Space Model
%   y_t = c + Z_t*alpha_t + eps_t, eps_t ~ N(0, H_t)
%   alpha_t+1 = T_t*alpha_t + R_t*eta_t+1, eta_t+1 ~ N(0, Q_t)
%   alpha_t ~ N(a_t, P_t)

%% Inputs
% Z: matrix transforming states to observations
% H0: observation noise covariance matrix
% T: matrix transforming prior states to current states
% R: linear transform of state noise
% Q0: state process noise covariance matrix
% y: observed data
% a0: initial state vector
% P0: initial state covariance matrix

% Set options
options = optimset('Display', 'off', 'TolX', 1e-15, 'TolFun', 1e-15);
%% Create storage objects
dim_y = size(y);
v = nan(dim_y(1), dim_y(2));
F = nan(dim_y(2), dim_y(2), dim_y(1));
K = nan(m, dim_y(2), dim_y(1));
L = nan(m, m, dim_y(1));
a = nan(dim_y(1)+1, m);
P = nan(m, m, dim_y(1)+1);
%% Loop over all observations
% Set starting a and P for loop
a(1,:) = a0;
%P0 = 1e6*eye(m);
func_p = @(p) p - T*p*T' - R*Q0*R';
[P0,func_p_sol,exitfflag] = fsolve(func_p, eye(m), options);
if exitfflag <= 0
    P0 = 1e6*eye(m);
end
    
P(:,:,1) = P0;

% Kalman Filter Loop
for i=1:dim_y(1)
    P_t = P(:,:,i);
    %P_t = T*P(:,:, i)*T' + R*Q0*R';
    v(i,:) = (y(i,:)' - c' - Z*a(i,:)')'; % Prediction error
    F(:,:,i) = Z*P_t*Z' + H0;         
    inv_F_t = (F(:,:,i))^-1;
    %K(:,:,i) = P_t*Z'*inv_F_t;        % Kalman Gain
    K(:,:,i) = T*P_t*Z'*inv_F_t;     % Kalman Gain
    L(:,:,i) = T - K(:,:,i)*Z;
    a(i+1,:) = (T*a(i,:)' + K(:,:,i)*v(i,:)')';     % Next period state
    P(:,:,i+1) = T*P_t*(T-K(:,:,i)*Z)' + R*Q0*R';   % Next period state noise covariance matrix
    %P(:,:,i+1) = (T - K(:,:,i)*Z)*P(:,:,i);
end
% Remove initialization from chain
a = a(2:dim_y(1)+1,:);
P = P(:,:,2:dim_y(1)+1);
end

