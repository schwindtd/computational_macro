%% Parameters Script
% Script to set parameters for Neoclassical Growth model with:
% 1. Labor-Leisure decision
% 2. Technological growth
% 3. Population growth
% 4. Calibration as in Cooley & Prescott, but sigma=20 and sigma_eps=0.035
% Author: Daniel Schwindt
% Date: 10/11/2022

%% 1. Collocation Parameters
n = 7;
m = 7;
%% 2. Define Preference and Production Parameters
% Technology
theta  = 0.4;
delta  = 0.012;
rho    = 0.95;
sigma_eps  = 0.035;
gamma  = (1.0156^0.25)-1;

% Preferences
beta   = 0.987;
sigma = 20;
alpha  = 0.64;
eta    = (1.012^0.25)-1;

% Parameters taken from Dynare results
% Steady State
css = oo_.dr.ys(1);
kss = oo_.dr.ys(2);
hss = oo_.dr.ys(3);
lss = 1-hss;
zss = oo_.dr.ys(5);
zscal = ((1/beta) - 1 + delta)/theta;

% Linear Decision Rule Coeffs
coef_ck = oo_.dr.ghx(3,1);
coef_cz = oo_.dr.ghx(3,2);
coef_ce = oo_.dr.ghu(3);
coef_kk = oo_.dr.ghx(1,1);
coef_kz = oo_.dr.ghx(1,2);
coef_ke = oo_.dr.ghu(1);
coef_hk = oo_.dr.ghx(4,1);
coef_hz = oo_.dr.ghx(4,2);
coef_he = oo_.dr.ghu(4);

% State variable spaces
k_max = 1.5*kss;
k_min = 0.5*kss;
%% 3. Fsolve Options
options_use = optimset('display','iter','MaxFunEvals',100000,'MaxIter',1000,'TolFun',10E-15,'FunValCheck','On');
%% 4. Basis Functions
x_n = cos(pi*(2*[n:-1:1]'-1)/(2*n));  % Roots of the nth order Chebysev polynomial
k = (k_max-k_min)*(x_n+1)/2+k_min;
[z,P] = tauchen(rho,sigma_eps,m);
z = z';
ar = autoreg(0, 0.035,0.95, 100000); % simulate 100,000 times
Z_max = 0.5; % Result of simulating 100000 times, sigma_eps = 0.035
Z_min = -0.5; % Result of simulating 100000 times, sigma_eps = 0.035
x_m = 2 * (z - Z_min) / (Z_max - Z_min) - 1;
psi_k = cos([0:n-1]'*acos(x_n')); % basis functions for K
psi_Z = cos([0:m-1]'*acos(x_m')); % basis functions for Z
%% 5. Tensor Products
psi = kron(psi_k,psi_Z);

% Construct the vectors that will carry the collocation points 
% in the order given by the Tensor product. 

K = kron(k,ones(m,1));

Z = kron(ones(n,1),z);
%% 6. Define structure to store parameter values
par.theta = theta;
par.delta = delta;
par.rho = rho;
par.sigma_eps = sigma_eps;
par.gamma = gamma;
par.beta = beta;
par.sigma = sigma;
par.alpha = alpha;
par.eta = eta;
par.k_min = k_min;
par.k_max = k_max;
par.Z_min = Z_min;
par.Z_max = Z_max;
par.n = n;
par.m = m;
par.zscal = zscal;