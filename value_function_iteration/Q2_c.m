%% PS2 - (c) Do #b but use a multi-grid approach with n1 = 50, n2 = 150, and then n3 = 300.
% Use linear interpolation when moving between grids.
%
% Daniel Schwindt
% 9/17/2022

%% Set-up environment
% check whether script is run directly or not. If directly, clear all etc.
stack = dbstack;
if length(stack) == 1
    clear all
    close all
    clc
end
stoch_growth_setup
mVF50 = ones(ngrid(1), m).*(1-bbeta)*log(consumSS);

%% Compute True Solution (Guess & Verify Method)
% V(k,z) = B + Cln(k) + Dz
% Normalized
ggamma = (1-bbeta);
C = aalpha*ggamma/(1-aalpha*bbeta);
D = (1-bbeta + bbeta*C)/(1-rrho);
B = log(ggamma) + bbeta*C/ggamma*log(bbeta*C) ...
    -(1-bbeta+bbeta*C)/(1-bbeta)*log(1-bbeta + bbeta*C) ...
    + (1-bbeta + bbeta*C)/(1-bbeta)*log(A);

mVF_true = B + C*log(vGridCapital') + D*prod;
mPF_true = bbeta*C/(1- bbeta + bbeta*C).*A.*exp(prod).*(vGridCapital').^aalpha;

%% VFI Computation
maxDiff = 10.0;
iter = 0;
time1 = tic;
% Call stoch_gr_vfi function to run value function iteration
[mVF50, mPF50, iter1] = stoch_gr_vfi(vGridCapital50, mVF50, P, mY50, bbeta, maxDiff, tol);

%% Multi-Grid VFI: n1=50 to n2=150 grid points
% Interpolate Value Function grid from n1=50 to n2=150
mVF150 = linear_interp(vGridCapital50, vGridCapital150, mVF50);

% Call stoch_gr_vfi function to run value function iteration
[mVF150, mPF150, iter2] = stoch_gr_vfi(vGridCapital150, mVF150, P, mY150, bbeta, maxDiff, tol);

%% Multi-Grid VFI: n2=150 to n3=300 grid points
% Interpolate Value Function grid from n1=50 to n2=150
mVF300 = linear_interp(vGridCapital150, vGridCapital300, mVF150);

% Call stoch_gr_vfi function to run value function iteration
[mVF300, mPF300, iter3] = stoch_gr_vfi(vGridCapital300, mVF300, P, mY300, bbeta, maxDiff, tol);

%% Multi-Grid VFI: n3=300 to n=1000 grid points
% Interpolate Value Function grid from n1=50 to n2=150
mVF = linear_interp(vGridCapital300, vGridCapital, mVF300);

% Call stoch_gr_vfi function to run value function iteration
[mVF, mPF, iter4] = stoch_gr_vfi(vGridCapital, mVF, P, mY, bbeta, maxDiff, tol);

time0 = toc(time1);
iter = iter1 + iter2 + iter3 + iter4;

% Compute the maximum absolute error (Value Function & Policy Function)
err = (mVF - mVF_true)./mVF_true;
max_err0 = max(max(abs(err)));
err_pf = (mPF - mPF_true)./mPF_true;
max_err_pf0 = max(max(abs(err_pf)));

%% Reset environment
stoch_growth_setup
mVF50 = ones(ngrid(1), m).*(1-bbeta)*log(consumSS);

%% Multi-Grid VFI from n1=50 to n=1000
time = tic;
% Call stoch_gr_vfi function to run value function iteration for first grid
[mVF50, mPF50, iter5] = stoch_gr_vfi(vGridCapital50, mVF50, P, mY50, bbeta, maxDiff, tol);

% Interpolate Value Function grid from n1=50 to n2=150
mVF = linear_interp(vGridCapital50, vGridCapital, mVF50);

% Call stoch_gr_vfi function to run value function iteration
[mVF, mPF, iter6] = stoch_gr_vfi(vGridCapital, mVF, P, mY, bbeta, maxDiff, tol);
time1 = toc(time);

% Compute the maximum absolute error (Value Function & Policy Function)
err = (mVF - mVF_true)./mVF_true;
max_err1 = max(max(abs(err)));
err_pf = (mPF - mPF_true)./mPF_true;
max_err_pf1 = max(max(abs(err_pf)));

%% Reset environment
stoch_growth_setup
mVF150 = ones(ngrid(2), m).*(1-bbeta)*log(consumSS);

%% Multi-Grid VFI from n1=150 to n=1000
time = tic;
% Call stoch_gr_vfi function to run value function iteration for first grid
[mVF150, mPF150, iter7] = stoch_gr_vfi(vGridCapital150, mVF150, P, mY150, bbeta, maxDiff, tol);

% Interpolate Value Function grid
mVF = linear_interp(vGridCapital150, vGridCapital, mVF150);

% Call stoch_gr_vfi function to run value function iteration
[mVF, mPF, iter8] = stoch_gr_vfi(vGridCapital, mVF, P, mY, bbeta, maxDiff, tol);
time2 = toc(time);

% Compute the maximum absolute error (Value Function & Policy Function)
err = (mVF - mVF_true)./mVF_true;
max_err2 = max(max(abs(err)));
err_pf = (mPF - mPF_true)./mPF_true;
max_err_pf2 = max(max(abs(err_pf)));

%% Reset environment
stoch_growth_setup
mVF300 = ones(ngrid(3), m).*(1-bbeta)*log(consumSS);

%% Multi-Grid VFI from n1=300 to n=1000
time = tic;
% Call stoch_gr_vfi function to run value function iteration for first grid
[mVF300, mPF300, iter9] = stoch_gr_vfi(vGridCapital300, mVF300, P, mY300, bbeta, maxDiff, tol);

% Interpolate Value Function grid
mVF = linear_interp(vGridCapital300, vGridCapital, mVF300);

% Call stoch_gr_vfi function to run value function iteration
[mVF, mPF, iter10] = stoch_gr_vfi(vGridCapital, mVF, P, mY, bbeta, maxDiff, tol);
time3 = toc(time);

% Compute the maximum absolute error (Value Function & Policy Function)
err = (mVF - mVF_true)./mVF_true;
max_err3 = max(max(abs(err)));
err_pf = (mPF - mPF_true)./mPF_true;
max_err_pf3 = max(max(abs(err_pf)));