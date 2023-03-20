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
params.ppsi = 0.005;             


%% Compute Asset Supply and Demand for vector of interest rates
params.xxi = xxi(1);                % Set borrowing limit to 0

% 0. Unpack Parameters
% Consumer
bbeta       = params.bbeta;
ssigma      = params.ssigma;
n_eff       = params.n_eff;       % Number of efficiency units possible
rrho        = params.rrho;
ssigma_u    = params.ssigma_u;
% Production
aalpha      = params.aalpha;
ddelta      = params.ddelta;
Z           = params.Z;
% Assets
n_asset     = params.n_asset;
n_asset_c   = params.n_asset_c;
xxi         = params.xxi;                % Borrowing limits parameter
a_upper_scale = params.a_upper_scale;
% Efficiency Units grid
tauchen_sd  = params.tauchen_sd;
% Iteration
tol         = params.tol;
max_iter_A    = params.max_iter_A;
max_iter_mu    = params.max_iter_mu;
ppsi        = params.ppsi;                % Parameter to control weighted average interest rate iterations

% 1. Productivity Grid
[v_log_eff, m_trans, pi] = tauchen(rrho, ssigma_u, ...
                                    n_eff, tauchen_sd);
v_eff = exp(v_log_eff);
N = v_eff*pi;                   % Total labor efficiency units

% 2. Steady State (Complete Markets)
r_ss = (1/bbeta - 1);
K_ss = N*((r_ss + ddelta)/(aalpha*Z))^(1/(aalpha-1));
w_ss = (1-aalpha)*Z*(K_ss/N)^aalpha;
a_ss = K_ss/N;
e_ss = v_eff(ceil(n_eff/2));
c_ss = r_ss*a_ss + w_ss*e_ss;
util_ss = (1-bbeta)*((c_ss)^(1-ssigma))/(1-ssigma);     % normalize utility

%% 3. Initialization
% Bound interest rate
r_lower = 0;
r_upper = 1/bbeta - 1;
% Define r grid; given r, z, N, define K and w
vr = linspace(r_lower+0.001, r_upper, 50);
vK = ((vr + ddelta)./(aalpha*Z*N^(1-aalpha))).^(1/(aalpha-1));
vw = (1-aalpha)*Z*vK.^aalpha*N^-aalpha;
a_upper = a_upper_scale * a_ss;

%% 3. Compute Asset Supply
A = nan(1,50);
start_time = tic;
for i = 1:length(vr)
    % Initialize asset grid
    r = vr(i);
    w = vw(i);
    nat_lim = -w*min(v_eff)/r; % natural borrowing limit
    a_lower = xxi*nat_lim; % borrowin constraints (taking into account xxi)
    %% Compute decision rules
    pack_mgm_input;
    output = mgm(mgm_input);
%     pack_egm_input;
%     output = egm(egm_input);
    
    %% Update mu iteratively using initial guess and decision rules
    % Decision rule index
    [~,m_pol_idx] = ismember(output.m_pol_fun, output.v_asset);
    mu = ones(n_asset,n_eff)./(n_asset*n_eff); % Initialize to uniform distribution
    [mu, iter, diff] = statdist(mu, m_trans, m_pol_idx, tol, max_iter_mu);
    
    %% Check Market Clearing
    A(i) = sum(output.v_asset'*mu);    % compute Aggregate assets
end

%% Plot
figure(1)
plot(vK, vr, 'b', 'LineWidth', 2)
hold on
plot(A, vr, 'r', 'LineWidth', 2)
legend('Capital Demand', 'Capital Supply', 'Location', 'Southeast')
xlabel('Assets/Capital')
ylabel('r')
saveas(gcf, '../output/cap_demand_supply_xxi0', 'epsc');


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
params.ppsi = 0.005;             


%% Compute Asset Supply and Demand for vector of interest rates
params.xxi = xxi(2);                % Set borrowing limit to 0

% 0. Unpack Parameters
% Consumer
bbeta       = params.bbeta;
ssigma      = params.ssigma;
n_eff       = params.n_eff;       % Number of efficiency units possible
rrho        = params.rrho;
ssigma_u    = params.ssigma_u;
% Production
aalpha      = params.aalpha;
ddelta      = params.ddelta;
Z           = params.Z;
% Assets
n_asset     = params.n_asset;
n_asset_c   = params.n_asset_c;
xxi         = params.xxi;                % Borrowing limits parameter
a_upper_scale = params.a_upper_scale;
% Efficiency Units grid
tauchen_sd  = params.tauchen_sd;
% Iteration
tol         = params.tol;
max_iter_A    = params.max_iter_A;
max_iter_mu    = params.max_iter_mu;
ppsi        = params.ppsi;                % Parameter to control weighted average interest rate iterations

% 1. Productivity Grid
[v_log_eff, m_trans, pi] = tauchen(rrho, ssigma_u, ...
                                    n_eff, tauchen_sd);
v_eff = exp(v_log_eff);
N = v_eff*pi;                   % Total labor efficiency units

% 2. Steady State (Complete Markets)
r_ss = (1/bbeta - 1);
K_ss = N*((r_ss + ddelta)/(aalpha*Z))^(1/(aalpha-1));
w_ss = (1-aalpha)*Z*(K_ss/N)^aalpha;
a_ss = K_ss/N;
e_ss = v_eff(ceil(n_eff/2));
c_ss = r_ss*a_ss + w_ss*e_ss;
util_ss = (1-bbeta)*((c_ss)^(1-ssigma))/(1-ssigma);     % normalize utility

%% 3. Initialization
% Bound interest rate
r_lower = 0;
r_upper = 1/bbeta - 1;
% Define r grid; given r, z, N, define K and w
vr = linspace(r_lower+0.001, r_upper, 50);
vK = ((vr + ddelta)./(aalpha*Z*N^(1-aalpha))).^(1/(aalpha-1));
vw = (1-aalpha)*Z*vK.^aalpha*N^-aalpha;
a_upper = a_upper_scale * a_ss;
%% 3. Compute Asset Supply
A = nan(1,50);
start_time = tic;
for i = 1:length(vr)
    % Initialize asset grid
    r = vr(i);
    w = vw(i);
    nat_lim = -w*min(v_eff)/r; % natural borrowing limit
    a_lower = xxi*nat_lim; % borrowin constraints (taking into account xxi)
    %% Compute decision rules
    pack_mgm_input;
    output = mgm(mgm_input);
%     pack_egm_input;
%     output = egm(egm_input);
    
    %% Update mu iteratively using initial guess and decision rules
    % Decision rule index
    [~,m_pol_idx] = ismember(output.m_pol_fun, output.v_asset);
    mu = ones(n_asset,n_eff)./(n_asset*n_eff); % Initialize to uniform distribution
    [mu, iter, diff] = statdist(mu, m_trans, m_pol_idx, tol, max_iter_mu);
    
    %% Check Market Clearing
    A(i) = sum(output.v_asset'*mu);    % compute Aggregate assets
end

%% Plot
figure(1)
plot(vK, vr, 'b', 'LineWidth', 2)
hold on
plot(A, vr, 'r', 'LineWidth', 2)
legend('Capital Demand', 'Capital Supply', 'Location', 'Southeast')
xlabel('Assets/Capital')
ylabel('r')
xlim([-20,100])
saveas(gcf, '../output/cap_demand_supply_xxi05', 'epsc');