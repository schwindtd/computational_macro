function [r, w, K, A, N, mu, v_asset, v_eff, m_pol_fun, m_V_fun] = srce(params)
%Wrapper function to compute SRCE

%% 0. Unpack Parameters
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
max_iter_A  = params.max_iter_A;
max_iter_mu    = params.max_iter_mu;
ppsi        = params.ppsi;                % Parameter to control weighted average interest rate iterations

%% 1. Productivity Grid
[v_log_eff, m_trans, pi] = tauchen(rrho, ssigma_u, ...
                                    n_eff, tauchen_sd);
v_eff = exp(v_log_eff);
N = v_eff*pi;                   % Total labor efficiency units

%% 2. Steady State (Complete Markets)
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
% Initialize random r; given r, z, N, initialize k; compute w
r = r_lower + (r_upper - r_lower)*rand();
K = ((r + ddelta)/(aalpha*Z*N^(1-aalpha)))^(1/(aalpha-1));
w = (1-aalpha)*Z*K^aalpha*N^-aalpha;
% Initialize asset grid
a_upper = a_upper_scale * a_ss;
nat_lim = -w*min(v_eff)/r; % natural borrowing limit
a_lower = xxi*nat_lim; % borrowin constraints (taking into account xxi)

%% 3. Loop over r to find r that clears capital market
diff_A = 100;
r_iter = 0;
start_time = tic;
while (diff_A > tol) && (r_iter < max_iter_A)
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
    A = sum(output.v_asset'*mu);    % compute Aggregate assets
    diff_A = abs(A-K);               % gap between assets and capital demand
    % Restrict total assets to be above 0 to prevent complex values
    if A < 0
        A = 0.01;
    end
    % Compute interest rate that would clear A assets
    r_clear = aalpha*Z*(N/A)^(1-aalpha)-ddelta;
    r_new = ppsi*r_clear + (1-ppsi)*r;    % define new r as weighted avg of guess & market clearing rate
    
    % Ensure new rate guess is within bounds
    if r_new < r_lower
        r_new = r_lower;
    elseif r_new > r_upper
        r_new = r_upper;
    end
    
    r_iter = r_iter + 1;

    if (mod(r_iter,10)==0 || r_iter ==1)
        fprintf('-------------------------------------------- \n')
        fprintf('r guess: %d, \n iterations: %d, \n diff_A: %e \n', r, r_iter, diff_A)
        fprintf('-------------------------------------------- \n')
    end
        
    %% Re-initialize for next iteration
    r = r_new;
    K = ((r + ddelta)/(aalpha*Z*N^(1-aalpha)))^(1/(aalpha-1));
    w = (1-aalpha)*Z*K^aalpha*N^-aalpha;
    a_upper = a_upper_scale * a_ss;
    nat_lim = -w*min(v_eff)/r; % natural borrowing limit
    a_lower = xxi*nat_lim; % borrowin constraints (taking into account xxi)
end
v_asset = output.v_asset;
m_pol_fun = output.m_pol_fun;
m_V_fun = output.m_V_fun;
end

