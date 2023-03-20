function output = mgm(input)
% Multi-grid method exploiting monotonicity
%   Computes policy function using multi-grid iteration method
%% Unpack input
 n_asset_c = input.c_grid_size;
 n_asset = input.f_grid_size;
 a_upper = input.a_upper;
 a_lower = input.a_lower;
 v_eff = input.v_eff;
 m_trans = input.m_trans;
 util_ss = input.util_ss;
 bbeta = input.bbeta;
 ssigma = input.ssigma;
 ddelta = input.ddelta;
 n_eff = input.n_eff;
 r = input.r;
 w = input.w;
 tol = input.tol;
 max_iter = input.max_iter_A;

 %% Initialization
times = nan(length(n_asset_c), 1);
mIterations = nan(length(n_asset_c), 2);

for nMultigrid = 1:length(n_asset_c)
    start_time = tic;
    %% Initialization of coarse grid
    
    % Coarse capital grid
    nAsset = n_asset_c(nMultigrid);
    v_asset = linspace(a_lower, a_upper, nAsset)';
    
    % initialize function values
    mV_fun = util_ss.*ones(nAsset, n_eff);
    mV_fun_new = zeros(nAsset, n_eff);
    mPol_fun = zeros(nAsset, n_eff);
    mV_fun_exp = zeros(nAsset, n_eff);
    
    
    %% Recursive solution on coarse grid
    diff = 100;
    
    iter = 0;
    
    while (diff > tol) && (iter < max_iter)
    
        % Calculate expected value function (only depends on mV_fun)
        mV_fun_exp = mV_fun*m_trans'; 
        
        % Loop over productivity 
        for i_eff = 1:n_eff
            % Fetch productivity (for readability)
            e = v_eff(i_eff);

            % initialise lower bound for k' deduced from monotonicity
            iPrevPolicy = 1;
            
            % Loop over capital 
            for iAsset = 1:nAsset
                % Calculate vector of consumption (for k' above threshold)
                cons = (1+r)*v_asset(iAsset) + w*e - v_asset(iPrevPolicy:end);
                % Set negative consumption equal to small number
                cons(cons <= 0) = 1e-10;

                % find new values of functions at (a, e)
                [mV_fun_new(iAsset, i_eff), index] = ...
                    max((1-bbeta).*(cons.^(1-ssigma))./(1-ssigma) + bbeta.*mV_fun_exp(iPrevPolicy:end,i_eff));
                mPol_fun(iAsset, i_eff) = v_asset(index+iPrevPolicy-1);

                % reset lower bound on k'
                iPrevPolicy = index+iPrevPolicy-1;

            end
           
         end
        
        diff = max(max(abs(mV_fun_new - mV_fun)));
        mV_fun = mV_fun_new;
    
        iter = iter + 1;
    %     if (mod(iter,10)==0 || iter ==1)
    %         fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, diff); 
    %     end
    
    end

%     fprintf(' Coarse grid: \n'); 
%     fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, diff); 
%     fprintf('\n')
    
    %% Interpolation 

    % save coarse capital grid
    pre_v_asset = v_asset;

    % create dense capital grid
    v_asset = linspace(a_lower, a_upper, n_asset)';

    % interpolate value function at each productivity level
    mV_fun = zeros(n_asset, n_eff);
    for i_eff = 1:n_eff
        mV_fun(:, i_eff) = interp1(pre_v_asset, mV_fun_new(:,i_eff), v_asset);
    end
    
    % initialise other function values
    mV_fun_new = zeros(n_asset, n_eff);
    mPol_fun = zeros(n_asset, n_eff);
    mV_fun_exp = zeros(n_asset, n_eff);


    %% Recursive solution on dense grid
    diff = 100;
    
    iter = 0;
    
    while (diff > tol) && (iter < max_iter)
    
        % Calculate expected value function (only depends on mV_fun)
        mV_fun_exp = mV_fun*m_trans'; 
        
        % Loop over productivity 
        for i_eff = 1:n_eff
            % Fetch productivity (for readability)
            e = v_eff(i_eff);

            % initialise lower bound for k' deduced from monotonicity
            iPrevPolicy = 1;
            
            % Loop over capital 
            for iAsset = 1:n_asset
                % Calculate vector of consumption (for k' above threshold)
                cons = (1+r)*v_asset(iAsset) + w*e - v_asset(iPrevPolicy:end);
                % Set negative consumption equal to small number
                cons(cons <= 0) = 1e-10;

                % find new values of functions at (a, e)
                [mV_fun_new(iAsset, i_eff), index] = ...
                    max((1-bbeta).*(cons.^(1-ssigma))./(1-ssigma) + bbeta.*mV_fun_exp(iPrevPolicy:end,i_eff));
                mPol_fun(iAsset, i_eff) = v_asset(index+iPrevPolicy-1);

                % reset lower bound on k'
                iPrevPolicy = index+iPrevPolicy-1;

            end

        end
        
        diff = max(max(abs(mV_fun_new - mV_fun)));
        mV_fun = mV_fun_new;
    
        iter = iter + 1;
%         if (mod(iter,10)==0 || iter ==1)
%             fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, diff); 
%         end
    
    end
    

%     fprintf(' Dense grid: \n'); 
%     fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, diff); 
%     fprintf('\n')
%     
    times(nMultigrid) = toc(start_time);
    mIterations(nMultigrid, 2) = iter;

end
%% Output 
% Initialize output structure
output = struct;
output.m_pol_fun = mPol_fun;
output.m_V_fun = mV_fun;
output.v_asset = v_asset;
end

