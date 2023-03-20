function [b0, b1, K, A, C, m_k_nxt, m_c_star, iter_b, max_diff_b, iter_k, max_diff_k] = solve_ks_im(params, z, v_eff, lbar, b0, b1, P_11, P_12, P_21, P_22,v_k, v_kbar, m_k, m_kbar, shocks_z, shocks_z_idx, shocks_e, shocks_e_idx)

% While loop to solve for b0 and b1 parameters
max_diff_b = 100;
iter_b = 0;
while (max_diff_b > params.tol_b && iter_b < params.max_iter)
    [m_k_nxt, m_c_star, iter_k, max_diff_k] = pfi_im(params, z, v_eff, lbar,b0, b1, P_11, P_12, P_21, P_22, v_k, v_kbar, m_k, m_kbar);
    %% Simulation
    z_idx = linspace(1,params.n_agg,params.n_agg)';
    e_idx = linspace(1,params.n_eff, params.n_eff)';
    A = nan(params.T,params.I);
    C = nan(params.T, params.I);
    A(1,:) = params.kbar_min*ones(1,params.I); % Start economy at min of kbar
    K = nan(params.T,1);
    K(1) = sum(A(1,:))/params.I;
    for i=2:params.T
        % Doing interpolation for both agg and ind shocks at same time is VERY
        % SLOW, break up into two steps instead
    %     A_tmp = interpn(v_k, v_kbar, z_idx, e_idx, m_k_nxt, A(i-1,:), K(i-1),shocks_z_idx(i), shocks_e_idx(i,:), 'spline');
    %     A_tmp2 = squeeze(A_tmp); % Convert 4D array to matrix
    %     A(i,:) = A_tmp2(1,:);
        % Interpolate for aggregate shock since K might not be on the true grid
        % Sample pts: vector of asset values, kbar from last period, and the
        % aggregate shock today
        A_tmp = interpn(v_k, v_kbar, z_idx, e_idx, m_k_nxt, v_k, K(i-1),shocks_z_idx(i), e_idx); 
        A_tmp = squeeze(A_tmp);
        C_tmp = interpn(v_k, v_kbar, z_idx, e_idx, m_c_star, v_k, K(i-1), shocks_z_idx(i), e_idx);
        C_tmp = squeeze(C_tmp);
        % Interpolate for individual shock since A_tmp might not be on true
        % grid
        % Sample pts: asset choice last period & idiosyncratic shock today
        A_tmp2 = interpn(v_k, z_idx, A_tmp, A(i-1,:), shocks_e_idx(i,:));
        C_tmp2 = interpn(v_k, z_idx, C_tmp, A(i-1,:), shocks_e_idx(i,:));
        % Check bounds for each person to make sure do not go outside bounds
        A_tmp2(A_tmp2 < params.k_min) = params.k_min;
        A_tmp2(A_tmp2 > params.k_max) = params.k_max;
        A(i,:) = A_tmp2;        
        % Build aggregate capital
        K(i) = sum(A(i,:))/params.I;
        if K(i) < params.kbar_min
            K(i) = params.kbar_min;
        elseif K(i) > params.kbar_max
            K(i) = params.kbar_max;
        end
    end

    %% Regression
    % Drop first BURN periods
    K_f = K(params.burn+1:end);
    shocks_z_idxf = shocks_z_idx(params.burn+1:end);
    % Create regression variables
    lK = log(K_f);
    s_z_dummy = (shocks_z_idxf>1)*0 + (shocks_z_idxf<=1)*1; % shocks dummy (=1 if good shock)
    interact = s_z_dummy.*lK;   % Interaction term
    y = lK(2:end);
    X = [ones(length(y),1),s_z_dummy(1:(end-1)),lK(1:(end-1)),interact(1:(end-1))]; 
    % Compute coefficients and save to new vectors
    B_reg = (X'*X)^-1*X'*y;
    b0_new = [B_reg(1) + B_reg(2), B_reg(1)];
    b1_new = [B_reg(3) + B_reg(4), B_reg(3)];
    i_b = 0; i_g = 0;
    x_b = 0; x_g = 0;
    y_b = 0; y_g = 0;
    for i=(params.burn+1):(params.T-1)
        if shocks_z_idx(i) == 1
            i_b = i_b + 1;
            x_b(i_b,1) = log(K(i));
            y_b(i_b,1) = log(K(i+1));
        else
            i_g = i_g + 1;
            x_g(i_g,1) = log(K(i));
            y_g(i_g,1) = log(K(i+1));
        end
    end
    % Bad states regression
    X_b = [ones(i_b,1) x_b];
    %b_b = (X_b'*X_b)^-1*X_b'*y_b;
    % Good states regression
    X_g = [ones(i_g,1) x_g];
    %b_g = (X_g'*X_g)^-1*X_g'*y_g;
    % store in new vectors
    %b0_new = [b_b(1) b_g(1)];
    %b1_new = [b_b(2) b_g(2)];
    % Compute R-squared
    %r2_b = rsq(y_b, X_b, b_b);
    %r2_g = rsq(y_g, X_g, b_g);
    r2_b = rsq(y_b, X_b, [b0_new(1) b1_new(1)]');
    r2_g = rsq(y_g, X_g, [b0_new(2) b1_new(2)]');
    max_diff_b = max(max(abs(b0_new - b0)),max(abs(b1_new - b1))); % compute coeff difference bw iterations
    % Update coefficients for next iteration as weighted average per Maliar,
    % Maliar, Valli (2010)
    b0 = params.psi_b*b0_new + (1-params.psi_b)*b0;
    b1 = params.psi_b*b1_new + (1-params.psi_b)*b1;

    iter_b = iter_b + 1;

       if (mod(iter_b,1)==0 || iter_b ==1)
          fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter_b, max_diff_b); 
          fprintf(' Bad State Coeffs: %2.4f %2.4f, Rsq = %2.4f\n', b0(1), b1(1), r2_b);
          fprintf(' Good State Coefs: %2.4f %2.4f, Rsq = %2.4f\n', b0(2), b1(2), r2_g);
          disp('----------------------------------------');
       end

    end
end

