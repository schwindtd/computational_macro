function ll = model_like(vp, data)
% Computes likelihood for model using parameters and data as inputs
%   Steps: 
%   1. Solves log-linearized model to obtain decision rules
%   2. Defines matrix objects for Kalman filter using outputs from (1)
%   3. Runs Kalman filter using outputs from (2)
%   4. Computes log-likelihood based on Kalman filter outputs from (3)
%% 1. Solve log-linearized model
    global M_ oo_ options_
    % Re-set parameter structure values
    params = struct;
    params.tau =      vp(1);
    params.kappa =    vp(2);
    params.psi1 =     vp(3);
    params.psi2 =     vp(4);
    params.rho_R =    vp(5);
    params.rho_G =    vp(6);
    params.rho_Z =    vp(7);
    params.r_A =      vp(8);
    params.dp_A =     vp(9);
    params.gamma_Q =  vp(10);
    params.sig_R =    vp(11);
    params.sig_G =    vp(12);
    params.sig_Z =    vp(13);
    set_param_value('tau', params.tau);
    set_param_value('kappa', params.kappa);
    set_param_value('psi1', params.psi1);
    set_param_value('psi2', params.psi2);
    set_param_value('dp_A', params.dp_A);
    set_param_value('gamma_Q', params.gamma_Q);
    set_param_value('r_A', params.r_A);
    set_param_value('rho_G', params.rho_G);
    set_param_value('rho_R', params.rho_R);
    set_param_value('rho_Z', params.rho_Z);
    set_param_value('sig_G', params.sig_G);
    set_param_value('sig_R', params.sig_R);
    set_param_value('sig_Z', params.sig_Z);
    % Re-solve model using re-set parameters
    [dr, info, M_, oo_] = resol(0,M_, options_, oo_);
    %% Compute model likelihood
    if info == 0
        %% Create State Space
        state_space_create;
        %% 3. Run Kalman Filter
        [v,F, K, L, a, P] = kalmanfilter(ssystem.Z, ssystem.H0, ssystem.T, ...
                                        ssystem.R, ssystem.Q0, data, ...
                                        ssystem.a0, ssystem.m, ssystem.cons);

        %% 4. Compute log-likelihood from Kalman Filter
        ll = kloglike(F, v);
    else
        ll= -Inf;
    end

end

