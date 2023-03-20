function [P_z, P_ze, P_1, P_2, P_11, P_12, P_21, P_22] = trans_prob(params)

    % Transition probabilities
    % Aggregate
    P_z = [1-1/params.d_b 1/params.d_b; 1/params.d_g 1-1/params.d_g];

    % Idiosyncratic
    % Col/Row (ses'e')
    P_ze = nan(params.n_agg, params.n_agg, params.n_eff, params.n_eff);
%     P_ze(2,1,2,1) = (params.dur_g - 1)/params.dur_g;
%     P_ze(1,1,1,1) = (params.dur_b - 1)/params.dur_b;
%     P_ze(2,1,2,2) = P_z(2,2) - P_ze(2,1,2,1);
%     P_ze(1,1,1,2) = P_z(1,1) - P_ze(1,1,1,1);
%     P_ze(2,1,1,1) = 1.25*P_z(2,1)/P_z(1,1)*P_ze(1,1,1,1);
%     P_ze(1,1,2,1) = 0.75*P_z(1,2)/P_z(2,2)*P_ze(2,1,2,1);
%     P_ze(2,1,1,2) = P_z(2,1) - P_ze(2,1,1,1);
%     P_ze(1,1,2,2) = P_z(1,2) - P_ze(1,1,2,1);
%     P_ze(1,2,1,1) = params.u_b*(P_z(1,1) - P_ze(1,1,1,1))/(1-params.u_b);
%     P_ze(1,2,2,1) = (params.u_g*P_z(1,2) - params.u_b*P_ze(1,1,2,1))/(1-params.u_b);
%     P_ze(2,2,1,1) = (params.u_b*P_ze(2,1) - params.u_g*P_ze(2,1,1,1))/(1-params.u_g);
%     P_ze(2,2,2,1) = params.u_g*(P_z(2,2) - P_ze(2,1,2,1))/(1-params.u_g);
%     P_ze(1,2,1,2) = P_z(1,1) - P_ze(1,2,1,1);
%     P_ze(1,2,2,2) = P_z(1,2) - P_ze(1,2,2,1);
%     P_ze(2,2,1,2) = P_z(2,1) - P_ze(2,2,1,1);
%     P_ze(2,2,2,2) = P_z(2,2) - P_ze(2,2,2,1);
    P_ze(1,1,1,1) = (params.dur_b - 1)/params.dur_b;
    P_ze(1,1,1,2) = 1 - P_ze(1,1,1,1);
    P_ze(1,2,1,1) = (params.u_b - params.u_b*P_ze(1,1,1,1))/(1-params.u_b);
    P_ze(1,2,1,2) = 1 - P_ze(1,2,1,1);
    
    P_ze(2,1,1,1) = 1.25*P_ze(1,1,1,1);
    P_ze(2,1,1,2) = 1 - P_ze(2,1,1,1);
    P_ze(2,2,1,1) = (params.u_b - params.u_g*P_ze(2,1,1,1))/(1-params.u_g);
    P_ze(2,2,1,2) = 1 - P_ze(2,2,1,1);
    
    P_ze(2,1,2,1) = (params.dur_g - 1)/params.dur_g;
    P_ze(2,1,2,2) = 1 - P_ze(2,1,2,1);
    P_ze(2,2,2,1) = (params.u_g - params.u_g*P_ze(2,1,2,1))/(1-params.u_g);
    P_ze(2,2,2,2) = 1 - P_ze(2,2,2,1);
    
    P_ze(1,1,2,1) = 0.75*P_ze(2,1,2,1);
    P_ze(1,1,2,2) = 1 - P_ze(1,1,2,1);
    P_ze(1,2,2,1) = (params.u_g - params.u_b*P_ze(1,1,2,1))/(1-params.u_b);
    P_ze(1,2,2,2) = 1 - P_ze(1,2,2,1);

    %% Convert into 4D format for (n_k, n_kbar, n_agg, n_eff)
    P_11 = nan(params.n_k, params.n_kbar, params.n_agg, params.n_eff);
    P_12 = nan(params.n_k, params.n_kbar, params.n_agg, params.n_eff);
    P_21 = nan(params.n_k, params.n_kbar, params.n_agg, params.n_eff);
    P_22 = nan(params.n_k, params.n_kbar, params.n_agg, params.n_eff);
    for i=1:params.n_agg
        for j=1:params.n_eff
%             P_11(:,:,i,j) = P_ze(i,j,1,1).*ones(params.n_k, params.n_kbar);
%             P_12(:,:,i,j) = P_ze(i,j,1,2).*ones(params.n_k, params.n_kbar);
%             P_21(:,:,i,j) = P_ze(i,j,2,1).*ones(params.n_k, params.n_kbar);
%             P_22(:,:,i,j) = P_ze(i,j,2,2).*ones(params.n_k, params.n_kbar);
            P_11(:,:,i,j) = P_z(i,1).*P_ze(i,j,1,1).*ones(params.n_k, params.n_kbar);
            P_12(:,:,i,j) = P_z(i,1).*P_ze(i,j,1,2).*ones(params.n_k, params.n_kbar);
            P_21(:,:,i,j) = P_z(i,2).*P_ze(i,j,2,1).*ones(params.n_k, params.n_kbar);
            P_22(:,:,i,j) = P_z(i,2).*P_ze(i,j,2,2).*ones(params.n_k, params.n_kbar);
        end
    end

    %& Convert P_z into 3D format for complete markets PFI
    P_1 = nan(params.n_k, params.n_kbar, params.n_agg);
    P_2 = nan(params.n_k, params.n_kbar, params.n_agg);

    for i=1:params.n_agg
        P_1(:,:,i) = P_z(i,1).*ones(params.n_k, params.n_kbar);
        P_2(:,:,i) = P_z(i,2).*ones(params.n_k, params.n_kbar);
    end
end

