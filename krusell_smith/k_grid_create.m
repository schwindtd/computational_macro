function [v_k, v_kbar, m_k,m_kbar] = k_grid_create(params, type)
% Create ND capital and capital bar grids
x = linspace(0,1,params.n_k);
y=x.^params.theta/max(x.^params.theta);
v_k = params.k_min + (params.k_max - params.k_min)*y;
v_kbar = linspace(params.kbar_min,params.kbar_max, params.n_kbar);  % linear per Krusell-Smith
if type == 1
    m_k = nan(params.n_k, params.n_kbar, params.n_agg, params.n_eff); 
    for i=1:params.n_agg
        for j=1:params.n_eff
            m_k(:,:,i,j) = v_k'*ones(1,params.n_kbar);
        end
    end

    m_kbar = nan(params.n_k, params.n_kbar, params.n_agg, params.n_eff);
    for i=1:params.n_agg
        for j=1:params.n_eff
            m_kbar(:,:,i,j) = ones(params.n_k,1)*v_kbar;
        end
    end
else
    m_k = nan(params.n_k, params.n_kbar, params.n_agg); 
    for i=1:params.n_agg
        m_k(:,:,i) = v_k'*ones(1,params.n_kbar);
    end

    m_kbar = nan(params.n_k, params.n_kbar, params.n_agg);
    for i=1:params.n_agg
        m_kbar(:,:,i) = ones(params.n_k,1)*v_kbar;
    end    
end
end

