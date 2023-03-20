%% Pack input

% grids and matrixes
mgm_input.c_grid_size = n_asset_c;
mgm_input.f_grid_size = n_asset;
mgm_input.a_upper = a_upper;
mgm_input.a_lower = a_lower;
mgm_input.v_eff = v_eff;
mgm_input.m_trans = m_trans;
mgm_input.util_ss = util_ss;

% Problem parameters
mgm_input.bbeta = bbeta;
mgm_input.ssigma = ssigma;
mgm_input.ddelta = ddelta;
mgm_input.n_eff = n_eff;
mgm_input.r = r;
mgm_input.w = w;

% Iteration parameters
mgm_input.tol = tol;
mgm_input.max_iter_mu = max_iter_mu;
mgm_input.max_iter_A = max_iter_A;