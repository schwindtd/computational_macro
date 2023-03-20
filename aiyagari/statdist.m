function [mu, iter, diff] = statdist(mu0, pi, g, tol, max_iter)
%Compute stationary distribution for Aiyagari
%   Use intitial distribution guess, decision rules, and markov chain
%   transition matrix to compute stationary distribution of assets
diff    = 100;
iter    = 0;
dim     = size(g);
n_asset = dim(1);
n_eff   = dim(2);

while (diff > tol) && (iter < max_iter)
    % Loop over endogenous state, exogenous state, and next period
    % exogenous state
    mu_new = zeros(n_asset, n_eff);
        for j=1:n_eff
          for i=1:n_asset
            for j_nxt = 1:n_eff
                mu_new(g(i,j),j_nxt) = mu_new(g(i,j),j_nxt) + ...
                    pi(j,j_nxt)*mu0(i,j);
            end
        end
    end

    diff = max(max(abs(mu_new - mu0)));
    mu0 = mu_new;
    iter = iter + 1;

end
    mu = mu0;
end