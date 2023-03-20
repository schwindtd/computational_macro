function y = draw_prior(n,dist,p)
    % Convert dist string to numeric
    dist_num = nan(1,1);
    dist_num(strcmp(dist,'normal')) = 1;
    dist_num(strcmp(dist,'gamma')) = 2;
    dist_num(strcmp(dist,'t')) = 3;
    dist_num(strcmp(dist,'beta')) = 4;
    dist_num(strcmp(dist,'invgamma')) = 5;
    dist_num(strcmp(dist, 'uniform')) = 6;
    
    % If statements to use built-in rand functions
    if dist_num == 1
        y = p(1) + p(2)*randn(n,1);
    elseif dist_num == 2
        y = gamrnd(repmat(p(1),[n,1]),repmat(p(2), [n,1]));
    elseif dist_num == 3
        y = trnd(p(1),1,n);
    elseif dist_num == 4
        y = betarnd(repmat(p(1),[n,1]),repmat(p(2), [n,1]));
    elseif dist_num == 5
        x = chi2rnd(p(2), n,1);
        z = p(1)^2*p(2)*(1./x);
        y = sqrt(z)./100;
    elseif dist_num == 6
        y = p(1) + (p(2) - p(1))*rand(n,1);
    else
        error('Specified unsupported distribution type.');
    end
end

