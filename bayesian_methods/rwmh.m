function [draws, acc_rate] = rwmh(init,N, tgt_mu, tgt_sigma, pr_mu, pr_sigma, scale)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
accept = 0;
reject = 0;
x = [init; zeros(N,2)];
eta = mvnrnd(pr_mu,scale*pr_sigma,N);
for i=2:N+1
    theta = x(i-1,:) + eta(i-1, :);
    p_c = mvnpdf(theta, tgt_mu, tgt_sigma);
    p_b = mvnpdf(x(i-1,:), tgt_mu, tgt_sigma);
%       p_c = mvnpdf(theta, x(i-1,:), c*sigma_eta);
%       p_b = mvnpdf(x(i-1,:), x(i-1,:), c*sigma_eta);
    p = min(1, p_c/p_b);
    if p == 1
        x(i,:) = theta;
        accept = accept + 1;
    else
        u = unifrnd(0,1);
        if p < u
            x(i,:) = x(i-1,:);
            reject = reject + 1;
        else
            x(i,:) = theta;
            accept = accept + 1;
        end
    end
    rej_rate = reject/(i-1);
    acc_rate = accept/(i-1);
end
draws = x;

end

