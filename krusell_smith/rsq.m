function r2 = rsq(y, X, b)
%Compute R-sq from linear regression
    ybar = mean(y);
    ssr = sum((y-X*b).^2);
    sst = sum((y-ybar).^2);
    r2 = 1 - ssr/sst;
end

