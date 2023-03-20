function yt = autoreg(alph, sigma ,rho, n)
xt = normrnd(0,sigma);
for i=2:n
    xt(i,1) = alph + rho*xt(i-1,1) + normrnd(0,sigma);
end
yt=xt;
end

