function [y] = foc(k, vect1, vect2, matr, a, pow, n, m, disc, dep)
%First Order Condition Function
%   Takes a k (captl tomorrow) and the grids/parameters from stochastic
%   growth model and computes the value of the FOC
    eval = interp1(vect2, matr(:,m), k, 'linear', 'extrap');
    y = (1-disc)/(a*exp(vect1(m))*vect2(n)^pow + (1-dep)*vect2(n) - k) - ...
       disc*eval;
end