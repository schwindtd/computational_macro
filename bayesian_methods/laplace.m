function y  = laplace(m, n, mu, b)

%   Generate i.i.d. laplacian random number drawn from laplacian distribution
%   with mean mu and standard deviation sigma. 
%   mu      : mean
%   sigma   : standard deviation
%   [m, n]  : the dimension of y.

%Check inputs
if nargin < 4
    error('All four inputs are required');
end

% Variance of Laplace distribution is = 2b^2
u = rand(m, n)-0.5;    % Bound uniform between (-0.5 , 0.5)
y = mu - b * sign(u).* log(1- 2* abs(u));   % Random number (matrix) generated
