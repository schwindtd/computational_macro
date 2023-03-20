function [lp, prior_pdf] = prior(v_params, dist_type, scal, hyper_mat)
%Function to compute joint prior density
%   INPUTS:     1. Vector of parameter values (double)
%               2. Cell array of distribution names (string): normal, gamma,
%               invgamma, t, beta
%               3. Scaling factors for each parameter
%               4. Matrix of distribution hyperparameters (double)
%   OUTPUTS:    1. log joint prior density
%               2. prior density for each individual parameter
%% 5. Prior Density
prior_pdf = nan(length(v_params),1);
% idxg = [1 2 3 4 8 9];
% idxb = [5 6 7];
% idxn = [10];
% idxivg = [11 12 13];
% Set numeric codes for distribution types
% Normal = 1, Gamma = 2, Student t = 3, Beta = 4, invgamma = 5
dist_type_num = nan(length(v_params),1);
dist_type_num(strcmp(dist_type,'normal')) = 1;
dist_type_num(strcmp(dist_type,'gamma')) = 2;
dist_type_num(strcmp(dist_type,'t')) = 3;
dist_type_num(strcmp(dist_type,'beta')) = 4;
dist_type_num(strcmp(dist_type,'invgamma')) = 5;
dist_type_num(strcmp(dist_type, 'uniform')) = 6;

% Define inverse gamma PDF
ivgpdf = @(y,s,v) ((2*((v*s^2)/2)^(v/2))/gamma(v/2))*100.*y.^(-1*(v+1)).*exp((-1*v*s^2)./(2*y.^2));
%ivgpdf = @(y,s,v) y.^(-1*(v+1)).*exp((-1*v*s^2)./(2*y.^2));

for i=1:length(v_params)
    if dist_type_num(i) == 1
        prior_pdf(i) = normpdf(scal(i)*v_params(i), hyper_mat(i,1), hyper_mat(i,2));
    elseif dist_type_num(i) == 2
        prior_pdf(i) = gampdf(scal(i)*v_params(i), hyper_mat(i,1), hyper_mat(i,2)); 
    elseif dist_type_num(i) == 3
        prior_pdf(i) = tpdf(scal(i)*v_params(i), hyper_mat(i,1));
    elseif dist_type_num(i) == 4
        prior_pdf(i) = betapdf(scal(i)*v_params(i), hyper_mat(i,1), hyper_mat(i,2));
    elseif dist_type_num(i) == 5
        prior_pdf(i) = ivgpdf(scal(i)*v_params(i), hyper_mat(i,1), hyper_mat(i,2));
    else
        prior_pdf(i) = unifpdf(v_params(i), hyper_mat(i,1), hyper_mat(i,2));
    end
end

% Joint prior likelihood
lp = sum(log(prior_pdf));
end

