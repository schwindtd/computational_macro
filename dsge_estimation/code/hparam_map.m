%% Hyperparameter Mapping
% Maps hyperparameters from An & Schorfheide (2007) to hyperparameters
% Matlab can use in built-in pdf functions
options = optimset('Display','off');
%% Load hyperparameters from An & Schorfheide (2007)
hyper_as = xlsread('../../data/prior_hyperparams.xlsx');
%hyper_as = xlsread('../../data/prior_hyperparams_test.xlsx');
%% Gamma Distribution parameters
f_gamma = @(x, c) [x(1)*x(2) - c(1) x(1)*x(2)^2 - c(2)^2];
dim_hyper = size(hyper_as);
hyper_mat = nan(13,2);
for i=1:4
    f = @(x) f_gamma(x, hyper_as(i,:));
    hyper_mat(i, :) = fsolve(f, [0.1, 0.1], options);
end
for i=8:9
    f = @(x) f_gamma(x, hyper_as(i,:));
    hyper_mat(i, :) = fsolve(f, [0.1, 0.1], options);
end

%% Beta Distribution parameters
f_beta = @(x,c) [x(1)/(x(1) + x(2)) - c(1) (x(1)*x(2))/((x(1) + x(2))^2*(x(1) + x(2) + 1)) - c(2)^2];
for i=5:7
    f = @(x) f_beta(x, hyper_as(i,:));
    hyper_mat(i, :) = fsolve(f, [0.1, 0.1], options);
end

%% Normal Distribution parameters
hyper_mat(10,:) = hyper_as(10,:);

%% Inverse gamma distribution parameters
hyper_mat(11:13,:) = hyper_as(11:13,:); 
% Keep same as we don't need to transform these like we would if we were using dynare
