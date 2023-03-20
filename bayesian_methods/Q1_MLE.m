%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem Set #4
% Part 1: Monte Carlo with Maximum Likelihood
% ECON630 FALL 2022
% Authors: Daniel Schwindt, Giuliano Simoncelli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We are interested in the linear model  y_t = alpha + beta*x_t + u_t   %
% We let alpha = 1, beta = 2                                            %
% We consider three alternatives for u_t:                               %
% 1) N(0,1)                                                             %
% 2) t distribution with 1 degree of freedom                            %
% 3) Laplace distribution with mu=0 and b=1                             %
%                                                                       %
% We run 1000 replications and let the sample size be both 20 and 100   %
% We will estimate alpha_MLE and beta_MLE                               %
%                                                                       %
% We will run a total of 18 estimations:                                %
% Three per true error distribution for T = 20                          %
% Three per true error distribution for T = 100                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 0.- Housekeeping

rng(123456);

aalpha    = 1;            % Intercept of linear model
bbeta     = 2;            % Slope coefficient of interest
N         = 1000;         % Number of replications
T         = [20,100];     % Sample size
nnu       = 1;            % Degrees of freedom of t-distribution
mmu       = 0;            % Mean of Laplace distribution
b         = 1;            % Scale parameter of Laplace distribution


%% 1.- Preallocating cell arrays
x_var   = {rand(T(1),N), rand(T(2),N)};                                         % Generate matrix for x
u1      = {randn(T(1),N), randn(T(2),N)};                                       % Generate alternative 1 for errors
u2      = {trnd(nnu,T(1),N), trnd(nnu,T(2),N)};                                 % Generate alternative 2 for errors
u3      = {laplace(T(1),N,mmu,b), laplace(T(2),N,mmu,b)};                       % Generate alternative 3 for errors

y1_var  = {aalpha + bbeta*x_var{1} + u1{1}, aalpha + bbeta*x_var{2} + u1{2}};   % Use alternative 1 to determine y_t
y2_var  = {aalpha + bbeta*x_var{1} + u2{1}, aalpha + bbeta*x_var{2} + u2{2}};   % Use alternative 2 to determine y_t
y3_var  = {aalpha + bbeta*x_var{1} + u3{1}, aalpha + bbeta*x_var{2} + u3{2}};   % Use alternative 3 to determine y_t

tic
%% 2.- Computing MLE estimates when true distribution is standard normal

% Preallocating matrices and cell arrays
sol_T20_normal    = {zeros(N,length(T)), zeros(N,length(T)), zeros(N,length(T))};
sol_T100_normal   = {zeros(N,length(T)), zeros(N,length(T)), zeros(N,length(T))};
numerator         = {zeros(T(1),1), zeros(T(1),1), zeros(T(1),1)  ;  zeros(T(2),1), zeros(T(2),1), zeros(T(2),1)};    % Array is 2x3
denominator       = {zeros(T(1),1), zeros(T(1),1), zeros(T(1),1)  ;  zeros(T(2),1), zeros(T(2),1), zeros(T(2),1)};    % Array is 2x3

% Main loop
for n = 1:N   % For each simulation

    for t1 = 1:T(1)  % For the case when sample = 20
        numerator{1,1}(t1,1)   = (x_var{1}(t1,n)- mean(x_var{1}(:,n)))*(y1_var{1}(t1,n)- mean(y1_var{1}(:,n)));
        numerator{1,2}(t1,1)   = (x_var{1}(t1,n)- mean(x_var{1}(:,n)))*(y2_var{1}(t1,n)- mean(y2_var{1}(:,n)));
        numerator{1,3}(t1,1)   = (x_var{1}(t1,n)- mean(x_var{1}(:,n)))*(y3_var{1}(t1,n)- mean(y3_var{1}(:,n)));
        
        denominator{1,1}(t1,1) = (x_var{1}(t1,n)- mean(x_var{1}(:,n)))^2;
        denominator{1,2}(t1,1) = (x_var{1}(t1,n)- mean(x_var{1}(:,n)))^2;
        denominator{1,3}(t1,1) = (x_var{1}(t1,n)- mean(x_var{1}(:,n)))^2;
    end

    for t2 = 1:T(2)  % For the case when sample = 100
        numerator{2,1}(t2,1)   = (x_var{2}(t2,n)- mean(x_var{2}(:,n)))*(y1_var{2}(t2,n)- mean(y1_var{2}(:,n)));
        numerator{2,2}(t2,1)   = (x_var{2}(t2,n)- mean(x_var{2}(:,n)))*(y2_var{2}(t2,n)- mean(y2_var{2}(:,n)));
        numerator{2,3}(t2,1)   = (x_var{2}(t2,n)- mean(x_var{2}(:,n)))*(y3_var{2}(t2,n)- mean(y3_var{2}(:,n)));

        denominator{2,1}(t2,1) = (x_var{2}(t2,n)- mean(x_var{2}(:,n)))^2;
        denominator{2,2}(t2,1) = (x_var{2}(t2,n)- mean(x_var{2}(:,n)))^2;
        denominator{2,3}(t2,1) = (x_var{2}(t2,n)- mean(x_var{2}(:,n)))^2;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For consistency we store every estimate in the same way        %
    %                                                                %
    % It is a 1x3 cell array, where:                                 %
    % The first matrix uses normal errors                            %
    % The second matrix uses t-distributed errors                    %
    % The third matrix uses Laplacian errors                         %
    %                                                                %
    % Inside every matrix:                                           %
    % The first column contains alpha_MLE for every replication      %
    % The second column contains beta_MLE for every replication      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % We are computing alpha and beta that arise from minimizing the
    % log-likelihood function (sum of squared residuals). Same as OLS when
    % (true) errors are N(0,1)
    sol_T20_normal{1}(n,:)      = [mean(y1_var{1}(:,n)) - sum(numerator{1,1})/sum(denominator{1,1})*mean(x_var{1}(:,n)) , sum(numerator{1,1})/sum(denominator{1,1})];
    sol_T20_normal{2}(n,:)      = [mean(y2_var{1}(:,n)) - sum(numerator{1,2})/sum(denominator{1,2})*mean(x_var{1}(:,n)) , sum(numerator{1,2})/sum(denominator{1,2})];
    sol_T20_normal{3}(n,:)      = [mean(y3_var{1}(:,n)) - sum(numerator{1,3})/sum(denominator{1,3})*mean(x_var{1}(:,n)) , sum(numerator{1,3})/sum(denominator{1,3})];

    sol_T100_normal{1}(n,:)     = [mean(y1_var{2}(:,n)) - sum(numerator{2,1})/sum(denominator{2,1})*mean(x_var{2}(:,n)) , sum(numerator{2,1})/sum(denominator{2,1})];
    sol_T100_normal{2}(n,:)     = [mean(y2_var{2}(:,n)) - sum(numerator{2,2})/sum(denominator{2,2})*mean(x_var{2}(:,n)) , sum(numerator{2,2})/sum(denominator{2,2})];
    sol_T100_normal{3}(n,:)     = [mean(y3_var{2}(:,n)) - sum(numerator{2,3})/sum(denominator{2,3})*mean(x_var{2}(:,n)) , sum(numerator{2,3})/sum(denominator{2,3})];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_normal = toc;

%% 3.- Computing MLE estimates when error has t distribution

% When the errors are t-distributed we need to estimate with a different
% method. See, for example 
% https://stats.stackexchange.com/questions/350550/mle-for-linear-regression-student-t-distributed-error

% Initialize solution arrays
sol_T20_tdist   = {zeros(N,length(T)), zeros(N,length(T)), zeros(N,length(T))};
sol_T100_tdist  = {zeros(N,length(T)), zeros(N,length(T)), zeros(N,length(T))};

% Main loop
% We will be solving the two score functions with fsolve
% We need to build the score functions
for n = 1:N

    % Initialize cell array of anonymous functions
    % We must keep these inside the loop so that I start computing the sums
    % from scratch for every different replication
    eqn_T20_tdist   = {{@(bet) 0, @(bet) 0},{@(bet) 0, @(bet) 0},{@(bet) 0, @(bet) 0}};
    eqn_T100_tdist  = {{@(bet) 0, @(bet) 0},{@(bet) 0, @(bet) 0},{@(bet) 0, @(bet) 0}};

   
    eqn_T20_tdist{1}  = {@(bet) sum( 2 * (y1_var{1}(:,n) - bet(1) - bet(2)*x_var{1}(:,n)) ./ ...
                            (nnu + (y1_var{1}(:,n) - bet(1) - bet(2)*x_var{1}(:,n)).^2) ) , ... 
                         @(bet) sum( 2 * x_var{1}(:,n) .* (y1_var{1}(:,n) - bet(1) - bet(2)*x_var{1}(:,n)) ./ ... 
                            (nnu + (y1_var{1}(:,n) - bet(1) - bet(2)*x_var{1}(:,n)).^2) ) };

    eqn_T20_tdist{2}  = {@(bet) sum( 2 * (y2_var{1}(:,n) - bet(1) - bet(2)*x_var{1}(:,n)) ./ ...
                            (nnu + (y2_var{1}(:,n) - bet(1) - bet(2)*x_var{1}(:,n)).^2) ) , ... 
                         @(bet) sum( 2 * x_var{1}(:,n) .* (y2_var{1}(:,n) - bet(1) - bet(2)*x_var{1}(:,n)) ./ ... 
                            (nnu + (y2_var{1}(:,n) - bet(1) - bet(2)*x_var{1}(:,n)).^2) ) };

    eqn_T20_tdist{3}  = {@(bet) sum( 2 * (y3_var{1}(:,n) - bet(1) - bet(2)*x_var{1}(:,n)) ./ ...
                            (nnu + (y3_var{1}(:,n) - bet(1) - bet(2)*x_var{1}(:,n)).^2) ) , ... 
                         @(bet) sum( 2 * x_var{1}(:,n) .* (y3_var{1}(:,n) - bet(1) - bet(2)*x_var{1}(:,n)) ./ ... 
                            (nnu + (y3_var{1}(:,n) - bet(1) - bet(2)*x_var{1}(:,n)).^2) ) };
        
    eqn_T100_tdist{1} = {@(bet) sum(  2 * (y1_var{2}(:,n) - bet(1) - bet(2)*x_var{2}(:,n)) ./ ...
                            (nnu + (y1_var{2}(:,n) - bet(1) - bet(2)*x_var{2}(:,n)).^2) ) , ...
                         @(bet) sum(  2 * x_var{2}(:,n) .* (y1_var{2}(:,n) - bet(1) - bet(2)*x_var{2}(:,n)) ./ ... 
                            (nnu + (y1_var{2}(:,n) - bet(1) - bet(2)*x_var{2}(:,n)).^2) ) };

    eqn_T100_tdist{2} = {@(bet) sum( 2 * (y2_var{2}(:,n) - bet(1) - bet(2)*x_var{2}(:,n)) ./ ...
                            (nnu + (y2_var{2}(:,n) - bet(1) - bet(2)*x_var{2}(:,n)).^2) ) , ...
                         @(bet) sum( 2 * x_var{2}(:,n) .* (y2_var{2}(:,n) - bet(1) - bet(2)*x_var{2}(:,n)) ./ ... 
                            (nnu + (y2_var{2}(:,n) - bet(1) - bet(2)*x_var{2}(:,n)).^2) ) };

    eqn_T100_tdist{3} = {@(bet) sum( 2 * (y3_var{2}(:,n) - bet(1) - bet(2)*x_var{2}(:,n)) ./ ...
                            (nnu + (y3_var{2}(:,n) - bet(1) - bet(2)*x_var{2}(:,n)).^2) ) , ...
                         @(bet) sum( 2 * x_var{2}(:,n) .* (y3_var{2}(:,n) - bet(1) - bet(2)*x_var{2}(:,n)) ./ ... 
                            (nnu + (y3_var{2}(:,n) - bet(1) - bet(2)*x_var{2}(:,n)).^2) ) };
    

    % Turn display off for fsolve
    options = optimset('Display','off','Algorithm','levenberg-marquardt');


    % Compute alpha_MLE and beta_MLE (2 equations, 2 unknowns)
    sol_T20_tdist{1}(n,:)   = fsolve(eqn_T20_tdist{1}, [1,2],options);   % Nx2 matrix where the first column is alpha_MLE and second is beta_MLE, for T=20
    clc; % This for now, I don't know how to deal with the Jacobian warning
    sol_T20_tdist{2}(n,:)   = fsolve(eqn_T20_tdist{2}, [1,2],options);   % Nx2 matrix where the first column is alpha_MLE and second is beta_MLE, for T=20
    clc; % This for now, I don't know how to deal with the Jacobian warning
    sol_T20_tdist{3}(n,:)   = fsolve(eqn_T20_tdist{3}, [1,2],options);   % Nx2 matrix where the first column is alpha_MLE and second is beta_MLE, for T=20
    clc; % This for now, I don't know how to deal with the Jacobian warning  
    sol_T100_tdist{1}(n,:)  = fsolve(eqn_T100_tdist{1}, [1,2],options);  % Nx2 matrix where the first column is alpha_MLE and second is beta_MLE, for T=100
    clc; % This for now, I don't know how to deal with the Jacobian warning
    sol_T100_tdist{2}(n,:)  = fsolve(eqn_T100_tdist{2}, [1,2],options);  % Nx2 matrix where the first column is alpha_MLE and second is beta_MLE, for T=100
    clc; % This for now, I don't know how to deal with the Jacobian warning
    sol_T100_tdist{3}(n,:)  = fsolve(eqn_T100_tdist{3}, [1,2],options);  % Nx2 matrix where the first column is alpha_MLE and second is beta_MLE, for T=100
    clc; % This for now, I don't know how to deal with the Jacobian warning

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_tdist = toc - time_normal;

%% 4.- Computing MLE estimates when error has Laplace distribution

% Modified version of 
% https://stats.stackexchange.com/questions/174881/linear-regression-with-laplace-errors
% (Just added the intercept)

% Initialize solution cell array
sol_T20_laplace   = {zeros(N,length(T)),zeros(N,length(T)),zeros(N,length(T))};
sol_T100_laplace  = {zeros(N,length(T)),zeros(N,length(T)),zeros(N,length(T))};

% Main loop
% We minimize the log-likelihood function
for n = 1:N

    % Turn display off
    options = optimset('Display','off');
    % Compute alpha_MLE and beta_MLE
    sol_T20_laplace{1}(n,:)   = fminunc(@(bet) 1/b*sum(abs(y1_var{1}(:,n)-bet(1)-bet(2)*x_var{1}(:,n))), [1,2], options);   % Nx2 matrix where the first column is alpha_MLE and second is beta_MLE, for T=20
    sol_T100_laplace{1}(n,:)  = fminunc(@(bet) 1/b*sum(abs(y1_var{2}(:,n)-bet(1)-bet(2)*x_var{2}(:,n))), [1,2], options);   % Nx2 matrix where the first column is alpha_MLE and second is beta_MLE, for T=100
    sol_T20_laplace{2}(n,:)   = fminunc(@(bet) 1/b*sum(abs(y2_var{1}(:,n)-bet(1)-bet(2)*x_var{1}(:,n))), [1,2], options);   % Nx2 matrix where the first column is alpha_MLE and second is beta_MLE, for T=20
    sol_T100_laplace{2}(n,:)  = fminunc(@(bet) 1/b*sum(abs(y2_var{2}(:,n)-bet(1)-bet(2)*x_var{2}(:,n))), [1,2], options);   % Nx2 matrix where the first column is alpha_MLE and second is beta_MLE, for T=100
    sol_T20_laplace{3}(n,:)   = fminunc(@(bet) 1/b*sum(abs(y3_var{1}(:,n)-bet(1)-bet(2)*x_var{1}(:,n))), [1,2], options);   % Nx2 matrix where the first column is alpha_MLE and second is beta_MLE, for T=20
    sol_T100_laplace{3}(n,:)  = fminunc(@(bet) 1/b*sum(abs(y3_var{2}(:,n)-bet(1)-bet(2)*x_var{2}(:,n))), [1,2], options);   % Nx2 matrix where the first column is alpha_MLE and second is beta_MLE, for T=100
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_laplace = toc-time_tdist;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% END OF ESTIMATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% 5.- Display information

disp('Finished computing maximum likelihood estimators. \n');
fprintf('Time running: Normal errors = %2.2f seconds, t-distribution errors = %2.2f seconds, Laplace errors = %2.2f seconds \n \n',time_normal,time_tdist,time_laplace)
disp('Preparing plots...');

%% 6.- Histograms for alpha_MLE
%figure(1)
%subplot(3,1,1)
%histogram(sol_T100_normal{1}(:,1))
%subplot(3,1,2)
%histogram(sol_T100_normal{2}(:,1))
%subplot(3,1,3)
%histogram(sol_T100_normal{3}(:,1))

[N1,edges1] = histcounts(sol_T100_normal{1}(:,1), 'Normalization','pdf');
edges1 = edges1(2:end) - (edges1(2)-edges1(1))/2;
[N2,edges2] = histcounts(sol_T100_normal{2}(:,1),1500, 'Normalization','pdf');
edges2 = edges2(2:end) - (edges2(2)-edges2(1))/2;
[N3,edges3] = histcounts(sol_T100_normal{3}(:,1), 'Normalization','pdf');
edges3 = edges3(2:end) - (edges3(2)-edges3(1))/2;

figure(1)
plot(edges1, N1,'blue--');
hold on
plot(edges2, N2,'Color','red');
plot(edges3, N3,'Color','green');
hold off
legend('N(0,1)','t-distribution (k = 1)','Laplace distribution (0,1)')
xlim([-1,3])
xlabel('Dotted line is the correct error specification')

%figure(2)
%subplot(3,1,1)
%histogram(sol_T100_tdist{1}(:,1))
%subplot(3,1,2)
%histogram(sol_T100_tdist{2}(:,1))
%subplot(3,1,3)
%histogram(sol_T100_tdist{3}(:,1))

[N1,edges1] = histcounts(sol_T100_tdist{1}(:,1), 'Normalization','pdf');
edges1 = edges1(2:end) - (edges1(2)-edges1(1))/2;
[N2,edges2] = histcounts(sol_T100_tdist{2}(:,1), 'Normalization','pdf');
edges2 = edges2(2:end) - (edges2(2)-edges2(1))/2;
[N3,edges3] = histcounts(sol_T100_tdist{3}(:,1), 'Normalization','pdf');
edges3 = edges3(2:end) - (edges3(2)-edges3(1))/2;

figure(2)
plot(edges1, N1,'Color','blue');
hold on
plot(edges2, N2,'red--');
plot(edges3, N3,'Color','green');
hold off
legend('N(0,1)','t-distribution (k = 1)','Laplace distribution (0,1)')
xlim([0.5,1.5])
xlabel('Dotted line is the correct error specification')

%figure(3)
%subplot(3,1,1)
%histogram(sol_T100_laplace{1}(:,1))
%subplot(3,1,2)
%histogram(sol_T100_laplace{2}(:,1))
%subplot(3,1,3)
%histogram(sol_T100_laplace{3}(:,1))

[N1,edges1] = histcounts(sol_T100_laplace{1}(:,1), 'Normalization','pdf');
edges1 = edges1(2:end) - (edges1(2)-edges1(1))/2;
[N2,edges2] = histcounts(sol_T100_laplace{2}(:,1), 'Normalization','pdf');
edges2 = edges2(2:end) - (edges2(2)-edges2(1))/2;
[N3,edges3] = histcounts(sol_T100_laplace{3}(:,1), 'Normalization','pdf');
edges3 = edges3(2:end) - (edges3(2)-edges3(1))/2;

figure(3)
plot(edges1, N1,'Color','blue');
hold on
plot(edges2, N2,'Color','red');
plot(edges3, N3,'green--');
hold off
legend('N(0,1)','t-distribution (k = 1)','Laplace distribution (0,1)')
xlim([-0.2,2.2])
xlabel('Dotted line is the correct error specification')

%% 7.- Histograms for beta_MLE
%figure(4)
%subplot(3,1,1)
%histogram(sol_T100_normal{1}(:,2))
%subplot(3,1,2)
%histogram(sol_T100_normal{2}(:,2))
%subplot(3,1,3)
%histogram(sol_T100_normal{3}(:,2))

[N1,edges1] = histcounts(sol_T100_normal{1}(:,2), 'Normalization','pdf');
edges1 = edges1(2:end) - (edges1(2)-edges1(1))/2;
[N2,edges2] = histcounts(sol_T100_normal{2}(:,2),1500, 'Normalization','pdf');
edges2 = edges2(2:end) - (edges2(2)-edges2(1))/2;
[N3,edges3] = histcounts(sol_T100_normal{3}(:,2), 'Normalization','pdf');
edges3 = edges3(2:end) - (edges3(2)-edges3(1))/2;

figure(4)
plot(edges1, N1,'blue--');
hold on
plot(edges2, N2,'Color','red');
plot(edges3, N3,'Color','green');
hold off
legend('N(0,1)','t-distribution (k = 1)','Laplace distribution (0,1)')
xlim([0.5,3.5])
xlabel('Dotted line is the correct error specification')

%figure(5)
%subplot(3,1,1)
%histogram(sol_T100_tdist{1}(:,2))
%subplot(3,1,2)
%histogram(sol_T100_tdist{2}(:,2))
%subplot(3,1,3)
%histogram(sol_T100_tdist{3}(:,2))

[N1,edges1] = histcounts(sol_T100_tdist{1}(:,2), 'Normalization','pdf');
edges1 = edges1(2:end) - (edges1(2)-edges1(1))/2;
[N2,edges2] = histcounts(sol_T100_tdist{2}(:,2), 'Normalization','pdf');
edges2 = edges2(2:end) - (edges2(2)-edges2(1))/2;
[N3,edges3] = histcounts(sol_T100_tdist{3}(:,2), 'Normalization','pdf');
edges3 = edges3(2:end) - (edges3(2)-edges3(1))/2;

figure(5)
plot(edges1, N1,'Color','blue');
hold on
plot(edges2, N2,'red--');
plot(edges3, N3,'Color','green');
hold off
legend('N(0,1)','t-distribution (k = 1)','Laplace distribution (0,1)')
xlim([1.8,2.2])
xlabel('Dotted line is the correct error specification')

%figure(6)
%subplot(3,1,1)
%histogram(sol_T100_laplace{1}(:,2))
%subplot(3,1,2)
%histogram(sol_T100_laplace{2}(:,2))
%subplot(3,1,3)
%histogram(sol_T100_laplace{3}(:,2))

[N1,edges1] = histcounts(sol_T100_laplace{1}(:,2), 'Normalization','pdf');
edges1 = edges1(2:end) - (edges1(2)-edges1(1))/2;
[N2,edges2] = histcounts(sol_T100_laplace{2}(:,2), 'Normalization','pdf');
edges2 = edges2(2:end) - (edges2(2)-edges2(1))/2;
[N3,edges3] = histcounts(sol_T100_laplace{3}(:,2), 'Normalization','pdf');
edges3 = edges3(2:end) - (edges3(2)-edges3(1))/2;

figure(6)
plot(edges1, N1,'Color','blue');
hold on
plot(edges2, N2,'Color','red');
plot(edges3, N3,'green--');
hold off
legend('N(0,1)','t-distribution (k = 1)','Laplace distribution (0,1)')
xlim([0.5,3.5])
xlabel('Dotted line is the correct error specification')

%% 8.- Storing estimates in a more efficient way

% Building a 1x9 cell array for the results
sol_T20 = {sol_T20_normal{1},sol_T20_normal{2},sol_T20_normal{3},sol_T20_tdist{1},sol_T20_tdist{2},sol_T20_tdist{3},sol_T20_laplace{1},sol_T20_laplace{2},sol_T20_laplace{3}};
sol_T100 = {sol_T100_normal{1},sol_T100_normal{2},sol_T100_normal{3},sol_T100_tdist{1},sol_T100_tdist{2},sol_T100_tdist{3},sol_T100_laplace{1},sol_T100_laplace{2},sol_T100_laplace{3}};

% T=100
% Summary Table for Alpha
s_alpha100 = zeros(2, 9);
for j=1:9
    s_alpha100(1,j) = mean(sol_T100{j}(:,1));
    s_alpha100(2,j) = std(sol_T100{j}(:,1));
end
% Summary Table for Beta
s_beta100 = zeros(2, 9);
for j=1:9
    s_beta100(1,j) = mean(sol_T100{j}(:,2));
    s_beta100(2,j) = std(sol_T100{j}(:,2));
end
% T=20
% Summary Table for Alpha
s_alpha20 = zeros(2, 9);
for j=1:9
    s_alpha20(1,j) = mean(sol_T20{j}(:,1));
    s_alpha20(2,j) = std(sol_T20{j}(:,1));
end
% Summary Table for Beta
s_beta20 = zeros(2, 9);
for j=1:9
    s_beta20(1,j) = mean(sol_T20{j}(:,2));
    s_beta20(2,j) = std(sol_T20{j}(:,2));
end

%% 9.- Boxplots

% Normal true - beta
figure(7)
boxplot([sol_T100_normal{1}(:,2),sol_T100_normal{2}(:,2),sol_T100_normal{3}(:,2)])
ylim([-2,6])
title('Estimation of beta when true distribution is N(0,1)')
ax = gca;
ax.XTick = [0,1,2,3,4];
ax.XTickLabel = {' ','N(0,1) errors','t errors','Laplacian errors',''};
line = [2,2,2,2,2];
hold on
plot(ax.XTick,line,'Color','Black')
hold off
xlim([0.5,3.5])

% t true - beta
figure(8)
boxplot([sol_T100_tdist{1}(:,2),sol_T100_tdist{2}(:,2),sol_T100_tdist{3}(:,2)])
ylim([1.7,2.3])
title('Estimation of beta when true distribution is t')
ax = gca;
ax.XTick = [0,1,2,3,4];
ax.XTickLabel = {' ','N(0,1) errors','t errors','Laplacian errors',''};
line = [2,2,2,2,2];
hold on
plot(ax.XTick,line,'Color','Black')
hold off
xlim([0.5,3.5])

% Laplace true - beta
figure(9)
boxplot([sol_T100_laplace{1}(:,2),sol_T100_laplace{2}(:,2),sol_T100_laplace{3}(:,2)])
ylim([0,4])
title('Estimation of beta when true distribution is Laplace')
ax = gca;
ax.XTick = [0,1,2,3,4];
ax.XTickLabel = {' ','N(0,1) errors','t errors','Laplacian errors',''};
line = [2,2,2,2,2];
hold on
plot(ax.XTick,line,'Color','Black')
hold off
xlim([0.5,3.5])

% Normal true - alpha
figure(10)
boxplot([sol_T100_normal{1}(:,1),sol_T100_normal{2}(:,1),sol_T100_normal{3}(:,1)])
ylim([-1,3])
title('Estimation of alpha when true distribution is N(0,1)')
ax = gca;
ax.XTick = [0,1,2,3,4];
ax.XTickLabel = {' ','N(0,1) errors','t errors','Laplacian errors',''};
line = [1,1,1,1,1];
hold on
plot(ax.XTick,line,'Color','Black')
hold off
xlim([0.5,3.5])

% t true - alpha
figure(11)
boxplot([sol_T100_tdist{1}(:,1),sol_T100_tdist{2}(:,1),sol_T100_tdist{3}(:,1)])
ylim([0.5,1.5])
title('Estimation of alpha when true distribution is t')
ax = gca;
ax.XTick = [0,1,2,3,4];
ax.XTickLabel = {' ','N(0,1) errors','t errors','Laplacian errors',''};
line = [1,1,1,1,1];
hold on
plot(ax.XTick,line,'Color','Black')
hold off
xlim([0.5,3.5])

% Laplace true - alpha
figure(12)
boxplot([sol_T100_laplace{1}(:,1),sol_T100_laplace{2}(:,1),sol_T100_laplace{3}(:,1)])
ylim([0,2])
title('Estimation of alpha when true distribution is Laplace')
ax = gca;
ax.XTick = [0,1,2,3,4];
ax.XTickLabel = {' ','N(0,1) errors','t errors','Laplacian errors',''};
line = [1,1,1,1,1];
hold on
plot(ax.XTick,line,'Color','Black')
hold off
xlim([0.5,3.5])

disp('END');