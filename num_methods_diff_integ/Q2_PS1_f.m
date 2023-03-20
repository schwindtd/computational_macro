%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem Set #1
% Question 2: Numerical Integration
% ECON630 FALL 2022
% Authors: Giuliano Simoncelli, Daniel Schwindt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc

% Set parameters
mu = [0 5];
sigma = [4, -1; -1, 1];
R = chol(sigma);
n = [10 100 1000 10000];
d = [1, 2, 3, 6, 14];

syms c s t

rng(10) % set seed for consistent results

%% Analytical Results
ares = zeros(2,6);
mgf = @(t, c, s)exp(c*t + (s/2)*t.^2); % MGF for generic univariate normal
for i=1:length(mu)
    for j=1:length(d)
        df = matlabFunction(diff(mgf, t, d(j)));
        ares(i,j) = df(mu(i), sigma(i,i),0);
    end
end
% Result for E[x^2*y^3]
ares(1,6) = 590; % hard-coded

%% Method 1: Monte-Carlo Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over to create the rvs
mcme = zeros(4, 11);

for i=1:length(n)
    % Create iid standard normal rvs
    x = mvnrnd(zeros(1,2), eye(2), n(i));
    % Transform x into RV with mean = mu and var-covar = sigma
    y = x*R + mu;
    for j=1:length(d)
        % Compute expectations for x, y components
        s = sum(y.^(d(j)), 1);
        et = (1/n(i))*s;
        if j==1
            mcme(i, j:j+1) = et;
        else
            mcme(i, 2*j-1:2*j) = et;
        end
    end
    % Compute expectation for 5th order correlation-type measure
    x2 = y(:,1).^2;
    y3 = y(:,2).^3;
    ec5t = (1/n(i))*sum(x2.*y3,1);
    mcme(i,length(mcme)) = ec5t;
end


%% Method 2: Gauss-Hermite Quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ghqe = zeros(3,11);

% Read in node/weight data for n=3, 7, 15
n3 = readtable('./ghq_xw_n3.csv');
n7 = readtable('./ghq_xw_n7.csv');
n15 = readtable('./ghq_xw_n15.csv');
n3_2d = readtable('./ghq_2d_n3.csv');
n7_2d = readtable('./ghq_2d_n7.csv');
n15_2d = readtable('./ghq_2d_n15.csv');

quads = [3, 7, 15];
% Convert data into vectors with dynamic naming
for i=1:length(quads)
    eval(['v' num2str(quads(i)) ' = table2array(n' num2str(quads(i)) '(:,"x_i"));' ])
    eval(['w' num2str(quads(i)) ' = table2array(n' num2str(quads(i)) '(:,"w_i"));' ])
    eval(['v' num2str(quads(i)) '_2d = table2array(n' num2str(quads(i)) '_2d(:,2:3));' ])
    eval(['w' num2str(quads(i)) '_2d = table2array(n' num2str(quads(i)) '_2d(:,"w_i"));' ])
    
end

% Compute expectations for each marginal variable using univariate
% quadrature
for i=1:length(quads)
    for j=1:length(d)
        etx = eval(['w' num2str(quads(i)) '.*(sqrt(2).*sqrt(sigma(1,1)).*v' num2str(quads(i)) ' + mu(1)).^d(j);']);
        etsx = (pi)^(-1/2)*sum(etx);
        ety = eval(['w' num2str(quads(i)) '.*(sqrt(2).*sqrt(sigma(2,2)).*v' num2str(quads(i)) ' + mu(2)).^d(j);']);
        etsy = (pi)^(-1/2)*sum(ety);
        ets = [etsx, etsy];
        if j==1
            ghqe(i, j:j+1) = ets;
        else
            ghqe(i, 2*j-1:2*j) = ets;
        end
    end
end
    
% Compute E[x^2*y^3]: quadrature of 3
ghv3 = [ reshape(meshgrid(v3,v3),[9,1]), repmat(v3,[3 1]) ];
ghw3 =  reshape(meshgrid(w3,w3),[9,1]).* repmat(w3,[3,1]); 
Vt3 = sqrt(2)*R*ghv3' + mu';
W3 = ghw3;
et3 = W3.*((pi)^-1.*(Vt3(1,:)'.^2).*(Vt3(2,:)'.^3));
ghqe(1,length(ghqe)) = sum(et3);

% Compute E[x^2*y^3]: quadrature of 7
ghv7 = [ reshape(meshgrid(v7,v7),[49,1]), repmat(v7,[7 1]) ];
ghw7 =  reshape(meshgrid(w7,w7),[49,1]).* repmat(w7,[7,1]); 
Vt7 = sqrt(2)*R*ghv7' + mu';
W7 = ghw7;
et7 = W7.*((pi)^-1.*(Vt7(1,:)'.^2).*(Vt7(2,:)'.^3));
ghqe(2,length(ghqe)) = sum(et7);

% Compute E[x^2*y^3]: quadrature of 15
ghv15 = [ reshape(meshgrid(v15,v15),[225,1]), repmat(v15,[15 1]) ];
ghw15 =  reshape(meshgrid(w15,w15),[225,1]).* repmat(w15,[15,1]); 
Vt15 = sqrt(2)*R*ghv15' + mu';
W15 = ghw15;
et15 = W15.*((pi)^-1.*(Vt15(1,:)'.^2).*(Vt15(2,:)'.^3));
ghqe(3,length(ghqe)) = sum(et15);

%% Method 3: Tauchen (1986)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute expectations for marginals of rvs
% Create point grid
m = 3; % constant
N = 10; % number of grid points
xN = m*R(1,1);
yN = m*1;
x1 = -xN;
y1 = -yN;
xspace = linspace(x1, xN, N);
yspace = linspace(y1, yN, N);
w_x = (xN - x1)/(N-1);
w_y = (yN - y1)/(N-1);

% Compute probabilities
p = zeros(N,2);
for i=1:N
    if i == 1
        pt_x = normcdf((xspace(i) + w_x/2)/R(1,1));
        pt_y = normcdf((yspace(i) + w_y/2)/1);
    elseif i < N
        pt_x = normcdf((xspace(i) + w_x/2)/R(1,1)) - normcdf((xspace(i) - w_x/2)/R(1,1));
        pt_y = normcdf((yspace(i) + w_y/2)/1) - normcdf((yspace(i) - w_y/2)/1);
    else
        pt_x = 1 - normcdf((xspace(i) - w_x/2)/R(1,1));
        pt_y = 1 - normcdf((yspace(i) - w_y/2)/1);
    end
    p(i,:) = [pt_x, pt_y];
end

% Compute expectations
tauche = zeros(1,11);
for i=1:length(d)
    etx = p(:,1).*((xspace + mu(1)).^d(i))';
    ety = p(:,2).*((yspace + mu(2)).^d(i))';
    ets = [sum(etx), sum(ety)];
    if i==1
            tauche(i:i+1) = ets;
    else
            tauche(2*i-1:2*i) = ets;
    end
end

% Computation for E[x^2*y^3]
x = linspace(-m, m, N);
y = linspace(-m, m, N);
[X, Y] = meshgrid(x, y);
Z = [X(:) Y(:)]';
w_x = (2*m/(N-1));
w_y = (2*m/(N-1));

p2 = zeros(N^2,1);
for i=1:size(Z,2)
    if i==1
        ptt_x = normcdf(Z(1,i)+ w_x/2);
        ptt_y = normcdf(Z(2,i)+ w_y/2);
    elseif i < N^2
        ptt_x = normcdf(Z(1,i) + w_x/2) - normcdf(Z(1,i) - w_x/2);
        ptt_y = normcdf(Z(2,i)+ w_y/2) - normcdf(Z(2,i) - w_y/2);
    else
        ptt_x = 1 - normcdf(Z(1,i) - w_x/2);
        ptt_y = 1 - normcdf(Z(2,i) - w_y/2);
    end
    ptt = ptt_x*ptt_y;
    p2(i) = ptt;
end
z = Z'*R + mu;
ett = p2.*(z(:,1)).^2.*(z(:,2)).^3;
tauche(1,length(tauche)) = sum(ett);

%% Plot results
ares_w = [ares(:,1)' ares(:,2)' ares(:,3)' ares(:,4)', ares(:,5)' ares(1,6)'];
% Compute errors
mcm_err = mcme - repmat(ares_w,4,1);
ghq_err = ghqe - repmat(ares_w,3,1);
tauch_err = tauche - ares_w;

%% Plot 1st Moments by variable and method
figure()
subplot(2,2,1)
hold on
plot(log(n), mcme(:,1),'.', 'MarkerSize',16)
ylim([-0.5, 0.5])
yline(0, 'k--')
xticks(log(n))
xticklabels(n)
title('(a) Monte Carlo Method')
xlabel('N')
ylabel('E[X]')
hold off

subplot(2,2,2)
hold on
plot(quads, ghqe(:,1),'.', 'MarkerSize',16)
ylim([-0.05, 0.05])
yline(0, 'k--')
title('(b) Gauss-Hermite Quadrature')
xlabel('Quadrature Points')
hold off

subplot(2,2,3)
hold on
plot(log(n), mcme(:,2),'.', 'MarkerSize',16)
ylim([4, 8])
yline(5, 'k--')
xticks(log(n))
xticklabels(n)
title('(c) Monte Carlo Method')
xlabel('N')
ylabel('E[Y]')
hold off

subplot(2,2,4)
hold on
plot(quads, ghqe(:,2),'.', 'MarkerSize',16)
ylim([4, 8])
yline(5, 'k--')
title('(d) Gauss-Hermite Quadrature')
xlabel('Quadrature Points')
hold off

saveas(gcf, 'm1_chart','epsc');

%% Plot 2nd Moments by method
figure()
subplot(2,2,1)
hold on
plot(log(n), mcme(:,3),'.', 'MarkerSize',16)
ylim([0, 8])
yline(4, 'k--')
xticks(log(n))
xticklabels(n)
title('(a) Monte Carlo Method')
xlabel('N')
ylabel('E[X^2]')
hold off

subplot(2,2,2)
hold on
plot(quads, ghqe(:,3),'.', 'MarkerSize',16)
ylim([3, 5])
yline(4, 'k--')
title('(b) Gauss-Hermite Quadrature')
xlabel('Quadrature Points')
hold off

subplot(2,2,3)
hold on
plot(log(n), mcme(:,4),'.', 'MarkerSize',16)
ylim([10, 40])
yline(26, 'k--')
xticks(log(n))
xticklabels(n)
title('(c) Monte Carlo Method')
xlabel('N')
ylabel('E[Y^2]')
hold off

subplot(2,2,4)
hold on
plot(quads, ghqe(:,4),'.', 'MarkerSize',16)
ylim([24, 28])
yline(26, 'k--')
title('(d) Gauss-Hermite Quadrature')
xlabel('Quadrature Points')
hold off

saveas(gcf, 'm2_chart','epsc');

%% Plot 3rd Moment
figure()
subplot(2,2,1)
hold on
plot(log(n), mcme(:,5),'.', 'MarkerSize',16)
ylim([-1,1])
yline(0, 'k--')
xticks(log(n))
xticklabels(n)
title('(a) Monte Carlo Method')
xlabel('N')
ylabel('E[X^3]')
hold off

subplot(2,2,2)
hold on
plot(quads, ghqe(:,5),'.', 'MarkerSize',16)
ylim([-1, 1])
yline(0, 'k--')
title('(b) Gauss-Hermite Quadrature')
xlabel('Quadrature Points')
hold off

subplot(2,2,3)
hold on
plot(log(n), mcme(:,6),'.', 'MarkerSize',16)
ylim([100, 200])
yline(140, 'k--')
xticks(log(n))
xticklabels(n)
title('(c) Monte Carlo Method')
xlabel('N')
ylabel('E[Y^3]')
hold off

subplot(2,2,4)
hold on
plot(quads, ghqe(:,6),'.', 'MarkerSize',16)
ylim([120, 160])
yline(140, 'k--')
title('(d) Gauss-Hermite Quadrature')
xlabel('Quadrature Points')
hold off

saveas(gcf, 'm3_chart','epsc');

%% Plot 6th Moment
figure()
subplot(2,2,1)
hold on
plot(log(n), mcme(:,7),'.', 'MarkerSize',16)
%ylim([0, 8])
yline(ares(1,4), 'k--')
xticks(log(n))
xticklabels(n)
title('(a) Monte Carlo Method')
xlabel('N')
ylabel('E[X^6]')
hold off

subplot(2,2,2)
hold on
plot(quads, ghqe(:,7),'.', 'MarkerSize',16)
%ylim([3, 5])
yline(ares(1,4), 'k--')
title('(b) Gauss-Hermite Quadrature')
xlabel('Quadrature Points')
hold off

subplot(2,2,3)
hold on
plot(log(n), mcme(:,8),'.', 'MarkerSize',16)
%ylim([10, 40])
yline(ares(2,4), 'k--')
xticks(log(n))
xticklabels(n)
title({'(c) Monte Carlo Method';' '})
xlabel('N')
ylabel('E[Y^6]')
hold off

subplot(2,2,4)
hold on
plot(quads, ghqe(:,8),'.', 'MarkerSize',16)
%ylim([24, 28])
yline(ares(2,4), 'k--')
title({'(d) Gauss-Hermite Quadrature'; ' '})
xlabel('Quadrature Points')
hold off

saveas(gcf, 'm6_chart','epsc');

%% Plot 14th Moment
figure()
subplot(2,2,1)
hold on
plot(log(n), mcme(:,9),'.', 'MarkerSize',16)
%ylim([0, 8])
yline(ares(1,5), 'k--')
xticks(log(n))
xticklabels(n)
title({'(a) Monte Carlo Method';' '})
xlabel('N')
ylabel('E[X^{14}]')
hold off

subplot(2,2,2)
hold on
plot(quads, ghqe(:,9),'.', 'MarkerSize',16)
%ylim([3, 5])
yline(ares(1,5), 'k--')
title({'(b) Gauss-Hermite Quadrature';' '})
xlabel('Quadrature Points')
hold off

subplot(2,2,3)
hold on
plot(log(n), mcme(:,10),'.', 'MarkerSize',16)
%ylim([10, 40])
yline(ares(2,5), 'k--')
xticks(log(n))
xticklabels(n)
title({'(c) Monte Carlo Method';' '})
xlabel('N')
ylabel('E[Y^{14}]')
hold off

subplot(2,2,4)
hold on
plot(quads, ghqe(:,10),'.', 'MarkerSize',16)
%ylim([24, 28])
yline(ares(2,5), 'k--')
title({'(d) Gauss-Hermite Quadrature'; ' '})
xlabel('Quadrature Points')
hold off

saveas(gcf, 'm14_chart','epsc');

%% Plot 5th Order Cross Moment
figure()
subplot(2,1,1)
hold on
plot(log(n), mcme(:,11),'.', 'MarkerSize',16)
%ylim([0, 8])
yline(590, 'k--')
xticks(log(n))
xticklabels(n)
title('(a) Monte Carlo Method')
xlabel('N')
ylabel('E[X^2Y^3]')
hold off

subplot(2,1,2)
hold on
plot(quads, ghqe(:,11),'.', 'MarkerSize',16)
ylim([580, 600])
yline(590, 'k--')
title('(b) Gauss-Hermite Quadrature')
xlabel('Quadrature Points')
hold off

saveas(gcf, 'm5x_chart','epsc');