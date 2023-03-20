%% Stochastic Growth Model Set-up File
%
% Daniel Schwindt
% 9/17/2022
% Creates calibration parameters and other needed vectors and matrices 
% for each of the different VFIs

%%  1. Calibration
A = 1;
rrho = 0.95;
ssigma = 0.007;
aalpha = 0.3;
ddelta = 1;
bbeta = 0.96;
tol = 1e-7;
n = 1000;
m=5; % # of productivity values
ngrid = [50, 150, 300];
%% 2. Steady State
capitalSS = (aalpha*A/(1/bbeta + ddelta - 1))^(1/(1-aalpha));
outputSS = A*capitalSS^aalpha;
consumSS = outputSS - capitalSS;
%% 3. Capital grid
vGridCapital = linspace(0.5*capitalSS, 1.5*capitalSS, n);
% Loop to create coarser grids
for i=1:length(ngrid)
    eval(['vGridCapital' num2str(ngrid(i)) ' = linspace(0.5*capitalSS, 1.5*capitalSS,' num2str(ngrid(i)) ');' ])
end
% Scalars for size of capital grid and productivity grid
nGridCap = length(vGridCapital);

%% 4. Productivity Grid
[prod, P, pi] = tauchen(rrho, ssigma, m, 3);

%% 4. Pre-Built Matrices/Vectors
% For n=1000 case
mY = A*(vGridCapital'.^aalpha)*exp(prod); % Pre-build possible outputs
mVF    = zeros(nGridCap,m);
mVFNew = zeros(nGridCap,m);
mPF   = zeros(nGridCap,m);
eVF = zeros(nGridCap,m);

% Loop to create coarser grids
for i=1:length(ngrid)
    ngrid_str = num2str(ngrid(i));
    eval(['mY' ngrid_str ' = A* (vGridCapital' ngrid_str '''.^aalpha)*exp(prod);' ])
    eval(['mVF' ngrid_str ' = zeros(' ngrid_str ',m);'])
    eval(['mVFNew' ngrid_str ' = zeros(' ngrid_str ',m);'])
    eval(['mPF' ngrid_str ' = zeros(' ngrid_str ',m);'])
    eval(['eVF' ngrid_str ' = zeros(' ngrid_str ',m);'])
end