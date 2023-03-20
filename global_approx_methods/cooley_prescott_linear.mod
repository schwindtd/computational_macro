% Cooley & Prescott RBC Model with Original Utility Form
%

% Edited by Daniel Schwindt

%----------------------------------------------------------------
% 0. Housekeeping
%----------------------------------------------------------------

close all;    % Do not use clear all -- Dynare automatically does that and if you include this command, it would create problems. 

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var c k h l z;
varexo e;

parameters theta delta rho sigma_eps gamma beta alpha sigma eta;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

% Technology
theta  = 0.4;
delta  = 0.012;
rho    = 0.95;
sigma_eps  = 0.035;
gamma  = (1.0156^0.25)-1;


% Preferences
beta   = 0.987;
sigma = 20;
alpha  = 0.64;
eta    = (1.012^0.25)-1;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model;
  alpha/(1-alpha)*c/l = (1-theta)*exp(z)*k(-1)^theta*h^-theta;
  (1+gamma)*(c^(1-alpha)*l^alpha)^-sigma*c^-alpha*l^alpha = beta*((c(+1)^(1-alpha)*l(+1)^alpha)^-sigma*c(+1)^-alpha*l(+1)^alpha)*(theta*exp(z(+1))*k^(theta-1)*h(+1)^(1-theta) + 1 - delta);
  c + (1+eta)*(1+gamma)*k = exp(z)*k(-1)^theta*h^(1-theta) + (1-delta)*k(-1);
  l + h = 1;
  z = rho*z(-1)+e;
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
  k = 10;
  c = 2;
  h = 0.31;
  z = 0;
  l = 0.69;
end;

shocks;
var e = sigma_eps^2;
end;

steady;

check;

stoch_simul(order = 1, nograph, nocorr, nomoments, nofunctions);