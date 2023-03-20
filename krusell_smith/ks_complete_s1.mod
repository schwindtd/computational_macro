% Krusell & Smith (1998) Complete Markets
%

% Edited by Daniel Schwindt

%----------------------------------------------------------------
% 0. Housekeeping
%----------------------------------------------------------------

close all;    % Do not use clear all -- Dynare automatically does that and if you include this command, it would create problems. 

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y c k r w z l;
varexo eps_z eps_l;

parameters bbeta ssigma aalpha ddelta rrho ssigma_eps_z;

%----------------------------------------------------------------
% 2. Set Parameters
%----------------------------------------------------------------

% Technology
aalpha = 0.36;
ddelta  = 0.025;
rrho    = 0.7622;
ssigma_eps_z  = 0.01;
ssigma_eps_l  = 0.03;

% Preferences
bbeta   = 0.99;
ssigma = 1;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model;
  # uc = c^-ssigma;
  # ucp = c(+1)^-ssigma;
  uc = bbeta*ucp*(1-ddelta+r(+1));
  c + k = r*k(-1) + w + (1-ddelta)*k(-1);
  r = aalpha*z*k^(aalpha-1)*l^(1-aalpha);
  w = (1-aalpha)*z*k^aalpha*l^-aalpha;
  y = z*k^aalpha*l^(1-aalpha);
  z = (1-rrho)*1 + rrho*z + eps_z;
  l = (1-rrho)*0.93 + rrho*l + eps_l;
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------
initval;
z = 1.01;
l = 0.96;
r = 1/bbeta - 1 + ddelta;
k = ((1/bbeta - 1 + ddelta)/(aalpha*l^(1-aalpha)))^(1/(aalpha-1));
w = (1-aalpha)*k^aalpha*l^-aalpha;
c = (r-ddelta)*k + w;
y = k^aalpha*l^(1-aalpha);
end;

shocks;
var eps_z = ssigma_eps_z^2;
var eps_l = ssigma_eps_l^2;
end;

steady;

check;

stoch_simul(order = 1, nograph);