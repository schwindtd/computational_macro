%%  New Keynesian DSGE Model (An & Schorfheide '07)

%--------------------------------------------------------------------------
% 1. Declare variables
%--------------------------------------------------------------------------

var c R g z y dp;
varexo eps_R eps_G eps_Z;

parameters tau kappa psi1 psi2 r_A dp_A gamma_Q; 
parameters rho_R rho_G rho_Z sig_R sig_G sig_Z;

%--------------------------------------------------------------------------
% 2. Initialize Parameters
%--------------------------------------------------------------------------
verbatim;
load('params.mat');
set_param_value('tau', params.tau);
set_param_value('kappa', params.kappa);
set_param_value('psi1', params.psi1);
set_param_value('psi2', params.psi2);
set_param_value('dp_A', params.dp_A);
set_param_value('gamma_Q', params.gamma_Q);
set_param_value('r_A', params.r_A);
set_param_value('rho_G', params.rho_G);
set_param_value('rho_R', params.rho_R);
set_param_value('rho_Z', params.rho_Z);
set_param_value('sig_G', params.sig_G);
set_param_value('sig_R', params.sig_R);
set_param_value('sig_Z', params.sig_Z);
end;
%--------------------------------------------------------------------------
% 3. Model
%--------------------------------------------------------------------------

model(linear);

# betta = 1/(1+r_A/400);        % Parameter dependence

y = y(+1) + g - g(+1) - (R - dp(+1) - z(+1)) / tau;
dp = betta * dp(+1) + kappa * (y - g);
y = c + g;

//R = rho_R * R(-1) + (1 - rho_R) * psi1 * dp + (1 - rho_R) * psi2 * (y - g) + sig_R * eps_R;
//g = rho_G * g(-1) + sig_G * eps_G;
//z = rho_Z * z(-1) + sig_Z * eps_Z;
R = rho_R * R(-1) + (1 - rho_R) * psi1 * dp + (1 - rho_R) * psi2 * (y - g) + eps_R;
g = rho_G * g(-1) + eps_G;
z = rho_Z * z(-1) + eps_Z;

end;

%--------------------------------------------------------------------------
% 4. Computation
%--------------------------------------------------------------------------

//initval;
//end;
//steady;
//check;

shocks;
    var eps_R = sig_R^2;
    var eps_G = sig_G^2;
    var eps_Z = sig_Z^2;
end;

stoch_simul(order=1,irf=12,periods=0, nograph, noprint);
