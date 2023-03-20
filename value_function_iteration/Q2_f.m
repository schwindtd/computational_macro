%% PS2 - (f) Use FOC for finding optimal k'
%
% Daniel Schwindt
% 9/17/2022

%% Set-up environment
% check whether script is run directly or not. If directly, clear all etc.
stack = dbstack;
if length(stack) == 1
    clear all
    close all
    clc
end
stoch_growth_setup
mVF = ones(nGridCap, m).*(1-bbeta)*log(consumSS);
% modify initial mVF guess so there is some curvature
fac = linspace(1.5, 0.5, n);
mVF = mVF.*fac';
options = optimoptions('fsolve','Display','off');

%% Compute True Solution (Guess & Verify Method)
% V(k,z) = B + Cln(k) + Dz
% Normalized
ggamma = (1-bbeta);
C = aalpha*ggamma/(1-aalpha*bbeta);
D = (1-bbeta + bbeta*C)/(1-rrho);
B = log(ggamma) + bbeta*C/ggamma*log(bbeta*C) ...
    -(1-bbeta+bbeta*C)/(1-bbeta)*log(1-bbeta + bbeta*C) ...
    + (1-bbeta + bbeta*C)/(1-bbeta)*log(A);

mVF_true = B + C*log(vGridCapital') + D*prod;
mPF_true = bbeta*C/(1- bbeta + bbeta*C).*A.*exp(prod).*(vGridCapital').^aalpha;

%% Normalized VFI with initial guess using Steady State & Using FOC
maxDiff = 10.0;
iter = 0;
time1 = tic;
while (maxDiff > tol)
    eVF = mVF*P';
    deriv = numderiv(mVF, vGridCapital);
    ederiv = deriv*P';
   
   % Loop over Productivity vector
   for nProd = 1:length(prod)
       % Loop over Capital Grid (today)
       for nCap = 1:nGridCap
           % Use numerical FOC to find k'
           capChoice = fsolve(@(k)foc(k, prod, vGridCapital, ...
               ederiv, A, aalpha, nCap, nProd, bbeta, ddelta), 0.123, options);
           cons = mY(nCap, nProd) - capChoice;
           eVFinterp = interp1(vGridCapital, eVF(:,nProd),capChoice, 'linear', 'extrap');
           vtmp = (1-bbeta)*log(cons) + bbeta*eVFinterp;
           
           % Store value and capital choice in Value Function and policy
           % matrices
           mVFNew(nCap, nProd) = vtmp;
           mPF(nCap, nProd) = capChoice;
           
       end
       
   end
   
   % Compute maximum difference & set value function equal to new value
   % function for next while loop iteration
   maxDiff = max(max(abs(mVFNew - mVF)));
   mVF = mVFNew;
   
   % Store and print iteration information
   iter = iter+1;
   if (mod(iter,10)==0 || iter ==1)
      fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, maxDiff); 
   end
end
time=toc(time1); % compute & store algorithm run time

%% Compute the maximum absolute error (Value Function & Policy Function)
err = (mVF - mVF_true)./mVF_true;
max_err = max(max(abs(err)));
err_pf = (mPF - mPF_true)./mPF_true;
max_err_pf = max(max(abs(err_pf)));