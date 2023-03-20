%% Stochastic Growth Model
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

%% Compute True Solution (Guess & Verify Method)
% V(k,z) = B + Cln(k) + Dz
% NOT Normalized
ggamma = 1-bbeta;
C = aalpha/(1-aalpha*bbeta);
D = (1+bbeta*C)/(1-rrho);
B = bbeta*C/ggamma*log(bbeta*C) ...
    - (1+bbeta*C)/ggamma*log(1+bbeta*C) ...
    +(1+bbeta*C)/ggamma*log(A);

mVF_true = B + C*log(vGridCapital') + D*prod;
mPF_true = bbeta*C/(1+bbeta*C).*A.*exp(prod).*(vGridCapital').^aalpha;

%% Naive VFI
maxDiff = 10.0;
iter = 0;
time1 = tic;
while (maxDiff > tol)
   eVF = mVF*P'; % create expected value function
   
   % Loop over Productivity vector
   for nProd = 1:length(prod)
       % Loop over Capital Grid (today)
       for nCap = 1:nGridCap
           vHigh = -1000.0;
           % Loop over Capital Grid (tomorrow)
           for nCapNext = 1:nGridCap
                cons = mY(nCap,nProd) - vGridCapital(nCapNext);
                vtmp = (cons.^-1)./(-1) + bbeta*eVF(nCapNext,nProd);
                % Select max value and accompanying capital choice
                if (vtmp > vHigh)
                    vHigh = vtmp;
                    capChoice = vGridCapital(nCapNext);
                else
                    break; % End grid search if maximum found
                end
           end
           
           % Store value and capital choice in Value Function and policy
           % matrices
           mVFNew(nCap, nProd) = vHigh;
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
naiveVFItime=toc(time1); % compute & store algorithm run time

%% Compute the maximum absolute error (Value Function & Policy Function)
err = (mVF - mVF_true)./mVF_true;
max_err = max(max(abs(err)));
err_pf = (mPF - mPF_true)./mPF_true;
max_err_pf = max(max(abs(err_pf)));
