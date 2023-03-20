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

%% VFI with monotonicity for n1=50
maxDiff = 10.0;
iter = 0;
time1 = tic;
while (maxDiff > tol)
   eVF50 = mVF50*P'; % create expected value function
   
   % Loop over Productivity vector
   for nProd = 1:length(prod)
       nCapNextStart = 1;
       % Loop over Capital Grid (today)
       for nCap = 1:ngrid(1)
           vHigh = -1000.0;
           % Loop over Capital Grid (tomorrow)
           for nCapNext = nCapNextStart:ngrid(1)
                cons = mY50(nCap,nProd) - vGridCapital50(nCapNext);
                vtmp = log(cons) + bbeta*eVF50(nCapNext,nProd);
                % Select max value and accompanying capital choice
                if (vtmp > vHigh)
                    vHigh = vtmp;
                    capChoice = vGridCapital50(nCapNext);
                    nCapNextStart = nCapNext;
                else
                    break; % End grid search if maximum found
                end
           end
           
           % Store value and capital choice in Value Function and policy
           % matrices
           mVF50New(nCap, nProd) = vHigh;
           mPF50(nCap, nProd) = capChoice;
           
       end
       
   end
   
   % Compute maximum difference & set value function equal to new value
   % function for next while loop iteration
   maxDiff = max(max(abs(mVF50New - mVF50)));
   mVF50 = mVF50New;
   
   % Store and print iteration information
   iter = iter+1;
   if (mod(iter,10)==0 || iter ==1)
      fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, maxDiff); 
   end
end
%% Interpolate 
% Interpolate Value Function grid from n1=50 to n=1000
mVF = linear_interp(vGridCapital50, vGridCapital, mVF50);

%% VFI for n=1000 using interpolated grid
maxDiff = 10.0;
while (maxDiff > tol)
   eVF = mVF*P'; % create expected value function
   
   % Loop over Productivity vector
   for nProd = 1:length(prod)
       nCapNextStart = 1;
       % Loop over Capital Grid (today)
       for nCap = 1:nGridCap
           vHigh = -1000.0;
           % Loop over Capital Grid (tomorrow)
           for nCapNext = nCapNextStart:nGridCap
                cons = mY(nCap,nProd) - vGridCapital(nCapNext);
                vtmp = log(cons) + bbeta*eVF(nCapNext,nProd);
                % Select max value and accompanying capital choice
                if (vtmp > vHigh)
                    vHigh = vtmp;
                    capChoice = vGridCapital(nCapNext);
                    nCapNextStart = nCapNext;
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
time=toc(time1); % compute & store algorithm run time

%% Compute the maximum absolute error (Value Function & Policy Function)
err = (mVF - mVF_true)./mVF_true;
max_err = max(max(abs(err)));
err_pf = (mPF - mPF_true)./mPF_true;
max_err_pf = max(max(abs(err_pf)));