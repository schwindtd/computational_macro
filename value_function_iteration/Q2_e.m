%% PS2 - (e) Do #b with monotonicity and include concavity with binary search.
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

%% Normalized VFI with initial guess using Steady State
maxDiff = 10.0;
iter = 0;
time1 = tic;
while (maxDiff > tol)
   eVF = mVF*P'; % create expected value function
   
   % Loop over Productivity vector
   for nProd = 1:length(prod)
       % Loop over Capital Grid (today)
       nCapStart = 1;
       for nCap = 1:nGridCap
           
           % Binary Search
           % Initialize index min and max
           i_min = nCapStart;
           i_max = nGridCap;
           check = i_max - i_min;
           
           % While Loop until check=2
           while (check > 2)
               % bisect remaining grid
               il = floor((i_min + i_max)/2);
               iu = il + 1;
               % compute consumption and Value function values
               cons_iu = (mY(nCap,nProd) - vGridCapital(iu));
               cons_il = (mY(nCap,nProd) - vGridCapital(il));
               V_iu = (1-bbeta)*log(cons_iu) + bbeta*eVF(iu, nProd);
               V_il = (1-bbeta)*log(cons_il) + bbeta*eVF(il, nProd);
               
               % discared range with decreasing slope
               if (V_iu > V_il)
                  i_min = il;
               else
                  i_max = iu;
               end

               check = i_max - i_min;
           end
            
            % evaluate at remaining three grid points
            if iu == i_max
                V_iMax = V_iu;
                V_iMin1 = V_il;
                V_iMin = (1-bbeta)*log(mY(nCap, nProd) ...
                    - vGridCapital(i_min)) + bbeta.*eVF(i_min, nProd);
            else
                V_iMax = (1-bbeta)*log(mY(nCap, nProd) ...
                    - vGridCapital(i_max)) + bbeta.*eVF(i_max, nProd);
                V_iMin1 = V_iu;
                V_iMin = V_il;
            end
            
            % Find maximizer among remaining three values
            [mVFNew(nCap, nProd), index] = max([V_iMin; V_iMin1; V_iMax]);
            mPF(nCap, nProd) = vGridCapital(index - 1 + i_min);
            
            % reset lower bound on k'
            nCapStart = index - 1 + i_min;
           
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