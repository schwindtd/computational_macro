%% Endogenous Grid Method
%
% Daniel Schwindt
% 9/21/2022

%% Set-up environment
% check whether script is run directly or not. If directly, clear all etc.
stack = dbstack;
if length(stack) == 1
    clear all
    close all
    clc
end
stoch_growth_setup
mVF = zeros(nGridCap, m);
mV  = ones(nGridCap, m);
mVNew = zeros(nGridCap, m);
% modify initial mV guess so there is some curvature
fac = linspace(0.5, 1, n);
mV = mV.*fac';
mYs = zeros(nGridCap, m);
options = optimoptions('fsolve','Display','off');

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

%% Endogenous Grid Method
% Note:
% vGridCapital = capital grid for tomorrow (k')
% mY = resources for tomorrow at each k' and z'
% mVF = V~

maxDiff = 10.0;
iter = 0;
time1=tic;
while (maxDiff > tol)
    mVk = numderiv(mV, vGridCapital);
    for nCap = 1:nGridCap
        for nProd = 1:length(prod)
            cons = mVk(nCap, nProd)^-1;
            mYs(nCap, nProd) = cons + vGridCapital(nCap);
            mVFNew(nCap, nProd) = log(cons) + mV(nCap, nProd);
        end
    end
    
    for nProd = 1:length(prod)
        mVF(:,nProd) = interp1(mYs(:,nProd), mVFNew(:,nProd), ...
        mY(:,nProd),'linear', 'extrap');
    end
    
    mVNew = bbeta*mVF*P';
    maxDiff = max(max(abs(mVNew - mV)));
    mV = mVNew;
    % Store and print iteration information
   iter = iter+1;
   if (mod(iter,10)==0 || iter ==1)
      fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, maxDiff); 
   end
end
time = toc(time1);
%% Extract k grid
vGridCapitalToday = zeros(nGridCap,length(prod));
fs = @(k, n, m) mYs(n, m) - exp(prod(m)).*k.^aalpha + (1-ddelta).*k;
for nCap = 1:nGridCap
    for nProd = 1:length(prod)
        vGridCapitalToday(nCap, nProd) = fsolve(@(k) fs(k, nCap, nProd), 0, options);
    end
end
time = toc(time1);

%% Compute the maximum absolute error (Value Function & Policy Function)
err = (mVF - mVF_true)./mVF_true;
max_err = max(max(abs(err)));
err_pf = (vGridCapitalToday - mPF_true)./mPF_true;
max_err_pf = max(max(abs(err_pf)));