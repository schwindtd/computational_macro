function [mValFun, mPF, iter] = stoch_gr_vfi(vGridCap,mValFun, p, mOutput,  disc, maxdiff, toler)
%Value Function Iteration function
%   Performs value function iteration using Steady State & Normalized
%   utility
mValFunNew = zeros(size(mValFun));
mPF = zeros(size(mValFun));
iter = 0;
while (maxdiff > toler)
   eVF = mValFun*p'; % create expected value function
   
   % Loop over Productivity vector
   for i = 1:length(p)
       % Loop over Capital Grid (today)
       for j = 1:length(vGridCap)
           vHigh = -1000.0;
           % Loop over Capital Grid (tomorrow)
           for k = 1:length(vGridCap)
                cons = mOutput(j,i) - vGridCap(k);
                vtmp = (1-disc)*log(cons) + disc*eVF(k,i);
                % Select max value and accompanying capital choice
                if (vtmp > vHigh)
                    vHigh = vtmp;
                    capChoice = vGridCap(k);
                else
                    break; % End grid search if maximum found
                end
           end
           
           % Store value and capital choice in Value Function and policy
           % matrices
           mValFunNew(j, i) = vHigh;
           mPF(j, i) = capChoice;
           
       end
       
   end
   
   % Compute maximum difference & set value function equal to new value
   % function for next while loop iteration
   maxdiff = max(max(abs(mValFunNew - mValFun)));
   mValFun = mValFunNew;
   
   % Store and print iteration information
   iter = iter+1;
   if (mod(iter,10)==0 || iter ==1)
      fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, maxdiff); 
   end
end
end

