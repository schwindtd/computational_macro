function [m_k_nxt, m_c_star, iter, max_diff_k] = pfi_im(params, z, v_eff, lbar,b0, b1, P_11, P_12, P_21, P_22, v_k, v_kbar, m_k, m_kbar)
%% Policy Function Iteration
% Based on Maliar, Maliar, Valli (2010) implementation
% Grids
n_k = params.n_k;
n_kbar = params.n_kbar;
n_agg = params.n_agg;
n_eff = params.n_eff;
k_min = params.k_min;
k_max = params.k_max;
kbar_min = params.kbar_min;
kbar_max = params.kbar_max;
% Utility & Production
ssigma = params.ssigma;
bbeta = params.bbeta;
aalpha = params.aalpha;
ddelta = params.ddelta;
% Tolerance & Convergence
tol = params.tol;
max_iter = params.max_iter;
psi_k = params.psi_k;

%% Create 3D Grids for computations: (k, kbar, z)
% Capital Grids
m_k_nxt = 0.9*m_k; % Initial policy function guess (based on what Maliar, Maliar, Valli use in their paper)

% labor, agg productivity, and employment status 
m_lbar=nan(n_k,n_kbar,n_agg,n_eff);
m_z=nan(n_k,n_kbar,n_agg,n_eff);
m_eff=nan(n_k,n_kbar,n_agg,n_eff);
for i=1:n_agg
    m_lbar(:,:,i,:)=lbar(i)*ones(n_k,n_kbar,n_eff);
    m_z(:,:,i,:) = z(i)*ones(n_k,n_kbar,n_eff);
end

for j=1:n_eff
    m_eff(:,:,:,j) = v_eff(j)*ones(n_k,n_kbar,n_agg);
end

% Current period interest rates and wages by states                            
r = aalpha*(m_z.*(m_kbar./m_lbar).^(aalpha-1));                             
w = (1-aalpha)*(m_z.*(m_kbar./m_lbar).^aalpha);                               
income=(1 + r - ddelta).*m_k+(w.*m_eff);

% Next period aggregate capital based on LoM
m_kbar_nxt = nan(n_k,n_kbar,n_agg,n_eff);
for i=1:n_agg
    for j=1:n_eff
        m_kbar_nxt(:,:,i,j) = exp(b0(i)*ones(n_k, n_kbar) + b1(i)*log(m_kbar(:,:,i,j)));
    end
end
m_kbar_nxt(m_kbar_nxt < kbar_min) = kbar_min;
m_kbar_nxt(m_kbar_nxt > kbar_max) = kbar_max;

% Next period interest rate and wage conditional on current states
for i=1:n_agg
    eval(['r_nxt_' num2str(i) '=aalpha*z(' num2str(i) ')*((m_kbar_nxt./(lbar(' num2str(i) '))).^(aalpha-1));']);
    eval(['w_nxt_' num2str(i) '=(1-aalpha)*z(' num2str(i) ')*((m_kbar_nxt./(lbar(' num2str(i) '))).^(aalpha));']);
end

%% Loop for solving household problem   
iter = 0;
max_diff_k=100;
   while (max_diff_k > tol && iter < max_iter)
   
     % 1. Input k' guess into k' to compute k''
     %      NOTE: used spline because linear seemed to give odd results
     % 2. Given k'', compute c' using next period budget constraint
     % 3. Set any c'<=0 to a tiny positive number (to avoid negatives and
     % Inf like in normal VFI)
     % 4. Compute MU(c') based on utility function [will use to back out c]
     for i=1:n_agg
         for j=1:n_eff
             eval(['m_k_nxt2_' num2str(i) num2str(j) '=interpn(v_k'', v_kbar'', m_k_nxt(:,:,' num2str(i) ',' num2str(j) '),m_k_nxt, m_kbar_nxt, ''spline'');']);
             eval(['m_c_nxt_' num2str(i) num2str(j) '=(1+ r_nxt_' num2str(i) '- ddelta).*m_k_nxt + w_nxt_' num2str(i) '.*(v_eff(' num2str(j) ')) - m_k_nxt2_' num2str(i) num2str(j) ';']);
             eval(['m_c_nxt_' num2str(i) num2str(j) '(m_c_nxt_' num2str(i) num2str(j) '<=0) = 1e-6;']);
             eval(['m_mu_nxt_' num2str(i) num2str(j) '=m_c_nxt_' num2str(i) num2str(j) '.^-ssigma;']);
         end
     end

    % 1. MU(c) = B*E[MU(c')]
    % 2. Therefore, c = MU^-1(B*E[MU(c')])
    % 3. Given c, compute k' based on budget constraint
    % 4. Check for if k' is in the bounds set
    m_emu_nxt=m_mu_nxt_11.*(1-ddelta+r_nxt_1).*P_11 + m_mu_nxt_12.*(1-ddelta+r_nxt_1).*P_12+...
    m_mu_nxt_21.*(1-ddelta+r_nxt_2).*P_21 + m_mu_nxt_22.*(1-ddelta+r_nxt_2).*P_22; % From intertemporal Euler equation
    m_c_star=(bbeta*m_emu_nxt).^(-1/ssigma);
    m_k_nxt_new=income-m_c_star;
    m_k_nxt_new(m_k_nxt_new > k_max) = k_max;
    m_k_nxt_new(m_k_nxt_new < k_min) = k_min;
    
    % Sup Difference
    max_diff_k=max(max(max(max(abs(m_k_nxt_new-m_k_nxt)))));
    
    m_k_nxt=psi_k*m_k_nxt_new+(1-psi_k)*m_k_nxt; % Update as weighted avg based on Maliar, Maliar, Valli
   
    % Store and print iteration information
   iter = iter+1;
%    if (mod(iter,10)==0 || iter ==1)
%       fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, max_diff_k); 
%    end   
   end
disp('--------------------------------------------');
fprintf(' PFI Completed. Iteration = %d, Sup Diff = %2.8f\n', iter, max_diff_k);
end

