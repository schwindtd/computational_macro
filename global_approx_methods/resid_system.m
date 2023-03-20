function G = residual(theta_in, params, psi_in, K, A, transition)

% Residual function for the stochastic growth model with labor 

% Author: Daniel Schwindt 
% Date: 10/11/2022

% Compute implied decision rules

H = psi_in' * theta_in;

row = 1;

% This loop computes the residual for each pair (k,a)

for i = 1:params.n  
    for j = 1:params.m
        
        % States
        a = A(j,:);
        k = K(i,:);
        
        % Decision rules
        h = H(row,:);
        
        % Error Check        
        if h < 0
            h = 10e-6;
            disp('Warning : At EE1, h < 0 encountered!!!')
        elseif h > 1
            h = 1-10e-6;
            disp('Warning: At EE1, h > 1 encountered!!!')
        end
        
        % Create consumption and leisure choice vars
        c = ((1-params.alpha)/(params.alpha))*(1-params.theta)* ...
            exp(a)*k^params.theta*h^(-params.theta)*(1-h);
        l = 1-h;
        y = exp(a) * (k ^ params.theta)*(h^(1-params.theta));
        u = (1-params.alpha)*c^-params.alpha * l^params.alpha* ...
            (c^(1-params.alpha) * l^params.alpha)^(-params.sigma);
        
        % Set capital next period
        k_prime = (y + (1 - params.delta) * k - c)/((1+params.eta)*(1+params.gamma));
        
        % Restrict k_prime to be within [k_min, k_max]
        if k_prime < params.k_min
            k_prime = 1.1*params.k_min;
        elseif k_prime > params.k_max
            k_prime = 0.9*params.k_max;
        end
        
        x_prime_k = 2 * (k_prime - params.k_min) / (params.k_max - params.k_min) -1; % This is k_prime projected into [-1,1]
        psi_k_prime = cos([0:(params.n-1)]'*acos(x_prime_k));  % These are n*1         
        
        % For EE
        
        for count = 1:params.m
    
            a_prime = A(count);    
           
            % Construct psi_a_prime
            
            x_prime_a = 2 * (a_prime - params.Z_min) / (params.Z_max - params.Z_min) -1; % This is a_prime projected into [-1,1]
    
            psi_a_prime = cos([0:(params.m-1)]'*acos(x_prime_a)); % These are m*1

            psi_prime = kron(psi_k_prime,psi_a_prime);
            
            h_prime = psi_prime' * theta_in;
            
            % Error Check        
            if h_prime < 0
                h = 10e-6;
                disp('Warning : At EE1, h_prime < 0 encountered!!!')
            elseif h_prime > 1
                h = 1-10e-6;
                disp('Warning: At EE1, h_prime > 1 encountered!!!')
            end
            
            % Create consumption and leisure choice vars
            c_prime = ((1-params.alpha)/(params.alpha))*(1-params.theta)* ...
                exp(a_prime)*k_prime^params.theta*h_prime^(-params.theta)*(1-h_prime);
            l_prime = 1-h_prime;
            
            % Define u' and f' for EE
            u_prime = (1-params.alpha)*(c_prime^(1-params.alpha) * ...
                l_prime ^ params.alpha)^(-params.sigma)* ...
                c_prime^-params.alpha * l_prime^params.alpha;            
            f_prime = params.theta * exp(a_prime) * (k_prime)^(params.theta-1)* ...
                (h_prime)^(1-params.theta);
            
            integrand(:,count) = u_prime * (f_prime + 1 - params.delta);        
                
        end

        disc_sol_integ = params.beta * integrand * transition(j,:)';
        lhs = (1+params.gamma)*u;
        
        res(row,:) = 1 - disc_sol_integ/lhs;
        
        row = row + 1;
        
    end
end


G = res;
end