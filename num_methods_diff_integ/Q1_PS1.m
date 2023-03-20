%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem Set #1
% Question 1: Numerical Optimization
% ECON630 FALL 2022
% Authors: Giuliano Simoncelli, Daniel Schwindt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc

%% 1- Initialise variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [6;-5];  %  USER: Set initial guess  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-7;   % Set the tolerance level
maxits = 50000;  % Number of iterations
syms x y


%% 2- USER: Manually insert function to be minimized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fun = @(x,y)((y.^2 - x.^2).^2 + (x - 1).^2);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% 3- Gradient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL: Input Jacobian manually  
grad_input = [ ; ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check to see if user wrote down Jacobian manually
if ~isempty(grad_input) 
    grad_fun = matlabFunction(grad_input);
else 
    grad_fun = matlabFunction([diff(fun,x) diff(fun,y)]);
end




%% 4- Hessian matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL: Input Hessian manually  
hessian_input = [ , ; , ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Check to see if user wrote down Hessian manually
if ~isempty(hessian_input)
    hessian_fun = matlabFunction(hessian_input);
else
    hessian_fun = matlabFunction([diff(diff(fun,x),x) diff(diff(fun,x),y);
    diff(diff(fun,x),y) diff(diff(fun,y),y)]);
end



%% 5- Optional tweaks to Newton's method in optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL: Modify the step size
option = NaN;   % 1 = Constant step size
scalar = NaN;   % Input a scalar to change the step size. Default is 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL: Use numerical derivative? 1 = YES
numerical_derivative = NaN;    % 1 = YES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Print options utilized
if isnan(option) && isnan(scalar)
    fprintf("Using constant step size? No \n\n");
    fprintf("Step size is scaled to a factor of 1 \n \n");
elseif ~isnan(option) && ~isnan(scalar)
    fprintf("Using constant step size? Yes \n\n");
    fprintf("Step size is scaled to a factor of %2.4f \n \n",scalar);
elseif ~isnan(option) && isnan(scalar)
    fprintf("Using constant step size? Yes \n\n");
    fprintf("Step size is scaled to a factor of 1 \n \n");
else
    fprintf("Using constant step size? No \n\n");
    fprintf("Step size is scaled to a factor of %2.4f \n \n",scalar);
end

if numerical_derivative == 1
    fprintf("Using numerical derivatives? Yes  \n\n")
else
    fprintf("Using numerical derivatives? No  \n\n")
end

%% 6- Newton's Method in Optimization
[minima, iters, flag] = newton_min(fun,grad_fun,hessian_fun,x0,tol,maxits,option,scalar,numerical_derivative);


%% 7- Text output with relevant information.
fprintf ("Newton's minimization method\n\n");
switch (flag)
    case 0
        fprintf ("There was a convergence to a nonoptimal point\n\n");
        fprintf ("Restart with new x0 \n\n");        
        fprintf("It took %d iterations.\n\n",iters);
    case 1
        fprintf ("There was a convergence on x\n\n");
        fprintf("The minima found is: \n");
        disp(minima);
        fprintf("It took %d iterations.\n\n",iters);
    otherwise
        fprintf ("There was no convergence\n\n");
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% 8- FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%----------------------------------------------------------------------%
%                       FUNCTION: Newton's Method                      %
%----------------------------------------------------------------------%


function [minima, iters, flag] = newton_min(funct,df,ddf,x0,maxtol,maxits,opt,scal,numderiv)
% The function returns minima, number of iterations, and corresponding text output.

% It takes as inputs the function itself, its Jacobian and Hessian, the initial
% guess, the tolerance allowed, and maximum number of iterations. 

% It also allows for the option to use an identity matrix instead of the 
% hessian to update x, to use a scaling parameter, and to use numerical
% derivatives.
    


    % Initialise values
    x = x0(1); y = x0(2);
    minima = NaN;
    flag = -1;
    
    % Main loop
    for iters = 1 : maxits
    
      
        x_old = [x;y];
        
        % Using numerical derivatives?
        if numderiv == 1
            df = @(x,y)(numgrad(funct,x,y));
            ddf = @(x,y)(numhess(funct,x,y));
        end 
        
        
        % Update value of x allowing for different step sizes
        if opt == 1 && ~isnan(scal)
            x_new = x_old - scal * eye(2) * df(x_old(1),x_old(2))';
        elseif opt == 1 && isnan(scal)
            scal = 1;
            x_new = x_old - scal * eye(2) * df(x_old(1),x_old(2))';
        elseif opt ~= 1 && ~isnan(scal)
            % Check if Hessian is singular
            if rank(ddf(x_old(1),x_old(2)))==2
                x_new = x_old - scal * inv(ddf(x_old(1),x_old(2))) * df(x_old(1),x_old(2))';
            else
                flag=-1;
            end
        else
            scal = 1;
            % Check if Hessian is singular
            if rank(ddf(x_old(1),x_old(2)))==2
                x_new = x_old - scal * inv(ddf(x_old(1),x_old(2))) * df(x_old(1),x_old(2))';
            else
                flag=-1;
            end
        end
        
        
        % Check if we found the minima
        if norm(x_new - x_old) < (eps + eps*norm(x_old))
            
            flag = 1;
            minima = x_new';
            
              
                if norm(df(x,y)) > maxtol*(1+norm(funct(x,y))) || (any(eig(ddf(x_old(1),x_old(2)))<=0)) % Second logical assures we are at global minima
            
                flag = 0;
                minima = x_new';
                return
                end
               
            
            return
        end    
            
        %Report values       
        x = x_new(1);      
        y = x_new(2);
    
      
    
    end
end





%%
%----------------------------------------------------------------------%
%                   FUNCTION: Numerical Derivatives                    %
%----------------------------------------------------------------------%

function grad = numgrad(fun,x1,x2)
h = 1e-7;
dx = (fun(x1+h,x2)-fun(x1,x2))/h;    
dy = (fun(x1,x2+h)-fun(x1,x2))/h;
grad = [dx dy];
end

function hess = numhess(fun,x1,x2)
h = 1e-7;
ddx = (fun(x1+h,x2)-2*fun(x1,x2) + fun(x1-h,x2))/h^2;
ddy = (fun(x1,x2+h)-2*fun(x1,x2) + fun(x1,x2-h))/h^2;
dxdy = (fun(x1+h,x2+h)-fun(x1+h,x2)-fun(x1,x2+h)+fun(x1,x2))/h^2;
hess = [ddx dxdy; dxdy ddy];
end