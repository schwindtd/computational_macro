%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem Set #1
% Question 1 Addendum: Heatmap
% ECON630 FALL 2022
% Authors: Giuliano Simoncelli, Daniel Schwindt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc

%% 1- Initialise variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-7;   % Set the tolerance level 
maxits = 5000;  % Number of iterations
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
hessian_input = [  ;  ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Check to see if user wrote down Hessian manually
if ~isempty(hessian_input)
    hessian_fun = matlabFunction(hessian_input);
else
    hessian_fun = matlabFunction([diff(diff(fun,x),x) diff(diff(fun,x),y);
    diff(diff(fun,x),y) diff(diff(fun,y),y)]);
end

%% 5- Creating grid of initial guesses
x = linspace(-6,6,66)'; % Number of gridpoints in heatmap (can edit)
y = x;
x0 = nan(1,2);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% 6- Loops to create Heatmap %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 6.a - Regular analytical (case a)
  z1 = nan(length(x),length(x));
    for i=1:length(x)
      for j=1:length(x)
          x0(1) = x(i);
          x0(2) = y(j);
          [minima1, iters1, flag1] = newton_min(fun,grad_fun,hessian_fun,x0,tol,maxits,NaN,NaN,0);
      if flag1 == 1
          z1(length(x)-j+1,i) = 1;     
      end
      end
   end
  
x = linspace(-6,6,66)';
y = x;
x0 = nan(1,2);

% 6.b - Regular numerical (case (b)
  z2 = nan(length(x),length(x));
    for i=1:length(x)
      for j=1:length(x)
          x0(1) = x(i);
          x0(2) = y(j);
          [minima1, iters1, flag1] = newton_min(fun,grad_fun,hessian_fun,x0,tol,maxits,NaN,NaN,1);
      if flag1 == 1
          z2(length(x)-j+1,i) = 1;     
      end
      end
   end
  
  
x = linspace(-6,6,66)';
y = x;
x0 = nan(1,2);  
  
  
% 6.c - Identity, analytical (case c)
  z3 = nan(length(x),length(x));
    for i=1:length(x)
      for j=1:length(x)
          x0(1) = x(i);
          x0(2) = y(j);
          [minima1, iters1, flag1] = newton_min(fun,grad_fun,hessian_fun,x0,tol,maxits,1,0.01,0);
      if flag1 == 1
          z3(length(x)-j+1,i) = 1;     
      end
      end
   end
  
x = linspace(-6,6,66)';
y = x;
x0 = nan(1,2);  
  
  
  
% 6.d - Scaled hessian, analytical (case d)
  z4 = nan(length(x),length(x));
    for i=1:length(x)
      for j=1:length(x)
          x0(1) = x(i);
          x0(2) = y(j);
          [minima1, iters1, flag1] = newton_min(fun,grad_fun,hessian_fun,x0,tol,maxits,NaN,0.01,0);
      if flag1 == 1
          z4(length(x)-j+1,i) = 1;     
      end
      end
    end
 
 
    
 % 6.e - Plot heatmaps  
 
 x = linspace(-6,6,66)';
 y = x;   
 % Subplots   
 subplot(2,2,1)   
 heatmap(z1)
 colorbar off
 subplot(2,2,2)
 heatmap(z2)
 colorbar off
 subplot(2,2,3)
 heatmap(z3) 
 colorbar off
 subplot(2,2,4)
 heatmap(z4) 
 colorbar off
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% 7- FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%
%----------------------------------------------------------------------%
%                                FUNCTION                              %
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
        elseif isnan(opt) && ~isnan(scal)
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