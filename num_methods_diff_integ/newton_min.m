%%
%----------------------------------------------------------------------%
%                                FUNCTION                              %
%----------------------------------------------------------------------%

% The function returns minima, iterations, and corresponding text.
% Takes as inputs the function itself, its Jacobian and Hessian, the initial
% guess, the tolerance and maximum number of iterations. It also allows for
% the option to use an identity matrix instead of the hessian to update x.
function [minima, iters, flag] = newton_min(funct,df,ddf,x0,maxtol,maxits,opt,scal,numderiv)
    
    


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
            x_new = x_old - scal * inv(ddf(x_old(1),x_old(2))) * df(x_old(1),x_old(2))';
        else
            scal = 1;
            x_new = x_old - scal * inv(ddf(x_old(1),x_old(2))) * df(x_old(1),x_old(2))';
        end
        
        
        % Check if we found the minima
        if norm(x_new - x_old) < (eps + eps*norm(x_old))
            
            flag = 1;
            minima = x_new;
            
            if isnan(opt)  % No need if we are using identity, I think
                if norm(df(x,y)) > maxtol*(1+norm(funct(x,y))) || (any(eig(ddf(x_old(1),x_old(2)))<=0)) % Second logical assures we are at global minima
            
                flag = 0;
                minima = x_new;
                return
                end
            end   
            
            return
        end    
            
        %Report values       
        x = x_new(1);      
        y = x_new(2);
    
      
    
    end
end





function grad = numgrad(fun,x1,x2)
h = 1e-7;
dx = (fun(x1+h,x2)-fun(x1,x2))/h;    
dy = (fun(x1,x2+h)-fun(x1,x2))/h;
grad = [dx dy];
end

function hess = numhess(fun,x1,x2)
h = 1e-7;
ddx = (fun(x1+2*h,x2)-2*fun(x1+h,x2) + fun(x1,x2))/h^2;
ddy = (fun(x1,x2+2*h)-2*fun(x1,x2+h) + fun(x1,x2))/h^2;
dxdy = (fun(x1+h,x2+h)-fun(x1+h,x2)-fun(x1,x2+h)+fun(x1,x2))/h^2;
hess = [ddx dxdy; dxdy ddy];
end
    

