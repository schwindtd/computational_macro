function grad = numgrad(fun,x1,x2)
h = 1e-7;
dx = (fun(x1+h,x2)-fun(x1,x2))/h;    
dy = (fun(x1,x2+h)-fun(x1,x2))/h;
grad = [dx dy];
end