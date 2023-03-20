function hess = numhess(fun,x1,x2)
h = 1e-7;
ddx = (fun(x1+2*h,x2)-2*fun(x1+h,x2) + fun(x1,x2))/h^2;
ddy = (fun(x1,x2+2*h)-2*fun(x1,x2+h) + fun(x1,x2))/h^2;
dxdy = (fun(x1+h,x2+h)-fun(x1+h,x2)-fun(x1,x2+h)+fun(x1,x2))/h^2;
hess = [ddx dxdy; dxdy ddy];
end
    


