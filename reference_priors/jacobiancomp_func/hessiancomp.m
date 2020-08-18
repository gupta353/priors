% Numerical computation of second derivative of a function

function hess = hessiancomp(fun,a,h)
      
    % computation of second derivatie of a one-dimenisonal function(three point-method)
    % der2 = (fun(x+h) + fun(x-h) - 2*fun(x))/h^2;
    
    % computation of second derivatie of a one-dimensional function (four-point method)
%     hess = (-fun(a-2*h) + 16*fun(a-h) - 30*fun(a) + 16*fun(a+h) - fun(a+2*h))/12/h^2;
    
    % computation of second derivative matrix (Hessian) of a two-dimensional function (three point method for pure derivative and central difference method for mixed derivatives)
    hess(1,1) = (fun(a) - 2*fun(a+[h,0]) + fun(a+[2*h,0]))/h^2;
    hess(2,2) = (fun(a) - 2*fun(a+[0,h]) + fun(a+[0,2*h]))/h^2;
    hess(1,2) = (fun(a+[h,h]) + fun(a+[-h,-h]) - fun(a+[-h,h]) - fun(a+[h,-h]))/4/h^2;
    hess(2,1) = hess(1,2);
    
end

