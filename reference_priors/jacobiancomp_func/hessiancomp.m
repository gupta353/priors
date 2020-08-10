% Numerical computation of second derivative of a function

function hess = hessiancomp(fun,a,h)
      
    % computation of second derivatie (three point-method)
    % der2 = (fun(x+h) + fun(x-h) - 2*fun(x))/h^2;
    
    % computation of second derivatie (four-point method)
    hess = (-fun(a-2*h) + 16*fun(a-h) - 30*fun(a) + 16*fun(a+h) - fun(a+2*h))/12/h^2;
end