% Iterative solution to Green-Ampt Equation
% Inputs: K=Hydraulic conductivity (in cm s^-1)
%         t=time at which infiltration is to be computed (in seconds)
%         psi=suction head (cm)
%         delta_theta=porosity-intial mositure content
% Outputs: Ft=cumulative infiltration at time t (in cm)

function Ft=Green_Ampt_solution(K,psi,delta_theta,t)

tol=0.000001;
a=psi*delta_theta;
Kt=K*t;

if length(t)==1 % if t is a scalar
    
    Ft_old=-inf;
    Ft_updated=K*t;

    while abs(Ft_updated-Ft_old)>tol
        Ft_old=Ft_updated;
        Ft_updated=Kt+a*log(1+Ft_old/a);
    end

    Ft=Ft_updated;
    
else  % it t is a vector
    
    Ft_old=-inf*length(t);
    Ft_updated=K*t;
    
    while max(abs(Ft_updated-Ft_old))>tol
        Ft_old=Ft_updated;
        Ft_updated=Kt+a*log(1+Ft_old./a);

    end
    Ft=Ft_updated;
end

end