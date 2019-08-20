% Iterative solution to Green-Ampt Equation
% Inputs: K=Hydraulic conductivity (in cm s^-1)
%         t=time at which infiltration is to be computed (in seconds)
%         psi=suction head (cm)
%         delta_theta=porosity-intial mositure content
%         H0=intial water level height from the surface (cm)
% Outputs: Ft=cumulative infiltration at time t (cm)
% reference: Schmidt, J. (2002). A model for transient flow to a subsurface
%            tile drain. Chapter 2.

function Ft=falling_head_Green_Ampt_solution(K,psi,delta_theta,H0,t)

tol=0.001;   % convregence-error tolerance 
a=1-delta_theta;
b=(H0+psi)*delta_theta;
Kt=a*K*t;
a_by_b=a/b;

if length(t)==1 % if t is a scalar
    
    Ft_old=-inf;
    Ft_updated=sqrt(2*Kt/a_by_b);;

    while abs(Ft_updated-Ft_old)>tol
        Ft_old=Ft_updated;
        Ft_updated=Kt+log(1+a_by_b*Ft_old)/a_by_b;
    end

    Ft=Ft_updated;
    
else  % it t is a vector
    
    Ft_old=-inf*length(t);
    Ft_updated=sqrt(2*Kt/a_by_b);
    
    while max(abs(Ft_updated-Ft_old))>tol
        Ft_old=Ft_updated;
        Ft_updated=Kt+log(1+a_by_b*Ft_old)/a_by_b;
    end
    Ft=Ft_updated;
end

end