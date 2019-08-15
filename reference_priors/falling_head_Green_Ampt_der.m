% Computation of derivative of falling-head Green-Ampt model w/r/t hydraulic conductivity
% Inputs: K=Hydraulic conductivity (in cm s^-1)
%         t=time at which infiltration is to be computed (in seconds)
%         psi=suction head (cm)
%         delta_theta=porosity-intial mositure content
% Outputs: der=Derivative of cumulative-infiltration w/r/t hydraulic
%          conductivity at input K

function der=falling_head_Green_Ampt_der(K,psi,delta_theta,H0,t)

a=1-delta_theta;
b=(H0+psi)*delta_theta;

Ft=falling_head_Green_Ampt_solution(K,psi,delta_theta,H0,t);
der=t.*(b./Ft+a);

end