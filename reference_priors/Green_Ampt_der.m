% Computation of derivative of Green-Ampt model w/r/t hydraulic conductivity
% Inputs: K=Hydraulic conductivity (in cm s^-1)
%         t=time at which infiltration is to be computed (in seconds)
%         psi=suction head (cm)
%         delta_theta=porosity-intial mositure content
% Outputs: der=Derivative of cumulative-infiltration w/r/t hydraulic
%          conductivity at input K

function der=Green_Ampt_der(K,psi,delta_theta,t)

Ft=Green_Ampt_solution(K,psi,delta_theta,t);
der=t.*(psi*delta_theta/Ft+1);

end