% This routine computes the log-likelihood of the infiltration data given
% the error term is zero-mean Gaussian with homoscedastic variance and the
% infiltration is modeled by Green-Ampt equation (Chow et al., 1988; Chapter
% 4)
% inputs: theta=an array such that theta(1) is the value of hydraulic
% conductivity, and theta(2) is the standard deviation
% output: L=log-likelihood value
% References: Chow, V. T., Maidment, D. R., & Mays, L. W. (1988). Chapter 4
% Subsurface flow. Applied Hydrology. New York: McGraw-Hill, 99-174.

function L=Gaussloglikeli(theta)

global GLOBAL_DATA

y0=GLOBAL_DATA.y0;
N=length(y0);

% other known parameters of Green-Ampt equation
psi=50;                % (in cm)
delta_theta=0.004787;  % change in moisture content
t=1148;                % time at which infiltration is computed (in s)
H0=13.7;

kh=theta(1);        % hydraulic conductivity in cm s^-1
sig=theta(2);       % standard-deviation

mu=falling_head_Green_Ampt_solution(kh,psi,delta_theta,H0,t);           % computation of cumulative infiltration at kh
err=sum((y0-mu).^2);                                                    % sum of suquare error
L=-0.5*N*log(2*pi)-N*log(sig)-err/2/sig^2;                              % computation of log-likelihood

end