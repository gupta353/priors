% This routine computes the log-likelihood of the infiltration data given
% the error term is zero-mean Gaussian with homoscedastic variance and the
% infiltration is modeled by Green-Ampt equation (Chow et al., 1988; Chapter
% 4)
% References: Chow, V. T., Maidment, D. R., & Mays, L. W. (1988). Chapter 4
% Subsurface flow. Applied Hydrology. New York: McGraw-Hill, 99-174.

function L=Gaussloglikeli(theta)

y0=10;              % observed cumulative infiltration value
kh=theta(1);        % hydraulic conductivity in cm s^-1
sig=theta(2);       % standard-deviation

% other known parameters of Green-Ampt equation
psi=16.68;          % (in cm)
delta_theta=0.340;  % change in moisture content
t=3600;   % time at which infiltration is computed (in s)
g=@(x)Green_Ampt_solution(x,psi,delta_theta,t);


mu=g(kh);           % computation of infiltration at kh

L=-0.5*log(2*pi)-log(sig)-(y0-mu)^2/2/sig^2; % computation of log-likelihood

end