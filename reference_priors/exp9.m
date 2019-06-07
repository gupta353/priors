% numerical implementation of Beranrdo's algorithm to determine hydraulic-condutivity 
% assuming Green-Ampt infiltration model to be a true model
% \alpha=\frac_{\sum_{i=1}{k}(x_i-\mu)^2/}{2*\sigma^2}
% initial prior used: pi_star=\partialg/\partial\theta
% model: x=g(\theta)+\epsilon, \epsilon~N(0,\sigma^2), g is the Green-Ampt
% model
% references: Berger et al. (2009), The formal definition of reference
% priors
% Chow et al. (1988), Applied Hydrology

clear all
close all
clc

k=1000;
m=10000;

% other known parameters of Green-Ampt equation
psi=16.68;          % (in cm)
delta_theta=0.340;  % change in moisture content
t=3600;   % time at which infiltration is computed (in s)
g=@(x)Green_Ampt_solution(x,psi,delta_theta,t);

kh=(1:0.5:30)/3600;          % hydraulic conductivity values at which prior is to be evaluated
sig2=[1,10];                  % variance values at which prior is to be evaluated

for i=1:length(kh)
    for ii=1:length(sig2)
        
        for j=1:m
            K=kh(i);
            mu_tmp=g(K);                                    % mean of the distribution
            sig2_tmp=sig2(ii);
            samps=normrnd(mu_tmp,sqrt(sig2_tmp),[k,1]);     % k samples drawn from the distributiion
            avg=sum(samps)/k;                               % sample average of drawn samples
            s2=sum((samps-avg).^2)/k;                       % sample variance of drawn samples
            alpha=sum((samps-mu_tmp).^2)/2;                 % defined in the header
            der=Green_Ampt_der(K,psi,delta_theta,t);
            
            log_asymp_post(j)=(k-2)*log(s2)/2+(k-1)/2*log(k)-...        % log of asymptotic posterior at the given point in parameter space
                0.5*log(2*pi)-k*log(sig2_tmp)/2+(k-4)/2*log(2)-...
                gammaln(k/2-1)-alpha/sig2_tmp+log(der);
        end
        PI(i,ii)=exp(sum(log_asymp_post)/m);                            % evaluation of the prior at the given point in parameter space
    end
end
