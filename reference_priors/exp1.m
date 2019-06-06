% Numerical computation of referefence priors for Normal distribution scale
% (standard deviation)
% reference: Berger et al. (2009)
% alpha=frac_{\sum_{i=1}{k}(x_i-\mu)}{2}
% initial prior, \pi^*=1
% 

clear all
close all
clc

k=1000;          % number of samples sets to be drawn
m=10000;         % number of samples in each set
mu=10;           % known mean of the normal distribution

pi_star=1;      % inital prior

sig2=1:0.1:10;      % vriance values at which prior has to be computed

for i=1:length(sig2)
    
    for j=1:m
        sig2_tmp=sig2(i);
        samps=normrnd(mu,sqrt(sig2_tmp),[k,1]);     % samples from normal distribution
        alpha=sum((samps-mu).^2)/2;                 % alpha defined in header
        log_asymp_post(j)=(k/2-1/2)*log(alpha)-...  % log of asymptotic posterior
            k*log(sig2_tmp)/2-gammaln(k/2-1/2)-alpha/sig2_tmp;
        
    end
    PI(i)=exp(mean(log_asymp_post));    % prior distribution function
    
end