% referefence priors for Normal distribution location and scale
% reference: Berger et al. (2009), The formal definition of reference
% priors
% intial prior used: pi_star=1;


clear all
close all
clc

k=1000;         % number of samples in each set
m=10000;        % number of sets of samples
pi_star=1;      % initial prior

mu=-10:10;      % mean values at which prior is to be evaluated
sig2=1:10;      % variance values at which prior is to be evaluated

for i=1:length(mu)
    
    for ii=1:length(sig2)
        for j=1:m
            mu_tmp=mu(i);
            sig2_tmp=sig2(ii);
            samps=normrnd(mu_tmp,sqrt(sig2_tmp),[k,1]);             % k samples from normal distribution
            avg=sum(samps)/k;                                       % sample average of drawn samples
            alpha=sum((samps-mu_tmp).^2)/2;                         % alpha as defined in the header
            s2=sum((samps-avg).^2)/k;                               % sample variance of drawn samples
            log_asymp_post(j)=(k-1)*log(s2)/2+(k/2-1/2)*log(k)-...  % log of asymptotic posterior evaulated at given points in parameter space
                1/2*log(2*pi)-(k/2-2)*log(2)-k*log(sig2_tmp)/2-...
                gammaln(k/2-1/2)-alpha/sig2_tmp;
        end
        PI(i,ii)=exp(sum(log_asymp_post)/m);                        % prior density at given point in parameter space
    end
    
end