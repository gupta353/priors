% numerical implementation of Bernardo's reference prior algorithm
% reference: Berger et al. (2009), The formal definition of reference
% priors
% intial prior used:
% pi_star=\beta^((\alpha-3)/2)*\gamma^0.5/(2*\pi)^0.5*\Gamma((\alpha-3)/2)*...
%          (1/\sigma^\alpha)*exp(-\beta/\sigma^2)*exp(-\gamma*(\mu-\mu_0)^2/2/\sigma^2),
% \alpha, \beta, \gamma, \mu_0 are the parameters of the intial prior
% distribution, \sigma denotes standard deviation and \mu denotes mean
% \beta_d=\beta+k*s^2/2+k*\gamma/(k+\gamma)*(mu_0-x_bar).^2,
% k denotes number of samples

clear all
close all
clc

k=1000;          % number of samples to be drawn in each set
m=1000;          % number of sets of samples to be drawn

mu=-10:10;       % mean values at which prior is to be evaluated
sig2=1:10;       % variance values at which prior is to be evaluated 

% parameters of iniitial prior as mentioned in header
gamma=0.01;     
alpha=5;
beta=1;
mu0=0;

for i=1:length(mu)
    for ii=1:length(sig2)
        
        for j=1:m
            mu_tmp=mu(i);
            sig2_tmp=sig2(ii);
            samps=normrnd(mu_tmp,sqrt(sig2_tmp),[k,1]);                 % k samples drawn from the distribution
            avg=sum(samps)/k;                                           % samples average of drawn sample
            s2=sum((samps-avg).^2)/k;                                   % sample variance of drawn samples
            beta_d=beta+k*s2/2+k*gamma/2/(k+gamma)*(mu0-avg)^2;         % explained in header
            
            log_asymp_post(j)=log(2)+0.5*log(gamma+k)+...               % log of asymptotic posterior distribution at the given point in parameter space
                ((k+alpha-2)/2)*log(beta_d)-0.5*log(2*pi)-...
                gammaln((alpha+k-2)/2)-(k+alpha)*log(sig2_tmp)/2-...
                k*s2/2/sig2_tmp-beta-k*(mu_tmp-avg)^2/2/sig2_tmp-...
                gamma*(mu_tmp-mu0)^2/2/sig2_tmp;
            
        end
        PI(i,ii)=exp(sum(log_asymp_post)/m);                            % evaluation of prior at the given point in parametric space
    end
end