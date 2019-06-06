% reference priors for Normal distribution location
% reference: Berger et al. (2009), The formal definition of reference
% priors

clear all
close all
clc

k=1000;          % number of samples in each set
m=100000;        % number of sets of samples
sig2=1;         % known variance 

pi_star=1;      % intial prior

mu=-10:0.1:10;      % values of mu at which prior function is to be evaluated

for i=1:length(mu)
    
    for j=1:m
        mu_tmp=mu(i);
        samps=normrnd(mu_tmp,sqrt(sig2),[k,1]); % k samples drawn from normal distribution
        avg=sum(samps)/k;                       % sample average of drawn samples
        log_asymp_post(j)=0.5*log(k)-0.5*log(2*pi)...
            -log(sig2)/2-k*(mu_tmp-avg)^2/2/sig2;
        
        % integration using trapz method
%         x=-20:0.1:20;
%         for l=1:length(x)
%             y(l)=func(x(l));
%         end
%         c(j)=trapz(x,y);
%         expr=prod(normpdf(samps,mu_tmp,sig2))*pi_star/c(j);
%         r(j)=log(expr);
    end
    
    PI(i)=exp(sum(log_asymp_post)/m);
end