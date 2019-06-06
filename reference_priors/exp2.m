%% referefence priors for Normal distribution location

clear all
close all
clc

k=100;
m=10000;
sig2=1;

pi_star=1;

mu=-10:10;

for i=1:length(mu)
    for j=1:m
        mu_tmp=mu(i);
        samps=normrnd(mu_tmp,sig2,[k,1]);
        func=@(x)prod(normpdf(samps,x,sig2));
%         c(j)=integral(func,1,10);
        % integration using trapz method
        x=-20:0.1:20;
        for l=1:length(x)
            y(l)=func(x(l));
        end
        c(j)=trapz(x,y);
        expr=prod(normpdf(samps,mu_tmp,sig2))*pi_star/c(j);
        r(j)=log(expr);
    end
    PI(i)=exp(sum(r)/m);
end