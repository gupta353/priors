%% referefence priors for Normal distribution scale

clear all
close all
clc

k=100;
m=1000;
mu=0;

pi_star=1;

sig2=1:10;

for i=1:length(sig2)
    for j=1:m
        sig2_tmp=sig2(i);
        samps=normrnd(mu,sig2_tmp,[k,1]);
        func=@(x)prod(integral_normpdf(samps,ones(length(samps),1)*mu,x));
        c(j)=integral(func,1,10);
        % integration using trapz method
%         x=0.1:0.1:20;
%         for l=1:length(x)
%             y(l)=func(x(l));
%         end
%         c(j)=trapz(x,y);
        expr=prod(normpdf(samps,mu,sig2_tmp))*pi_star/c(j);
        r(j)=log(expr);
    end
    PI(i)=exp(sum(r)/m);
end