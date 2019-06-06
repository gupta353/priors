% parameters of a binomial distribution

clear all
close all
clc

k=100;
m=1000;
n=100;

pi_star=1;

theta=0:0.1:1;

for i=1:length(theta)
    for j=1:m
        theta_tmp=theta(i);
        samps=binornd(n,theta_tmp,[k,1]);
        func=@(x)prod(binopdf(samps,n,x));
%         c(j)=integral(func,1,10);
        % integration using trapz method
        x=0:0.01:1;
        for l=1:length(x)
            y(l)=func(x(l));
        end
        c(j)=trapz(x,y);
        expr=prod(binopdf(samps,n,theta_tmp))*pi_star/c(j);
        r(j)=log(expr);
    end
    PI(i)=exp(sum(r)/m);
end