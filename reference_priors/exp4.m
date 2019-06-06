%% referefence priors for Normal distribution location and scale

clear all
close all
clc

k=100;
m=100;
pi_star=1;

mu=-10:10;
sig2=1:10;

for i=1:length(mu)
    for ii=1:length(sig2)
        for j=1:m
            mu_tmp=mu(i);
            sig2_tmp=sig2(ii);
            samps=normrnd(mu_tmp,sig2_tmp,[k,1]);
            func=@(x,y)prod(normpdf(samps,x,y));
    %         c(j)=integral(func,1,10);
            % integration using trapz method
            x=-20:20;
            y=1:0.1:10;
            for a=1:length(x)
                for b=1:length(y)
                    z(a,b)=func(x(a),y(b));
                end
            end

            c=trapz(x,z,1);
            c(j)=trapz(y,c);
            expr=prod(normpdf(samps,mu_tmp,sig2_tmp))*pi_star/c(j);
            r(j)=log(expr);
        end
        PI(i,ii)=exp(sum(r)/m);
    end
end