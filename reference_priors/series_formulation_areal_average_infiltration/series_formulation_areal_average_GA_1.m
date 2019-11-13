% This rootine implements an alternate scheme of computing areal-average-
% infiltration using series formulation: the numerical integration sceheme
% is used to compute the second term in Eq. (13) of Govindaraju et al. (2001)
% Ref: Govindaraju et al. (2001); Areal Infiltration modeling...

clear all
close all
clc

mu_ks=1;
sigma_ks=0.86;
psi=699.64;
delta_theta=0.122;
r=[4.4,13.6,2.80,1.2];
obs=[4.4,9.94,2.76,1.152];
t=0:0.5:2;
n_ks=1000;

mu_y=log(mu_ks/(1+(sigma_ks/mu_ks)^2)^0.5); % mean of log(Ks)
sigma_y=(log(1+(sigma_ks/mu_ks)^2))^0.5;    % standard deviation of log(Ks)

% critical hydraulic conductivity during each time-interval
r_avg=cumsum(r.*t(2:end))./cumsum(t(2:end));
Kc=r./(psi*delta_theta./(r_avg.*t(2:end))+1);

% draw samples from log-normal distribution of Ks
y_samps=normrnd(mu_y,sigma_y,n_ks,1);     % random samples of y=ln(x)
K_samps=exp(y_samps);                     % random samples of Ks


%
F=zeros(n_ks,1);

for t_ind=2:length(t)
    
    t_tmp=t(t_ind);
    r_tmp=r(t_ind-1);
    Kc_tmp=Kc(t_ind-1);
    mu_omega(t_ind-1)=1-0.5*erfc((log(Kc_tmp)-mu_y)/sigma_y/2^0.5);
    
    % compute cumulative infiltration at the end of each time-step for
    % different values of Ks
    for K_ind=1:length(K_samps)
        K_tmp=K_samps(K_ind);
        F_tmp=F(K_ind,t_ind-1)+r_tmp*(t_tmp-t(t_ind-1));
        if K_tmp<Kc_tmp
            
            Fp(K_ind,t_ind-1)=K_tmp*psi*delta_theta./(r_tmp-K_tmp);
            tp(K_ind,t_ind-1)=t(t_ind-1)+...
                (Fp(K_ind,t_ind-1)-F(K_ind,t_ind-1))/r_tmp;
            
%             if tp(K_ind,t_ind-1)<t(t_ind-1)
%                 tp(K_ind,t_ind-1)=t(t_ind-1);
%                 Fp(K_ind,t_ind-1)=F(K_ind,t_ind-1);
%             end
            
            F(K_ind,t_ind)=min(F_tmp,Fp(K_ind,t_ind-1)+...
                sqrt(2*K_tmp*psi*delta_theta)*(t_tmp^0.5-tp(K_ind,t_ind-1)^0.5)+...
                2/3*K_tmp*(t_tmp-tp(K_ind,t_ind-1))+...
                1/18*(2*K_tmp^3/psi/delta_theta)^0.5*(t_tmp^1.5-tp(K_ind,t_ind-1)^1.5));
            
        else
            F(K_ind,t_ind)=F_tmp;
        end
    end
    
    % compute the second term in Eq. (13) of Govindaraju et al. (2001) at
    % time-step t_ind, for each value of K
    T2(:,t_ind-1)=K_samps.*(1+psi*delta_theta./F(:,t_ind));
    ind=find(K_samps>Kc_tmp);
    T2(ind,t_ind-1)=0;
end

I=r.*(1-mu_omega)+mean(T2);
scatter(obs,I,'filled');
hold on
xlim([0 max([obs,I])+1]); ylim([0 max([obs,I])+1]);
plot([0 max([obs,I])+1],[0 max([obs,I]+1)],'color','black')