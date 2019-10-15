% This routine computes areal avearge infiltration for a constant and
% spatially homogenous amount of rainfall over a time-interval; series
% formulation presented in Govindaraju et al. (2006) has been implemented
% with minor modifications (evaluation of equation 26 is different in this routine)

clear all
close all
clc

psi=1116;               % wetting front suction head in mm
mu_ks=4.17;             % mean of Ks in mm/h
sigma_ks=4.17;           % Standard deviation of Ks in mm/h
delta_theta=0.1666;     % change in moisture content

r=10;                   % rainfall rate in mm/h

t=0:0.5:10;            % time at which areal average cumulative infiltration is to be computed in hours

% computation of central moments of y=log(Ks)
mu_y=log(mu_ks/(1+(sigma_ks/mu_ks)^2)^0.5);
sigma_y=(log(1+(sigma_ks/mu_ks)^2))^0.5;

% computation of critical hydraulic conductivity at different time-steps
Kc=r^2*t./(psi*delta_theta+r*t);

%% computation of expected infiltration rate
% mu_omega=logncdf(Kc,mu_ks,sigma_ks);       % check this expression
mu_omega=1-0.5*erfc((log(Kc)-mu_y)/sigma_y/2^0.5);

% draw samples from log-normal distribution of Ks
y_samps=normrnd(mu_y,sigma_y,100000,1);     % random samples of y=ln(x)
K_samps=exp(y_samps);                   % random samples of Ks
%
% T2(1)=0;
for t_ind=2:length(t)
    %
    % numerical implementation
    for K_ind=1:length(K_samps)
        
        Ks_temp=K_samps(K_ind);
        Kc_temp=Kc(t_ind);
        
        if Ks_temp<Kc_temp
            % compute cumulative infiltration each value of Ks in K_samps
            
            tp=Ks_temp*psi*delta_theta/r./(r-Ks_temp);
            t_temp=t(t_ind)-tp;
            F(K_ind,t_ind)=r*tp...
                +(2*Ks_temp*psi*delta_theta).^0.5*(t_temp^0.5)...
                +2/3*Ks_temp*t_temp...
                +1/18*(2*Ks_temp^3/psi/delta_theta)^0.5*(t_temp^1.5);
            T2(K_ind,t_ind)=Ks_temp*((psi*delta_theta/F(K_ind,t_ind))+1);
        else
            T2(K_ind,t_ind)=0;
        end
        
    end
    
    % semi-analytrical implementation (from Govindaraju et al., 2001)
%         Kc_temp=Kc(t_ind);
%         tp=Kc_temp/2*psi*delta_theta/r/(r-Kc_temp/2);
%         t_temp=t(t_ind)-tp;
%     
%         F(t_ind)=r*tp...
%                  +(Kc_temp*psi*delta_theta).^0.5*(t(t_ind)^0.5-tp^0.5)...
%                  +1/3*Kc_temp*t_temp...
%                  +1/18*(Kc_temp^3/4/psi/delta_theta)^0.5*(t(t_ind)^1.5-tp^1.5);
%     
%         T2(t_ind)=(1+psi*delta_theta/F(t_ind))*...
%             exp(mu_y+sigma_y^2/2)*...
%             (1-erfc((log(Kc_temp)-mu_y)/sigma_y/2^0.5-sigma_y/2^0.5)/2);
end

E_T2=mean(T2);

I=r*(1-mu_omega)+E_T2;