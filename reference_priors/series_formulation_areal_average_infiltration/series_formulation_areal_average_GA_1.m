% This rootine implements an alternate scheme of computing areal-average-
% infiltration using series formulation: the numerical integration sceheme
% is used to compute the second term in Eq. (13) of Govindaraju et al.
% (2001) Ref: Govindaraju et al. (2001); Areal Infiltration modeling...
% Distribution of hydraulic conductivity is assumed to be lognormal inputs:
% mu_ks = mean parameter of the lognormal distribution (in mm h^{-1})
%         sigma_ks = standard deviation of the lognormal distribution (in
%         mm h^{-1}) psi = water front suction head (in mm) delta_moisture
%         = change in moisture content as the water-front penetrates the
%         soil n_ks = number of samples of hydraulic conductivity to be
%         drawn for numerical extimation of expectation r = vector of
%         rainfall rate during each time interval (in mm h^{-1}) t = time
%         steps between which rainfall rate is available (in h)
% output: I=infiltration rate during each time-interval (in mm h^{-1})

function [I]=series_formulation_areal_average_GA_1(mu_ks,sigma_ks,psi,delta_theta,n_ks,r,t)
    
    persistent psi_delta_theta r_avg Kc normal_standard
    
    if isempty(psi_delta_theta)
        psi_delta_theta=psi*delta_theta;
        % critical hydraulic conductivity during each time-interval
        r_avg=cumsum(r.*t(2:end))./cumsum(t(2:end));
        Kc=r./(psi_delta_theta./(r_avg.*t(2:end))+1);
        normal_standard=normrnd(0,1,n_ks,1);        % random samples of standard normal distribution
    end
    
    mu_y=log(mu_ks/(1+(sigma_ks/mu_ks)^2)^0.5); % mean of log(Ks)
    sigma_y=(log(1+(sigma_ks/mu_ks)^2))^0.5;    % standard deviation of log(Ks)
    
    % fraction of cells above critical hydraulic conductivity at each time=
    mu_omega=1-0.5*erfc((log(Kc)-mu_y)./sigma_y/2^0.5);
    
    % draw samples from log-normal distribution of Ks
    % y_samps1=normrnd(mu_y,sigma_y,n_ks,1);     % random samples of y=ln(x)
    y_samps=mu_y+sigma_y*normal_standard;
    K_samps=exp(y_samps);                     % random samples of Ks
    K_samps_psi_delta_theta=(2*K_samps*psi_delta_theta).^0.5;
    K_samps_cu_by_delta_theta=1/18*(2*K_samps.^3/psi_delta_theta).^0.5;
    K_samps_term_2=2/3*K_samps;
    %
    F=zeros(n_ks,length(t)-1);            % cumulative infiltration at each value of hydraulic conductivity, at the end of each time-interval
    EI_np=zeros(n_ks,length(t)-1);           %
    Fp=zeros(n_ks,1);
    tp=zeros(n_ks,1);
    
    for t_ind=2:length(t)       % loop for each time-interval
        
        t_tmp=t(t_ind);
        r_tmp=r(t_ind-1);
        Kc_tmp=Kc(t_ind-1);
        
        rainfall_depth_tmp=r_tmp*(t_tmp-t(t_ind-1));
        
        ind=K_samps<Kc_tmp;
        F_tmp=F(:,t_ind-1)+rainfall_depth_tmp;
        Fp(ind)=K_samps(ind)*psi_delta_theta./(r_tmp-K_samps(ind));
        tp(ind)=t(t_ind-1)+(Fp(ind)-F(ind,t_ind-1))/r_tmp;
        
        F(ind,t_ind)=Fp(ind)+...
            K_samps_psi_delta_theta(ind).*(t_tmp^0.5-tp(ind).^0.5)+...
            K_samps_term_2(ind).*(t_tmp-tp(ind))+...
            K_samps_cu_by_delta_theta(ind).*(t_tmp^1.5-tp(ind).^1.5);
        
        ind_zero=K_samps>Kc_tmp;
        F(ind_zero)=F_tmp(ind_zero);
        
        EI_np(:,t_ind-1)=K_samps.*(1+psi_delta_theta./F(:,t_ind));
        ind_tmp=Fp>F_tmp;
        F(ind_tmp,t_ind)=F_tmp(ind_tmp);
        EI_np(ind_tmp,t_ind-1)=r_tmp;
        EI_np(ind_zero,t_ind-1)=0;
        
        
        % same algorithm implemented in for loop: for readability
        %{
     t_tmp=t(t_ind);
    r_tmp=r(t_ind-1);
    Kc_tmp=Kc(t_ind-1);
    indices_r=[];
    
    % compute cumulative infiltration at the end of each time-step for
    % different values of Ks
    for K_ind=1:length(K_samps)     % loop for each value of K in K_samps
        K_tmp=K_samps(K_ind);
        F_tmp=F(K_ind,t_ind-1)+r_tmp*(t_tmp-t(t_ind-1));
        if K_tmp<Kc_tmp
            
            % cumulative infiltration at ponding and time to ponding
            Fp(K_ind,t_ind-1)=K_tmp*psi_delta_theta/(r_tmp-K_tmp);
            tp(K_ind,t_ind-1)=t(t_ind-1)+...
                (Fp(K_ind,t_ind-1)-F(K_ind,t_ind-1))/r_tmp;
                       
            % actual cumulative infiltration at the end of each
            % time-interval
            F(K_ind,t_ind)=Fp(K_ind,t_ind-1)+...
                (2*K_tmp*psi_delta_theta)^0.5*(t_tmp^0.5-tp(K_ind,t_ind-1)^0.5)+...
                2/3*K_tmp*(t_tmp-tp(K_ind,t_ind-1))+...
                1/18*(2*K_tmp^3/psi_delta_theta)^0.5*(t_tmp^1.5-tp(K_ind,t_ind-1)^1.5);
            
            if F(K_ind,t_ind)>F_tmp
                F(K_ind,t_ind)=F_tmp;
                indcices_r=[indices_r;K_ind];  % check
            end
        else
            F(K_ind,t_ind)=F_tmp;
        end
    end
    
    % compute the second term in Eq. (13) of Govindaraju et al. (2001) at
    % time-step t_ind, for each value of K
    T2(:,t_ind-1)=K_samps.*(1+psi_delta_theta./F(:,t_ind));
    ind=K_samps>Kc_tmp;
    T2(ind,t_ind-1)=0;
    T2(indices_r,t_ind-1)=r_tmp;        % check
        %}
        
    end
    
    I=r.*(1-mu_omega)+sum(EI_np)/n_ks;     % Infiltration rate during each time-interval
end