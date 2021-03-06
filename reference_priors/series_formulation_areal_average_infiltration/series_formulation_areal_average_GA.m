% This routine computes areal avearge infiltration for a constant and
% spatially homogenous amount of rainfall over a time-interval; series
% formulation presented in Govindaraju et al. (2006) has been implemented
% with minor modifications (evaluation of equation 26 is different in this routine)
% inputs: psi = wetting front suction head in mm
% mu_ks = mean of Ks in mm/h
% sigma_ks = standard deviation of Ks in mm/h
% delta_theta = change in moisture content
% r = rainfall rate in mm/h
% t = time at which areal average cumulative infiltration is to be computed in hours

function [I,F]=series_formulation_areal_average_GA(mu_ks,sigma_ks,psi,delta_theta,r,t)
    
    r_backup=r;
    % check if r is a scalar, if yes convert it into a vector of length equal to
    % t with the same value of each element, otherwise compute average
    % rainfall rate upto each time-step
    if length(r)==1
        r=r*ones(length(t),1);
    else
        %     r=cumsum(r)./(1:length(r));
    end
    
    % computation of central moments of y=log(Ks)
    mu_y=log(mu_ks/(1+(sigma_ks/mu_ks)^2)^0.5); % mean
    sigma_y=(log(1+(sigma_ks/mu_ks)^2))^0.5;    % standard deviation
    
    % computation of critical hydraulic conductivity at different time-steps
    Kc=r./(psi*delta_theta./cumsum(r.*t)+1);
    
    %% computation of expected infiltration rate
    % mu_omega=logncdf(Kc,mu_ks,sigma_ks);       % check this expression
    mu_omega=1-0.5*erfc((log(Kc)-mu_y)/sigma_y/2^0.5);
    
    % draw samples from log-normal distribution of Ks
    % y_samps=normrnd(mu_y,sigma_y,100000,1);     % random samples of y=ln(x)
    % K_samps=exp(y_samps);                   % random samples of Ks
    %
    
    F(1)=0;
    for t_ind=1:length(t)
        %
        %     if t(t_ind)~=0
        % numerical implementation
        %     for K_ind=1:length(K_samps)
        %
        %         Ks_temp=K_samps(K_ind);
        %         Kc_temp=Kc(t_ind);
        %
        %         if Ks_temp<Kc_temp
        %             % compute cumulative infiltration each value of Ks in K_samps
        %
        %             tp=Ks_temp*psi*delta_theta/r./(r-Ks_temp);
        %             t_temp=t(t_ind)-tp;
        %             F(K_ind,t_ind)=r*tp...
        %                 +(2*Ks_temp*psi*delta_theta).^0.5*(t_temp^0.5)...
        %                 +2/3*Ks_temp*t_temp...
        %                 +1/18*(2*Ks_temp^3/psi/delta_theta)^0.5*(t_temp^1.5);
        %             T2(K_ind,t_ind)=Ks_temp*((psi*delta_theta/F(K_ind,t_ind))+1);
        %         else
        %             T2(K_ind,t_ind)=0;
        %         end
        %
        %     end
        
        % semi-analytrical implementation (from Govindaraju et al., 2001, p. 153)
        Kc_temp=Kc(t_ind);
        r_temp=r(t_ind);
        r_t=r_backup(t_ind);
        t_temp=t(t_ind);
        %         tp=(Kc_temp/2)*psi*delta_theta/r_temp/(r_t-Kc_temp/2);
        Fp=(Kc_temp/2)*psi*delta_theta/(r_t-Kc_temp/2);
        if t_ind>1
            tp=t(t_ind-1)-F(t_ind-1)/r_t+Fp/r_t;
        else
            tp=Fp/r_t;
        end
        t_temp=t(t_ind);
        
        %         if tp>0
        F(t_ind+1)=Fp...
            +(Kc_temp*psi*delta_theta).^0.5*(t_temp^0.5-tp^0.5)...
            +1/3*Kc_temp*(t_temp-tp)...
            +1/18*(Kc_temp^3/4/psi/delta_theta)^0.5*(t_temp^1.5-tp^1.5);
        
        T2(t_ind)=(1+psi*delta_theta/F(t_ind+1))*...
            exp(mu_y+sigma_y^2/2)*...
            log_normal_moment(Kc_temp,mu_y,sigma_y,1);
        %     else
        %         T2(t_ind)=0;
        %     end
        
        %     else
        %         T2(t_ind)=0;
        %     end
    end
    
    % E_T2=mean(T2);
    
    I=r.*(1-mu_omega)+T2;
end

function lnm=log_normal_moment(Kc,mu_y,sigma_y,order)
    
    lnm=1-0.5*erfc((log(Kc)-mu_y)/sigma_y/2^0.5-sigma_y*order/2^0.5);
    
end