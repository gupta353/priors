% This routine computes Bernardo prior for mean and standard deviation parameetrs of
% areal average infiltraion model (Govindaraju et al., 2001)
% Likelihood function is assumed to be multi-variate Gaussain with mean equal to the
% area-average-infiltration model output, and homoscedastic varaince and
% non-diagonal elements of covariance matrix equal to zero

clear all
close all
clc

n=100;          % number of sets to be drawn
k=1000;         % number of samples in each set
na=100;         % number of samples drawn for numerical approximation of the integral

% parameters of the model
psi=1116;               % wetting front suction head in mm
delta_theta=0.1663;     % change in moisture content
r=10;                   % rainfall rate in mm h^{-1}

% range of parameters
mu_ks_min=0;        % min values of mu_ks
mu_ks_max=10;           % max values of mu_ks
sigma_ks_min=0.01;      % min values of sigma_ks
sigma_ks_max=10;        % max values of sigma_ks
sigma_err_min=0.1;     % min values of sigma_err
sigma_err_max=1;       % max values of sigma_err
t=0:0.5:10;             % time at which cumulative infiltraion is available

% list of parameters at which prior is to be computed
mu_ks_list=0.01:0.1:10;          % list of values of mu_ks
sigma_ks_list=0.1:0.1:10;        % list of  values of sigma_ks
sigma_err_list=0.1:0.1:1;        % list of values of sigma_err

count=0;
for mu_ks_ind=1:length(mu_ks_list)                      % loop of mu_ks values
    for sigma_ks_ind=1%:length(sigma_ks_list)            % loop for sigma_ks values
        for sigma_err_ind=1%:length(sigma_err_list)      % loop for sigma_err values
            
            count=count+1;
            mu_ks_tmp=mu_ks_list(mu_ks_ind);
            sigma_ks_tmp=sigma_ks_list(sigma_ks_ind);
            sigma_err_tmp=sigma_err_list(sigma_err_ind);
            
            for i=1:n                                   % loop for each set of samples
                
                % drawn k samples from the likelihood function
                mu_tmp=series_formulation_areal_average_GA...
                    (mu_ks_tmp,sigma_ks_tmp,psi,delta_theta,r,t);
                samps=mvnrnd(mu_tmp,sigma_err_tmp^2*eye(length(t)),k);
                
                %%%%%%%%% numerical integration to evaluate normalizing constant %%%%%%%%%
                par_samps=[unifrnd(mu_ks_min,mu_ks_max,na,1),unifrnd(sigma_ks_min,sigma_ks_max,na,1)];
                for par_ind=1:size(par_samps,1)
                    
                    infil_rate=series_formulation_areal_average_GA...
                    (par_samps(par_ind,1),par_samps(par_ind,2),psi,delta_theta,r,t);
                    diff=bsxfun(@minus,samps,infil_rate);
                    l(par_ind)=sum(sum(diff.*diff));
                    
                end
                l_min=min(l);
                
                diff1=bsxfun(@minus,samps,mu_tmp);
                diff1=sum(sum(diff1.*diff1));
                
                num_int_sum=sum((l/l_min).^(-(length(t)*k-1)/2));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  
                % log of asymptotic distribution
                log_post(i)=-length(t)*k*log(sigma_err_tmp)-...
                    diff1/2/sigma_err_tmp^2+...
                    (length(t)*k-1)/2*log(l_min)-...
                    log(num_int_sum);
            end
            
            E_log_post(count)=mean(log_post);
            PI(count,:)=[mu_ks_tmp,sigma_ks_tmp,sigma_err_tmp,E_log_post(count)];
            log_post=[];
        end
    end
end
