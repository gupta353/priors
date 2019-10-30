% This routine computes the Bernardo prior over mean, scale, and shape
% parameters for one-dimensional generalized-Gaussian distribution (Gomez
% et al., 1998)

clear all
close all
clc

n=100;                  % number of sets to be drawn to be able to apply LLN
k=100;                  % number of samples in a set

% range of parameters
mu_range=-10:1:10;
phi_range=0.1:0.1:1;
beta_range=0.1:0.1:1;

% samples of mu and beta for numerical integration
mu_int=-10:0.1:10;              % values of mu at which the integrand will be evaulated
beta_int=0.1:0.01:1;            % values of beta at which the integrand will be evaulated
delta_mu=0.1;                   
delta_beta=0.01;

count=0;
for mu_ind=1%:length(mu_range)
    for phi_ind=1%:length(phi_range)
        for beta_ind=1:length(beta_range)
        
             count=count+1;
             mu_tmp=mu_range(mu_ind);
             phi_tmp=phi_range(phi_ind);
             beta_tmp=beta_range(beta_ind);
             
            for n_ind=1:n           % loop over number of sets to be able to apply LLN
                
                samps=gennorm(k,mu_tmp,phi_tmp,beta_tmp);
                
                %%%%%%%%%%%%%%%%%% numerical integration %%%%%%%%%%%%%%%%%
                for int_i=1:length(beta_int)            % loop over various values of beta to evaluate the integral
                    beta_int_tmp=beta_int(int_i);
                    for int_ii=1:length(mu_int)
                        
                        mu_int_tmp=mu_int(int_ii);
                        SS_beta(int_ii)=sum((abs(samps-mu_int_tmp)).^(2*beta_int_tmp));
                        
                    end
                    sum_scaled_SS_beta=(SS_beta/min(SS_beta)).^(-(k-1)/2/beta_int_tmp);
                    sum_scaled_SS_beta(1)=sum_scaled_SS_beta(1)/2; sum_scaled_SS_beta(end)=sum_scaled_SS_beta(end)/2;
                    sum_scaled_SS_beta=(sum(sum_scaled_SS_beta))^(-2*beta_int_tmp/(k-1));
                    SS_beta_min(int_i)=min(SS_beta)*sum_scaled_SS_beta;
                    b(int_i)=gammaln((k-1)/2/beta_int_tmp+2)-k*gammaln(1+1/2/beta_int_tmp)-...
                        1/2/beta_int_tmp*log(2)-log(beta_int_tmp)-...
                        ((k-1)/2/beta_int_tmp)*log(SS_beta_min(int_i));
                end
                b_max=max(b);
                exp_b=exp(b-b_max);
                exp_b(1)=exp_b(1)/2; exp_b(end)=exp_b(end)/2;
                sum_exp_b=sum(exp_b);
                log_I=log(delta_beta)+log(delta_mu)-(k+1)*log(2)+...
                    b_max+log(sum_exp_b); % log of integral
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % beta-sum-of-squares
                SS_beta=sum(((samps-mu_tmp).^2).^beta_tmp);
                
                log_asymp_post(n_ind)=-k*(log(phi_tmp)+gammaln(1+1/2/beta_tmp)+(1+1/2/beta_tmp)*log(2))-...
                    0.5*SS_beta/phi_tmp^(2*beta_tmp)-log_I;
            end
            
            E_log_asymp_post=sum(log_asymp_post)/length(log_asymp_post);
            PI(count,:)=[mu_tmp,phi_tmp,beta_tmp,E_log_asymp_post];
        end
    end
end