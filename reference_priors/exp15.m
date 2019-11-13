% This routine computes the Bernardo prior over mean, scale, and shape
% parameters for one-dimensional generalized-Gaussian distribution (Gomez
% et al., 1998)

clear all
close all
clc

save_dir=['D:\Research\Thesis_work\Non_informative_priors\'...
    'matlab_codes\reference_priors\results\generalized_Gaussian'];

n=100;                  % number of sets to be drawn to be able to apply LLN
k=100;                  % number of samples in a set

% range of parameters
mu_range=-5:1:5;
phi_range=0.1:0.1:10;
beta_range=0.5:0.001:1;
mu_interval=1;
phi_interval=0.1;
beta_interval=0.001;

% samples of mu and beta for numerical integration
mu_int=-5:0.1:5;              % values of mu at which the integrand will be evaulated
beta_int=0.1:0.01:1;            % values of beta at which the integrand will be evaulated
phi_int=0.1:0.01:1;
delta_mu=0.1;                   
delta_beta=0.01;
delta_phi=0.01;

count=0;
E_log_asymp_post=zeros(length(mu_range),length(phi_range),...
    length(beta_range));
% PI=zeros(length(mu_range)*length(phi_range)*length(beta_range),4);
for mu_ind=6%:length(mu_range)
    for phi_ind=100%:length(phi_range)
        for beta_ind=1:length(beta_range)
        
             count=count+1;
             mu_tmp=mu_range(mu_ind);
             phi_tmp=phi_range(phi_ind);
             beta_tmp=beta_range(beta_ind);
             
             log_asymp_post=zeros(n,1);
            for n_ind=1:n           % loop over number of sets to be able to apply LLN
                
                samps=gennorm(k,mu_tmp,phi_tmp,beta_tmp);
                
                %%%%%%%%%%%%%%%%%% numerical integration %%%%%%%%%%%%%%%%%
                % first method:needs to checked
                %{
                SS_beta_min=zeros(length(beta_int),1);
                b=zeros(length(beta_int),1);
                for int_i=1:length(beta_int)            % loop over various values of beta to evaluate the integral
                    beta_int_tmp=beta_int(int_i);
                    SS_beta=zeros(length(mu_int),1);
                    for int_ii=1:length(mu_int)
                        
                        mu_int_tmp=mu_int(int_ii);
                        SS_beta(int_ii)=sum((abs(samps-mu_int_tmp)).^(2*beta_int_tmp));
                        
                    end
                    sum_scaled_SS_beta=(SS_beta/min(SS_beta)).^(-(k-1)/2/beta_int_tmp);
                    sum_scaled_SS_beta(1)=sum_scaled_SS_beta(1)/2; sum_scaled_SS_beta(end)=sum_scaled_SS_beta(end)/2;
                    sum_scaled_SS_beta=(sum(sum_scaled_SS_beta))^(-2*beta_int_tmp/(k-1));
                    SS_beta_min(int_i)=min(SS_beta)*sum_scaled_SS_beta;
                    b(int_i)=gammaln((k-1)/2/beta_int_tmp)-k*gammaln(1+1/2/beta_int_tmp)-...
                        log(2)/2/beta_int_tmp-log(beta_int_tmp)-...
                        ((k-1)/2/beta_int_tmp)*log(SS_beta_min(int_i));
                end
                b_max=max(b);
                exp_b=exp(b-b_max);
                exp_b(1)=exp_b(1)/2; exp_b(end)=exp_b(end)/2;
                sum_exp_b=sum(exp_b);
                log_I(n_ind)=log(delta_beta)+log(delta_mu)-(k-1)*log(2)+...
                    b_max+log(sum_exp_b); % log of integral
                %}
                
                % second method : needes to checked
                %{
                for int_i=1:length(beta_int)
                    for int_ii=1:length(mu_int)
                        
                        beta_int_tmp=beta_int(int_i);
                        mu_int_tmp=mu_int(int_ii);
                        SS_beta_int(int_i,int_ii)=sum((abs(samps-mu_int_tmp)).^(2*beta_int_tmp));
                        b(int_i,int_ii)=gammaln((k-1)/2/beta_int_tmp)-...
                            k*log(1+1/2/beta_int_tmp)-...
                            (k+1+1/2/beta_int_tmp)*log(2)-...
                            ((k-1)/2/beta_int_tmp)*log(SS_beta_int(int_i,int_ii))-...
                            log(beta_int_tmp);
                        
                    end
                end
                b_max=max(b(:));
                exp_b=exp(b-b_max);
                exp_b(1,:)=0.5*exp_b(1,:);
                exp_b(end,:)=0.5*exp_b(end,:);
                exp_b(:,1)=0.5*exp_b(:,1);
                exp_b(:,end)=0.5*exp_b(:,end);
                sum_exp_b=sum(exp_b(:));
                log_I(n_ind)=log(delta_mu)+log(delta_beta)+b_max+log(sum_exp_b);
                %}
                
                % third method:needs to be checked
                %{
                for int_i=1:length(beta_int)
                    for int_ii=1:length(mu_int)
                        for int_iii=1:length(phi_int)
                            
                            beta_int_tmp=beta_int(int_i);
                            mu_int_tmp=mu_int(int_ii);
                            phi_int_tmp=phi_int(int_iii);
                            
                            SS_beta_int=sum((abs(samps-mu_int_tmp)).^(2*beta_int_tmp));
                            b(int_i,int_ii)=-k*(log(phi_int_tmp)+gammaln(1+1/2/beta_int_tmp)+...
                                (1+1/2/beta_int_tmp)*log(2))-...
                                0.5*SS_beta_int/phi_int_tmp^(2*beta_int_tmp);
                        end
                    end
                end
                b_max=max(b(:));
                exp_b=exp(b-b_max);
                exp_b(1,:,:)=0.5*exp_b(1,:,:); exp_b(end,:,:)=0.5*exp_b(end,:,:);
                exp_b(:,1,:)=0.5*exp_b(:,1,:); exp_b(:,end,:)=0.5*exp_b(:,end,:);
                exp_b(:,:,1)=0.5*exp_b(:,:,1); exp_b(:,:,end)=0.5*exp_b(:,:,end);
                sum_exp_b=sum(exp_b(:));
                log_I(n_ind)=log(delta_mu)+log(delta_phi)+log(delta_beta)+b_max+log(sum_exp_b);
                %}
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % fourth method
                % integration is to be carried out only with respect to
                % beta
                for int_i=1:length(beta_int)
                    beta_int_tmp=beta_int(int_i);
                    SS_beta_int(beta_ind,int_i)=sum((abs(samps-mu_tmp)).^(2*beta_int_tmp));
                    b(beta_ind,int_i)=-k*log(phi_tmp)-k*gammaln(1+1/2/beta_int_tmp)-...
                        k*(1+1/2/beta_int_tmp)*log(2)-0.5*SS_beta_int(beta_ind,int_i)/phi_tmp^(2*beta_int_tmp);
                end
                %}
                b_max=max(b(beta_ind,:));
                exp_b=exp(b(beta_ind,:)-b_max);
                exp_b(1)=0.5*exp_b(1);  exp_b(end)=0.5*exp_b(end);
                sum_exp_b=sum(exp_b);
                log_I(beta_ind,n_ind)=b_max+log(sum_exp_b);
%                 
                % beta-sum-of-squares
                SS_beta(beta_ind,n_ind)=sum((abs(samps-mu_tmp)).^(2*beta_tmp));
                
                log_asymp_post(n_ind)=-k*(log(phi_tmp)+gammaln(1+1/2/beta_tmp)+(1+1/2/beta_tmp)*log(2))-...
                    0.5*(SS_beta(beta_ind,n_ind)/phi_tmp^(2*beta_tmp))-log_I(beta_ind,n_ind);
            end
%             
            E_log_asymp_post(mu_ind,phi_ind,beta_ind)=sum(log_asymp_post)/length(log_asymp_post);
            PI(count,:)=[mu_tmp,phi_tmp,beta_tmp,E_log_asymp_post(mu_ind,phi_ind,beta_ind)];
        end
    end
end

%% normalization of bernardo pdf
% if the prior over all the parameters is estimated simultaneously
%{
E_log_asymp_post=E_log_asymp_post-mean(E_log_asymp_post(:));
% %%%%%%%%%%% computation of normalization constant%%%%%%%%%%%%%%%%%%%%%%%%
exp_E_log_asymp_post=exp(E_log_asymp_post);
% multiply each fase of teh cube by 0.5
exp_E_log_asymp_post(:,:,1)=0.5*exp_E_log_asymp_post(:,:,1);
exp_E_log_asymp_post(:,:,end)=0.5*exp_E_log_asymp_post(:,:,end);
exp_E_log_asymp_post(:,1,:)=0.5*exp_E_log_asymp_post(:,1,:);
exp_E_log_asymp_post(:,end,:)=0.5*exp_E_log_asymp_post(:,end,:);
exp_E_log_asymp_post(1,:,:)=0.5*exp_E_log_asymp_post(1,:,:);
exp_E_log_asymp_post(end,:,:)=0.5*exp_E_log_asymp_post(end,:,:);

A=mu_interval*phi_interval*beta_interval*sum(exp_E_log_asymp_post(:));          % computation of total area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PI(:,4)=exp(PI(:,4)-mean(PI(:,4)))/A;
%}

% if the prior is computed only over beta 
PI(:,4)=PI(:,4)-mean(PI(:,4));
PI(:,4)=exp(PI(:,4));
A=trapz(PI(:,3),PI(:,4));
PI(:,4)=PI(:,4)/A;
plot(PI(:,3),PI(:,4));

% save data to textfile
fname='GG_prior_beta_mu=0_phi=10_11_06_2019';
save_filename=fullfile(save_dir,fname);
dlmwrite(save_filename,PI);