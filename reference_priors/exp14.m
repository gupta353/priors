% This routine computes Bernardo prior for mean and standard deviation parameetrs of
% areal average infiltraion model (Govindaraju et al., 2001)
% Likelihood function is assumed to be multi-variate Gaussain with mean equal to the
% area-average-infiltration model output, and homoscedastic varaince and
% non-diagonal elements of covariance matrix equal to zero

clear all
close all
clc

n=2000;          % number of sets to be drawn
k=2000;         % number of samples in each set
na=1000;         % number of samples drawn for numerical approximation of the integral


% extract observed data from etxtfiles
direc=['D:/Research/Thesis_work/Non_informative_priors/'...
    'matlab_codes/reference_priors/data/experiment_data'];

theta_sat=0.36;                         % saturated moisture content
fname='val_8.txt';
filename=fullfile(direc,fname);
fid=fopen(filename,'r');
data=textscan(fid,'%s%s%s','delimiter','\t');
fclose(fid);
% extract different parameters and variables
data1=data{1};
theta_ini=strsplit(data1{2},'='); psi=strsplit(data1{3},'=');
theta_ini=str2double(theta_ini{2});     % intital moisture content
psi=str2double(psi{2});                 % wetting front suction head in mm
delta_theta=theta_sat-theta_ini;        % change in moisture content

t=data1(5:end);                         % time-steps between rainfall rate is available
t_interval=length(t)-1;
rainfall=data{2}(5:end);
runoff=data{3}(5:end);

t=cellfun(@str2num,t);
runoff=cellfun(@str2num,runoff);
delta_t=t(2:end)-t(1:end-1);
r=cellfun(@str2num,rainfall);
infil=r-runoff;
r=r(2:end)./delta_t;                    % rainfall rate in mm h^{-1}
obs=infil(2:end)./delta_t;              % observed infiltration rate in mm h^{-1}

r=r'; obs=obs'; t=t';

n_ks=100;                              % number of samples of hydraulic conductivity to be drawn
%                                        for numerical extimation of expectation
%{
profile on;
[I]=series_formulation_areal_average_GA_1(1,1,psi,delta_theta,n_ks,r,t);
scatter(obs,I,'filled');
xlim([0 max([obs,I])]); ylim([0 max([obs,I])]);
hold on
plot([0 max([I1,I])],[0 max([I1,I])],'color','black');
profile off; profile report
%}

% range of parameters
mu_ks_min=0;        % min values of mu_ks
mu_ks_max=2;           % max values of mu_ks
sigma_ks_min=0.1;      % min values of sigma_ks
sigma_ks_max=1;        % max values of sigma_ks
sigma_err_min=0.1;     % min values of sigma_err
sigma_err_max=5;       % max values of sigma_err
% t=0:0.5:10;             % time at which cumulative infiltraion is available


% list of parameters at which prior is to be computed
mu_ks_list=0.1:0.1:2;          % list of values of mu_ks
sigma_ks_list=0.1:0.1:1;        % list of  values of sigma_ks
sigma_err_list=0.1:0.1:5;        % list of values of sigma_err

E_log_post=zeros(length(mu_ks_list),length(sigma_ks_list),length(sigma_err_list));
PI=zeros(length(mu_ks_list)*length(sigma_ks_list)*length(sigma_err_list),4);
count=0;
for mu_ks_ind=10%:length(mu_ks_list)                      % loop of mu_ks values
    for sigma_ks_ind=1:length(sigma_ks_list)            % loop for sigma_ks values
        for sigma_err_ind=30%:length(sigma_err_list)      % loop for sigma_err values
            
            count=count+1;
            mu_ks_tmp=mu_ks_list(mu_ks_ind);
            sigma_ks_tmp=sigma_ks_list(sigma_ks_ind);
            sigma_err_tmp=sigma_err_list(sigma_err_ind);
            
            log_post=zeros(n,1);
            for i=1:n                                   % loop for each set of samples
                
                % drawn k samples from the likelihood function
                mu_tmp=series_formulation_areal_average_GA_1(mu_ks_tmp,sigma_ks_tmp,psi,delta_theta,n_ks,r,t);
                samps=mvnrnd(mu_tmp,sigma_err_tmp^2*eye(t_interval),k);
                
                %%%%%%%%% numerical integration to evaluate normalizing constant %%%%%%%%%
                
                % first method
                %
                par_samps=[unifrnd(mu_ks_min,mu_ks_max,na,1),unifrnd(sigma_ks_min,sigma_ks_max,na,1)];
                l=zeros(size(par_samps,1),1);
                parfor par_ind=1:size(par_samps,1)
                    
                    infil_rate=series_formulation_areal_average_GA_1(par_samps(...
                        par_ind,1),par_samps(par_ind,2),psi,delta_theta,n_ks,r,t);
                    diff=bsxfun(@minus,samps,infil_rate);
                    diff_square=diff.*diff;
                    l(par_ind)=sum(diff_square(:));
                    
                end
                l_min=min(l);
                
                diff1=bsxfun(@minus,samps,mu_tmp);
                diff1=sum(sum(diff1.*diff1));
                
                num_int_sum=sum((l/l_min).^(-(length(t)*k-1)/2));
                %}
                % second method
                %{
                mu_ks_int_list=mu_ks_min:(mu_ks_max-mu_ks_min)/na:mu_ks_max;
                sigma_ks_int_list=sigma_ks_min:(sigma_ks_max-sigma_ks_min)/na:sigma_ks_max;
                SS=zeros(na^2,1);
                for mu_ks_int_ind=1:length(mu_ks_int_list)
                    for sigma_ks_int_ind=1:length(sigma_ks_int_list)
                        infil_rate=series_formulation_areal_average_GA_1(mu_ks_tmp,sigma_ks_tmp,psi,delta_theta,n_ks,r,t);
                        diff=bsxfun(@minus,samps,infil_rate);
                        diff_square=diff.*diff;
                        SS(par_ind)=sum(diff_square(:));
                    end
                end
                %}
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % log of asymptotic distribution
                log_post(i)=-length(t)*k*log(sigma_err_tmp)-...
                    diff1/2/sigma_err_tmp^2+...
                    (length(t)*k-1)/2*log(l_min)-...
                    log(num_int_sum);
            end
            
            E_log_post(mu_ks_ind,sigma_ks_ind,sigma_err_ind)=sum(log_post)/n;
            PI(count,:)=[mu_ks_tmp,sigma_ks_tmp,sigma_err_tmp,E_log_post(mu_ks_ind,sigma_ks_ind,sigma_err_ind)];
            log_post=[];
        end
    end
end
