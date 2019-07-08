% numerical implementation of Beranrdo's algorithm to determine hydraulic-condutivity 
% assuming Green-Ampt infiltration model to be a true model
% \alpha=\frac_{\sum_{i=1}{k}(x_i-\mu)^2/}{2*\sigma^2}
% initial prior used: pi_star=\partialg/\partial\theta
% model: x=g(\theta)+\epsilon, \epsilon~N(0,\sigma^2), g is the Green-Ampt
% model
% references: Berger et al. (2009), The formal definition of reference
% priors
% Chow et al. (1988), Applied Hydrology

clear all
close all
clc

save_dir=['D:/Research/Thesis_work/Non_informative_priors',...
    '/matlab_codes/reference_priors/plots'];

k=1000;
m=10000;

% other known parameters of Green-Ampt equation
psi=16.68;          % (in cm)
delta_theta=0.340;  % change in moisture content
t=3600;   % time at which infiltration is computed (in s)
g=@(x)Green_Ampt_solution(x,psi,delta_theta,t);

kh=(0.01:0.01:1)/3600;          % hydraulic conductivity values at which prior is to be evaluated
sig2=0.01:0.01:1;                  % variance values at which prior is to be evaluated

for i=1:length(kh)
    for ii=1:length(sig2)
        
        for j=1:m
            K=kh(i);
            mu_tmp=g(K);                                    % mean of the distribution
            sig2_tmp=sig2(ii);
            samps=normrnd(mu_tmp,sqrt(sig2_tmp),[k,1]);     % k samples drawn from the distributiion
            avg=sum(samps)/k;                               % sample average of drawn samples
            s2=sum((samps-avg).^2)/k;                       % sample variance of drawn samples
            alpha=sum((samps-mu_tmp).^2)/2;                 % defined in the header
            der=Green_Ampt_der(K,psi,delta_theta,t);
            
            log_asymp_post(j)=(k-2)*log(s2)/2+(k-1)/2*log(k)-...        % log of asymptotic posterior at the given point in parameter space
                0.5*log(2*pi)-k*log(sig2_tmp)/2+(k-4)/2*log(2)-...
                gammaln(k/2-1)-alpha/sig2_tmp+log(der);
        end
        E_log_asymp_post(i,ii)=sum(log_asymp_post)/m;                   % expectation of log of asymptotic posterior at the given point in parameter space
    end
end

% Normalization of the prior
% a contant equal to mean of teh matrix E_log_asymp_post is sbtracted from
% E_log_asymp_post before taking exponential
PI=exp(E_log_asymp_post-sum(sum(E_log_asymp_post))/length(kh)/length(sig2));   % unnormalized density (ith row denotes ith 'kh' value and jth column denotes jth 'sig2' values)
for i=1:length(sig2)
    A1(i)=trapz(kh,PI(:,i));
end
A=trapz(sqrt(sig2),A1); % total area under the curve


PI=PI/A; % normalized density


% plots of prior density
%
% density of hydraulic conductivity for each value of standard
% deviation
plot(kh*3600,PI(:,1:floor(end/5)-1:end),'linewidth',2)
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',12,box)
xlabel('hydraulic conductivity (K_h, cm h^{-1})',...
    'fontname','arial','fontsize',12);
ylabel('prior density','fontname','arial','fontsize',12);

sname='GA_prior_K_intital_prior_1_07_07_2019';
save_filename=fullfile(save_dir,sname);
print(save_filename,'-r300','-djpeg');
clear box

% density of standard deviation for each value of 
% hydraulic conductivity
plot(sqrt(sig2),PI(1:floor(end/2)-1:end,:),'linewidth',2)
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',12,box)
xlabel('standard deviation (\sigma, cm) ',...
    'fontname','arial','fontsize',12);
ylabel('prior density','fontname','arial','fontsize',12);

sname='GA_prior_sigma_intital_prior_1_07_07_2019';
save_filename=fullfile(save_dir,sname);
print(save_filename,'-r300','-djpeg');
clear box

% joint density of hydraulic conductivity and standard deviation
[X,Y]=meshgrid(sqrt(sig2),3600*kh');
surf(X,Y,PI,'LineStyle','none','facealpha',0.8)
colorbar;
xlabel('standard deviation ( \sigma cm)',...
    'fontname','arial','fontsize',12);
ylabel(' K_h (cm h^{-1})',...
    'fontname','arial','fontsize',12);
zlabel('prior density','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12)

sname='GA_prior_joint_intital_prior_1_07_07_2019';
save_filename=fullfile(save_dir,sname);
print(save_filename,'-r300','-djpeg');
%}


% save the data
% reformat the data such that every row contains a point in multi-dimensional
% space and corresponding Bernardo's density
count=0;
for i=1:length(kh)
    for ii=1:length(sig2)
        count=count+1;
        save_data(count,:)=[kh(i),sqrt(sig2(ii)),PI(i,ii)];
    end
end
fname='prior_density_data_initial_prior_1_07_07_2019';
save_filename=fullfile(save_dir,fname);
dlmwrite(save_filename,save_data,',')





