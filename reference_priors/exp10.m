% numerical implementation of Bernardo's reference prior algorithm to
% determine hydraulic conductivity assuming Green-Ampt to be true model
% reference: Berger et al. (2009), The formal definition of reference
% priors
% Chow et al. (1988), Applied Hydrology
% intial prior used:
% pi_star=\beta^((\alpha-3)/2)*\gamma^0.5/(2*\pi)^0.5*\Gamma((\alpha-3)/2)*...
%          (1/\sigma^\alpha)*exp(-\beta/\sigma^2)*exp(-\gamma*(g(\theta)-\mu_0)^2/2/\sigma^2)*...
%           \partialg(\theta)/\partial\theta,
% \alpha, \beta, \gamma, \mu_0 are the parameters of the intial prior
% distribution, \sigma denotes standard deviation and \mu denotes mean
% \beta_d=\beta+k*s^2/2+k*\gamma/(k+\gamma)*(mu_0-x_bar).^2,
% k denotes number of samples

clear all
close all
clc

save_dir=['D:\Research\Thesis_work\Non_informative_priors',...
    '\matlab_codes\reference_priors\plots'];

k=1000;          % number of samples to be drawn in each set
m=10000;          % number of sets of samples to be drawn

% other known parameters of Green-Ampt equation
psi=16.68;          % (in cm)
delta_theta=0.340;  % change in moisture content
t=3600;   % time at which infiltration is computed (in s)
g=@(x)Green_Ampt_solution(x,psi,delta_theta,t);


kh=(1:0.5:30)/3600;      % hydraulic conductivity values at which prior is to be evaluated
sig2=1:100;       % variance values at which prior is to be evaluated 

% parameters of iniitial prior as mentioned in header
gamma=0.01;     
alpha=5;
beta=0.0001;
mu0=10;

for i=1:length(kh)
    for ii=1:length(sig2)
        
        for j=1:m
            mu_tmp=g(kh(i));
            sig2_tmp=sig2(ii);
            samps=normrnd(mu_tmp,sqrt(sig2_tmp),[k,1]);                 % k samples drawn from the distribution
            avg=sum(samps)/k;                                           % samples average of drawn sample
            s2=sum((samps-avg).^2)/k;                                   % sample variance of drawn samples
            beta_d=beta+k*s2/2+k*gamma/2/(k+gamma)*(mu0-avg)^2;         % explained in header
            der=Green_Ampt_der(kh(i),psi,delta_theta,t);
            
            log_asymp_post(j)=log(2)+0.5*log(gamma+k)+...               % log of asymptotic posterior distribution at the given point in parameter space
                ((k+alpha-2)/2)*log(beta_d)-0.5*log(2*pi)-...
                gammaln((alpha+k-2)/2)-(k+alpha)*log(sig2_tmp)/2-...
                k*s2/2/sig2_tmp-beta-k*(mu_tmp-avg)^2/2/sig2_tmp-...
                gamma*(mu_tmp-mu0)^2/2/sig2_tmp+log(der);
            
        end
        E_log_asymp_post(i,ii)=sum(log_asymp_post)/m;                            % expectation of log of asymptotic posterior at the given point in parametric space
    end
end

% unnormalized density (ith row denotes ith 'kh' value and jth column denotes jth 'sig2' values)
PI=exp(E_log_asymp_post-sum(sum(E_log_asymp_post))/length(sig2)/length(kh));     

for i=1:length(sig2)
    A1(i)=trapz(kh,PI(:,i));                                                            
end
A=trapz(sqrt(sig2),A1);                                                          % total area under the unnormalized prior curve

PI=PI/A;                                                                         % normalized prior

% plots of prior density
%
% density of hydraulic conductivity for each value of standard
% deviation
plot(kh*3600,PI(:,1:floor(end/5)-1:end),'linewidth',2)
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',12,box)
xlabel('hydraulic conductivity (K, cm h^{-1})',...
    'fontname','arial','fontsize',12);
ylabel('prior density','fontname','arial','fontsize',12);
clear box
%
sname='GA_prior_K_intital_prior_2';
save_filename=fullfile(save_dir,sname);
print(save_filename,'-r300','-djpeg');
%}
% density of standard deviation for each value of 
% hydraulic conductivity
plot(sqrt(sig2),PI(1:floor(end/2)-1:end,:),'linewidth',2)
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',12,box)
xlabel('standard deviation (\sigma)',...
    'fontname','arial','fontsize',12);
ylabel('prior density','fontname','arial','fontsize',12);
%
sname='GA_prior_sigma_intital_prior_2';
save_filename=fullfile(save_dir,sname);
print(save_filename,'-r300','-djpeg');

%}
% joint density of hydraulic conductivity and standard deviation
[X,Y]=meshgrid(sqrt(sig2),3600*kh');
surf(X,Y,PI,'LineStyle','none','facealpha',0.8)
colorbar;
xlabel('standard deviation ( \sigma)',...
    'fontname','arial','fontsize',12);
ylabel(' K (cm h^{-1})',...
    'fontname','arial','fontsize',12);
zlabel('prior density','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12)
clear box
%
sname='GA_prior_joint_intital_prior_2';
save_filename=fullfile(save_dir,sname);
print(save_filename,'-r300','-djpeg');
%}

% save the data
% reformat the data such that every row contains point in multi-dimensional
% space and corresponding Bernardo's density
count=0;
for i=1:length(kh)
    for ii=1:length(sig2)
        count=count+1;
        save_data(count,:)=[kh(i),sqrt(sig2(ii)),PI(i,ii)];
    end
end
fname='prior_density_data_initial_prior_2';
save_filename=fullfile(save_dir,fname);
dlmwrite(save_filename,save_data,',')