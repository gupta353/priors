% This routine computes induced prior over cumulative infiltration given a
% prior over hydraulic conductivity

clear all
close all
clc

save_dir=['D:/Research/Thesis_work/Non_informative_priors/'...
    'matlab_codes/reference_priors/'];

%% induced prior over cumulative infiltration in Green-Ampt model
%
%% draw samples of hydraulic conductivity from the prior
% samples from Bernardo prior
filename=['D:/Research/Thesis_work/Non_informative_priors/matlab_codes/'...
    'reference_priors/results/onwards_august_2019/'...
    'prior_density_data_initial_prior_2_sig2=1_08_15_2019'];
fid=fopen(filename,'r');
data=textscan(fid,'%f%f%f','delimiter',',');
fclose(fid);

kh=data{1};
kh_prior_pdf=data{3};

% samples of hydraulic conductivity
kh_samps=drsampsdens_ddim(kh,kh_prior_pdf,10000);

%% compute prior over cumulative infiltration obtained by Green-AMpt model
psi=16.68;          % (in cm)
delta_theta=0.340;  % change in moisture content
t=3600;   % time at which infiltration is computed (in s)
g=@(x)Green_Ampt_solution(x,psi,delta_theta,t);

for kh_ind=1:length(kh_samps)
    F(kh_ind)=g(kh_samps(kh_ind));
end

% histogram of induced prior
hist(F)
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',12,...
    'xlim',[0.4 4],box)
xlabel('F(t=1 hour) (cm)',...
    'fontname','arial','fontsize',12);
ylabel('Number of samples in the bin','fontname','arial','fontsize',12);

sname='induced_prior_cum_infiltration_10_30_2019.svg';
save_filename=fullfile(save_dir,'Manuscript_plots',sname);
% print(save_filename,'-r300','-djpeg');
plot2svg(save_filename);
clear box
%}

%% induced prior over output realizations of generalized Gaussian distribution
%{
mu=0;phi=1;
n_beta_uniform=10000;
n_beta_bernardo=10000;
% draw samples of beta from uniform distribtion
beta_samps=unifrnd(0.5,1,[n_beta_uniform,1]);
y_samps_uniform=[];
for beta_ind=1:length(beta_samps)
    y_samps_uniform=[y_samps_uniform;gennorm(100,mu,phi,beta_samps(beta_ind))];
end

% draw samples of beta from Bernardo prior
fname='GG_prior_beta_mu=0_phi=1_11_06_2019';
filename=fullfile(save_dir,'results','generalized_Gaussian',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%f%f%f%f','delimiter',',');
fclose(fid);
bernardo_prior=[data{3},data{4}];
beta_samps=drsampsdens_ddim(bernardo_prior(:,1),bernardo_prior(:,2),n_beta_bernardo);

y_samps_bernardo=[];
for beta_ind=1:length(beta_samps)
    y_samps_bernardo=[y_samps_bernardo;gennorm(100,mu,phi,beta_samps(beta_ind))];
end
[f1,x1]=ksdensity(y_samps_uniform);
[f2,x2]=ksdensity(y_samps_bernardo);

plot(x1,f1,'linewidth',2);
hold on
plot(x2,f2,'--','linewidth',2,'color','green');
xlabel('sample from generalized Gaussian',...
    'fontname','arial','fontsize',12);
ylabel('prior density',...
    'fontname','arial','fontsize',12);
legend({'Uniform prior on \beta','Bernardo prior on \beta'},...
    'fontname','arial','fontsize',12);
legend('boxoff');
box('on'); box.linewidth=2;
set(gca,'fontname','arial','fontsize',12,box);
ylim([-0.1 0.7])
clear box

sname='induced_prior_on_samples_from_generalized_Gaussian';
save_filename=fullfile(save_dir,sname);
print(save_filename,'-r300','-djpeg');
%}