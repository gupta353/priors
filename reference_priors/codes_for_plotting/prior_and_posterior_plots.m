% This routine plots prior and posterior distribution for various
% cases-studies conducted in thi study
clear all
close all
clc

dir=['D:/Research/Thesis_work/Non_informative_priors'...
    '/matlab_codes/reference_priors'];

%% prior plots
%{
fname='FSI_bernardo_prior_02_03_2020';
filename=fullfile(dir,'results','field_scale_infiltration',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%f%f%f%f','delimiter',',');
fclose(fid);
bernardo_prior=[data{1},data{2},data{3},data{4}];              % factor of 3600 used to convert the K values from cm s^-1 to cm h^-1
mu_ks_range=unique(bernardo_prior(:,1));
sigma_ks_range=unique(bernardo_prior(:,2));
sigma_err_range=unique(bernardo_prior(:,3));
% beta_range=unique(bernardo_prior(:,3));
% sig_range=unique(bernardo_prior(:,2));

% density of hydraulic conductivity for each value of standard
% deviation
%{
sig_values=[0.774597000000000,...
    1,...
    1.140175000000000];
sig_values_legend=round(sig_values*1000)/1000;
colorspec={'b'};
line_pattern={'-','--',':'};
figure;
hold on
for i=1:length(sig_values)
    ind=find(bernardo_prior(:,2)==sig_values(i));
    kh=bernardo_prior(ind,1);
    dens=bernardo_prior(ind,3);
    plot(kh,dens,line_pattern{i},'linewidth',2,'color','black')
    
end

box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',14,'xlim',[-2 52],box)
xlabel('hydraulic conductivity (K_h, cm h^{-1})',...
    'fontname','arial','fontsize',14);
ylabel('prior density','fontname','arial','fontsize',14);
legend({['\sigma=',num2str(sig_values_legend(1)),char(32),'cm'],...
    ['\sigma=',num2str(sig_values_legend(2)),char(32),'cm'],...
    ['\sigma=',num2str(sig_values_legend(3)),char(32),'cm']},...
    'fontname','arial','fontsize',14,'Location','northeast');
legend('boxoff');
clear box

%
sname='FHGA_prior_kh_sig_09_17_2019';
save_filename=fullfile(dir,'plots','onwards_august_2019',sname);
print(save_filename,'-r300','-djpeg');
%}
% density of standard deviation for each value of 
% hydraulic conductivity
%{
kh_values=[5.000400000000000,...
    30.002399999999998,...
    50.000400000000000];          % in cm h^-1
kh_values_legend=round(kh_values*100)/100;
colorspec={'r','g','b'};
line_pattern={'-','--',':'};
figure;
hold on
for i=1:length(kh_values)
    ind=find(bernardo_prior(:,1)==kh_values(i));           % factor of 10^4 is used becuase of precision errors in the values of kh in bernardo_prior
    sig=bernardo_prior(ind,2);
    dens=bernardo_prior(ind,3);
    plot(sig,dens,line_pattern{i},'linewidth',2,'color','black')
    
end
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',14,'xlim',[0.30 1.5],box)
xlabel('standard deviation (\sigma, cm)',...
    'fontname','arial','fontsize',14);
ylabel('prior density','fontname','arial','fontsize',14);
legend({['K_h=',num2str(kh_values_legend(1)),char(32),'cm h^{-1}'],...
    ['K_h=',num2str(kh_values_legend(2)),char(32),'cm h^{-1}'],...
    ['K_h=',num2str(kh_values_legend(3)),char(32),'cm h^{-1}']},...
    'fontname','arial','fontsize',14);
legend('boxoff');
clear box

sname='FHGA_prior_kh_sig_2_09_17_2019';
save_filename=fullfile(dir,'plots','onwards_august_2019',sname);
print(save_filename,'-r300','-djpeg');
%}
% joint density of hydraulic conductivity and standard deviation
%{
bernardo_prior=sortrows(bernardo_prior,2);
[X,Y]=meshgrid(sig_range,kh_range);
Z=reshape(bernardo_prior(:,3)',length(kh_range),length(sig_range));
figure;
surf(X,Y,Z,'LineStyle','none','facealpha',0.8)
colorbar;
xlabel('standard deviation ( \sigma, cm)',...
    'fontname','arial','fontsize',12);
ylabel(' K (cm h^{-1})',...
    'fontname','arial','fontsize',12);
zlabel('prior density','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12)
clear box
%
sname='FHGA_prior_joint_prior_09_17_2019';
save_filename=fullfile(dir,'plots','onwards_august_2019',sname);
print(save_filename,'-r300','-djpeg');
%}

% density of shape parameter beta for fixed mu (mean) and phi (scale parameter)
%{
mu=bernardo_prior(1,1);
phi=bernardo_prior(1,2);
plot(bernardo_prior(:,3),bernardo_prior(:,4),'linewidth',2);
xlabel('shape parameter (\beta)',...
    'fontname','arial','fontsize',12);
ylabel('Bernardo prior density',...
    'fontname','arial','fontsize',12);
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',12,box)
title(['\mu=',num2str(mu),',',char(32),'\phi=',num2str(phi)]);
ylim([-5 40]);
hold on
plot([min(bernardo_prior(:,3)),max(bernardo_prior(:,3))],[0,0],...
    '--','color','black','linewidth',2);

sname='GG_prior_beta_mu=0_phi=1_11_06_2019';
save_filename=fullfile(dir,sname);
print(save_filename,'-r300','-djpeg');

%}

% joint-density of mu_ks and sigma_ks
%

Z=reshape(bernardo_prior(:,4),[length(sigma_ks_range),length(mu_ks_range)]);
[X,Y]=meshgrid(mu_ks_range,sigma_ks_range);
surf(X,Y,Z,'edgecolor','none','facealpha',0.7);

% find maximum density and corresponding mu_ks for each sigma_ks
%{
for sigma_ks_ind=1:size(Y,1)
    ind=find(Z(sigma_ks_ind,:)==max(Z(sigma_ks_ind,:)));
    X_max(sigma_ks_ind,1)=X(1,ind);
    Z_max(sigma_ks_ind,1)=Z(sigma_ks_ind,ind);
end
hold on; plot3(X_max,Y(:,1),Z_max','--','color','black','linewidth',2);
%}
xline=reshape(X(4:11,1:3),numel(X(4:11,1:3)),1);
yline=reshape(Y(4:11,1:3),numel(Y(4:11,1:3)),1);
zline=reshape(Z(4:11,1:3),numel(Z(4:11,1:3)),1);
k = convhull(xline,yline,zline);
k=convhull([xline,yline,zline],'Simplify',true);
hold on; trisurf(k,xline,yline,zline,'facealpha',0.2,'facecolor','none');

% hold on; plot3(k(:,1),k(:,2),k(:,3),'o','color','black','linewidth',2);

colormap summer
break
colorbar
xlabel('\mu_K (mm h^{-1})','fontname','arial','fontsize',14);
ylabel('\sigma_K (mm h^{-1})','fontname','arial','fontsize',14);
zlabel('prior','fontname','arial','fontsize',14);

set(gca,'fontname','arial','fontsize',14,'xlim',[0 0.3],...
    'ylim',[0.7 2.1]);

save_fname='FSI_bernardo_prior_cal4_02_03_2020_musigma_limited';
save_filename=fullfile(dir,'results','field_scale_infiltration',save_fname);
print(save_filename,'-r600','-djpeg');
%}

% marginal prior over mu_ks
%{
sigma_ks_values=[0.9,...
    1.3,...
    1.9];
colorspec={[1,0,0],[0,1,0],[0,0,1]};
line_pattern={'-','--',':'};
figure;
hold on
for i=1:length(sigma_ks_values)
    ind=find(bernardo_prior(:,2)==sigma_ks_values(i));
    mu_ks=bernardo_prior(ind,1);
    dens=bernardo_prior(ind,4);
    plot(mu_ks,dens,line_pattern{i},'linewidth',2,'color',colorspec{i})
    
end

box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',14,'xlim',[min(mu_ks_range) max(mu_ks_range)],box)
xlabel('\mu_K, (mm h^{-1})',...
    'fontname','arial','fontsize',14);
ylabel('prior density','fontname','arial','fontsize',14);
legend({['\sigma_K=',num2str(sigma_ks_values(1)),char(32),'mm h^{-1}'],...
    ['\sigma_K=',num2str(sigma_ks_values(2)),char(32),'mm h^{-1}'],...
    ['\sigma_K=',num2str(sigma_ks_values(3)),char(32),'mm h^{-1}']},...
    'fontname','arial','fontsize',14,'Location','northeast');
legend('boxoff');
clear box

%
sname='FSI_prior_mu_ks_02_03_2020';
save_filename=fullfile(dir,'results','field_scale_infiltration',sname);
print(save_filename,'-r600','-djpeg');
%}

% marginal prior over sigma_ks
%{
mu_ks_values=[0.9,...
    1.3,...
    1.9];
colorspec={[1,0,0],[0,1,0],[0,0,1]};
line_pattern={'-','--',':'};
figure;
hold on
for i=1:length(mu_ks_values)
    
    ind=find(bernardo_prior(:,1)==mu_ks_values(i));
    sigma_ks=bernardo_prior(ind,2);
    dens=bernardo_prior(ind,4);
    plot(sigma_ks,dens,line_pattern{i},'linewidth',2,'color',colorspec{i})
    
end

box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',14,'xlim',[min(mu_ks_range) max(mu_ks_range)],box)
xlabel('\sigma_K, (mm h^{-1})',...
    'fontname','arial','fontsize',14);
ylabel('prior density','fontname','arial','fontsize',14);
legend({['\mu_K=',num2str(mu_ks_values(1)),char(32),'mm h^{-1}'],...
    ['\mu_K=',num2str(mu_ks_values(2)),char(32),'mm h^{-1}'],...
    ['\mu_K=',num2str(mu_ks_values(3)),char(32),'mm h^{-1}']},...
    'fontname','arial','fontsize',14,'Location','northwest');
legend('boxoff');
clear box

%
sname='FSI_prior_sigma_ks_02_03_2020';
save_filename=fullfile(dir,'results','field_scale_infiltration',sname);
print(save_filename,'-r300','-djpeg');
%}
%% chain plots
%{
sig2='2';
y0='1.2';
prior='bernardo';
fname=strcat(['kh_sig2=',sig2,'_y0=',y0,...
    '_Posterior_mhsample_',prior,'_prior_07_06_2019']);
filename=fullfile(dir,'results',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%f');
fclose(fid);
kh=3600*data{1};            % in cm h^-1

plot(kh,'color','black')
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',14,box)
xlabel('step-number','fontname','arial','fontsize',14);
ylabel('K_h (cm h^{-1})','fontname','arial','fontsize',14);
save_name=strcat(['chain_plots_kh_sig2=',sig2,'_y0=',y0,...
    prior,'_prior_07_06_2019.jpg']);
save_filename=fullfile(dir,'plots',save_name);
print(save_filename,'-r300','-djpeg')
%}

%% Posterior histogram
% hydraulic conductivity (K) and shape parameter (beta)
%{
burn_in=2000;
sig2='2';
y0='3.17';
fname=strcat(['kh_sigma=',sig2,'_y0=',y0,...
    '_Posterior_mhsample_bernardo_prior_08_15_2019']);
filename=fullfile(dir,'results','onwards_august_2019',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%f');
fclose(fid);
kh_ber=3600*data{1}(burn_in+1:end);            % in cm h^-1

fname=strcat(['kh_sigma=',sig2,'_y0=',y0,...
    '_Posterior_mhsample_uniform_prior_08_15_2019']);
filename=fullfile(dir,'results','onwards_august_2019',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%f');
fclose(fid);
kh_unif=3600*data{1}(burn_in+1:end);            % in cm h^-1

% kernel density
[f_ber,x_ber]=ksdensity(kh_ber,'support',[0,50.0001]);
[f_unif,x_unif]=ksdensity(kh_unif,'support',[0,50.0001]);

% hist(kh_ber,20);
% hold on
% hist(kh_unif,20)
% h=findobj(gca,'Type','patch');
% set(h(1),'facecolor',[0.7,1,0.7],'edgecolor','w','facealpha',0.75)
% set(h(2),'facecolor',[1,0.4,0.4],'edgecolor','w','facealpha',0.75)

plot(x_ber,f_ber,x_unif,f_unif,'linewidth',2);
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',24,'xlim',[0,1],box)
legend({'Bernardo prior','Uniform prior'},'Location','best',...
    'fontname','arial','fontsize',24);
legend('boxoff');
xlabel('K_h (cm h^{-1})','fontname','arial','fontsize',24);

fname=strcat(['GA_kh_sigma=',sig2,'_y0=',y0,...
    '_Posterior_mhsample_bernardo_prior_08_15_2019.jpeg']);
save_filename=fullfile(dir,'AGU_poster_plots',fname);
print(save_filename,'-r300','-djpeg');
%}

% mu_ks and sigma_ks
%
burn_in=2000;
sig2=0.01;
fname=strcat('FSI_sig2=',num2str(sig2),'_Posterior_mhsample_bernardo_prior_cal4_02_03_2020');
filename=fullfile(dir,'results','field_scale_infiltration',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%f%f','delimiter',',');
fclose(fid);
mu_ks_bernardo=data{1}(burn_in+1:end); sigma_ks_bernardo=data{2}(burn_in+1:end);

fname=strcat('FSI_sig2=',num2str(sig2),'_Posterior_mhsample_uniform_prior_cal4_02_03_2020');
filename=fullfile(dir,'results','field_scale_infiltration',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%f%f','delimiter',',');
fclose(fid);
mu_ks_uniform=data{1}(burn_in+1:end); sigma_ks_uniform=data{2}(burn_in+1:end);

% kernel density of mu_ks
[f_ber,x_ber]=ksdensity(mu_ks_bernardo,'support',[0,2.11]);
[f_unif,x_unif]=ksdensity(mu_ks_uniform,'support',[0,2.11]);
plot(x_ber,f_ber,x_unif,f_unif,'-.','linewidth',2);
xlabel('\mu_K (mm h^{-1})','fontname','arial','fontsize',14);
ylabel('posterior density','fontname','arial','fontsize',14);
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',14,'xlim',[0.55 0.60],box);

legend('using bernardo prior','using uniform prior',...
    'fontname','arial','fontsize',14,'location','northwest');
legend('boxoff')
clear box;

sname=strcat('FSI_posterior_mu_ks_sig2=',num2str(sig2),'_cal4_02_03_2020.jpg');
save_filename=fullfile(dir,'results','field_scale_infiltration',sname);
print(save_filename,'-r300','-djpeg');
pause;
% kernel density of sigma_ks
[f_ber,x_ber]=ksdensity(sigma_ks_bernardo,'support',[0,2.11]);
[f_unif,x_unif]=ksdensity(sigma_ks_uniform,'support',[0,2.11]);
plot(x_ber,f_ber,x_unif,f_unif,'-.','linewidth',2);
xlabel('\sigma_K (mm h^{-1})','fontname','arial','fontsize',14);
ylabel('posterior density','fontname','arial','fontsize',14);
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',14,'xlim',[0 2.1],box);

legend('using bernardo prior','using uniform prior',...
    'fontname','arial','fontsize',14,'Location','northwest');
legend('boxoff')
clear box;

sname=strcat('FSI_posterior_sigma_ks_sig2=',num2str(sig2),'_cal4_02_03_2020.jpg');
save_filename=fullfile(dir,'results','field_scale_infiltration',sname);
print(save_filename,'-r300','-djpeg');
%}