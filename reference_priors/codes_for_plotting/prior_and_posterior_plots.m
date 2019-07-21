% This routine plots prior and posterior distribution for various
% cases-studies conducted in thi study
clear all
close all
clc

dir=['D:/Research/Thesis_work/Non_informative_priors'...
    '/matlab_codes/reference_priors'];

%% prior plots
%
fname='prior_density_data_initial_prior_1_07_07_2019';
filename=fullfile(dir,'results',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%f%f%f','delimiter',',');
fclose(fid);
bernardo_prior=[3600*data{1},data{2},data{3}];              % factor of 3600 used to convert the K values from cm s^-1 to cm h^-1
kh_range=unique(bernardo_prior(:,1));
sig_range=unique(bernardo_prior(:,2));

% density of hydraulic conductivity for each value of standard
% deviation
%{
sig_values=[0.1,0.5,1];
colorspec={'r','g','b'};
figure;
hold on
for i=1:length(sig_values)
    ind=find(bernardo_prior(:,2)==sig_values(i));
    kh=bernardo_prior(ind,1);
    dens=bernardo_prior(ind,3);
    plot(kh,dens,'linewidth',2,'color',colorspec{i})
    
end
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',14,'xlim',[0 1],box)
xlabel('hydraulic conductivity (K_h, cm h^{-1})',...
    'fontname','arial','fontsize',14);
ylabel('prior density','fontname','arial','fontsize',14);
legend({['\sigma=',num2str(sig_values(1)),char(32),'cm'],...
    ['\sigma=',num2str(sig_values(2)),char(32),'cm'],...
    ['\sigma=',num2str(sig_values(3)),char(32),'cm']},...
    'fontname','arial','fontsize',14,'Location','northwest');
legend('boxoff');
clear box

%
sname='GA_prior_K_intital_prior_1_07_07_2019';
save_filename=fullfile(dir,'plots',sname);
print(save_filename,'-r300','-djpeg');
%}
% density of standard deviation for each value of 
% hydraulic conductivity
%{
kh_values=[0.1,0.5,1];          % in cm h^-1
colorspec={'r','g','b'};
figure;
hold on
for i=1:length(kh_values)
    ind=find(floor(bernardo_prior(:,1)*10^4)/10^4==kh_values(i));           % factor of 10^4 is used becuase of precision errors in the values of kh in bernardo_prior
    sig=bernardo_prior(ind,2);
    dens=bernardo_prior(ind,3);
    plot(sig,dens,'linewidth',2,'color',colorspec{i})
    
end
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',14,'xlim',[0 1],box)
xlabel('standard deviation (\sigma, cm)',...
    'fontname','arial','fontsize',14);
ylabel('prior density','fontname','arial','fontsize',14);
legend({['K_h=',num2str(kh_values(1)),char(32),'cm h^{-1}'],...
    ['K_h=',num2str(kh_values(2)),char(32),'cm h^{-1}'],...
    ['K_h=',num2str(kh_values(3)),char(32),'cm h^{-1}']},...
    'fontname','arial','fontsize',14);
legend('boxoff');
clear box

sname='GA_prior_sigma_intital_prior_1_07_07_2019';
save_filename=fullfile(dir,'plots',sname);
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
sname='GA_prior_joint_intital_prior_1_07_07_2019';
save_filename=fullfile(dir,'plots',sname);
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
%
burn_in=2000;
sig2='1';
y0='3.17';
fname=strcat(['kh_sig2=',sig2,'_y0=',y0,...
    '_Posterior_mhsample_bernardo_prior_07_06_2019']);
filename=fullfile(dir,'results',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%f');
fclose(fid);
kh_ber=3600*data{1}(burn_in+1:end);            % in cm h^-1

fname=strcat(['kh_sig2=',sig2,'_y0=',y0,...
    '_Posterior_mhsample_uniform_prior_07_06_2019']);
filename=fullfile(dir,'results',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%f');
fclose(fid);
kh_unif=3600*data{1}(burn_in+1:end);            % in cm h^-1

% kernel density
[f_ber,x_ber]=ksdensity(kh_ber,'support',[0,1.0001]);
[f_unif,x_unif]=ksdensity(kh_unif,'support',[0,1.0001]);

% hist(kh_ber,20);
% hold on
% hist(kh_unif,20)
% h=findobj(gca,'Type','patch');
% set(h(1),'facecolor',[0.7,1,0.7],'edgecolor','w','facealpha',0.75)
% set(h(2),'facecolor',[1,0.4,0.4],'edgecolor','w','facealpha',0.75)

plot(x_ber,f_ber,x_unif,f_unif,'linewidth',2);
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',14,'xlim',[0,1],box)
legend({'Bernardo prior','Uniform prior'},'Location','northwest',...
    'fontname','arial','fontsize',14);
legend('boxoff');
xlabel('K_h (cm h^{-1})','fontname','arial','fontsize',14);
ylabel('count','fontname','arial','fontsize',14);

fname=strcat(['kh_sig2=',sig2,'_y0=',y0,...
    '_ksdensity','_07_06_2019.jpeg']);
save_filename=fullfile(dir,'plots',fname);
print(save_filename,'-r300','-djpeg');
%}