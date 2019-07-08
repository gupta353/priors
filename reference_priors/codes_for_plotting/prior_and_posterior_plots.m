% This routine plots prior and posterior distribution for various
% cases-studies conducted in thi study
clear all
close all
clc

dir=['D:/Research/Thesis_work/Non_informative_priors'...
    '/matlab_codes/reference_priors'];


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
burn_in=2000;
sig2='2';
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

hist(kh_ber,20);
hold on
hist(kh_unif,20)
h=findobj(gca,'Type','patch');
set(h(1),'facecolor',[0.7,1,0.7],'edgecolor','w','facealpha',0.75)
set(h(2),'facecolor',[1,0.4,0.4],'edgecolor','w','facealpha',0.75)
box('on');
box.linewidth=2;
set(gca,'fontname','arial','fontsize',14,'xlim',[0,1],box)
legend({'Bernardo prior','Uniform prior'},'Location','northwest',...
    'fontname','arial','fontsize',14);
legend('boxoff');
xlabel('K_h (cm h^{-1})','fontname','arial','fontsize',14);
ylabel('count','fontname','arial','fontsize',14);

fname=strcat(['kh_sig2=',sig2,'_y0=',y0,...
    '_hist','_07_06_2019.jpeg']);
save_filename=fullfile(dir,'plots',fname);
print(save_filename,'-r300','-djpeg');