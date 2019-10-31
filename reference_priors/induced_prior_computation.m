% This routine computes induced prior over cumulative infiltration given a
% prior over hydraulic conductivity

clear all
close all
clc

save_dir=['D:/Research/Thesis_work/Non_informative_priors/'...
    'matlab_codes/reference_priors/plots'];
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
set(gca,'fontname','arial','fontsize',12,box)
xlabel('F(t=1 hour) (cm)',...
    'fontname','arial','fontsize',12);
ylabel('count','fontname','arial','fontsize',12);

sname='induced_prior_cum_infiltration_10_30_2019';
save_filename=fullfile(save_dir,sname);
print(save_filename,'-r300','-djpeg');
clear box