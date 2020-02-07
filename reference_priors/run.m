% This routine is the main script to compute posterior distribution of
% parameters of a model, using the Bernardo's prior as an objective prior

% clear all
% close all
% clc

global GLOBAL_DATA

save_dir=['D:/Research/Thesis_work/Non_informative_priors'...
    '/matlab_codes/reference_priors'];


d=3; % dimension of the parameter space
% sig2=0.1;
% read Bernardo's prior data
fname='FSI_bernardo_prior_cal4_02_03_2020';
filename=fullfile(save_dir,'results','field_scale_infiltration',fname);
formatspec=repmat('%f',1,d+1);
fid=fopen(filename,'r');
data=textscan(fid,formatspec,'delimiter',',');
fclose(fid);
Bernardo_pdf=[data{1},data{2},data{3},data{4}];

% observed cumulative infiltration
%{
fname='infiltrometer_A_outer_ring.txt';
filename=fullfile(save_dir,'data','Infiltration_data',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%f%f','delimiter','\t','headerlines',1);
fclose(fid);
y0=data{2};
t=data{1};
 % remove the time-step zero from the data
y0(1)=[];      
t(1)=[];

GLOBAL_DATA.t=t;
GLOBAL_DATA.y0=y0;
kh_range=unique(Bernardo_pdf(:,1));
sigma_range=unique(Bernardo_pdf(:,2));

GLOBAL_DATA.Bernardo_pdf=Bernardo_pdf;
GLOBAL_DATA.kh_range=kh_range;
GLOBAL_DATA.sigma_range=sigma_range;
%}

% observed field scale infiltration rate
%
% extract observed data from etxtfiles
direc=['D:/Research/Thesis_work/Non_informative_priors/'...
    'matlab_codes/reference_priors/data/experiment_data'];

theta_sat=0.36;                         % saturated moisture content
fname='cal_4.txt';                      % experiment name
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
r(end)=[]; obs(end)=[]; t(end)=[]; infil(end)=[];
r=r'; obs=obs'; t=t';

GLOBAL_DATA.t=t;
GLOBAL_DATA.y0=obs;
GLOBAL_DATA.r=r;
GLOBAL_DATA.psi=psi;
GLOBAL_DATA.delta_theta=delta_theta;

mu_ks_range=unique(Bernardo_pdf(:,1));
sigma_ks_range=unique(Bernardo_pdf(:,2));

GLOBAL_DATA.Bernardo_pdf=Bernardo_pdf;
GLOBAL_DATA.mu_ks_range=mu_ks_range;
GLOBAL_DATA.sigma_ks_range=sigma_ks_range;
%}
%% computation of marginal pdf of hydraulic conductivity K_h
%{
if length(sigma_range)>1            % if sigma is assumed to be unknowna and prior over sigma is also computed
    pdf_kh=sortrows(Bernardo_pdf,1);
    for i=1:length(kh_range)
        kh=kh_range(i);
        ind=find(pdf_kh(:,1)==kh);
        sigma=pdf_kh(ind,2);
        base_area=sigma(2)-sigma(1);
        pdf_temp=pdf_kh(ind,3);
        pdf_temp(1)=pdf_temp(1)/2;
        pdf_temp(end)=pdf_temp(end)/2;
        marginal_pdf_values_kh(i)=sum(pdf_temp)*base_area;
    end
   GLOBAL_DATA.marginal_pdf_values_kh=marginal_pdf_values_kh;
else                                % is sigma is assumed to known
    marginal_pdf_values_kh=Bernardo_pdf(:,3);
    GLOBAL_DATA.marginal_pdf_values_kh=marginal_pdf_values_kh;
end
%}
%% computation of marginal pdf of sigma
%{
if length(sigma_range)>1
    pdf_sigma=sortrows(Bernardo_pdf,2);
    % compute marginal pdf
    for i=1:length(sigma_range)
        sigma=sigma_range(i);
        ind=find(pdf_sigma(:,2)==sigma);
        kh=pdf_sigma(ind,1);
        base_area=kh(2)-kh(1);
        pdf_temp=pdf_sigma(ind,3);
        pdf_temp(1)=pdf_temp(1)/2;
        pdf_temp(end)=pdf_temp(end)/2;
        marginal_pdf_values_sigma(i)=base_area*sum(pdf_temp);
    end
GLOBAL_DATA.marginal_pdf_values_sigma=marginal_pdf_values_sigma;
end
%}
%% computation of posterior using DREAM
%{

global DREAM_dir EXAMPLE_dir CONV_dir
DREAM_dir=['D:/Research/Thesis_work/Non_informative_priors'...
    '/matlab_codes/reference_priors/Dream'];
CONV_dir=['D:/Research/Thesis_work/Non_informative_priors'...
    '/matlab_codes/reference_priors/Dream/diagnostics'];
EXAMPLE_dir = ['D:/Research/Thesis_work/Non_informative_priors'...
    '/matlab_codes/reference_priors'];
addpath(DREAM_dir)
addpath(CONV_dir)
addpath(EXAMPLE_dir)

Func_name='Gaussloglikeli';                 % target pdf
DREAMPar.d=2;                               % dimensionality of the problem
DREAMPar.N=8;                               % Number of chains
DREAMPar.T=5000;                             % Number of generations
DREAMPar.lik=2;                             % Type of likelihood (2 indicates log-likelihood)

Par_info.initial='prior';                 % type of prior 
Par_info.min=[0.0001,0.0001];                   
Par_info.max=[0.0083333,10];          % upper bound
Par_info.boundhandling='fold';              % type of boundhandling
Par_info.prior=@(x)berpdf(x);

[chain,output,fx]=DREAM(Func_name,...       % main function to Run Dream
    DREAMPar,Par_info);
%}

%% computation of posterior using MATLAB mhsample
%
sig2=0.01;
known_sigma=sqrt(sig2);
start=[1,1];
nsamples=10000;
min_value=[0.01,0.01];                   
max_value=[2.10,2.10];
pdf=@(x)exp(Gaussloglikeli_fieldscaleGA(x,known_sigma))*bernardo_pdf([x,2]);%*prod(unifpdf(x,min_value(1),max_value(1)));;       % target distribution
proppdf=@(x,y)prod(unifpdf(x,min_value,max_value));                                % proposal distribution
proprnd=@(x)unifrnd(min_value,max_value);                                          % random number geenrator from proposal distribution
smpl = mhsample(start,nsamples,'pdf',pdf,'proppdf',proppdf, 'proprnd',proprnd);

% save the smpl file
%
save_fname=strcat('FSI_sig2=',num2str(sig2),'_Posterior_mhsample_bernardo_prior_cal4_02_03_2020');
save_filename=fullfile(save_dir,'results','field_scale_infiltration',save_fname);
dlmwrite(save_filename,smpl);
%}