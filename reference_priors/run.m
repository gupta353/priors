% This routine is the main script to compute posterior distribution of
% parameters of a model, using the Bernardo's prior as an objective prior

clear all
close all
clc

global GLOBAL_DATA

save_dir=['D:/Research/Thesis_work/Non_informative_priors'...
    '/matlab_codes/reference_priors'];


d=2; % dimension of the data

% read Bernardo's prior data
fname='prior_density_data_initial_prior_2_sig2=1_07_06_2019';
filename=fullfile(save_dir,'results',fname);
formatspec=repmat('%f',1,d+1);
fid=fopen(filename,'r');
data=textscan(fid,formatspec,'delimiter',',');
fclose(fid);
Bernardo_pdf=[data{1},data{2},data{3}];

% observed cumulative infiltration
y0=1.2;            % in cm
GLOBAL_DATA.y0=y0;
kh_range=unique(Bernardo_pdf(:,1));
sigma_range=unique(Bernardo_pdf(:,2));

GLOBAL_DATA.Bernardo_pdf=Bernardo_pdf;
GLOBAL_DATA.kh_range=kh_range;
GLOBAL_DATA.sigma_range=sigma_range;

%% computation of marginal pdf of hydraulic conductivity K_h
if length(sigma_range)>1
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
else
    marginal_pdf_values_kh=Bernardo_pdf(:,3);
    GLOBAL_DATA.marginal_pdf_values_kh=marginal_pdf_values_kh;
end
%% computation of marginal pdf of sigma
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
known_sigma=1;
start=[0.04/3600];
nsamples=10000;
min_value=0.01/3600;                   
max_value=1/3600; 
pdf=@(x)exp(Gaussloglikeli([x,known_sigma]))*marginal_bernardo_pdf_kh(x);
proppdf=@(x,y)unifpdf(x,min_value,max_value);
proprnd=@(x)unifrnd(min_value,max_value);
smpl = mhsample(start,nsamples,'pdf',pdf,'proppdf',proppdf, 'proprnd',proprnd);

% save the smpl file
%
save_fname='kh_sigma=1_y0=1.2_Posterior_mhsample_bernardo_prior_07_06_2019';
save_filename=fullfile(save_dir,save_fname);
dlmwrite(save_filename,smpl);
%}