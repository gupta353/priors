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
fname='prior_density_data_initial_prior_2';
filename=fullfile(save_dir,'results',fname);
formatspec=repmat('%f',1,d+1);
fid=fopen(filename,'r');
data=textscan(fid,formatspec,'delimiter',',');
fclose(fid);
Bernardo_pdf=[data{1},data{2},data{3}];

% observed cumulative infiltration
y0=1.27;            % in cm

% GLOBAL_DATA
GLOBAL_DATA.Bernardo_pdf=Bernardo_pdf;
GLOBAL_DATA.y0=y0;

% computation of posterior using DREAM
%

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
DREAMPar.T=500;                             % Number of generations
DREAMPar.lik=2;                             % Type of likelihood (2 indicates log-likelihood)

Par_info.initial='uniform';                 % type of prior 
Par_info.min=[0.0001,0.0001];                   
Par_info.max=[0.0083333,100];          % upper bound
Par_info.boundhandling='fold';              % type of boundhandling

[chain,output,fx]=DREAM(Func_name,...       % main function to Run Dream
    DREAMPar,Par_info);
%}