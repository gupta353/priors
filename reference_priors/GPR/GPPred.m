% this script computes streamflow using a calibrated treed-Gaussian process


clear all
close all
clc

direc = 'D:/Research/Thesis_work/Non_informative_priors/matlab_codes/reference_priors';

% optimal parameter values
fname = 'treedGP.mat';
filename = fullfile(direc,'results/pdm_giuh',fname);
load(filename);

% train-test data
fname = 'pdm_giuh_runs_11000_test.mat';
filename = fullfile(direc,'results/pdm_giuh',fname);
load(filename);

% convert ks into log_ks space
XTEST(:,3) = log(XTEST(:,3))/log(10);
ind = find(XTEST(:,2)<=0.04);
XTEST(ind,:) = [];
YTEST(:,ind) = [];


sig2 = 0.0001;          % a small observation variance is added to keep the covariance matrix positive definite

region_def = GP.region_def;
profile on;
for samp = 1:size(XTEST,1)
    xtest = XTEST(samp,:);
    %% determine he region in which xtest lies
   
    reg = region_def(xtest);
    
    ind = find(reg==1);
    
    theta = GP.params(ind,:);
    sigf2 = theta(6);
    l = 1./theta(1:5);
    M = diag(l);
    
    Xtrain = GP.Xtrain{ind};
    ytrain = GP.ytrain{ind};
    Ktrain = GP.Ktrain{ind};
    
    % normalize the test sample
    xtest = bsxfun(@minus,xtest,GP.m(ind,:));
    xtest = bsxfun(@rdivide,xtest,GP.stand_dev(ind,:));
    
    % computation of covariance matrix (of response) between training and test
    % points
    Ktrt=KerComp_pred(Xtrain,xtest,GP.S1{ind},M,sigf2);
    
    alpha = GP.alpha{ind};
    % computation of covariance matrix (of response) between test data
%     Ktest=KerComp_iden(xtest,xtest,M,sigf2);
    
    ftest(:,samp) = alpha'*Ktrt;
end
profile off
