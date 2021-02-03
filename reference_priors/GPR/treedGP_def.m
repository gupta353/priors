% this script defines a treed Gaussian process to compute streamflow for a
% given parameter set


clear all
close all
clc

direc = 'D:/Research/Thesis_work/Non_informative_priors/matlab_codes/reference_priors';

% optimal parameter values
fname = 'theta_opt_region_wise.mat';
filename = fullfile(direc,'results/pdm_giuh',fname);
load(filename);

% train-test data
fname = 'GP_train_test_data_updated.mat';
filename = fullfile(direc,'results/pdm_giuh',fname);
load(filename);

% convert ks into log_ks space
Xtrain(:,3) = log(Xtrain(:,3))/log(10);
Xtest(:,3) = log(Xtest(:,3))/log(10);

sig2 = 0.0001;          % a small observation variance is added to keep the covariance matrix positive definite

%% region 1
ind = find(Xtrain(:,3)>-5 & Xtrain(:,2)>0.04 & Xtrain(:,1)>50);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{1} = Xtrain_r;
ytrain_final{1} = ytrain_r;
m_final(1,:) = m;
stand_dev_final(1,:) = stand_dev;

theta = theta_opt_region_wise(1,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{1} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{1},'lower');         % cholesky decomposition of covariance matrix
alpha{1} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{1} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 2
ind = find(Xtrain(:,3)>-5 & Xtrain(:,2)>0.04 & Xtrain(:,1)<=50);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{2} = Xtrain_r;
ytrain_final{2} = ytrain_r;
m_final(2,:) = m;
stand_dev_final(2,:) = stand_dev;

theta = theta_opt_region_wise(2,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{2} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{2},'lower');         % cholesky decomposition of covariance matrix
alpha{2} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{2} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 3
ind = find(Xtrain(:,3)<=-5 & Xtrain(:,2)>0.04 & Xtrain(:,2)>1 & Xtrain(:,4)>1);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{3} = Xtrain_r;
ytrain_final{3} = ytrain_r;
m_final(3,:) = m;
stand_dev_final(3,:) = stand_dev;

theta = theta_opt_region_wise(3,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{3} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{3},'lower');         % cholesky decomposition of covariance matrix
alpha{3} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{3} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 4
ind = find(Xtrain(:,3)<=-5 & Xtrain(:,2)>0.04 & Xtrain(:,2)>1 & Xtrain(:,4)<=1 & Xtrain(:,4)>0.2);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{4} = Xtrain_r;
ytrain_final{4} = ytrain_r;
m_final(4,:) = m;
stand_dev_final(4,:) = stand_dev;

theta = theta_opt_region_wise(4,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{4} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{4},'lower');         % cholesky decomposition of covariance matrix
alpha{4} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{4} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 5
ind = find(Xtrain(:,3)<=-5 & Xtrain(:,2)>0.04 & Xtrain(:,2)>1 & Xtrain(:,4)<=0.2 & Xtrain(:,4)>0.04);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{5} = Xtrain_r;
ytrain_final{5} = ytrain_r;
m_final(5,:) = m;
stand_dev_final(5,:) = stand_dev;

theta = theta_opt_region_wise(5,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{5} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{5},'lower');         % cholesky decomposition of covariance matrix
alpha{5} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{5} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 6
ind = find(Xtrain(:,3)<=-5 & Xtrain(:,2)>0.04 & Xtrain(:,2)>1 & Xtrain(:,4)<=0.04);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{6} = Xtrain_r;
ytrain_final{6} = ytrain_r;
m_final(6,:) = m;
stand_dev_final(6,:) = stand_dev;

theta = theta_opt_region_wise(6,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{6} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{6},'lower');         % cholesky decomposition of covariance matrix
alpha{6} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{6} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 7
ind = find(Xtrain(:,3)<=-5 & Xtrain(:,2)>0.04 & Xtrain(:,2)<=1 & Xtrain(:,2)>0.2 & Xtrain(:,4)>2);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{7} = Xtrain_r;
ytrain_final{7} = ytrain_r;
m_final(7,:) = m;
stand_dev_final(7,:) = stand_dev;

theta = theta_opt_region_wise(7,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{7} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{7},'lower');         % cholesky decomposition of covariance matrix
alpha{7} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{7} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 8
ind = find(Xtrain(:,3)<=-5 & Xtrain(:,2)>0.04 & Xtrain(:,2)<=1 & Xtrain(:,2)>0.2 & Xtrain(:,4)<=2 & Xtrain(:,4)>0.2);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{8} = Xtrain_r;
ytrain_final{8} = ytrain_r;
m_final(8,:) = m;
stand_dev_final(8,:) = stand_dev;

theta = theta_opt_region_wise(8,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{8} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{8},'lower');         % cholesky decomposition of covariance matrix
alpha{8} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{8} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 9
ind = find(Xtrain(:,3)<=-5 & Xtrain(:,2)>0.04 & Xtrain(:,2)<=1 & Xtrain(:,2)>0.2 & Xtrain(:,4)<=0.2);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{9} = Xtrain_r;
ytrain_final{9} = ytrain_r;
m_final(9,:) = m;
stand_dev_final(9,:) = stand_dev;

theta = theta_opt_region_wise(9,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{9} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{9},'lower');         % cholesky decomposition of covariance matrix
alpha{9} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{9} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 10
ind = find(Xtrain(:,3)<=-5 & Xtrain(:,2)>0.04 & Xtrain(:,2)<=0.2 & Xtrain(:,2)>0.12 & Xtrain(:,4)>2);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{10} = Xtrain_r;
ytrain_final{10} = ytrain_r;
m_final(10,:) = m;
stand_dev_final(10,:) = stand_dev;

theta = theta_opt_region_wise(10,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{10} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{10},'lower');         % cholesky decomposition of covariance matrix
alpha{10} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{10} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 11
ind = find(Xtrain(:,3)<=-5 & Xtrain(:,2)>0.04 & Xtrain(:,2)<=0.2 & Xtrain(:,2)>0.12 & Xtrain(:,4)<=2 & Xtrain(:,4)>0.2);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{11} = Xtrain_r;
ytrain_final{11} = ytrain_r;
m_final(11,:) = m;
stand_dev_final(11,:) = stand_dev;

theta = theta_opt_region_wise(11,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{11} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{11},'lower');         % cholesky decomposition of covariance matrix
alpha{11} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{11} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 12
ind = find(Xtrain(:,3)<=-5 & Xtrain(:,2)>0.04 & Xtrain(:,2)<=0.2 & Xtrain(:,2)>0.12 & Xtrain(:,4)<=0.2);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{12} = Xtrain_r;
ytrain_final{12} = ytrain_r;
m_final(12,:) = m;
stand_dev_final(12,:) = stand_dev;

theta = theta_opt_region_wise(12,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{12} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{12},'lower');         % cholesky decomposition of covariance matrix
alpha{12} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{12} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 13
ind = find(Xtrain(:,3)<=-5 & Xtrain(:,3)>-8.5 & Xtrain(:,2)>0.04 & Xtrain(:,2)<=0.12 & Xtrain(:,4)>2);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{13} = Xtrain_r;
ytrain_final{13} = ytrain_r;
m_final(13,:) = m;
stand_dev_final(13,:) = stand_dev;

theta = theta_opt_region_wise(13,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{13} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{13},'lower');         % cholesky decomposition of covariance matrix
alpha{13} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{13} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 14
ind = find(Xtrain(:,3)<=-8.5 & Xtrain(:,2)>0.04 & Xtrain(:,2)<=0.12 & Xtrain(:,4)>2);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{14} = Xtrain_r;
ytrain_final{14} = ytrain_r;
m_final(14,:) = m;
stand_dev_final(14,:) = stand_dev;

theta = theta_opt_region_wise(14,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{14} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{14},'lower');         % cholesky decomposition of covariance matrix
alpha{14} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{14} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 15
ind = find(Xtrain(:,3)<=-5 & Xtrain(:,2)>0.04 & Xtrain(:,2)<=0.12 & Xtrain(:,4)<=2 & Xtrain(:,4)>0.2);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{15} = Xtrain_r;
ytrain_final{15} = ytrain_r;
m_final(15,:) = m;
stand_dev_final(15,:) = stand_dev;

theta = theta_opt_region_wise(15,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{15} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{15},'lower');         % cholesky decomposition of covariance matrix
alpha{15} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{15} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region 16
ind = find(Xtrain(:,3)<=-5 & Xtrain(:,2)>0.04 & Xtrain(:,2)<=0.12 & Xtrain(:,4)<=0.2);
Xtrain_r = Xtrain(ind,:);
ytrain_r = ytrain(ind,:);
stand_dev = std(Xtrain_r);
m = mean(Xtrain_r);
Xtrain_r = bsxfun(@minus,Xtrain_r,m);
Xtrain_r = bsxfun(@rdivide,Xtrain_r,stand_dev);

Xtrain_final{16} = Xtrain_r;
ytrain_final{16} = ytrain_r;
m_final(16,:) = m;
stand_dev_final(16,:) = stand_dev;

theta = theta_opt_region_wise(16,:);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);
Ktrain{16} = KerComp_iden(Xtrain_r,Xtrain_r,M,sigf2)+sig2*eye(size(Xtrain_r,1));

L = chol(Ktrain{15},'lower');         % cholesky decomposition of covariance matrix
alpha{16} = L'\(L\ytrain_r);

S = Xtrain_r*M*Xtrain_r';
S1{16} = diag(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define regions
region_def = @(x)[(x(:,3)>-5 & x(:,2)>0.04 & x(:,1)>50),...                             % region 1
    (x(:,3)>-5 & x(:,2)>0.04 & x(:,1)<=50),...                                          % region 2
    (x(:,3)<=-5 & x(:,2)>0.04 & x(:,2)>1 & x(:,4)>1),...                                % region 3
    (x(:,3)<=-5 & x(:,2)>0.04 & x(:,2)>1 & x(:,4)<=1 & x(:,4)>0.2),...                  % region 4
    (x(:,3)<=-5 & x(:,2)>0.04 & x(:,2)>1 & x(:,4)<=0.2 & x(:,4)>0.04),...               % region 5
    (x(:,3)<=-5 & x(:,2)>0.04 & x(:,2)>1 & x(:,4)<=0.04),...                            % region 6
    (x(:,3)<=-5 & x(:,2)>0.04 & x(:,2)<=1 & x(:,2)>0.2 & x(:,4)>2),...                  % region 7
    (x(:,3)<=-5 & x(:,2)>0.04 & x(:,2)<=1 & x(:,2)>0.2 & x(:,4)<=2 & x(:,4)>0.2),...    % region 8
    (x(:,3)<=-5 & x(:,2)>0.04 & x(:,2)<=1 & x(:,2)>0.2 & x(:,4)<=0.2),...               % region 9
    (x(:,3)<=-5 & x(:,2)>0.04 & x(:,2)<=0.2 & x(:,2)>0.12 & x(:,4)>2),...               % region 10
    (x(:,3)<=-5 & x(:,2)>0.04 & x(:,2)<=0.2 & x(:,2)>0.12 & x(:,4)<=2 & x(:,4)>0.2),... % region 11
    (x(:,3)<=-5 & x(:,2)>0.04 & x(:,2)<=0.2 & x(:,2)>0.12 & x(:,4)<=0.2),...            % region 12
    (x(:,3)<=-5 & x(:,3)>-8.5 & x(:,2)>0.04 & x(:,2)<=0.12 & x(:,4)>2),...              % region 13
    (x(:,3)<=-8.5 & x(:,2)>0.04 & x(:,2)<=0.12 & x(:,4)>2),...                          % region 14
    (x(:,3)<=-5 & x(:,2)>0.04 & x(:,2)<=0.12 & x(:,4)<=2 & x(:,4)>0.2),...              % region 15
    (x(:,3)<=-5 & x(:,2)>0.04 & x(:,2)<=0.12 & x(:,4)<=0.2)];                           % region 16
%% define a structure containing GP information
GP.Xtrain = Xtrain_final;
GP.ytrain = ytrain_final;
GP.m = m_final;
GP.stand_dev = stand_dev_final;
GP.Ktrain = Ktrain;
GP.region_def = region_def;
GP.params = theta_opt_region_wise;
GP.alpha = alpha;
GP.S1 = S1;

sname = 'treedGP.mat';
filename = fullfile(direc,'results/pdm_giuh',sname);
save(filename,'GP')