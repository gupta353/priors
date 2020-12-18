% this routine creates dataset for Gaussian emulation of PDM-GIUH after 1st
% calibration; after ist calibration, samples in a particular region of
% parameter space were found to be poorly predicted; this routine drawn
% samples from that region such that the drawn samples are at maximum
% distance from already drawn samples

clear all
close all
clc


direc = 'D:/Research/Thesis_work/Non_informative_priors/matlab_codes/reference_priors';
sname = 'GP_optim_data.mat';
filename = fullfile(direc,'results/pdm_giuh',sname);
load(filename);

Xtrain(:,3) = log(Xtrain(:,3))/log(10);

stand_dev = std(Xtrain);
m = mean(Xtrain);
Xtrain = bsxfun(@minus,Xtrain,m);
Xtrain = bsxfun(@rdivide,Xtrain,stand_dev);

param_samps = [];
%% draw samples from 1st region
n =10000;
par_range = [1, 2000;...                      % storage capacity (in mm)
    0.01, 10;...                              % storage distribution parameter
    -8, -5.23;...                               % logarithm of baseflow reservior constant (in s^-1)
    0.01, 10;...                              % in-stream-velocity (in m s^-1)
    0.01, 10];                                % hill slope velocity (in m s^-1)

X1 = LHsampAG(n,size(par_range,1),par_range);
X = bsxfun(@minus,X1,m);
X = bsxfun(@rdivide,X,stand_dev);

% select 1000 samples which are at maximum distance from the samples in the
% training data
for pind = 1:n
    
    Xt =  X(pind,:);
    d = bsxfun(@minus,Xtrain,Xt);
    d = sum(d.^2,2);
    d_min(pind) = min(d);
    
end

X1 = [d_min',X1];
X1 = sortrows(X1);
X1(:,1) = [];
X1 = X1(n-1000+1:n,:);
param_samps = [param_samps;X1];
%% draw samples from 2nd region
n =10000;
par_range = [1, 2000;...                      % storage capacity (in mm)
    0.01, 10;...                              % storage distribution parameter
    -3.067, -2;...                               % logarithm of baseflow reservior constant (in s^-1)
    0.01, 10;...                              % in-stream-velocity (in m s^-1)
    0.01, 10];                                % hill slope velocity (in m s^-1)

X1 = LHsampAG(n,size(par_range,1),par_range);
X = bsxfun(@minus,X1,m);
X = bsxfun(@rdivide,X,stand_dev);

% select 1000 samples which are at maximum distance from the samples in the
% training data
for pind = 1:n
    
    Xt =  X(pind,:);
    d = bsxfun(@minus,Xtrain,Xt);
    d = sum(d.^2,2);
    d_min(pind) = min(d);
    
end

X1 = [d_min',X1];
X1 = sortrows(X1);
X1(:,1) = [];
X1 = X1(n-1000+1:n,:);
param_samps = [param_samps;X1];
%% draw samples from 3rd region
%{
n =10000;
par_range = [1, 2000;...                      % storage capacity (in mm)
    0.01, 10;...                              % storage distribution parameter
    -15, -2;...                               % logarithm of baseflow reservior constant (in s^-1)
    0.01, 0.5;...                              % in-stream-velocity (in m s^-1)
    0.01, 10];                                % hill slope velocity (in m s^-1)

X1 = LHsampAG(n,size(par_range,1),par_range);
X = bsxfun(@minus,X1,m);
X = bsxfun(@rdivide,X,stand_dev);

% select 1000 samples which are at maximum distance from the samples in the
% training data
for pind = 1:n
    
    Xt =  X(pind,:);
    d = bsxfun(@minus,Xtrain,Xt);
    d = sum(d.^2,2);
    d_min(pind) = min(d);
    
end

X1 = [d_min',X1];
X1 = sortrows(X1);
X1(:,1) = [];
X1 = X1(n-500+1:n,:);
param_samps = [param_samps;X1];

%% draw samples from 4th region
n =10000;
par_range = [1, 2000;...                      % storage capacity (in mm)
    0.01, 2;...                              % storage distribution parameter
    -15, -2;...                               % logarithm of baseflow reservior constant (in s^-1)
    0.01, 10;...                              % in-stream-velocity (in m s^-1)
    0.01, 10];                                % hill slope velocity (in m s^-1)

X1 = LHsampAG(n,size(par_range,1),par_range);
X = bsxfun(@minus,X1,m);
X = bsxfun(@rdivide,X,stand_dev);

% select 1000 samples which are at maximum distance from the samples in the
% training data
for pind = 1:n
    
    Xt =  X(pind,:);
    d = bsxfun(@minus,Xtrain,Xt);
    d = sum(d.^2,2);
    d_min(pind) = min(d);
    
end

X1 = [d_min',X1];
X1 = sortrows(X1);
X1(:,1) = [];
X1= X1(n-500+1:n,:);
param_samps = [param_samps;X1];

param_samps(:,3) = 10.^param_samps(:,3);

save('GP_optim_data_additional','param_samps');
%}