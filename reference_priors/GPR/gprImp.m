%% implementation of GPR
% this function calls all the necessary functions to implement GPR
% observation noise is zero in this example

clear all
close all
clc

% read data
direc = 'D:/Research/Thesis_work/Non_informative_priors/matlab_codes/reference_priors';
fname = 'GP_optim_data.mat';
filename = fullfile(direc,'results/pdm_giuh',fname);
load(filename);

% remove training and test samples with values log_ks values less than
% -5.3
ind = find(log(Xtrain(:,3))/log(10)<-6);
Xtrain = Xtrain(ind,:);
ytrain = ytrain(ind,:);

ind = find(log(Xtest(:,3))/log(10)<-6);
Xtest = Xtest(ind,:);
ytest = ytest(ind,:);

% observation variance
sig2 = 0.0001;          % a small observation variance is added to keep the covariance matrix positive definite

% read first nr data points trainign data and rest as test data
% rng(1);
% test_ind = randsample(1:11000,1000);
% train_ind = setdiff(1:11000,test_ind);
% X = param_samps(train_ind,:);
% ytrain = streamflow_samps(train_ind,:);
X = Xtrain;
[ntr,d] = size(X);
dt = size(ytrain,2);

% read next nr data points points as test data
% Xtest = param_samps(test_ind,:);
% ytest = streamflow_samps(test_ind,:);

% normalize the predictor space
stand_dev = std(X);
m = mean(X);
Xtrain = bsxfun(@minus,X,m);
Xtrain = bsxfun(@rdivide,Xtrain,stand_dev);

Xtest = bsxfun(@minus,Xtest,m);
Xtest = bsxfun(@rdivide,Xtest,stand_dev);

%% GPR at optimal parameters
%{
% sigf2 = 122.6309;
% l = 1./[30.2480   87.4062    0.0145   44.3737   96.4016];
% M = diag(l);
theta = [10.9954892460445,42.9464729186181,0.000100000000315161,8.39768787181954,57.0322299756572,38.9699184188749,0.688245726769408,100.957704540310];
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);

% computation of covariance matrix (of response) between training data
profile on;
tic;
Ktrain=KerComp_iden(Xtrain,Xtrain,M,sigf2)+sig2*eye(ntr);

% computation of covariance matrix (of response) between training and test
% points
Ktrt=KerComp(Xtrain,Xtest,M,sigf2);

% computation of covariance matrix (of response) between test data
Ktest=KerComp_iden(Xtest,Xtest,M,sigf2);

% computation of covariance matrix in time dimension
% Kt = KerComp_iden((1:dt)',(1:dt)',Mt,sig2t);

% GPR implementation
[ftest,V]=GPR(Xtrain,ytrain,Ktrain,Ktrt,Ktest);
toc;
profile off;
% rmse=sqrt(mean((ytest-ftest).^2));
%}

% parameter optimization
%
loss=@(theta)GPRobj(theta,Xtrain,ytrain,Xtest,ytest,sig2);      % loss function
% simulated annealing
parent = [1,1,1,1,1,10,2,10];
% parent = [4.47607032455577,57.4743548804763,0.000100000000000110,2.79062368567944,28.0900342186303,55.1299559566544];
lb=[0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001];                 % lower bound
ub=[1000,1000,1000,1000,1000,1000,1000,1000];                                  % upper bound
options = optimset('TolFun',10^-8,'MaxFunEvals',10);
tic;
[theta_opt,fval] = simulannealbnd(loss,parent,lb,ub,options);
toc;
% save('optimal_GP_param_11000_2','theta_opt','fval');
%}

%% debug positive definiteness error
%{
lb=[0.0001,0.0001,0.0001,0.0001,0.0001,0.0001];           % lower bound
ub=[100,100,100,100,100,1000];                                  % upper bound

count = 0;
while count<inf
theta = unifrnd(lb,ub);
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);

% computation of covariance matrix (of response) between training data
tic;
Ktrain=KerComp(Xtrain,Xtrain,M,sigf2)+sig2*eye(ntr);

% computation of covariance matrix (of response) between training and test
% points
Ktrt=KerComp(Xtrain,Xtest,M,sigf2);

% computation of covariance matrix (of response) between test data
Ktest=KerComp(Xtest,Xtest,M,sigf2);

% GPR implementation
[ftest,V]=GPR(Xtrain,ytrain,Ktrain,Ktrt,Ktest);

count = count+1;
end
%}

% test the effectiveness of GPR on test samples
% read data
%{
fname = 'pdm_giuh_runs_1000_test.mat';
filename = fullfile(direc,'results/pdm_giuh',fname);
load(filename);

% remove training and test samples with values log_ks values less than
% -5.3
ind = find(log(XTEST(:,3))/log(10)>-5.3);
XTEST = XTEST(ind,:);
YTEST = YTEST(:,ind);

% normalize the XTEST
XTEST = bsxfun(@minus,XTEST,m);
XTEST = bsxfun(@rdivide,XTEST,stand_dev);

% parameters
%theta = [34.2597102882363,21.5131335625532,0.000100000000002855,16.9366278337041,15.8116876803111,69.6723723138751,0.684457347826273,76.7386368080396];
%theta = [61.0752102913007,57.8158266606964,0.0109558571044910,42.3504230297851,117.567889567609,59.7881438784855,80.1902988798573,67.4976467147845];
theta = [35.6344532939392,23.2478620583043,0.0145223340005544,28.8085373555612,27.5225839263258,4.26879082052562,7.13740873772817,11.5963839457466];
sigf2 = theta(6);
l = 1./theta(1:5);
M = diag(l);

tic;
% computation of covariance matrix (of response) between training data
Ktrain=KerComp_iden(Xtrain,Xtrain,M,sigf2)+sig2*eye(ntr);

% computation of covariance matrix (of response) between training and test
% points
Ktrt=KerComp(Xtrain,XTEST,M,sigf2);

% computation of covariance matrix (of response) between test data
Ktest=KerComp_iden(XTEST,XTEST,M,sigf2);

% GPR implementation
[ftest,V]=GPR(Xtrain,ytrain,Ktrain,Ktrt,Ktest);
toc;

% plot observed vs. predicted data

for ind = 1:size(XTEST,1)
    SS = sum((YTEST(:,ind)-ftest(:,ind)).^2);
    Var = sum((YTEST(:,ind)-mean(YTEST(:,ind))).^2);
    NSE(ind) = 1-SS/Var;
%     scatter(YTEST(:,ind),ftest(:,ind),'filled'); hold on
%     xlim([0 500]); ylim([0 500]); hold on
%     plot([0 500], [0 500],'color','black');
%     title(['NSE = ',num2str(NSE(ind))],'fontname','arial','fontsize',12);
%     pause(1); hold off
end

%}