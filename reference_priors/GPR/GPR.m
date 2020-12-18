% Implementation of Gaussian Process Regression
%%% Inputs:
% Xtrain=(nxd) data matrix with n examples and d features to be used for
% training
% ytrain=n-vector of noisy targets to be used for training
% Ktrain=covariance matrix (of response) of the data-matrix Xtrain
% Ktest=covariance matrix (of response) of test matrix Xtest
% Ktrt=matrix of covariance (of response) between test points and training-points
% note: size of ktrt should be (# of training samples x # of test samples)
% sig2=variance of random noise
%%% Outputs:
% ftest=mean of the predictive distribution
% V=variance of the predictive distribution
%%% Reference:
% Williams and Rasmussen, (2006). Gaussian Process for Regression

function [ftest,V]=GPR(Xtrain,ytrain,Ktrain,Ktrt,Ktest)
    
    % for scalar valued response
%     [~,d]=size(Xtrain);
%     L=chol(Ktrain,'lower');         % cholesky decomposition of covariance matrix
%     alpha=L'\(L\ytrain);
%     ftest=Ktrt'*alpha;              % mean of the predictive distribution
%     v=L\Ktrt;
%     V=Ktest-v'*v;                   % predictive variance

    % for vector valued response
    % mean estimation
    [~,d]=size(Xtrain);
    L=chol(Ktrain,'lower');         % cholesky decomposition of covariance matrix
    alpha=L'\(L\ytrain);
    ftest = alpha'*Ktrt;
    
    % variance estimation
    v = L\Ktrt;
    V = Ktest - v'*v;
end
