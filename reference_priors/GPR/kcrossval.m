% Leave-one-out cross-validation routine for GPR
%%% Inputs:
% X=(nxd) data matrix with n examples and d features to be used for
% training
% y=an n column vector of noisy targets
% M=an inverse covariance matrix of features (to be used for covariance between samples)
% l=length-scale (scalar)
% sig2=variance of random noise (scalar)
% sigf2=signal-variance (scalar)
%%% Outputs:
% z=a k-vector of negative-log-likelihood of for different folds
% Note: Woodbury-Sherman-Morrison formula was used for compuatational
% efficiency (Williams and Rasmussen, Chapter:5)

function [L]=kcrossval(X,y,M,l,sig2,sigf2)
[n,~]=size(X);
K=KerComp(X,X,M,l,sigf2)+sig2*eye(n);              % computation of covariance matrix (of response)
C=chol(K,'lower');
Kinv=C'\(C\eye(n)); 
alpha=C'\(C\y); 
V1=(diag(diag(Kinv)));                             % diagonal matrix with (i,i) entry equal to 1/variance of ith response
V=V1^(-1);                                         % diagonal covariance matrix of response
mu=y-V*alpha;
Q=(y-mu)'*V1*(y-mu);
L=0.5*sum(log(diag(V)))+0.5*Q+0.5*n*log(2*pi);

end
