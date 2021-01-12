% computation of covariance matrix (of response) between points given in
% two (identical or non-identical) matrices
% Covariance function = squared-exponential kernel
%%% Inputs:
% X1= a (n1xd) matrix containing a data-point in each row
% X2= a (n2xd) matrix containing a data-point in each row
% M = an (dxd) inverse covariance matrix with (i,j) entry equal to the covariance between ith and jth feature
% sigf2=signal variance(scalar)
% note: (1) random-noise is not added in the implementation
%       (2) this script is written to compute covariance matrix in
%       prediction mode; therefore, the matrix S1 (for traning predictors)
%       is already computed which accelerates the prediction


function K=KerComp_pred(X1,X2,S1,M,sigf2)
    
    n1=size(X1,1);
    n2=size(X2,1);
       
    S1 = repmat(S1,1,n2);
    
    S2 = X2*M*X2';
    S2 = diag(S2);
    S2 = repmat(S2',n1,1);
    
    S3 = X1*M*X2';
    
    K = S1 + S2 - 2*S3;
    K = sigf2*exp(-K);
end