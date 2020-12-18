% % computation of covariance matrix (of response) between points given in
% two identical matrices
% Covariance function = squared-exponential kernel
%%% Inputs:
% X1= a (n1xd) matrix containing a data-point in each row
% X2= a (n2xd) matrix containing a data-point in each row
% M = an (dxd) inverse covariance matrix with (i,j) entry equal to the covariance between ith and jth feature
% sigf2=signal variance(scalar)
% note: random-noise is not added in the implementation

function K = KerComp_iden(X1,X2,M,sigf2)
    
    n1=size(X1,1);
    n2=size(X2,1);
    
    S = X1*M*X2';
    S1 = diag(S);
    S2 = repmat(S1,1,length(S1));
    K =  S2+ S2' - 2*S;
    K = sigf2*exp(-K);
    
end