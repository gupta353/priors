% computation of covariance matrix (of response) between points given in
% two (identical or non-identical) matrices
% Covariance function = squared-exponential kernel
%%% Inputs:
% X1= a (n1xd) matrix containing a data-point in each row
% X2= a (n2xd) matrix containing a data-point in each row
% M = an (dxd) inverse covariance matrix with (i,j) entry equal to the covariance between ith and jth feature
% sigf2=signal variance(scalar)
% note: random-noise is not added in the implementation


function K=KerComp(X1,X2,M,sigf2)
    
    n1=size(X1,1);
    n2=size(X2,1);
%     K=zeros(n1,n2);

% brute force method
%     for row=1:n1
%         for col=1:n2
%             x1=X1(row,:)';
%             x2=X2(col,:)';
%             K(row,col)=KernelExp(x1,x2,M,sigf2); % KernelExp computes squared exponential kernel
%         end
%     end

% first efficient method when X1 and X2 are identical
%     S = X1*M*X2';
%     
%     for row=1:n1
%         for col=1:n2
%             K(row,col) = S(row,row) + S(col,col) - S(row,col) - S(col,row);
%         end
%     end
%     K = sigf2*exp(-K);

% second efficient method when X1 and X2 are identical
%     S = X1*M*X2';
%     S1 = diag(S);
%     S2 = repmat(S1,1,length(S1));
%     K =  S2+ S2' - 2*S;
%     K = sigf2*exp(-K);

% first efficient method when X1 and X2 are non-identical
    S1 = X1*M*X1';
    S1 = diag(S1);
    S1 = repmat(S1,1,n2);
    
    S2 = X2*M*X2';
    S2 = diag(S2);
    S2 = repmat(S2',n1,1);
    
    S3 = X1*M*X2';
    
    K = S1 + S2 - 2*S3;
    K = sigf2*exp(-K);
end