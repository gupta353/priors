% covariance (of response) matrix computation for squared-exponential
% covariance function
%%% Inputs: x1 and x2= two column vectors (of same length)
%           M = a (dxd) inverse covariance matrix with (i,j) entry equal to the
%               covariance between ith and jth feature
%           sigf2 = variance (scalar)
%%% Output: K = covariance between response at given inputs
% Note: Observation variance has not been added to covariance
%%% References
% Williams and Rasmussen, (2006). Gaussian Process for Regression
function K=KernelExp(x1,x2,M,sigf2)
    
    d=x1-x2;
    K=sigf2*exp(-d'*M*d);
    
end



