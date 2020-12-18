% objective function defined as negative-log likelihood of the
% observations; the parameter of covariance matrix and the hydrologic model
% are to be estimated
% input: param: param(1)=Cmax (maximum storage capacity)
%               param(2)=b (spatial variation parameter of storage)
%               capacity distribution
%               param(3)=Kb (baseflow reservoir constant)
%               param(4)=vel_stream (in-stream velocity of water-drops)
%               param(5)=vel_hillslope (hillslope velocity of water-drops)
%               param(6)=correlation parameter in square exponential kernel
%               param(7)=variance of structural error
%               param(8)=variance of observation errors
% output: negative-log likelihood evaluated at given 'parameters' 
% reference: Kennedy and O'Hagan (2001).
function L = objectiveComputation(param,GLOBAL_DATA,GEOMORPH)


strmobs=GLOBAL_DATA.strmobs;

% parameters to be optimized
theta = param(1:5);
theta(3) = 10^theta(3);
alpha = param(6);
sigma1_sq = param(7);
sigma2_sq = param(8);

% computation of residuals
strmsim=int_pdm_giuh(theta,GLOBAL_DATA,GEOMORPH);
minlen=min(length(strmsim),length(strmobs));
strmsim=strmsim(1:minlen)';           % conversion from row to column vector
strmobs_temp=strmobs(1:minlen);
err=errcompute(strmobs_temp,strmsim);  %% computation of error vector
d = length(err);
%% sum-of square erros
%{
L=sum(err.^2);
%}
%% Gaussian process covariance matrix (Kennedy and O'Hagan, 2001)
t = 1:d;
diff_mat = repmat(t',1,length(t))-repmat(t,length(t),1);
diff_mat = diff_mat.^2;

cov_1 = sigma1_sq*exp(-alpha*diff_mat);     % covariance matrix of Gaussian process
cov_1 = cov_1 + sigma2_sq*(eye(length(t))); % add the error varaince to cov_1
L1 = chol(cov_1);                           % Cholesky decomposition
detlogcov = 2*sum(log(diag(L1)));           % log of determinant of covariance matrix

% compute sum of squares of errors normalized by covariance values
SS1 = L1\err;                               
SS = SS1'*SS1;

L = d/2*log(2*pi) + 0.5*detlogcov + 0.5*SS;        % negtaive of log-likelihood

end