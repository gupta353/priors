% this script computes log of multivariate normal distribtuion
% Inputs: x = point at which log of pdf is to be computed (row vector)
%         mean_vec = mean of the pdf (row vector)
%         cov_mat = covariance matrix of the pdf
% Outpts: logpdf = log of pdf

function logpdf = logmvnpdf(x,mean_vec,cov_mat)
    
    n = length(x);  % number of variables
    
    % log of determinant of covariance matrix
    L = chol(cov_mat);
    logdet = 2*sum(log(diag(L)));
    
    % log of pdf
    err = x-mean_vec;
    if err*err'<10^-5
        logpdf = -1/2*log(2*pi) - 1/2*logdet;
    else
        logpdf = -n/2*log(2*pi) - 1/2*logdet - 1/2*err*cov_mat^(-1)*err';
    end
    
end