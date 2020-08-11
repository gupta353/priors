% this script computes log of multivariate normal distribtuion
% Inputs: x = point at which log of pdf is to be computed (row vector)
%         mean_vec = mean of the pdf (row vector)
%         cov_mat = covariance matrix of the pdf
% Outpts: logpdf = log of pdf

function logpdf = logmvnpdf(x,mean_vec,cov_mat)
    
    % log of determinant of covariance matrix
    L = chol(cov_mat);
    logdet = 2*sum(log(diag(L)));
    
    % log of pdf
    logpdf = -1/2*log(2*pi) - 1/2*logdet - 1/2*(x-mean_vec)*cov_mat^(-1)*(x-mean_vec)';

end