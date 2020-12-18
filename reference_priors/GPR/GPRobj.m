% this routine computes log-likelihood for GP model
% inputs: theta = vector of paremeters
%         Xtrain = an (ntrain x v) matrix of containing predictors of training samples in each
%                  row, where ntrain = number of trainign samples and v = number of
%                  dimensions in predictor space
%         ytrain = an (ntrain x dt) matrix containing multivariate response
%                  of  training samples in each row, where dt = numbee of time-steps
%                  or number of responses
%         Xtest = an (ntest x v) matrix containing predictors of test
%                 samples
%         ytest = an (ntest x dt) ,matrix containing response of test
%                 samples
%         sig2 = observation-noise
% outputs: L = negative log-likelihood

function L = GPRobj(theta,Xtrain,ytrain,Xtest,ytest,sig2)
    
    % read parameters
    l = 1./theta(1:5);              % length scale of five normalized predictor variables
    M = diag(l);                    
    sigf2 = theta(6);               % signal variance
%     Mt = 1/theta(7);                % length scale of respnse in time dimension
%     sig2t = theta(8);               % signal variance in time-dimensions
    
    ntr = size(Xtrain,1);             % number of training samples
    ntest = size(Xtest,1);      % number of test samples
    dt = size(ytest,2);               % number of time-steps in response
    
    ytest = ytest';
    
    % computation of covariance matrix (of response) between training data
    Ktrain=KerComp_iden(Xtrain,Xtrain,M,sigf2)+sig2*eye(ntr);
    
    % computation of covariance matrix (of response) between training and test
    % points
    Ktrt=KerComp(Xtrain,Xtest,M,sigf2);
    
    % computation of covariance matrix (of response) between test data
    Ktest=KerComp_iden(Xtest,Xtest,M,sigf2);
    
    % computation of covariance matrix in time dimension
    %Kt = KerComp_iden((1:dt)',(1:dt)',Mt,sig2t) + 0.0001*eye(dt);
    
    % GPR implementation
    [ftest,V]=GPR(Xtrain,ytrain,Ktrain,Ktrt,Ktest); % ftest =  mean prediction; C = covariance matrix multiplier
    %V = V + 0.0001*eye(ntest);
    
    % sum of square error
    %
    sq = (ytest-ftest).^2;
    L = sum(sq(:));
    %}
    
    % negative-log-likelihood error
    %{
    % computation of sum-of-square
    Lkt = chol(Kt);  
    Lc = chol(V);     % Cholesky decompositions
    beta = Lkt\(ytest-ftest);
    SS = beta'*beta;
    SS = Lc'\(Lc\eye(ntest)).*SS;
    SS = sum(SS(:));
    
    % computation of log of determinant
    logdetKt = sum(log(diag(Lkt)));
    logdetV = sum(log(diag(Lc)));
    logdet = dt*logdetV + ntest*logdetKt;
    
    L = dt*ntest/2*log(2*pi) + 0.5*logdet + 0.5*SS; % negative log-likelihood
    
    if L > 4.689*10^8
        
        L = inf;
        
    end
    %}
end