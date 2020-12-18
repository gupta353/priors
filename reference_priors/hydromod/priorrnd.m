% this routine draws samples from prior distribution of model and error
% parameters
% inputs: cmax_l and cmax_u = lower and upper bounds over storage capacity (in mm)
%         b_l and b_u = lower and upper bounds over storage distribution parameter
%         cmax_l and cmax_u = lower and upper bounds over storage capacity (in mm)
%         b_l and b_u = lower and upper bounds over storage distribution parameter
%         log_k_l and log_k_u = lower and upper bounds over logarithm of baseflow reservior constant (in s^-1)
%         vs_l and vs_u = lower and upper bounds over in-stream-velocity (in m s^-1)
%         vh_l and vh_u = lower and upper bounds over hill slope velocity (in m s^-1)
%         alpha_l and alpha_u = lower and upper bounds over correlation parameter
%         log_sigma1_sq_l and log_sigma1_sq_u = lower and upper bounds over log structural uncertainty variance
%         log_sigma2_sq_l and log_sigma2_sq_u = lower and upper bounds over log observation uncertainty variance
% outputs: param = parameters in order of inputs described
function param = priorrnd(cmax_l,cmax_u,b_l,b_u,log_k_l,log_k_u,vs_l,vs_u,vh_l,vh_u,alpha_l,alpha_u,log_sigma1_sq_l,log_sigma1_sq_u,log_sigma2_sq_l,log_sigma2_sq_u)
    
    % cmax_l = 1; cmax_u = 2000;                      % storage capacity (in mm)
    % b_l = 0.01; b_u = 10;                           % storage distribution parameter
    % log_k_l = -15;   log_k_u = -2;                  % logarithm of baseflow reservior constant (in s^-1)
    % vs_l = 0.01; vs_u = 10;                         % in-stream-velocity (in m s^-1)
    % vh_l = 0.01; vh_u = 10;                         % hill slope velocity (in m s^-1)
    % alpha_l = 0.01; alpha_u = 7;                   % correlation parameter
    % log_sigma1_sq_l = 0.01; log_sigma1_sq_u = 7;    % log structural uncertainty variance
    % log_sigma2_sq_l = 0.01; log_sigma2_sq_u = 7;    % log observation uncertainty variance
    
    % model parameter
    cmax = unifrnd(cmax_l,cmax_u);
    b = unifrnd(b_l,b_u);
    log_k = unifrnd(log_k_l,log_k_u);
    vs = unifrnd(vs_l,vs_u);
    vh = unifrnd(vh_l,vh_u);
    
    % correlation parameter log-uniform prior
    alpha = exp(unifrnd(alpha_l,alpha_u));
    
    % variance of structural uncertainty Jeffery's prior (log-uniform prior)
    sigma1_sq = exp(unifrnd(log_sigma1_sq_l,log_sigma1_sq_u));
    
    % variance of observation uncertainty Jeffery's prior (log-uniform prior)
    sigma2_sq = exp(unifrnd(log_sigma2_sq_l,log_sigma2_sq_u));
    
    param = [cmax,b,log_k,vs,vh,alpha,sigma1_sq,sigma2_sq];
    
end