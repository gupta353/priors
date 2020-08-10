% numerical implementation of Beranrdo's algorithm to determine hydraulic-condutivity
% assuming Green-Ampt infiltration model to be a true model
% \alpha=\frac_{\sum_{i=1}{k}(x_i-\mu)^2/}{2*\sigma^2}
% initial prior used: pi_star=\partialg/\partial\theta
% model: x=g(\theta)+\epsilon, \epsilon~N(0,\sigma^2), g is the Green-Ampt
% model
% references: Berger et al. (2009), The formal definition of reference
% priors
% Chow et al. (1988), Applied Hydrology
% this scripts implements an approximate method different from that in
% Berger et al. (2009)

clear all
close all
clc

save_dir=['D:/Research/Thesis_work/Non_informative_priors',...
    '/matlab_codes/reference_priors/'];

k=10000;

% other known parameters of Green-Ampt equation
psi=16.68;          % (in cm)
delta_theta=0.340;  % change in moisture content
t=3600;   % time at which infiltration is computed (in s)
g=@(x)Green_Ampt_solution(x,psi,delta_theta,t);

kh=(0.01:0.01:1)/3600;          % hydraulic conductivity values at which prior is to be evaluated
sig2 = 1;
opt_seeds = (0.5:0.5:10)/3600;  % seeds for different optimization routines
opt_min = 0.0001/3600;          % lower bound over hydraulic conductivity (in cm s^-1)
opt_max = 11/3600;              % upper bound over hydraulic conductivity (in cm s^-1)

for count=1:length(kh)
    K = kh(count);
    sig2_tmp = sig2;                
    
    mu_tmp=g(K);                                    % mean of the distribution
    sig2_tmp=sig2;
    samps=normrnd(mu_tmp,sqrt(sig2_tmp),[k,1]);     % k samples drawn from the distribution
    
    for opt_ind = 1:20
        
        % optimization
        fun = @(x)(k/2*log(2*pi)+k/2*log(sig2_tmp)+1/2/sig2_tmp^2*sum((samps-g(x)).^2));
        Ksol(opt_ind) = fmincon(fun,opt_seeds(opt_ind),[],[],[],[],opt_min,opt_max);
        
        % computation of variance at Ksol
        fun1 = @(x)(-k/2*log(2*pi)-k/2*log(sig2_tmp)-1/2/sig2_tmp^2*sum((samps-g(x)).^2));
        hess = hessiancomp(fun1,Ksol(opt_ind),0.001/3600);
        varest(opt_ind) = -1*hess^(-1);
        
        % computation of density
        if abs(Ksol(opt_ind)-K)<10^-5
            log_post(count,opt_ind) = -1/2*log(2*pi)-1/2*log(varest(opt_ind));
        else
            log_post(count,opt_ind) = -1/2*log(2*pi)-1/2*log(varest(opt_ind))-1/2*(K-Ksol(opt_ind)).^2/varest(opt_ind);
        end
    end
end

%% normalize the density
h = kh(2)-kh(1);
log_post_max = max(log_post(:));
log_post = log_post - log_post_max;
exp_log_post = exp(log_post);
exp_log_post = mean(exp_log_post,2);
den = exp_log_post;
den(1) = den(1)/2;
den(end) = den(end)/2;
den = sum(den)*h;
p = exp_log_post/den;

%% write data to a text file
sname = 'prior_density_approx_method_sig2=1';
filename = fullfile(save_dir,'results/onwards_august_2019',sname);
fid = fopen(filename,'w');
fprintf(fid,'%s\t%s\n','Kh(cm/s)','normalized_density');
for wind = 1:length(kh)
    fprintf(fid,'%f\t%f\n',kh(wind),p(wind));
end
fclose(fid);