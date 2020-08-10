% Estimation of hydraulic conductivity using Green-Ampt model on actual
% dataset, using only one sample; change in water content \Delta\theta is
% assumed to be known
% Refs: Green-Ampt model: Chow et al. (1988)
%             Schmidt, J. (2002). A model for transient flow to a subsurface tile
%             drain. Chapter 2. Masters Thesis
%
%       Data: Stillman, J.S., N.W. Haws, R.S. Govindaraju, and P.S.C. Rao (2006).
%             “A model for transient flow to a subsurface tile drain under
%             macropore-dominated flow conditions.” Journal of Hydrology, 317,
%             49-62.
%             Schmidt, J. (2002). A model for transient flow to a subsurface tile
%             drain. Chapter 2. Masters Thesis.
% this scripts implements an approximate method different from that in
% Berger et al. (2009)

clear all
close all
clc

save_dir=['D:/Research/Thesis_work/Non_informative_priors',...
    '/matlab_codes/reference_priors'];

k=1000;          % number of samples to be drawn in each set
m=10000;          % number of sets of samples to be drawn

% other known parameters of Green-Ampt equation for outer ring data
psi=50;          % (in cm)
delta_theta=0.004787;  % change in moisture content
t=1148;             % time at which infiltration is computed (in s)
H0=13.7;            % Intital water surface level from ground (cm)

kh=(0.1:0.1:50)/3600;      % hydraulic conductivity values (in cm s^-1) at which the prior is to be evaluated
sig2=1;           % variance values (cumulatve infiltration) at which prior is to be evaluated
g = @(x)falling_head_Green_Ampt_solution(x,psi,delta_theta,H0,t);
options = optimset('tolFun',10^-12,'MaxIter',10000,'MaxFunEvals',10000);
opt_seeds = (0.001:5:50)/3600;      % seeds for optimization routines
opt_min = 0.0001/3600;          % lower bound over hydraulic conductivity (in cm s^-1)
opt_max = 51/3600;              % upper bound over hydraulic conductivity (in cm s^-1)

for count=1:length(kh)
    
    K = kh(count);
    sig2_tmp = sig2;
    mu_tmp=g(K);                                    % mean of the distribution
    samps=normrnd(mu_tmp,sqrt(sig2_tmp),[k,1]);     % k samples drawn from the distributiion
    
    for opt_ind = 1:length(opt_seeds)
        
        fun = @(x)(k/2*log(2*pi)+k/2*log(sig2_tmp)+1/2/sig2_tmp^2*sum((samps-g(x)).^2));
        Ksol(opt_ind) = fmincon(fun,opt_seeds(opt_ind),[],[],[],[],opt_min,opt_max,[],options);
        fun1 = @(x)(-k/2*log(2*pi)-k/2*log(sig2_tmp)-1/2/sig2_tmp^2*sum((samps-g(x)).^2));
        hess = hessiancomp(fun1,Ksol(opt_ind),0.0001/3600);
        varest(opt_ind) = -1*hess^(-1);
        
        if abs(K-Ksol(opt_ind))<10^-4
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
sname = 'FHGA_prior_density_approx_method_sig2=1';
filename = fullfile(save_dir,'results/onwards_august_2019',sname);
fid = fopen(filename,'w');
fprintf(fid,'%s\t%s\n','Kh(cm/s)','normalized_density');
for wind = 1:length(kh)
    fprintf(fid,'%f\t%f\n',kh(wind),p(wind));
end
fclose(fid);