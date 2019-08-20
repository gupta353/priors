% This routine computes Bernado prior on hydrualic conductivity and
% statndard deviation such that cumulative infiltration data available at
% multiple time-steps is utilized
% Refs: Green-Ampt model: Chow et al. (1988)
%             Schmidt, J. (2002). A model for transient flow to a subsurface tile
%             drain. Chapter 2. Masters Thesis
%
%       Data: Stillman, J.S., N.W. Haws, R.S. Govindaraju, and P.S.C. Rao (2006).
%             �A model for transient flow to a subsurface tile drain under
%             macropore-dominated flow conditions.� Journal of Hydrology, 317,
%             49-62.
%             Schmidt, J. (2002). A model for transient flow to a subsurface tile
%             drain. Chapter 2. Masters Thesis.

clear all
close all
clc

save_dir=['D:/Research/Thesis_work/Non_informative_priors',...
    '/matlab_codes/reference_priors/'];

k=1000;          % number of samples to be drawn in each set
m=1000;          % number of sets of samples to be drawn
unif_integration_samples=500;
                
% raed the data
fname='Infiltrometer_A_outer_ring.txt';
filename=fullfile(save_dir,'data','infiltration_data',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%f%f','delimiter','\t','headerlines',1);
fclose(fid);
t=data{1};              % time-steps (in secs) at which data is available

% other known parameters of Green-Ampt equation for outer ring data
psi=50;          % (in cm)
delta_theta=0.004787;  % change in moisture content
H0=13.7;            % Intital water surface level from ground (cm)
g=@(x)falling_head_Green_Ampt_solution(x,psi,...        % a fucntion handle for Green-Ampt model
                delta_theta,H0,t);

kh=(1)/3600;      % hydraulic conductivity values (in cm s^-1) at which the prior is to be evaluated
sig2=1;           % variance values (cumulatve infiltration) at which prior is to be evaluated

for i=1:length(kh)
    for ii=1:length(sig2)
        
         K=kh(i);
         mu_tmp=g(K);                         % mean of the distribution
         sig2_tmp=sig2(ii);                              % variance of the distribution
         sigma=diag(sqrt(sig2_tmp)*ones(length(t),1));   % element-wise square root of covariance matrix
            
        for j=1:m
           
            samps=mvnrnd(mu_tmp,sigma,k);               % k samples drawn from the distributiion  
            tmp_matrix=bsxfun(@minus,samps,mu_tmp');
            tmp_matrix=tmp_matrix.^2;
            T1=sum(tmp_matrix(:));
            T1=-1/2/sig2_tmp*T1;
            
            fun=@(theta)-1/2/sig2_tmp*...
                sum(sum(...
                (bsxfun(@minus,samps,g(theta)').^2)...
                ));
            unif_samp=unifrnd(min(kh),max(kh),[unif_integration_samples,1]);
            for fun_i=1:length(unif_samp)
                q(fun_i)=fun(unif_samp(fun_i));
            end
            T2=log(max(kh)-min(kh)/unif_integration_samples)+max(q);
            
            log_asymp_post(j)=T1+T2;        % log of asymptotic posterior at the given point in parameter space
        end
        E_log_asymp_post(i,ii)=sum(log_asymp_post)/m;                   % expectation of log of asymptotic posterior at the given point in parameter space
    end
end