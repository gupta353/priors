% This routine computes Bernado prior on hydrualic conductivity and
% statndard deviation such that cumulative infiltration data available at
% multiple time-steps is utilized
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

clear all
close all
clc

save_dir=['D:/Research/Thesis_work/Non_informative_priors',...
    '/matlab_codes/reference_priors/'];

k=1000;          % number of samples to be drawn in each set
m=3000;          % number of sets of samples to be drawn
unif_integration_samples=1000;

% raed the data
fname='Infiltrometer_A_outer_ring.txt';
filename=fullfile(save_dir,'data','infiltration_data',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%f%f','delimiter','\t','headerlines',1);
fclose(fid);
t=data{1};              % time-steps (in secs) at which data is available
t(1)=[];                % remove t=0

% other known parameters of Green-Ampt equation for outer ring data
psi=50;          % (in cm)
delta_theta=0.004787;  % change in moisture content
H0=13.7;            % Intital water surface level from ground (cm)
g=@(x)falling_head_Green_Ampt_solution(x,psi,...        % a fucntion handle for Green-Ampt model
    delta_theta,H0,t);

kh=(0.1:0.1:50)/3600;      % hydraulic conductivity values (in cm s^-1) at which the prior is to be evaluated
sig2=1;           % variance values (cumulatve infiltration) at which prior is to be evaluated
max_kh=50.1/3600;
min_kh=0.01/3600;
% Computation of Green-Ampt solutions for uniformly drawn samples of
% hydraulic conuctivity

unif_samp=min_kh:(max_kh-min_kh)/unif_integration_samples:max_kh;
infil=nan(unif_integration_samples,length(t));
for infil_i=1:unif_integration_samples
    infil(infil_i,:)=g(unif_samp(infil_i));
end

E_log_asymp_post=nan(length(kh),length(sig2));
for i=1:length(kh)
    for ii=1:length(sig2)
        
        K=kh(i);
        mu_tmp=g(K);                         % mean of the distribution
        sig2_tmp=sig2(ii);                              % variance of the distribution
        sigma=diag(sqrt(sig2_tmp)*ones(length(t),1));   % element-wise square root of covariance matrix
        
        log_asymp_post=nan(1,m);
        parfor j=1:m
            
            samps=mvnrnd(mu_tmp,sigma,k);               % k samples drawn from the distributiion
            tmp_matrix=bsxfun(@minus,samps,mu_tmp');
            tmp_matrix=tmp_matrix.^2;
            T1=sum(tmp_matrix(:));
            T1=-1/2/sig2_tmp*T1;
            
            q=nan(1,unif_integration_samples);
            for fun_i=1:unif_integration_samples
                tmp_matrix=bsxfun(@minus,samps,infil(fun_i,:));
                tmp_matrix=tmp_matrix.*tmp_matrix;
                q_tmp=(-1/2/sig2*sum(tmp_matrix(:)));
                q(fun_i)=exp(q_tmp-T1);
            end
            q(1)=q(1)/2; q(end)=q(end)/2;
            T2=T1+log(((max_kh-min_kh)/unif_integration_samples)*sum(q));
            
            log_asymp_post(j)=T1-T2;        % log of asymptotic posterior at the given point in parameter space
        end
        E_log_asymp_post(i,ii)=sum(log_asymp_post)/m;                   % expectation of log of asymptotic posterior at the given point in parameter space
    end
end

% unnormalized density
mean_E_log_asymp_post=sum(E_log_asymp_post(:))/numel(E_log_asymp_post);
PI=exp(E_log_asymp_post-mean_E_log_asymp_post);

%% normalized density 
% computation of normalizing contant using trapezoidal integration
if length(sig2)==1 % if the sigma is assumed to be known
    Y=PI;
    Y(1)=Y(1)/2; Y(end)=Y(end)/2;
    X=kh;
    A=sum(Y)*(kh(2)-kh(1));    % normalizing constant
else
    Y=PI;
    for kh_ind=1:length(sig2)
        tmp_Y=Y(:,1);
        tmp_Y(1)=tmp_Y(1)/2; tmp_Y(end)=tmp_Y(end)/2;
        A1(kh_ind)=(kh(2)-kh(1))*sum(tmp_Y);
    end
    A1(1)=A1(1)/2; A1(end)=A1(end)/2;
    A=sum(A1)*(sig2(2)-sig2(1));  % normalizing constant
end

% normalized density
PI=PI/A;            
plot(kh*3600,PI(:,1));

% write the data
wfname='FHGA_prior_kh_sig2=2_09_02_2019';
write_filename=fullfile(save_dir,'results/onwards_august_2019',wfname);
fid=fopen(write_filename,'wt');
for i=1:length(kh)
    for ii=1:length(sig2)
        fprintf(fid,'%f\t%f\t%f\n',kh(i),sqrt(sig2(ii)),PI(i,ii),'delimiter','\t');
    end
end
fclose(fid);

