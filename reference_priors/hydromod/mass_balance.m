% this routine computes different inflows and outflow of the watershed
% system to close the mass balance equation

global GLOBAL_DATA GEOMORPH

save_dir=GLOBAL_DATA.save_dir;
strmobs=GLOBAL_DATA.strmobs;
support=GLOBAL_DATA.support;
R_time_scale=GLOBAL_DATA.R_time_scale;
S_time_scale=GLOBAL_DATA.S_time_scale;
Ea_warmup=GLOBAL_DATA.Ea_warmup;
S0=GLOBAL_DATA.S0;
nodeobs=GLOBAL_DATA.nodeobs;
darea=GLOBAL_DATA.darea;

% computation of total rainfall volume
P=GLOBAL_DATA.R;                % rainfall at horly time-scale
R=reshape(P,24,4383);
daily_P=sum(R,1);               % rainfall at daily scale
daily_P(isnan(daily_P))=0;      % nans are replaced by zeros
cumulative_rainfall_volume=cumsum(daily_P)*10^(-3)*darea*10^4; % in m^3

% computation of total evaporation
Ea=GLOBAL_DATA.Ea;
daily_Ea=reshape(Ea,24,4383);
daily_Ea=sum(daily_Ea,1);
cumulative_evaporation_volume=cumsum(daily_Ea)*10^(-3)*darea*10^4; % in m^3

% Posterior parameters
load(['D:/Research/EPA_Project/input_uncertainty/'...
    '2019_01_16/huc_04100003/result/DREAM_06_04_2019.mat']);

sample=chain;
chain=1;                                         % chain number for plotting
loglikeli= sample(:,end,chain);                          % log-likehood of samples
sample(:,end-1:end,:)=[];
N=size(sample,3);                                   % number of chains
d=size(sample(:,:,:),2);                            % number of parameters
d_model=5;                                          % number of model-parameter
d_error=7;                                          % number of error parameters
nsamples=size(sample,1);                            % number of samples
burn_in_thresh=1.5;                                 % threshold for R-diagnostic

% Computation of diagnostic R-statistic and determination of burn-in
%
clear W B Rd
for i=1:N
    cum_sum(:,:,i)=cumsum(sample(:,:,i));         % cumulative sum of parameters
    cum_sum2(:,:,i)=cumsum(sample(:,:,i).^2);     % cumulative sum of parameters squared
end
cum_count=repmat((1:nsamples)',1,d,N);            % cumulative count
cum_mean=cum_sum./cum_count;                      % cumulative mean
cum_mean2=cum_sum2./cum_count;                    % cumulative mean of squared sum
cum_mean(1,:,:)=[]; cum_mean2(1,:,:)=[];
cum_count(1,:,:)=[];
cum_var=cum_count./(cum_count-1).*...
    (cum_mean2-cum_mean.^2);                      % cumulative variance
W=mean(cum_var,3);                                % cumulative within-chain variance of parameters
across_chain_mean=mean(cum_mean,3);
B=N/(N-1)*mean((cum_mean-repmat...
    (across_chain_mean,1,1,N)).^2,3);             % cumulative between-chain variance

sigma2_hat=B+(nsamples-1)/nsamples*W;             
Rd=sqrt((N+1)/N*sigma2_hat./W-1/N*(1-1./cum_count(:,:,1)));

%%% determination posterior after burn-in
ind=[];
for i=1:d
    ind=[ind ,find(Rd(:,i)>burn_in_thresh, 1, 'last')];
end
burn_in_samples=max(ind);
sample=sample(burn_in_samples+1:end,:,chain);         % remove burn-in from sample and select 1st chain
loglikeli=loglikeli(burn_in_samples+1:end);         % remove burn-in from loglikeli

% computation of total rainfall depth
R=reshape(P,24,4383);
daily_P=sum(R,1);               % rainfall at daily scale
daily_P(isnan(daily_P))=0;      % nans are replaced by zeros
cumulative_rainfall=cumsum(daily_P); % in mm

% computation of streamflow and storage time-series
theta=[300,1,10^-8,1,1];
[streamflow,surfaceflow,baseflow,S,V,C,gam,Ea_updated]=int_pdm_giuh(theta,GLOBAL_DATA,GEOMORPH);

% daily and cumulative daily updated-evaporation values
daily_Ea_updated=reshape(Ea_updated,24,4383);
daily_Ea_updated=sum(daily_Ea_updated,1);
cumulative_evaporation_volume=cumsum(daily_Ea_updated)*10^(-3)*darea*10^4; % in m^3

% compute soil mositure depth at daily scale
S_initial=S(1);
S(1)=[];
S_daily=reshape(S,24,4383);
S_daily=S_daily(24,:);  % storage at the end of each day (in mm)
cumulative_delta_S=S_daily-S_initial;
delta_S=S_daily-[S_initial,S_daily(1:end-1)];

% computation of total streamflow throutgh the nodeobs
streamflow_depth=24*3600*10^3*streamflow/darea/10^4;     % factor 10^4 converts area  from hectare to m2
cumulative_streamflow_depth=cumsum(streamflow_depth);

% computation of cumulative excess rainfall
reshape_V=reshape(V,24,length(V)/24);
daily_V=sum(reshape_V);
cumulative_V=cumsum(sum(reshape_V));


% computation of cumulative baseflow rainfall
baseflow_depth=24*3600*10^3*baseflow/darea/10^4;     % factor 10^4 converts area  from hectare to m2
cumulative_baseflow_depth=cumsum(baseflow_depth);