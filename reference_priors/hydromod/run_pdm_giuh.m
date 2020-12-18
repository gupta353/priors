% Integrated PDM and GIUH
% User-inputs: save_dir=directory that contains all the relevant data about
%                       the watershed
%              wt=watershed topology matrix (see mainScript.m for more details)
%              support=length of the unit-hydrograph to be computed (in seconds; recommended: 0 to 2000000)
%              R_time_scale=the time-interval at which rainfall data is
%                           avaialable (in mins)
%              S_time_scale=the time-interval at which streamflow data is
%              available (in mins)
%              S0=intitial soil moisture depth (in mm; strongly recomemnded value: 0)
%              nodeobs=node at which streamflow data is available
% Outputs:     surfaceflow=surfaceflow (in m^3 s^-1)
%              baseflow=baseflow (in m^3 s^-1)
%              streamflow=surafceflow+baseflow (in m^3 s^-1)
clear all
close all
clc

direc = 'D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';
addpath(direc);

global GEOMORPH GLOBAL_DATA
%% User-inputs
save_dir = fullfile(direc,'huc_04100003');
wt_matrix=fullfile(save_dir,'wt_sjrw');      % watershed topology matrix 
load(wt_matrix);
wt=wt_sjrw;
support=0:2000000;           % support for instantaneous unit-hydrograph computation (in seconds)
R_time_scale=60;             % time-scale (in mins) of rainfall data (use streamflow data at hourly scale)
S_time_scale=24*60;          % time-scale (in mins) of streamflow data
S0=0;                        % intitial moisture content
nodeobs=18;

warm_cal_data = 'cal_warmup_data.mat';              % name of the mat file containing warm-up and calibration data
area_length_fname='area_length.txt';                % name of the text-file containing direc-draiange-area and link-length of each subabsin
%% direct drainage area and link-length of each sub-watershed
arealen_filename=fullfile(save_dir,area_length_fname);
fid=fopen(arealen_filename,'r');
area_length=textscan(fid,'%f%f%f','headerlines',1);
fclose(fid);
ddrain_area=area_length{2};
length_link=area_length{3};

%% rainfall data
filename = fullfile(direc,'huc_04100003',warm_cal_data);
load(filename);
P_warmup = rain_warmup;           % rainfall series in warm-up period (in mm)
P = rain_cal;                      % rainfall series in calibration period (in mm)

% replace NaNs by average of previous and next step
% for warmup period
indnan = find(isnan(P_warmup));
for ind  = 1:length(indnan)
    ind_tmp = indnan(ind);
    P_warmup(ind_tmp) = (P_warmup(ind_tmp-1) + P_warmup(ind_tmp+1))/2;
end

% for calibration period
indnan = find(isnan(P));
for ind  = 1:length(indnan)
    ind_tmp = indnan(ind);
    P(ind_tmp) = (P(ind_tmp-1) + P(ind_tmp+1))/2;
end

% generate error message if the data still contains NaNs
if ~isempty(find(isnan(P_warmup)))
    error('Rainfall data in warm-up period contains NaNs');
end

if ~isempty(find(isnan(P)))
    error('Rainfall data in calibration period contains NaNs');
end

%% evaporation data
%%% for warm-up
% coversion of evaporation data from daily to hourly scale
Ea_hourly=evap_warmup'/24;
Ea_hourly=repmat(Ea_hourly,24,1);
Ea_warmup=reshape(Ea_hourly,[24*length(Ea_hourly),1]);

%%%% for calibration and posterior computation
% coversion of evaporation data from daily to hourly scale
Ea_hourly=evap_cal'/24;
Ea_hourly=repmat(Ea_hourly,24,1);
Ea=reshape(Ea_hourly,[24*length(Ea_hourly),1]);

%% streamflow data
strmobs = strm_cal;
%% Run the mainScript
GEOMORPH=mainScript(wt,ddrain_area,length_link,save_dir);

%% read the text-file containing the list of links draining into nodes
filename=fullfile(save_dir,'links_draining_into_each_node.txt');
fid=fopen(filename,'r');
node_link=textscan(fid,'%s%s');
links=node_link{2}{nodeobs+1};
links=strsplit(links,',');
links=cellfun(@str2num,links);

%% total drainage of user-supplied node
darea=sum(ddrain_area(links));          % drainage of nodeobs (in hectares)
%% data to pass
GLOBAL_DATA.save_dir=save_dir;
GLOBAL_DATA.R=P;
GLOBAL_DATA.R_warmup=P_warmup;
GLOBAL_DATA.darea=darea;
GLOBAL_DATA.nodeobs=nodeobs;
GLOBAL_DATA.strmobs=strmobs;
GLOBAL_DATA.support=support;
GLOBAL_DATA.R_time_scale=R_time_scale;
GLOBAL_DATA.S_time_scale=S_time_scale;
GLOBAL_DATA.Ea=Ea;
GLOBAL_DATA.Ea_warmup=Ea_warmup;
GLOBAL_DATA.S0=S0;

%% parameters and other important variables
%
% Cmax=param(1);        % maximum storage capacity (in mm)
% b=param(2);           % spatial variation parameter
% Kb=param(3);          % constant of baseflow reservoir (in s^-1)
% theta(1)=param(4);    % in-stream velocity (m s^-1)
% theta(2)=param(5);    % hillslope velocity (m s^-1)

param=[54.9070447714919,9.35621268850241,1.18071468426315e-08,3.06617953734208,2.98410560301773];
% profile on;
[streamflow,surfaceflow,baseflow,S,V,C,gam,Ea_updated]=int_pdm_giuh(param,GLOBAL_DATA,GEOMORPH);
% profile off; profile report
%}

%% define prior distribution parameters
%{
cmax_l = 1; cmax_u = 2000;                      % storage capacity (in mm)
b_l = 0.01; b_u = 10;                           % storage distribution parameter
log_k_l = -15;   log_k_u = -2;                  % logarithm of baseflow reservior constant (in s^-1)
vs_l = 0.01; vs_u = 10;                         % in-stream-velocity (in m s^-1)
vh_l = 0.01; vh_u = 10;                         % hill slope velocity (in m s^-1)
alpha_l = 0.01; alpha_u = 7;                   % correlation parameter
log_sigma1_sq_l = 0.01; log_sigma1_sq_u = 7;    % log structural uncertainty variance
log_sigma2_sq_l = 0.01; log_sigma2_sq_u = 7;    % log observation uncertainty variance
%}
%% parameter optimization
%{
loss=@(param)objectiveComputation(param,GLOBAL_DATA,GEOMORPH);     % loss function
nopts = 2;         % number of optimizations

for opt_ind = 1:nopts
    
    parent=priorrnd(cmax_l,cmax_u,b_l,b_u,log_k_l,log_k_u,vs_l,vs_u,vh_l,vh_u,alpha_l,alpha_u,log_sigma1_sq_l,log_sigma1_sq_u,log_sigma2_sq_l,log_sigma2_sq_u);
    % initial guess
    
    % sceua method
    %{
lb=[500,0.01,-14,0.1,0.1,0.01,0.01,0.01,0.01,0.01,0.5,0.1];           % lower bound
ub=[1500,15,-4,4,10,10,1,1,1,1,1,100];                               % upper bound
[bestx,bestf] = sceua(parent,lb,ub,100000,5,0.01,[],40,2000,1);
    %}
    % simulated annealing
    %
    lb=[100,0.01,-14,0.1,0.1,0.0001,0.00001,0.00001];           % lower bound
    ub=[2000,15,-4,4,10,100,100,100];                              % upper bound
    options = optimset('TolFun',10^-4,'MaxFunEvals',200);
    [param(opt_ind,:),fval(opt_ind)] = simulannealbnd(loss,parent,lb,ub,options);
    
    % computation of probability
    fval = -fval;
    fval = fval - max(fval);
    probs = exp(fval)/sum(fval);
    
end
%}

%% computation of uncertainty band using the parameters obtained by
% optimization
%{
% parameters to be optimized
param = minimum;
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

% compute covriance matrix
t = 1:d;
diff_mat = repmat(t',1,length(t))-repmat(t,length(t),1);
diff_mat = diff_mat.^2;

cov_1 = sigma1_sq*exp(-alpha*diff_mat);     % covariance matrix of Gaussian process
cov_1 = cov_1 + sigma2_sq*(eye(length(t))); % add the error varaince to cov_1

% draw samples from multivraite Gaussian with mean 'strmsim' and covariance 'cov_1'
samps = mvnrnd(strmsim,cov_1,100000);
strm_l = prctile(samps,0.05);
strm_u = prctile(samps,99.95);
plot(strmobs,'o'); hold on
plot(strmsim,'r');
plot(strm_l); plot(strm_u);
%}
%% plots of surface- and base-flow
%{
subplot(2,1,1)
plot(streamflow)
hold on
plot(strmobs,'r')
xlabel('time-step (days starting from 01/01/2005)',...
    'fontname','arial','fontsize',12)
ylabel('streamflow (m^3 s^{-1})',...
    'fontname','arial','fontsize',12)
set(gca,'fontname','arial','fontsize',12)

title(['C_{max}=',char(32),num2str(param(1)),char(32),'mm,',char(32),...
    'b=',char(32),num2str(param(2)),',',char(32),...
    'K_b=',char(32),num2str(param(3)),char(32),'s^{-1},',char(32),...
    'v_s=',char(32),num2str(param(4)),char(32),'m s^{-1},',char(32),...
    'v_h=',char(32),num2str(param(5)),char(32),'m s^{-1},'],...
    'fontname','arial','fontsize',12);

legend('simulated','observed','fontname','arial','fontsize',12);
legend('boxoff')

subplot(2,1,2)
plot(S)
xlabel('time-step (days starting from 01/01/2005)',...
    'fontname','arial','fontsize',12)
ylabel('storage depth (mm)',...
    'fontname','arial','fontsize',12)
set(gca,'fontname','arial','fontsize',12)
%}
%% mass balance
%{
% computation of total rainfall depth
R=reshape(P,24,length(strmobs));
daily_P=sum(R,1);               % rainfall at daily scale
cumulative_rainfall=cumsum(daily_P); % in mm

% evaporation
daily_Ea_updated=reshape(Ea_updated,24,length(strmobs));
daily_Ea_updated=sum(daily_Ea_updated,1);
cumulative_evaporation_depth=cumsum(daily_Ea_updated); % in m^3

% compute soil mositure depth at daily scale
S_initial=S(1);
S(1)=[];
S_daily=reshape(S,24,length(strmobs));
S_daily=S_daily(24,:);  % storage at the end of each day (in mm)
cumulative_delta_S=S_daily-S_initial;
delta_S=S_daily-[S_initial,S_daily(1:end-1)];

% computation of cumulative excess rainfall
reshape_V=reshape(V,24,length(V)/24);
daily_V=sum(reshape_V);
cumulative_V=cumsum(sum(reshape_V));

% computation of cumulative baseflow rainfall
baseflow_depth=24*3600*10^3*baseflow/darea/10^4;     % factor 10^4 converts area  from hectare to m2
cumulative_baseflow_depth=cumsum(baseflow_depth);
%}

%% computation of posterior using DREAM
%{
profile on
global DREAM_dir EXAMPLE_dir CONV_dir
DREAM_dir=['D:/Research/EPA_Project/input_'...
    'uncertainty/2019_01_16/Dream'];
CONV_dir=['D:/Research/EPA_Project/input_'...
    'uncertainty/2019_01_16/Dream/diagnostics'];
EXAMPLE_dir = ['D:/Research/EPA_Project/input_'...
    'uncertainty/2019_01_16/generalized_gaussian'];
addpath(DREAM_dir)
addpath(CONV_dir)
addpath(EXAMPLE_dir)

Func_name='GenGausspdf';                    % target pdf
DREAMPar.d=12;                               % dimensionality of the problem
DREAMPar.N=8;                               % Number of chains
DREAMPar.T=500;                             % Number of generations
DREAMPar.lik=2;                             % Type of likelihood (2 indicates log-likelihood)

% options.parallel='yes';
% options.modout='yes';

% Meas_info.Y=strmobs;

Par_info.initial='uniform';                 % type of prior 
Par_info.min=[100,0.001,-12,0.01,0.01,...           % lower bound on parameters
    0.001,0.001,0.001,0.001,0.001,0.5,0.00001];                   
Par_info.max=[1500,5,-4,3,1,...
    1,1,1,1,1,1,10];          % upper bound
Par_info.boundhandling='fold';              % type of boundhandling

[chain,output,fx]=DREAM(Func_name,...       % main function to Run Dream
    DREAMPar,Par_info);
profile off; profile report
%}