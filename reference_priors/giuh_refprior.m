% Routine for reference prior computation

clear all
close all
clc

global GLOBAL_DATA GEOMORPH

%% input data
save_dir=['D:/Research/Thesis_work/Structural_uncertainty'...
    '/MatLab_codes/20180222/huc_0512011115'];
wt_matrix=['D:/Research/Thesis_work/Structural_uncertainty/MatLab'...
     '_codes/20180222/huc_0512011115/wt_0512011115'];
% load(['D:/Research/Thesis_work/Structural_uncertainty/MatLab'...
%     '_codes/20180222/results/05_18_2018']);

load(wt_matrix);
wt=wt_0512011115;
nodeobs=1;
support=0:2000000;           % support for instantaneous unit-hydrograph computation
R_time_scale=60;             % time-scale (in mins) of rainfall data (use streamflow data at hourly scale)
S_time_scale=24*60;             % time-scale (in mins) of streamflow data   
% param=bestx;
%% read rainfall data
p_filename=fullfile(save_dir,'rainfall_excess.txt');
fid=fopen(p_filename,'r');
R=textscan(fid,'%s');
fclose(fid);
R=R{1};
R=cellfun(@str2num,R);     % daily rainfall data

%% read the text-file containing the list of links draining into nodes
filename=fullfile(save_dir,'links_draining_into_each_node.txt');
fid=fopen(filename,'r');
node_link=textscan(fid,'%s%s');
links=node_link{2}{nodeobs+1};
links=strsplit(links,',');
links=cellfun(@str2num,links);
%% read area-length data
filename=fullfile(save_dir,'area_length.txt');
fid=fopen(filename,'r');
area_length=textscan(fid,'%s%s%s');
fclose(fid);
ddrain_area=area_length{2}(2:end);
ddrain_area=cellfun(@str2num,ddrain_area);
darea=sum(ddrain_area(links));
length_link=area_length{3}(2:end);
length_link=cellfun(@str2num,length_link);
%% read baseflow separated surface-flow data (in m^3/s)
strm_filename=fullfile(save_dir,'streamflow_data',...
    strcat('surfaceflow_event_',num2str(nodeobs),'.txt'));
fid=fopen(strm_filename,'r');
strmobs=textscan(fid,'%s%s');
fclose(fid);
time_steps=strmobs{1}(2:end);
% time_steps_date=strmobs{1}(2:end);
% time_steps_time=strmobs{2}(2:end);
% wrapper=@(x1,x2)strcat(x1,32,x2);
% time_steps=cellfun(wrapper,time_steps_date,time_steps_time,...
%     'UniformOutput',false);
wrapper_1=@(x)datenum(x,'mm/dd/yyyy');
time_steps=cellfun(wrapper_1,time_steps);
strmobs=strmobs{2}(2:end);
strmobs=0.0283168*cellfun(@str2num,strmobs); % 0.0283168 for conversion of units from ft^3/s to m^3/s
% sampling of streamflow data at hourly-time scale (average streamflow of each is measured)
%{
scale=24*60*(time_steps(2)-time_steps(1));  % in mins
skip_tsteps=60/scale;
count=0;
for i=1:round(skip_tsteps):length(strmobs)
    count=count+1;
    if length(strmobs)>i+3
        strm(count)=mean(strmobs(i:i+3));
        tsteps(count)=time_steps(i);
    end
end
strmobs=strm';
time_steps=tsteps';
%}
%% data to pass
GLOBAL_DATA.save_dir=save_dir;
GLOBAL_DATA.R=R;
GLOBAL_DATA.darea=darea;
GLOBAL_DATA.nodeobs=nodeobs;
GLOBAL_DATA.strmobs=strmobs;
GLOBAL_DATA.support=support;
GLOBAL_DATA.R_time_scale=R_time_scale;
GLOBAL_DATA.S_time_scale=S_time_scale;
%% Run the mainScript
GEOMORPH=mainScript(wt,ddrain_area,length_link,save_dir);
%% Reference prior computations
k=100;
m=100;
fac=5;
% parameters 
% vs (1:5), vh(0.2:0.2:1), sig2 (2:2:10)
pi_star=1;
vs=(1:5)';
vh=(0.2:0.2:1)';
sig2=(2:2:10)';

% tot_rows=length(vs)*length(vh)*length(sig2);
% A=[repmat(vs,tot_rows/length(vs),1),repmat(vh,tot_rows/length(vh),1),...
%     repmat(sig2,tot_rows/length(sig2),1)];

vs_temp=vs(1);
vh_temp=vh(1);
sig2_temp=sig2(1);

for j=1:m
    % samples from the likelihood model
    streamflow=hydrograph([vs_temp,vh_temp],...
        GLOBAL_DATA,GEOMORPH);
    randoms=normrnd(0,sig2_temp,[k,length(streamflow)]);         % samples from zero-mean, sig2_temp variance Gaussian
    Y=randoms+repmat(streamflow,k,1);
   
    % marginal of drawn samples
    func=@(x)exp(sum(sum(log(fac*normpdf(Y,repmat(streamflow,k,1),...
        repmat(x,k,length(streamflow)))))));
    x=0.1:0.1:20;
    for l=1:length(x)
        y(l)=func(x(l));
    end
    c(j)=fac^(-m*length(streamflow))*trapz(x,y);
    break
    r(j)=sum(sum(log(normpdf(Y,repmat(streamflow,k,1),...
            repmat(sig2_temp,k,length(streamflow))))))+log(pi_star)-log(c(j));
end
