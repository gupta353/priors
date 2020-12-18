%%% Estimation of streamflow by convolution of rainfall and instantanous unit hydrograph
% inputs: theta;[in-stream velocity, hillslope velocity] 
%         V=excess rainfall (in mm)
%         GLOBAL_DATA=a structure array of user provided data (look at run_pdm_giuh.m for more details)
%         GEOMORPH=a structure array of geomorphology of the watershed (look at
%        the mainScript.m documentation for details)
% output: streamflow=average streamflow in m^3/s at scale provided in
%                    GLOBAL_data
%         uh=instantanous unit hydrograph at point of
%            observation/simulation (at every second)

function [streamflow,uh]=hydrograph(theta,V,GLOBAL_DATA,GEOMORPH)

R=V;
darea=GLOBAL_DATA.darea;
node=GLOBAL_DATA.nodeobs;
support=GLOBAL_DATA.support;
R_time_scale=GLOBAL_DATA.R_time_scale;
S_time_scale=GLOBAL_DATA.S_time_scale;

vel_stream=theta(1);
vel_hillslope=theta(2);
uh=unitHydrograph(node,vel_stream,vel_hillslope,support,GEOMORPH,GLOBAL_DATA);

% if strcmpi(time_scale,'hourly')
%     Tr=3600; % for hourly rainfall data 
% elseif strcmpi(time_scale,'daily')
%     Tr=24*3600; % for daily rainfall data
% else
%     error('time_scale should be either ''hourly'' or ''daily''');
% end
Tr=60*min(R_time_scale,S_time_scale);
% convsersion of uh-scale to minimum of rainfall time-scale and stteamflow scale
l_uh=length(uh);
addit=Tr-Tr*(l_uh/Tr-floor(l_uh/Tr));
uh=[uh;zeros(ceil(addit),1)];
uh=reshape(uh,Tr,length(uh)/Tr);
uh=sum(uh)/size(uh,1);
uh=[0;uh'];

% re-arranging rainfall data
% for i=1:length(R)
%     Rs(1+Tr*(i-1):i*Tr)=R(i)/Tr;
% end

streamflow_s=10*darea*conv(uh,R); % 1-mm x 1-hectare=10-m^3

%convsersion of streamflow from hourly to daily scale
%
if S_time_scale==60*24
    count=0;
    Ts=24;
    slength=floor(length(streamflow_s)/Ts); % total number of days in streamflow time-series
    streamflow=nan(1,slength);
    for ii=1:Ts:length(streamflow_s)
        count=count+1;
        ind=ii:ii+Ts-1;
        if ind(end)<length(streamflow_s)
            streamflow(count)=sum(streamflow_s(ind))/length(ind);
        end
    end
else
    streamflow=streamflow_s';
end
%}
end