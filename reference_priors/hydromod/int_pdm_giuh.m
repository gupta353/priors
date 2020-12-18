% computation of surface- and base-flow, for a given rainfall data and
% watershed geomorphology
% inputs: GLOBAL_DATA=A structure array (see run_pdm_giuh.m for details)
%         GEOMORPH=A structure array (see mainScript.m and run_pdm_giuh.m for details)
%         param= a vector of model-parameters
%                param(1)=Cmax: maximum storage capacity
%                param(2)=b: spatial variation parameter of storage
%                capacity distribution
%                param(3)=baseflow reservoir constant
%                param(4)=in-stream velocity of water-drops
%                param(5)=hillslope velocity of water-drops
% outputs: streamflow=streamflow (in m^3 s^-1)
%          surfaceflow=surface runoff (in m^3 s^-1)
%          baseflow=baseflow (in m^3 s^-1)
%          S=storage depth of water in soil column (in mm)
function [streamflow,surfaceflow,baseflow,S,V,C,gam,Ea_updated]=int_pdm_giuh(param,GLOBAL_DATA,GEOMORPH)


R=GLOBAL_DATA.R;         % total rainfall time-series
Ea=GLOBAL_DATA.Ea;       % actual evaporation time-series
S0=GLOBAL_DATA.S0;       % initial soil-moisture depth
darea=GLOBAL_DATA.darea; % drainage area in hectares


Cmax=param(1);
b=param(2);
Kb=param(3);
theta(1)=param(4);
theta(2)=param(5);


% determination of intital soil mositure depth (boundary condition)
R_warmup=GLOBAL_DATA.R_warmup;         % rainfall time-series for warm-up period
Ea_warmup=GLOBAL_DATA.Ea_warmup;       % evaporation time-series for warm-up period
[~,~,S]=pdm(R_warmup,Ea_warmup,S0,Cmax,b,Kb);
S0=S(end);

% computation of excess rainfall
[V,gam,S,C,Ea_updated]=pdm(R,Ea,S0,Cmax,b,Kb);
baseflow=gam*darea*10;

% conversion of base-flow from hourly to daily scale
baseflow=sum(reshape(baseflow,[24,length(baseflow)/24]))/24;

% computation of streamflow
[surfaceflow]=hydrograph(theta,V,GLOBAL_DATA,GEOMORPH);
minlength=min(length(surfaceflow),length(baseflow));
streamflow=surfaceflow(1:minlength)+baseflow(1:minlength);
% end


end