% Probabilistic Distributed Model (PDM);
% Equal-moisture depth is assumed in each (conceptual) storage unit;
% Distribution function of storage capacity:
% F(c)=1-(1-frac_{c}{c_max})^b,
% where, c_max=maximum storage capacity
% b=parameter describing spatial variation
% inputs: P=hourly total rainfall (in mm)
%         Ea=hourly actual evaporation (in mm)
%         S0=intital soil moisture depth (in mm)
%         Cmax=maximum storage capacity (in mm)
%         b=spatial variation parameter of storage distribution
%         Kb=constant of baseflow reservoir (in s^-1)
% output: V=excess rainfall during each time-interval (in mm)
%         gam=baseflow rate (in mm s^-1)
%         S=soil-moisture storage (in mm)
% Note: this rourtine updates input hourly evaporation series such that
% evaporation during each hour is proportional to soil stoarge at the end
% of previous hour
% Ref: Moore R.J. (1985), The probabilistic-distributed principle and
% runoff production at point and basin scales.


function [V,gam,S,C,Ea]=pdm(P,Ea,S0,Cmax,b,Kb)
    
    Ea_err_tol=0.01;    % error tolerance for evaporation update
    Ea_num_updates=10;  % maximum number of updates
    
    n=length(P);        % length of input time-series
    delta_t=3600;       % length of each time-step (in seconds)
    Smax=Cmax/(b+1);    % total storage-depth
    %% error-checks
    if length(P)~=length(Ea)
        error('length of rainfall and evaporation series must be the same');
    elseif S0>Smax
        error('Intital depth cannot be greater than the total-stoarage depth');
    end
    
    %%
    % Ea=zeros(n,1);      % Actual evaporation
    net_P=zeros(n,1);   % net rainfall rate
    V=zeros(n,1);       % excess rainfall (direct-runoff)
    C=zeros(n+1,1);     % critial capacity
    S=zeros(n+1,1);     % storage depth
    gam=zeros(n,1);     % baseflow
    
    C(1)=Cmax*(1-((1-S0/Smax)^(1/(b+1))));          % initial critial depth
    S(1)=S0;                                        % moisture storage before the first loop
    
    hour_count=0;
    update_count=0;
    step=1;
    vblock = [];
    while step<n+1
        
        step=step+1;
        hour_count=hour_count+1;
        ind=step-1;
        
        gam(ind)=Kb*S(ind);                        % drainage (mm/s) due to storage at previous time-step
        net_P(ind)=(P(ind)-Ea(ind))...
            /delta_t-gam(ind);                     % net rainfall rate (mm/s)
        C(step)=min(Cmax,...
            C(ind)+net_P(ind)*delta_t);            % critial depth (mm) at the end of this time-period
        C(step)=max(0,C(step));
        
        % Excess rainfall (direct runoff) and storage-depth computation
        if net_P(ind)>0 % if net rainfall rate is positive
            
            if C(step)<Cmax % all stores are not full at the end of this interval
                V(ind)=net_P(ind)*delta_t-Smax*...
                    ((1-C(ind)/Cmax)^(b+1)-...
                    (1-C(step)/Cmax)^(b+1));
                S(step)=S(ind)+net_P(ind)*delta_t-V(ind);  % soil moisture depth at the end of this interval
                
                norm_str = S(step)/Smax;
                
                if abs(norm_str-1)<10^-15
                    C(step)=Cmax; % computation of critical storge depth
                else
                    C(step)=Cmax*(1-((1-norm_str)^(1/(b+1)))); % computation of critical storge depth
                end
                
            elseif C(step)>=Cmax && C(step-1)<=Cmax % all stores are full at the end of this interval
                V(ind)=net_P(ind)*delta_t-Smax*...
                    (1-C(ind)/Cmax)^(b+1);
                S(step)=Smax;
                C(step)=Cmax;
                
            else % all stores were full before the begining of this interval
                V(ind)=net_P(ind)*delta_t;
                S(step)=Smax;
                C(step)=Cmax;
                
            end
            
        else  % if net rainfall rate is negative
            
            V(ind)=0;
            S(step)=S(ind)+net_P(ind)*delta_t;
            S(step)=max(0,S(step));
            
            % computation of critical storge depth (C)
            norm_str = S(step)/Smax;
            if abs(norm_str-1)<10^-15
                C(step)=Cmax; 
            else
                C(step)=Cmax*(1-((1-norm_str)^(1/(b+1)))); 
            end
            
        end
        
        % update hourly evaporation data
        if hour_count==24
            
            Ea_update=evaporation_redistribution(...
                Ea(ind-hour_count+1:ind),S(ind-hour_count+1:ind));
            err=abs(Ea(ind-hour_count+1:ind)-Ea_update);
            
            % if Ea_updated is different from Ea
            if max(err)>Ea_err_tol && update_count<Ea_num_updates
                update_count=update_count+1;
                Ea(ind-hour_count+1:ind)=Ea_update;
                step=step-hour_count;
                
            elseif  max(err)>Ea_err_tol && update_count>=Ea_num_updates
                update_count=0;
            end
            
            hour_count=0;
        end
        
    end
end

function Ea_updated=evaporation_redistribution(Ea,S) % distrubution of daily evaporation within a day evaporation
    
    % compute the ratio of total daily evaporation, and total stoarge at the
    % beginning of each time-interval
    ratio=sum(Ea)./sum(S);
    Ea_updated=S*ratio;
end
