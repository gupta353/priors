function Ea_updated=evaporation_redistribution(Ea,S) % distrubution of daily evaporation within a day evaporation 

% compute the ratio of total daily evaporation, and total stoarge at the
% beginning of each time-interval
ratio=sum(Ea)./sum(S);
Ea_updated=S*ratio;
end