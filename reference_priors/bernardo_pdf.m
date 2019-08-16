% computation of Bernardo pdf value of a given point in parameter space
% using nearest-neighbour method; essentially, any point is assigned a pdf
% equal to that of nearest point at which pdf value is available
% inputs: theta=row vector of parameters at which pdf is to be computed
% outputs: pdf=Bernardo pdf value at a given point
% global variables: GLOBAL_DATA=contains a metrix Bernardo_pdf
% each row of Bernardo_pdf conatins a point in d-dimensional parameter
% space in first d columns and the bernardo pdf of that point in (d+1)th
% column

function pdf=bernardo_pdf(theta)

global GLOBAL_DATA

Bernardo_pdf=GLOBAL_DATA.Bernardo_pdf;
avail_points=Bernardo_pdf(:,1:end-1);               % points at which Bernardo pdf is avaialable
[n,d]=size(avail_points);                           % n=number of points, d=dimensionality

% computation of distance, ds, of theta from each point 
diff=(avail_points-repmat(theta,n,1)).^2;
ds=sqrt(sum(diff,2));

ind= ds==min(ds);
pdf=Bernardo_pdf(ind,d+1);          % Bernardo pdf

% in case of more than one possible nearest-neighbours
if length(pdf)>1
    pdf=randsample(pdf,1);
end

end