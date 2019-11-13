% This routine draws samples from a d-dimesnional population given the probability
% density of the distribution at discrete point
% Nearest neighbour interpolation has been used
% input: x= Equi-spaced points at which pdfs are avaialble (as row vectors)
%           (each row contains a 1xd vector)
%        pdf= pdf value at equi-spaced points
%        n= number of samples to be drawn from the distribution
% output: y=samples drawn from the distribution

function y=drsampsdens_ddim(x,pdf,n)

               
[ns,d]=size(x);         % ns=number of points at which density is available; d=number of dimensions


sep=[]; % minimum separation between points in each dimension (not equal to zero)
bpoints=[]; % boundary points in each dimension, each row contains lower (in first column) and upper (in second column) boundary points in a dimension
for i=1:d
    unique_i=unique(x(:,i));
    min1=min(unique_i);
    unique_i(unique_i==min1)=[];
    min2=min(unique_i);
    max1=max(unique_i);
    sep=[sep,min2-min1];
    bpoints=[bpoints;min1,max1];  
end
base_area=prod(sep);

bpoints=reshape(bpoints,1,numel(bpoints));
diff=repmat(x,1,2)-repmat(bpoints,size(x,1),1);

% variable bpoints changed to contain different infotrmation: first d
% columns conatin '-1's and '0's, where '-1' in (i,j) element indictaes
% that ith sample is on the lower boundary in the jth dimension; similarly,
% second d columns contain '1's and '0's, where '1's in (i,d+j)element  
% indicates that ith sample is on the upper boundary of the jth dimension
bpoints=diff(:,1:d)==0; % replace zeros in in first d-columns of diff by 1 to indicate a lower boundary point 
bpoints=[bpoints,-1*(diff(:,d+1:2*d)==0)]; % replace zeros in in last d-columns of diff by -1 to indicate an upper boundary point

for i=1:size(diff,1)
    num_bpoints(i,1)=length(find(diff(i,:)==0));
end
base_area=repmat(base_area,ns,1)./2.^num_bpoints;
% calculate the areas under the density curve for each segment
A=pdf.*base_area;

% draw n intervals from the distribution represented by 'A'
cum_A=cumsum(A);                % cumulative probability of different intervals

count=0;
for i=1:n
    u1=unifrnd(0,1);
    ind=find(cum_A<u1,1,'last');
    if ~isempty(ind)
        count=count+1;
        segment(count)=ind;
    end
end

% draw sample from each segment
for i=1:length(segment)
    temp_seg=segment(i);
    if num_bpoints(temp_seg)==0
        u2=unifrnd(-0.5,0.5,d,1);
        y(i,:)=x(temp_seg,:)+sep.*u2';
    else
        ind_low=find(bpoints(temp_seg,:)==1);
        ind_up=find(bpoints(temp_seg,:)==-1);
        unif_low_lim=-0.5*ones(d,1);
        unif_low_lim(ind_low)=0;
        unif_up_lim=0.5*ones(d,1);
        unif_up_lim(ind_up-d)=0;
        
        u2=unifrnd(unif_low_lim,unif_up_lim);
        y(i,:)=x(temp_seg,:)+sep.*u2';
    end
end

end