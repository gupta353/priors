% % This routine draws samples from a d-dimesnional population given the probability
% density of the distribution at discrete point
% Nearest neighbour interpolation has been used

clear all
close all
clc

x1=-4:0.1:4;
x2=-4:0.1:4;
count=0;
for i=1:length(x1)
    for j=1:length(x2)
        count=count+1;
        x(count,:)=[x1(i),x2(j)];
    end
end
pdf=mvnpdf(x,[0,0],[1,0;0,1]);

n=10000;                  % number of samples to be drawn
[ns,d]=size(x);         % ns=number of points at which density is available; d=number of dimensions

sep=[0.1,0.1];
base_area=prod(sep);

bpoints=[-4,4;-4,4];
bpoints=reshape(bpoints,1,numel(bpoints));
diff=repmat(x,1,2)-repmat(bpoints,size(x,1),1);
% replace zeros in in first d-columns of diff by 1 to indicate a lower boundary point 
bpoints=diff(:,1:d)==0;
% replace zeros in in last d-columns of diff by -1 to indicate an upper boundary point
bpoints=[bpoints,-1*(diff(:,d+1:2*d)==0)];

for i=1:size(diff,1)
    num_bpoints(i,1)=length(find(diff(i,:)==0));
end
base_area=repmat(base_area,ns,1)./2.^num_bpoints;
A=pdf.*base_area;

% calculate areas of the intervals
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

hist(y)