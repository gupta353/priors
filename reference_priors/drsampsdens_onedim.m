% This routine draws samples from a one-dimesnional population given the probability
% density of the distribution at discrete point

clear all
close all
clc

x=[0:0.01:15]';
pdf=chi2pdf(x,4);

n=10000;

% find the areas of the trapezoids created by adjacent points
pdf1=[0;pdf];
T1=(pdf+pdf1(1:end-1))/2;       % first multiplication term in area of trapezoid

x1=[0;x];
T2=(x-x1(1:end-1));             % second multiplication term in the area of the trpezoid

A=T1.*T2;
A(1)=[];

% draw n intervals from the distribution represented by 'A'
cum_A=cumsum(A);                % cumulative probability of different intervals

count=0;
for i=1:n
    u1=unifrnd(0,1);
    ind=find(cum_A<u1,1,'last');
    if ~isempty(ind)
        count=count+1;
        interval(count)=ind;
    end
end

% draw a sample from each interval
slope=(pdf-pdf1(1:end-1))./(x-x1(1:end-1));     % calculation of the slope of density curve in each interval
slope(1)=[];

for i=1:length(interval)
    temp_inter=interval(i);
    u2=unifrnd(0,1);
    
    if slope(temp_inter)~=0
        y(i)=x(temp_inter)+(-pdf(temp_inter)+...                    % random sample from the distribution
        sqrt(pdf(temp_inter)^2+2*slope(temp_inter)*A(temp_inter)*u2))/slope(temp_inter);
    else
        y(i)=x(temp_inter)+A(temp_inter)*u2/pdf(temp_inter);
    end
end

hist(y)