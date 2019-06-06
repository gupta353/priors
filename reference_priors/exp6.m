% illustration of effect of prior on a simple modeling process
% modeling process:
% y=exp(-\theta*t)+\epsilon
% where,  t = time-step
%         \theta=parameter
%         \epsilon~N(0,\sigma^2)
clear all
close all
clc

%% generate the artifical dataset with Gaussian noise
t=[1]';
sig2=0.4;
theta_true=0;
% yobs=exp(-theta_true*t)+ normrnd(0,1,[length(t),1]);
yobs=0.1;
%% set the prior
likeli=@(theta)prod(normpdf(yobs,exp(-theta*t),sig2));
theta=0:0.1:1;
for i=1:length(theta)
    post(i)=(1/3)*likeli(theta(i));
end

plot(theta,post)