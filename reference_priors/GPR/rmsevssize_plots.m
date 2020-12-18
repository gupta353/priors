% script for plotting rmse vs number of training samples
clear all
close all
clc

% read rmse value for first parameter set
cd('Results\50_20_20')                  % change folder
list=dir('*mat');                       % list all the variables with .mat extension
list1={list.name};                      % access the 'name' field of the structure list in a cell array

for i=1:length(list1)
    load(list1{i})                      % load mat-file into workspace
    RMSE(i)=rmse;                       % store rmse value
    N(i)=length(ytrain);                % store number of training samples
end

X=[N',RMSE'];
X1=sortrows(X);

% read rmse values for second parameter set
cd ..
cd ..
cd('Results\81_50_50')                  % change folder
list=dir('*mat');                       % list all the variables with .mat extension
list1={list.name};                      % access the 'name' field of the structure list in a cell array

for i=1:length(list1)
    load(list1{i})                      % load mat-file into workspace
    RMSE(i)=rmse;                       % store rmse value
    N(i)=length(ytrain);                % store number of training samples
end

X=[N',RMSE'];
X2=sortrows(X);

plot(X1(:,1),X1(:,2),'-o',X2(:,1),X2(:,2),'-o','linewidth',2,'MarkerSize',10,'MarkerEdgeColor','r')
xlabel('Training-set size','fontsize',16)
ylabel('Root-mean-squared error','fontsize',16)
legend({'l=50,\sigma_f^2=20,\sigma^2=20','l=81,\sigma_f^2=50,\sigma^2=50'},'fontsize',16,'fontname','arial')
legend('boxoff')

