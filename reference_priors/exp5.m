% reference priors for a hypothetical process model-parameters
% model: y=exp(-1*\theta)+\epsilon, \epsilon~N(0,\sigma^2)
% 

clear all
close all
clc

model=@(theta)exp(-theta);


% computation of reference priors (Berger et al., 2010)
%{
k=100;
m=1000;
sig2=0.4;
input=1;

pi_star=1;

theta=0:0.1:1;

for i=1:length(theta)
    for j=1:m
        theta_tmp=theta(i);
        samps=normrnd(model(theta_tmp),sig2,[k,1]);
        func=@(x)prod(normpdf(samps,model(x),sig2));
%         c(j)=integral(func,1,10);
        % integration using trapz method
        x=-10:0.1:10;
        for l=1:length(x)
            y(l)=func(x(l));
        end
        c(j)=trapz(x,y);
        expr=prod(normpdf(samps,model(theta_tmp),sig2))*pi_star/c(j);
        r(j)=log(expr);
    end
    PI(i)=exp(sum(r)/m);
end
%}
%% Jeffery's prior for \epsilon~GN(0,\sigma^2,\beta) (ref: Gomez et al., 1998)
% elements of Fisher-Information matrix
theta_list=0:0.5:5;
s_list=0.0001:0.5:5;
beta_list=0.51:0.24:1;

count=0;
for theta_iter=1:length(theta_list)
    for s_iter=1:length(s_list)
        for beta_iter=1:length(beta_list)
            
            count=count+1;
            
            theta=theta_list(theta_iter);
            s=s_list(s_iter);
            beta=beta_list(beta_iter);
            
            k=1+0.5/beta;
            
            if k~=2
                i11=-(2*(2*beta-1)/s^2)*(gamma(2-k)/gamma(k))*exp(-2*theta);
            else
                i11=0;
            end
            i12=0;
            
            i13=0;
            
            i22=-2*beta/s^2;
            
            i23=(1/beta/s)*(1+log(2)+psi(k));
            
            i33=-(1/beta^3)*(log(2)+0.5*log(2)^2)-(1/beta^3)*(1+log(2))*psi(k)...
                -(0.5/beta^3)*(0.5/beta+1)*(psi(1,k)+psi(k)^2)+(0.25/beta^4)*psi(k)^2;
            
            % detrminant computation
            matrix{count}=-1*[i11,i12,i13;i12,i22,i23;i13,i23,i33];
            
            determ=det(matrix{count});
            

            
            check_data(count)=determ;
            
            pi_shi(count,:)=[theta,s,beta,determ^0.5];
        end
    end
end
%% computation of posterior for a given observation
%{
global DREAM_dir EXAMPLE_dir CONV_dir yobs
DREAM_dir=['D:/Research/Thesis_work/Structural_'...
    'uncertainty/MatLab_codes/20180222/Dream'];
CONV_dir=['D:/Research/Thesis_work/Structural_'...
    'uncertainty/MatLab_codes/20180222/Dream/diagnostics'];
EXAMPLE_dir = ['D:/Research/Thesis_work/Structural_'...
    'uncertainty/MatLab_codes/20180222/reference_priors'];
addpath(DREAM_dir)
addpath(CONV_dir)
addpath(EXAMPLE_dir)


yobs=0.5;                                   % hypothetical observed value

Func_name='exp5_likeli';                    % target pdf
DREAMPar.d=1;                               % dimensionality of the problem
DREAMPar.N=8;                               % Number of chains
DREAMPar.T=10000;                           % Number of generations
DREAMPar.lik=1;                             % Type of likelihood (2 indicates log-likelihood)

Par_info.initial='prior';                 % type of prior 
Par_info.prior={'exppdf(1)'};
Par_info.min=0;                             % lower bound on parameters                   
Par_info.max=1;                             % upper bound
Par_info.boundhandling='fold';              % type of boundhandling

[chain,output,fx]=DREAM(Func_name,...       % main function to Run Dream
    DREAMPar,Par_info);
%}