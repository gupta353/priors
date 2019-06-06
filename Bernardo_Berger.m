% Bernardo's prior for a genric process model and zero-mean Gaussian
% error-model

clear all
close all
clc

%% inputs
k=100000;       % number of samples to be drawn
n=100000;       % number of iterations for each parameter

% definition of process model (example: 4.3.1)
psi=16.68;          % (in cm)
delta_theta=0.340;  % change in moisture content
t=3600;   % time at which infiltration is computed (in s)
g=@(x)Green_Ampt_solution(x,psi,delta_theta,t);


% parameters of variance prior
alpha=5;                            
beta=10;
gam=1;

theta_list=1/3600*(0:1:10);       % range of values of hydraulic conductivity (in cm s^-1)
sigma_list=sqrt(1:0.5:10);                 

count=0;
for sigma=1
    for theta=theta_list         

        count=count+1;
        
        for samp=1:n
            Ft=g(theta);            % infiltration at time t
            
            random=normrnd(0,1,[k,1]);
            X=Ft*ones(k,1)+random;
            x_bar=mean(X);
            alpha_d=0.5*k+alpha-2;
            beta_d=0.5*sum((X-x_bar).^2)+beta;
            
            der=Green_Ampt_der(theta,psi,delta_theta,t);
            log_asymp_post(samp)=0.5*log(k)+gammaln(alpha_d-1)...
                +log(abs(der))-beta_d/sigma^2-k/2/sigma^2*...
                (Ft-x_bar)^2-(2*alpha_d+4)*log(sigma)-...
                (alpha_d-1)*log(beta_d)-0.5*log(2*pi);
        end
        
        factor=mean(abs(log_asymp_post))-...            % factor included for numerical stability
            mod(abs(mean(log_asymp_post)),1000);
        prior(count,:)=[theta,sigma,factor,...
            exp(factor+mean(log_asymp_post))];
    end  
end

save_filename='data';
save(save_filename,'prior')