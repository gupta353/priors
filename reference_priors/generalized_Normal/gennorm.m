% This routine geenrates random samples from generalized Gaussian distribution
% Inputs: n=number of samples to be drawn
%         mu=mean
%         phi=scale parameter
%         beta=shape parameter
% output: samps=n samples from generalized-Gaussian
% reference: Gomez et al. (1998), A multivariate generalization of the
% power exponential family of distributions

function samps=gennorm(n,mu,phi,beta)

% draw samples from gamma distribtuion
x=gamrnd(1/2/beta,1,[n,1]);
r=x.^(1/2/beta);

% draw +1 and -1s with 50-50 chance
unif=unifrnd(0,1,[n,1]);
sphe_uni=ones(n,1);
sphe_uni(unif<=0.5)=-1;

% draw samples from geeralized-Gaussian
samps=mu+phi*r.*sphe_uni;
end