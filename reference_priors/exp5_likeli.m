% likelihood function for exp5

function p=exp5_likeli(theta)

global yobs

sig2=0.2;
mu=exp(-theta);
p=normpdf(yobs,mu,sig2);

end