% A function that computes pdf values at given sample values with a given
% mean and different values of sigma

function p=integral_normpdf(x,mu,sigs)
x
mu
sigs
    for i=1:length(sigs)
        sigs(i)
        p(i)=prod(normpdf(x,mu,sigs(i)*ones(length(x),1)));
    end
end