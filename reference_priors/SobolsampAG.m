% this routine drawn samples from a parameter space using Sobolev sequences
% inputs: n = number of samples to be drawn
%         d = number of dimensions
%         par_range  = parameter range
%         ni = interval from which samples are drawn
% outputs: param_samps = parameter sets drawn using sobolev sequences



function param_samps = SobolsampAG(n,d,par_range,ni)
    
    % number of dimensions (each dimension corresponds to one parameter)
    sobol_samps = 	sobolset(d,'Skip',1000);
    sobol_samps =   sobol_samps(ni,:);
    
    param_samps = zeros(n,d);
    for pind = 1:size(sobol_samps,2)
        
        param_samps(:,pind) = par_range(pind,1) + sobol_samps(:,pind)*(par_range(pind,2)-par_range(pind,1));
        
    end
end