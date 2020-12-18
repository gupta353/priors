% this routine drawn samples from a parameter space using latin hypercube
% (LH) sampling
% inputs: n = number of samples to be drawn
%         d = number of dimensions
%         par_range  = parameter range
% outputs: param_samps = parameter sets drawn using LH sampling



function param_samps = LHsampAG(n,d,par_range)
    
    % number of dimensions (each dimension corresponds to one parameter)
    lh_samps = lhsdesign(n,d);
      
    param_samps = zeros(n,d);
    for pind = 1:size(lh_samps,2)
        
        param_samps(:,pind) = par_range(pind,1) + lh_samps(:,pind)*(par_range(pind,2)-par_range(pind,1));
        
    end
end