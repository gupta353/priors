% Kronecker delta function
%%% inputs:
% m=a scalar
% n=a scalar
%%% outputs
% 1 or 0

function f=kroneckerDelta(m,n)
if m==n
    f=1;
else
    f=0;
end
end