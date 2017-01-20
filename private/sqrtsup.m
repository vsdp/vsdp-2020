function res = sqrtsup(a)
% SQRTSUP  Calculates a verified upper bound for the square root.
%
%    setround(1) is assumed
%    Only used in case the sqrt function does not regard the rounding mode.

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

res = sqrt(a);
res = res + realmin * (res.*(-res) > -a);
end
