function res = sqrtsup(a)
%% SQRTSUP:  calculates a verified upper bound for the square root
%   -> setround(1) is assumed
% only used in case the sqrt function does not regard the rounding mode
res = sqrt(a);
res = res + realmin * (res.*(-res) > -a);
