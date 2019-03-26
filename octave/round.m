function y = round (x, n)
% ROUND  More Matlab compatible rounding function.
%
%   Example:
%
%       y = round (pi, 3)  % = 3.1420
%
%   See also round.
%

% Copyright 2016-2019 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

if (nargin == 1)
  n = 0;
endif
y = builtin ("round", x * 10^n) / 10^n;

endfunction
