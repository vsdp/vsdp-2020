function verifyEqual (~, a, b, varargin)
% VERIFYEQUAL  Verify equality of a and b.
%
%   This function is only intended to by used by GNU Octave (<= 4.2.0), as it
%   is not present there.  Basically it is a wrapper to call assert in a Matlab
%   compatible way.
%
%   Example:
%
%       verifyEqual ([], 1, 2)
%
%   See also demovsdp.

% Copyright 2016-2017 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

if (nargin == 3)
  assert(a, b);
else
  switch (varargin{1})
    case "AbsTol"
      assert(a, b, varargin{2});
    case "RelTol"
      assert(a, b, -varargin{2});
  endswitch
endif
endfunction
