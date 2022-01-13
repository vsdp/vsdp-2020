function verifyEqual (~, a, b, varargin)
% VERIFYEQUAL  Verify equality of a and b.
%
%   This MATLAB function is not present in GNU Octave (<= 5.1.0).
%
%   Example:
%
%       verifyEqual ([], 1, 2)
%
%   See also testVSDP, testVSDP2006.
%

% Copyright 2016-2020 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

% Octave's assert has problems with "intval" class.
if (isa (a, "intval") || isa (b, "intval"))
  if (~ isequal (a, b))
    error ("verifyEqual: intval test for a == b failed.");
  endif
  return;
endif

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
