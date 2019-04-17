function obj = filter (obj, varargin)
% FILTER  Index subset of the VSDP benchmark object 'obj'.
%
%   The whole set of benchmarks is in the variable 'obj.BENCHMARK'.
%
%   obj.filter ()  Leaves this variable untouched.  This is equivalent to:
%
%       obj.filter ([], [])
%       obj.filter ('.*', '.*')
%
%   obj.filter (bm_library, name)  Optionally, the benchmark can be reduced to
%     a subset of the data by applying filters, i.e. regular expressions
%     machted with the "regexp()" function.  To filter the benchmark library
%     provide a string for 'bm_library' and for the test case name set 'name'.
%
%   Example:
%
%     obj.filter ('SDPL.*', 'arch.*');
%
%   See also vsdp_benchmark.
%

% Copyright 2004-2019 Christian Jansson (jansson@tuhh.de)

idx = true (size (obj.BENCHMARK(:,1)));
if (isempty (idx))
  return,
end

% Filter benchmark librariess.
if ((nargin > 1) && ~isempty (varargin{1}))
  bm_library  = varargin{1};
  idx = idx & cellfun (@(x) ~isempty (regexp (x, bm_library, 'once')), ...
    obj.BENCHMARK(:,1));
else
  bm_library = '.*';
end

% Filter benchmark names.
if ((nargin > 2) && ~isempty (varargin{2}))
  name = varargin{2};
  idx = idx & cellfun (@(x) ~isempty (regexp (x, name, 'once')), ...
    obj.BENCHMARK(:,2));
else
  name = '.*';
end
if (~any (idx))
  error ('VSDP_BENCHMARK:filter:noMatch', ...
    ['filter: bm_library = ''%s'' and name = ''%s'' do not match ', ...
    'any benchmark test case.'], bm_library, name);
else
  obj.BENCHMARK = obj.BENCHMARK(idx,:);
end

end
