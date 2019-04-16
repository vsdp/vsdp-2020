function f = filter (obj, varargin)
% FILTER  Index subset of the VSDP benchmark object 'obj'.
%
%   The whole set of benchmarks is in the variable 'obj.BENCHMARK' and the
%   whole set of solvers in the variable 'obj.SOLVER(3:end).name'.
%
%   f = obj.filter ()  Return indices for all benchmarks and solvers in
%     'f.benchmark' and 'f.solver', respectively.  This is equivalent to:
%
%       obj.filter ([], [], [])
%       obj.filter ('.*', '.*', '.*')
%
%   f = obj.filter (bm_library, name, solver)  Optionally, the returned solver
%     and benchmark indices can be reduced to a subset of the data by applying
%     filters, i.e. regular expressions machted with the "regexp()" function.
%     To filter the benchmark library provide a string for 'bm_library', for
%     the test case name set 'name', and for the solver set 'solver'.
%
%   Example:
%
%     f = obj.filter ('SDPL.*', 'arch.*');
%     obj.BENCHMARK(f.benchmark).name
%
%   See also vsdp_benchmark.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

f.benchmark = 1:length(obj.BENCHMARK);
f.solver    = 3:length(obj.SOLVER);  % Ignore VSDP and INTLAB.

% Filter benchmark librariess.
if ((nargin > 1) && ~isempty (varargin{1}))
  bm_library  = varargin{1};
  f.benchmark = find (cellfun ( ...
    @(x) ~isempty (regexp (x, bm_library, 'once')), {obj.BENCHMARK.lib}));
else
  bm_library = '.*';
end

% Filter benchmark names.
if ((nargin > 2) && ~isempty (varargin{2}))
  name = varargin{2};
  idx = cellfun (@(x) ~isempty (regexp (x, name, 'once')), ...
    {obj.BENCHMARK(f.benchmark).name});
  f.benchmark = f.benchmark(idx);
else
  name = '.*';
end
if (isempty (f.benchmark))
  error ('VSDP_BENCHMARK:filter:noMatch', ...
    ['filter: bm_library = ''%s'' and name = ''%s'' do not match ', ...
    'any benchmark test case.'], bm_library, name);
end

% Filter solver.
if ((nargin > 3) && ~isempty (varargin{3}))
  solver = varargin{3};
  idx = cellfun (@(x) ~isempty (regexp (x, solver, 'once')), ...
    {obj.SOLVER(f.solver).name});
  f.solver = f.solver(idx);
  if (isempty (f.solver))
    error ('VSDP_BENCHMARK:filter:noMatch', ...
      'filter: solver = ''%s'' does not match any solver.', solver);
  end
end

end
