function obj = solve_sdpt3 (obj)
% SOLVE_SDPT3  Approximately solve conic problem instance with SDPT3.
%
%   For information about SDPT3, see:
%
%      http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
%      https://github.com/sqlp/sdpt3
%
%   See also vsdp.solve.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% In case of interval data solve midpoint problem.
A = mid (obj.At);
b = mid (obj.b);
c = mid (obj.c);

% Should initial solution guess be taken into account?
if ((obj.options.USE_STARTING_POINT) && (~isempty (obj.solution)))
  [x0, y0, z0] = deal (obj.solution.x, obj.solution.y, obj.solution.z);
else
  [x0, y0, z0] = deal ([], [], []);
end
% Clear previous solutions or initial points.
obj.solution = [];

if (~isempty (obj.options.SOLVER_OPTIONS))
  OPTIONS = obj.options.SOLVER_OPTIONS;
else
  OPTIONS = [];
end
if (~obj.options.VERBOSE_OUTPUT)
  if (exist ('sqlpmain.m', 'file') == 2) % if SDPT3-4.0
    OPTIONS.printlevel = 0; % default: 3
  else
    OPTIONS.printyes = 0;   % default: 1
  end
end
warning ('off', 'VSDP:svec:justScale');
A = mat2cell (vsdp.svec (obj, A, sqrt(2)), obj.K.dims, obj.m);
warning ('on', 'VSDP:svec:justScale');
c = mat2cell (c, obj.K.dims, 1);
for i = 1:length(obj.K.s)
  c{i} = vsdp.smat ([], c{i}, 1);
end
blk = obj.K.blk;
% SDPT3 terminology is "unconstrained" instead of "free" variables.
if (blk{1,1} == 'f')
  blk{1,1} = 'u';
end

% Call solver.
[~, x, y, z, INFO] = sqlp (blk, A, c, b, OPTIONS, x0, y0, z0);

% Store results.
x = vsdp.svec (obj, vsdp.cell2mat (x(:)), 2);
z = vsdp.svec (obj, vsdp.cell2mat (z(:)), 1);
f_objective = [obj.c'*x; obj.b'*y];
if (isstruct (INFO))
  info = INFO.termcode;  % SDPT3-4.0 output
else
  info = INFO(1);  % SDPT3-3.x output
end
obj.add_solution(x, y, z, f_objective, obj.options.SOLVER, info)

end
