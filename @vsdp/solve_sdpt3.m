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
if (obj.options.USE_STARTING_POINT)
  [x0, y0, z0] = deal (obj.x, obj.y, obj.z);
else
  [x0, y0, z0] = deal ([], [], []);
end
% Clear previous solutions or initial points.
[obj.x, obj.y, obj.z] = deal ([], [], []);

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

% Call solver.
[~, x, y, z, INFO] = sqlp (obj.K.blk, A, c, b, OPTIONS, x0, y0, z0);

% Store results.
obj.x = vsdp.svec (obj, vsdp.cell2mat (x(:)), 2);
obj.y = y;
obj.z = vsdp.svec (obj, vsdp.cell2mat (z(:)), 1);
if (isstruct (INFO))
  info = INFO.termcode;  % SDPT3-4.0 output
else
  info = INFO(1);  % SDPT3-3.x output
end

end
