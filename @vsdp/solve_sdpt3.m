function obj = solve_sdpt3 (obj, sol_type)
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

narginchk (1, 2);

% Check solver availability.
if (exist ('sqlp', 'file') ~= 2)
  error ('VSDP:solve_sdpt3:notAvailable', ...
    ['solve_sdpt3: SDPT3 does not seem to be ready.  ', ...
    'Did you run ''install_sdpt3()'' inside the solver directory?\n\n', ...
    'To select another solver, run:  %s.solve()'], inputname(1));
end

if (nargin == 1)
  sol_type = 'Approximate';
  [A, b, c] = deal (obj.At, obj.b, obj.c);
else
  [A, b, c] = obj.get_perturbed_midpoint_problem ();
end

% Should initial solution guess be taken into account?
if ((obj.options.USE_INITIAL_GUESS) && (~isempty (obj.solution('Initial'))))
  isol = obj.solution('Initial');
  [x0, y0, z0] = deal (isol.x, isol.y, isol.z);
else
  [x0, y0, z0] = deal ([], [], []);
end

% Should special solver options be taken into account?
if (~isempty (obj.options.SOLVER_OPTIONS))
  OPTIONS = obj.options.SOLVER_OPTIONS;
else
  OPTIONS = [];
end

% Adapt output verbosity.
if (~obj.options.VERBOSE_OUTPUT)
  if (exist ('sqlpmain.m', 'file') == 2) % if SDPT3-4.0
    OPTIONS.printlevel = 0; % default: 3
  else
    OPTIONS.printyes = 0;   % default: 1
  end
end

% Prepare data for solver.
blk = obj.K.blk;
warning ('off', 'VSDP:svec:justScale');
A = mat2cell (vsdp.svec (obj, A, sqrt(2)), obj.K.dims, obj.m);
warning ('on', 'VSDP:svec:justScale');
c = mat2cell (c, obj.K.dims, 1);
for i = find ([blk{:,1}] == 's', 1):size(blk, 1)
  c{i} = vsdp.smat ([], c{i}, 1);
end
% In SDPT3 terminology "unconstrained" are "free" variables.
if (blk{1,1} == 'f')
  blk{1,1} = 'u';
end

% Call solver.
tic;
[~, x, y, z, info] = sqlp (blk, A, c, b, OPTIONS, x0, y0, z0);
solver_info.elapsed_time = toc;

% Store solution.
x = vsdp.svec (obj, vsdp.cell2mat (x(:)), 2);
z = vsdp.svec (obj, vsdp.cell2mat (z(:)), 1);
f_objective = [obj.c'*x; obj.b'*y];
solver_info.name = 'sdpt3';
if (isstruct (info))
  termcode = info.termcode;  % SDPT3-4.0 output
else
  termcode = info(1);  % SDPT3-3.x output
end
switch (termcode)
  case 0
    solver_info.termination = 'Normal termination';
  case 1
    solver_info.termination = 'Primal infeasible';
  case 2
    solver_info.termination = 'Dual infeasible';
  case 3
    solver_info.termination = 'Primal and dual infeasibile';
  otherwise
    solver_info.termination = 'Unknown';
end

obj.add_solution (sol_type, x, y, z, f_objective, solver_info);

end
