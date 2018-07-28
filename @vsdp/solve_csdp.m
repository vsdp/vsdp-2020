function obj = solve_csdp (obj)
% SOLVE_CSDP  Approximately solve conic problem instance with CSDP.
%
%   For more information about CSDP, see:
%
%     https://projects.coin-or.org/Csdp/
%     https://github.com/coin-or/Csdp
%
%   See also vsdp.solve.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% Check solver availability.
if (exist ('csdp', 'file') ~= 2)
  error ('VSDP:solve_csdp:notAvailable', ...
    ['solve_csdp: CSDP does not seem to be ready.\n\n', ...
    'To select another solver, run:  %s.solve()'], inputname(1));
end

% Check cones.
if (sum (obj.K.q) > 0)
  error ('VSDP:solve_csdp:unsupportedCone', ...
    'solve_csdp: Second-order cones (K.q) are not supported by CSDP.');
end

% Should initial solution guess be taken into account?
if ((obj.options.USE_STARTING_POINT) && (~isempty (obj.solution)))
  [x0, y0, z0] = deal (obj.solution.x, obj.solution.y, obj.solution.z);
else
  [x0, y0, z0] = deal ([], [], []);
end

% Should special solver options be taken into account?
if (~isempty (obj.options.SOLVER_OPTIONS))
  pars = obj.options.SOLVER_OPTIONS;
else
  pars = [];
end

% Adapt output verbosity.
if (~obj.options.VERBOSE_OUTPUT)
  pars.printlevel = 0;
end

% Prepare data for solver.
[A, b, c] = obj.get_perturbed_midpoint_problem ();
[b, c, x0, y0, z0] = deal (full (b), full (c), full (x0), full (y0), full (z0));

% Convert to SeDuMi-Format (same as CSDP format).
A = vsdp.smat (obj, A, 1);
c = vsdp.smat (obj, c, 1);
K = obj.K;
if (K.f > 0)
  warning('VSDP:solve_csdp:unsupportedCone', ...
    ['solve_csdp: CSDP supports free variables (K.f) by converting them ' ...
    'to the difference of positive variables.  The resulting problem is ', ...
    'ill-posed.']);
  [A, b, c, K] = convertf (A, b, c, K);  % CSDP function.
end

args = {A, b, c, K, pars};
% CSDP requires all initial solutions to be not empty!
if (~isempty (x0) && ~isempty(y0) && ~isempty (z0))
  args = [args, {x0, y0, z0}];
end

% If using Windows, switch to the 'csdp.exe' directory.
if (ispc ())
  p = fileparts (which ('csdp.exe'));
  p = cd (p);
end

% Call solver.
tic;
[x, y, z, INFO] = csdp (args{:});
elapsed_time = toc;

if (ispc ())
  cd (p);
end

% Store solution.
x = vsdp.svec (obj, x, 2);
z = vsdp.svec (obj, z, 1);
f_objective = [obj.c'*x; obj.b'*y];
if (any (INFO == [0, 1, 2]))
  info = INFO;
end
obj.add_solution(x, y, z, f_objective, obj.options.SOLVER, info, elapsed_time);

end
