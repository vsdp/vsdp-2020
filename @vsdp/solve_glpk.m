function obj = solve_glpk (obj, sol_type)
% SOLVE_GLPK  Approximately solve conic problem instance with GLPK.
%
%   For more information about GLPK, see:
%
%     https://www.gnu.org/software/glpk/
%
%   See also vsdp.solve.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (1, 2);

% Check solver availability.
if (exist ('glpk', 'file') ~= 2)
  error ('VSDP:solve_glpk:notAvailable', ...
    ['solve_glpk: GLPK does not seem to be ready.\n\n', ...
    'To select another solver, run:  %s.solve()'], inputname(1));
end

% Check cones.
if ((sum (obj.K.q) > 0) || (sum (obj.K.s) > 0))
  error ('VSDP:solve_glpk:unsupportedCones', ...
    ['solve_glpk: Second-order cones (K.q) and semidefinite cones (K.s) ', ...
    'are not supported by GLPK.']);
end

if (nargin == 1)
  sol_type = 'Approximate';
end
[A, b, c] = obj.get_midpoint_problem_data (sol_type);

% Should initial solution guess be taken into account?
if (obj.options.USE_INITIAL_GUESS)
  warning ('VSDP:solve_glpk:ignoreInitialGuess', ...
    ['solve_glpk: GLPK does not support initial guesses (x0,y0,z0) ' ...
    'and proceeds without them.']);
end

% Should special solver options be taken into account?
if (~isempty (obj.options.SOLVER_OPTIONS))
  param = obj.options.SOLVER_OPTIONS;
else
  param = [];
end

% Adapt output verbosity.
if (~obj.options.VERBOSE_OUTPUT)
  param.msglev = 0;
end

% Prepare data for solver.
[A, b, c] = deal (A', full (b), full (c));

sense = 1;  % Minimization
lbound = [ ...
  -inf( obj.K.f, 1); ...
  zeros(obj.K.l, 1)];                 % lower bound
ubound = inf (length (c), 1);         % upper bound
vtype = repmat ('C', length (c), 1);  % variable   types: continuous
ctype = repmat ('S', length (b), 1);  % constraint types: equality

% Call solver.
tic;
[x, ~, errnum, extra] = glpk ...
  (c, A, b, lbound, ubound, ctype, vtype, sense, param);
solver_info.elapsed_time = toc;

% Store solution.
if (isfield (extra, 'lambda'))
  y = extra.lambda;
  z = obj.c - obj.At*y;
end
f_objective = [obj.c'*x; obj.b'*y];
solver_info.name = 'glpk';
switch (errnum)
  case 0
    solver_info.termination = 'Normal termination';
  case 10
    solver_info.termination = 'Primal infeasible';
  case 11
    solver_info.termination = 'Dual infeasible';
  case 15
    solver_info.termination = 'Primal and dual infeasible';
  otherwise
    solver_info.termination = 'Unknown';
end

obj.add_solution (sol_type, x, y, z, f_objective, solver_info);

end
