function obj = solve_linprog (obj, sol_type)
% SOLVE_LINPROG  Approximately solve conic problem instance with LINPROG.
%
%   For more information about LINPROG, see:
%
%     https://www.mathworks.com/help/optim/ug/linprog.html
%
%   See also vsdp.solve.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (1, 2);

% Check solver availability.
if (exist ('linprog', 'file') ~= 2)
  error ('VSDP:solve_linprog:notAvailable', ...
    ['solve_linprog: LINPROG does not seem to be ready.\n\n', ...
    'To select another solver, run:  %s.solve()'], inputname(1));
end

% Check cones.
if ((sum (obj.K.q) > 0) || (sum (obj.K.s) > 0))
  error ('VSDP:solve:unsupportedCones', ...
    ['solve_linprog: Second-order cones (K.q) and semidefinite cones ', ...
    '(K.s) are not supported by LINPROG.']);
end

if (nargin == 1)
  sol_type = 'Approximate solution';
  [A, b, c] = deal (obj.At, obj.b, obj.c);
else
  [A, b, c] = obj.get_perturbed_midpoint_problem ();
end

% Should initial solution guess be taken into account?
if ((obj.options.USE_INITIAL_GUESS) ...
    && (~isempty (obj.solution('Initial guess'))))
  x0 = full (obj.solution('Initial guess').x);
else
  x0 = [];
end

% Should special solver options be taken into account?
options = optimoptions ('linprog', 'Algorithm', 'interior-point-legacy');
if (~isempty (obj.options.SOLVER_OPTIONS))
  options = optimoptions (options, obj.options.SOLVER_OPTIONS);
end

% Adapt output verbosity.
if (~obj.options.VERBOSE_OUTPUT)
  options = optimoptions (options, 'Display', 'off');
end

% Prepare data for solver.
[A, b, c] = obj.get_perturbed_midpoint_problem ();
[A, b, c] = deal (A', full (b), full (c));
lbound = [ ...
  -inf(obj.K.f, 1); ...
  zeros(obj.K.l, 1)];       % lower bound
ubound = inf(length(c),1);  % upper bound

% Call solver.
tic;
[x, ~, flag, ~, lambda] = linprog ...
  (c, [], [], A, b, lbound, ubound, x0, options);
solver_info.elapsed_time = toc;

% Store solution.
if (isfield (lambda, 'eqlin'))
  y = -lambda.eqlin;
  z = c - A'*y;
end
f_objective = [obj.c'*x; obj.b'*y];
solver_info.name = 'linprog';
switch (flag)
  case 1
    solver_info.termination = 'Normal termination';
  case -2
    solver_info.termination = 'Primal infeasible';
  case -3
    solver_info.termination = 'Dual infeasible';
  case -5
    solver_info.termination = 'Primal and dual infeasibile';
  otherwise
    solver_info.termination = 'Unknown';
end

obj.add_solution (sol_type, x, y, z, f_objective, solver_info);

end
