function obj = solve_glpk (obj)
% SOLVE_GLPK  Approximately solve conic problem instance with GLPK.
%
%   For more information about GLPK, see:
%
%     https://www.gnu.org/software/glpk/
%
%   See also vsdp.solve.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

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
[A, b, c] = obj.get_perturbed_midpoint_problem ();
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
elapsed_time = toc;

% Store solution.
if (isfield (extra, 'lambda'))
  y = extra.lambda;
  z = c - A*obj.y;
end
f_objective = [obj.c'*x; obj.b'*y];
switch (errnum)
  case 0
    info = 0; % normal termination
  case 10
    info = 1; % primal infeasible
  case 11
    info = 2; % dual infeasible
  case 15
    info = 3; % primal and dual infeasible
  otherwise
    info = -1; % an error occured
end

obj.add_solution(x, y, z, f_objective, obj.options.SOLVER, info, elapsed_time);

end

