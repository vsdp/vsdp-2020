function obj = solve_lp_solve (obj, sol_type)
% SOLVE_LP_SOLVE  Approximately solve conic problem instance with LP_SOLVE.
%
%   For more information about LP_SOLVE, see:
%
%     http://lpsolve.sourceforge.net/5.5/index.htm
%
%   See also vsdp.solve.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (1, 2);

% Check solver availability.
if (exist ('lp_solve', 'file') ~= 2)
  error ('VSDP:solve_lp_solve:notAvailable', ...
    ['solve_lp_solve: LP_SOLVE does not seem to be ready.\n\n', ...
    'To select another solver, run:  %s.solve()'], inputname(1));
end

% Check cones.
if ((sum (obj.K.q) > 0) || (sum (obj.K.s) > 0))
  error ('VSDP:solve_lp_solve:unsupportedCones', ...
    ['solve_lp_solve: Second-order cones (K.q) and semidefinite cones ', ...
    '(K.s) are not supported by LP_SOLVE']);
end

if (nargin == 1)
  sol_type = 'Approximate';
  [A, b, c] = deal (mid (obj.At), mid (obj.b), mid (obj.c));
else
  [A, b, c] = obj.get_perturbed_midpoint_problem ();
end

% Should initial solution guess be taken into account?
if (obj.options.USE_INITIAL_GUESS)
  warning ('VSDP:solve_lp_solve:ignoreInitialGuess', ...
    ['solve_lp_solve: LP_SOLVE does not support initial guesses (x0,y0,z0) ' ...
    'and proceeds without them.']);
end

% Prepare data for solver.
[A, b, c] = deal (A', full (b), full (c));
lbound = [ ...
  -inf(obj.K.f, 1); ...
  zeros(obj.K.l, 1)];           % lower bounds
vtypes = zeros (length (b), 1); % variable types: 0 == equality

% Call solver.
tic;
[~, x, y, stat] = lp_solve (-c, A, b, vtypes, lbound);
solver_info.elapsed_time = toc;

% Store solution.
y = -y;
z = c - A'*y;
f_objective = [obj.c'*x; obj.b'*y];
solver_info.name = 'lp_solve';
switch (stat)
  case 0
    solver_info.termination = 'Normal termination';
  case 2
    solver_info.termination = 'Primal infeasible';
  case 3
    solver_info.termination = 'Dual infeasible';
  otherwise
    solver_info.termination = 'Unknown';
end

obj.add_solution (sol_type, x, y, z, f_objective, solver_info);

end
