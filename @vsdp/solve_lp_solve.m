function obj = solve_lp_solve (obj)
% SOLVE_LP_SOLVE  Approximately solve conic problem instance with LP_SOLVE.
%
%   For more information about LP_SOLVE, see:
%
%     http://lpsolve.sourceforge.net/5.5/index.htm
%
%   See also vsdp.solve.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

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

% Prepare data for solver.
[A, b, c] = obj.get_perturbed_midpoint_problem ();
[A, b, c] = deal (A', full (b), full (c));
lbound = [ ...
  -inf(obj.K.f, 1); ...
  zeros(obj.K.l, 1)];           % lower bounds
vtypes = zeros (length (b), 1); % variable types: 0 == equality

% Call solver.
tic;
[~, x, y, stat] = lp_solve (-c, A, b, vtypes, lbound);
elapsed_time = toc;

% Store solution.
y = -y;
z = c - A*obj.y;
switch (stat)
  case 0
    info = 0; % normal termination
  case 2
    info = 1; % primal infeasible
  case 3
    info = 2; % dual infeasible (primal unbounded)
  otherwise
    info = -1; % an error occured
end

obj.add_solution(x, y, z, f_objective, obj.options.SOLVER, info, elapsed_time);

end
