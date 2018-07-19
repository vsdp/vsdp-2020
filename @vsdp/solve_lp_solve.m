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

% Check cones.
if ((sum (obj.K.q) > 0) || (sum (obj.K.s) > 0))
  error ('VSDP:solve_lp_solve:unsupportedCones', ...
    ['solve_lp_solve: Second-order cones (K.q) and semidefinite cones ', ...
    '(K.s) are not supported by LP_SOLVE']);
end

% In case of interval data solve midpoint problem.
A = mid (obj.At);
b = full (mid (obj.b));
c = full (mid (obj.c));
lbound = [ ...
  -inf(obj.K.f, 1); ...
  zeros(obj.K.l, 1)];           % lower bounds
vtypes = zeros (length (b), 1); % variable types: 0 == equality

% Call solver.
[~, obj.x, obj.y, stat] = lp_solve (-c, A', b, vtypes, lbound);

% Store solution.
obj.y = -obj.y;
obj.z = c - A*obj.y;

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
end
