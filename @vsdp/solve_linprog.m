function obj = solve_linprog (obj)
% SOLVE_LINPROG  Approximately solve conic problem instance with LINPROG.
%
%   For more information about LINPROG, see:
%
%     https://www.mathworks.com/help/optim/ug/linprog.html
%
%   See also vsdp.solve.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% Check cones.
if ((sum (obj.K.q) > 0) || (sum (obj.K.s) > 0))
  error ('VSDP:solve:unsupportedCones', ...
    ['solve: Second-order cones (K.q) and semidefinite cones (K.s) ', ...
    'are not supported by ''%s'''], obj.options.SOLVER);
end

options = optimoptions ('linprog', 'Algorithm', 'interior-point-legacy');
if (~obj.options.VERBOSE_OUTPUT)
  options = optimoptions (options, 'Display', 'off');
end

% In case of interval data solve midpoint problem.
A = mid (obj.At);
b = full (mid (obj.b));
c = full (mid (obj.c));
lbound = [ ...
  -inf(obj.K.f, 1); ...
  zeros(obj.K.l, 1)];       % lower bound
ubound = inf(length(c),1);  % upper bound

% Should initial solution guess be taken into account?
if (obj.options.USE_STARTING_POINT)
  x0 = obj.x;
else
  x0 = [];
end

% Call solver
[obj.x, ~, flag, ~, lambda] = linprog ...
  (c, [], [], A', b, lbound, ubound, x0, options);

% Transform results to VSDP format
if (isfield (lambda, 'eqlin'))
  obj.y = -lambda.eqlin;
  obj.z = c - A*obj.y;
end

switch (flag)
  case 1
    info = 0; % normal termination
  case -2
    info = 1; % primal infeasible
  case -3
    info = 2; % dual infeasible (primal unbounded)
  case -5
    info = 3; % primal and dual infeasible
  otherwise
    info = -1; % an error occured
end
end
