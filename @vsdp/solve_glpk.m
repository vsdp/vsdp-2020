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

% Check cones.
if ((sum (obj.K.q) > 0) || (sum (obj.K.s) > 0))
  error ('VSDP:solve_glpk:unsupportedCones', ...
    ['solve_glpk: Second-order cones (K.q) and semidefinite cones (K.s) ', ...
    'are not supported by GLPK.']);
end

if (~isempty (obj.options.SOLVER_OPTIONS))
  param = obj.options.SOLVER_OPTIONS;
else
  param = [];
end
if (~obj.options.VERBOSE_OUTPUT)
  param.msglev = 0;
end

% In case of interval data solve midpoint problem.
A = mid (obj.At);
b = full (mid (obj.b));
c = full (mid (obj.c));

sense = 1;  % Minimization
lbound = [ ...
  -inf( obj.K.f, 1); ...
  zeros(obj.K.l, 1)];                 % lower bound
ubound = inf (length (c), 1);         % upper bound
vtype = repmat ('C', length (c), 1);  % variable   types: continuous
ctype = repmat ('S', length (b), 1);  % constraint types: equality

% Call solver
[obj.x, ~, errnum, extra] = glpk ...
  (c, A', b, lbound, ubound, ctype, vtype, sense, param);

% Transform results to VSDP format
if (isfield (extra, 'lambda'))
  obj.y = extra.lambda;
  obj.z = c - A*obj.y;
end
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