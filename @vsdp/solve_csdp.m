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

% Check cones.
if (sum (obj.K.q) > 0)
  error ('VSDP:solve_csdp:unsupportedCone', ...
    'solve_csdp: Second-order cones (K.q) are not supported by CSDP.');
end

if (~isempty (obj.options.SOLVER_OPTIONS))
  pars = obj.options.SOLVER_OPTIONS;
else
  pars = [];
end
if (~obj.options.VERBOSE_OUTPUT)
  pars.printlevel = 0;
end

% In case of interval data solve midpoint problem.
A = mid (obj.At);
b = full (mid (obj.b));
c = full (mid (obj.c));

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
if (obj.options.USE_STARTING_POINT)
  % CSDP requires all variables to be given!
  x0 = full (obj.x);
  y0 = full (obj.y);
  z0 = full (obj.z);
  if (~isempty (x0) && ~isempty(y0) && ~isempty (z0))
    args = [args, {x0, y0, z0}];
  end
end

% Call solver.
[obj.x, obj.y, obj.z, INFO] = csdp (args{:});

if (any (INFO == [0, 1, 2]))
  info = INFO;
end

end
