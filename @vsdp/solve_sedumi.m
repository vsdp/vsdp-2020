function obj = solve_sedumi (obj)
% SOLVE_SEDUMI  Approximately solve conic problem instance with SeDuMi.
%
%   For information about SeDuMi, see:
%
%      http://sedumi.ie.lehigh.edu/
%      https://github.com/sqlp/sedumi
%
%   See also vsdp.solve.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% In case of interval data solve midpoint problem.
A = mid (obj.At);
b = mid (obj.b);
c = mid (obj.c);

if (obj.options.USE_STARTING_POINT)
  warning ('VSDP:solve_sedumi:ignoreStartingPoint', ...
    ['solve_sedumi: SeDuMi does not support initial guesses (x0,y0,z0) ' ...
    'and proceeds without them.']);
end

if (~isempty (obj.options.SOLVER_OPTIONS))
  pars = obj.options.SOLVER_OPTIONS;
else
  pars = [];
end
if (~obj.options.VERBOSE_OUTPUT)
  pars.fid = 0;
end

A = vsdp.smat (obj, A, 1);
c = vsdp.smat (obj, c, 1);

% Call solver.
[x, y, info] = sedumi (A, b, c, obj.K, pars);

% Store results.
obj.x = vsdp.svec (obj, x, 2);
obj.y = y;
obj.z = vsdp.svec (obj, c - A*y, 2);

info = info.pinf + 2*info.dinf;

end
