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

% Check solver availability.
if (exist ('sedumi', 'file') ~= 2)
  error ('VSDP:solve_sedumi:notAvailable', ...
    ['solve_sedumi: SeDuMi does not seem to be ready.  ', ...
    'Did you run ''install_sedumi()'' inside the solver directory?\n\n', ...
    'To select another solver, run:  %s.solve()'], inputname(1));
end

% Should initial solution guess be taken into account?
if (obj.options.USE_STARTING_POINT)
  warning ('VSDP:solve_sedumi:ignoreStartingPoint', ...
    ['solve_sedumi: SeDuMi does not support initial guesses (x0,y0,z0) ' ...
    'and proceeds without them.']);
end

% Should special solver options be taken into account?
if (~isempty (obj.options.SOLVER_OPTIONS))
  pars = obj.options.SOLVER_OPTIONS;
else
  pars = [];
end

% Adapt output verbosity.
if (~obj.options.VERBOSE_OUTPUT)
  pars.fid = 0;
end

% Prepare data for solver.
[At, b, c] = obj.get_perturbed_midpoint_problem ();
A = vsdp.smat (obj, A, 1);
c = vsdp.smat (obj, c, 1);

% Call solver.
tic;
[x, y, info] = sedumi (A, b, c, obj.K, pars);
elapsed_time = toc;

% Store solution.
x = vsdp.svec (obj, x, 2);
z = vsdp.svec (obj, c - A*y, 1);
f_objective = [obj.c'*x; obj.b'*y];
info = info.pinf + 2*info.dinf;
obj.add_solution(x, y, z, f_objective, obj.options.SOLVER, info, elapsed_time);

end
