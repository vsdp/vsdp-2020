function obj = solve_mosek (obj, sol_type)
% SOLVE_MOSEK  Approximately solve conic problem instance with MOSEK.
%
%   For more information on the MOSEK format, see:
%
%     [1] https://docs.mosek.com/8.1/toolbox/data-types.html.
%
%   See also vsdp.solve.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (1, 2);

% Check solver availability.
if (exist ('mosekopt', 'file') ~= 3)
  error ('VSDP:solve_mosek:notAvailable', ...
    ['solve_mosek: MOSEK does not seem to be ready.  ', ...
    'To select another solver, run:  %s.solve()'], inputname(1));
end

if (nargin == 1)
  sol_type = 'Approximate';
end
[A, b, c] = obj.get_midpoint_problem_data (sol_type);

% Should initial solution guess be taken into account?
if (obj.options.USE_INITIAL_GUESS)
  warning ('VSDP:solve_sedumi:ignoreInitialGuess', ...
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
[r, res] = mosekopt('symbcon');

prob.bardim = obj.K.s;
for i = 1:length(obj.K.s)
  N = obj.K.s(i);
  cdim = N * (N + 1) / 2;
  % cone number
  prob.barc.subj = ones(cdim, 1) * i;
  % colums
  prob.barc.subl = zeros(cdim, 1);
  prob.barc.subl([1, cumsum(1:N-1) + 1],1) = 1;
  prob.barc.subl = cumsum(prob.barc.subl);
  % rows
  prob.barc.subk = ones(cdim, 1);
  prob.barc.subk(cumsum(1:N-1) + 1,1) = 0:-1:(-N+2);
  prob.barc.subk = cumsum(prob.barc.subk);
  prob.barc.val(c(1:K.s));
end

% Call solver.
tic;
[r,res] = mosekopt('minimize info',prob);
solver_info.elapsed_time = toc;

% Store solution.
x = vsdp.svec (obj, x, 2);
z = vsdp.svec (obj, c - A*y, 1);
f_objective = [obj.c'*x; obj.b'*y];
solver_info.name = 'mosek';
switch (info.pinf + 2*info.dinf)
  case 0
    solver_info.termination = 'Normal termination';
  case 1
    solver_info.termination = 'Primal infeasible';
  case 2
    solver_info.termination = 'Dual infeasible';
  case 3
    solver_info.termination = 'Primal and dual infeasibile';
  otherwise
    solver_info.termination = 'Unknown';
end

obj.add_solution (sol_type, x, y, z, f_objective, solver_info);

end
