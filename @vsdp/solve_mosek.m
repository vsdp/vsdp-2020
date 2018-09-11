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

% Split into individual SDP cones.
c = mat2cell(c, obj.K.dims, 1);
A = mat2cell(A, obj.K.dims, ones(1, obj.m));

% Compute the lower triangular matrix in each cell.
c = cellfun(@(x) tril (vsdp.smat([], x, 1)), c, 'UniformOutput', false);
A = cellfun(@(x) tril (vsdp.smat([], x, 1)), A, 'UniformOutput', false);

% Get the non-zero entries including the indices.
[prob.barc.subk, prob.barc.subl, prob.barc.val] = cellfun(@find, c, ...
  'UniformOutput', false);
% 
prob.barc.subj = cellfun(@(x, i) i * ones(size(x)), prob.barc.subk, ...
  num2cell(1:3)', 'UniformOutput', false);
prob.barc.subj = vertcat(prob.barc.subj{:})'; % cone   index
prob.barc.subk = vertcat(prob.barc.subk{:})'; % row    index
prob.barc.subl = vertcat(prob.barc.subl{:})'; % column index
prob.barc.val  = vertcat(prob.barc.val{:})';  % values
prob.bardim    = obj.K.s;

prob.a = [];
prob.c = [];

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
