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
  warning ('VSDP:solve_mosek:ignoreInitialGuess', ...
    ['solve_mosek: MOSEK does not support initial guesses (x0,y0,z0) ' ...
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

% Split into individual cones.
c = vsdp_indexable (c, obj);
A = vsdp_indexable (A, obj);

% Regard the free, linear, and second-order cone.
prob.a = sparse([A.f', A.l', A.q']);
prob.c = [c.f', c.l', c.q'];

% Define lower bounds for free, linear, and second-order cone variables.
prob.blx = [-inf * ones(1, obj.K.f), zeros(1, obj.K.l), ...
  -inf * ones(1, sum (obj.K.q))];

% If there are second-order cones.
if (sum (obj.K.q) > 0)
  [~, res] = mosekopt('symbcon');
  prob.cones.type   = repmat (res.symbcon.MSK_CT_QUAD, 1, length (obj.K.q));
  % The indices of all second-order cone variables.
  prob.cones.sub    = obj.K.idx.q(1,1):obj.K.idx.q(end,end);
  % The indices of all second-order cone starts.
  prob.cones.subptr = obj.K.idx.q(:,1)';
end

% If there are semidefinite cones.
if (sum (obj.K.s) > 0)
  prob.bardim = obj.K.s;
  sdp_dims = obj.K.dims(end-length(obj.K.s)+1:end);
  [~, prob.barc.subj, ...
    prob.barc.subk, prob.barc.subl, prob.barc.val] = to_mosek_fmt ( ...
    mat2cell(c.s, sdp_dims, 1));
  [prob.bara.subi, prob.bara.subj, ...
    prob.bara.subk, prob.bara.subl, prob.bara.val] = to_mosek_fmt ( ...
    mat2cell(A.s, sdp_dims, ones(1, obj.m)));
end

% As the VSDP format support equality contraints only, we have to set the
% lower and upper bound for 'A*x' equal to 'b'.
prob.blc = b';
prob.buc = b';

% Adapt output verbosity.
if (obj.options.VERBOSE_OUTPUT)
  echo_level = 3;
else
  echo_level = 0;
end

% Call solver.
tic;
[r, res] = mosekopt (sprintf ('minimize info echo(%d)', echo_level), prob);
solver_info.elapsed_time = toc;

% Store solution after normal termination.
if (isscalar(r) && isnumeric(r) && (r == 0))
  solver_info.name = 'mosek';
  solver_info.termination = 'Normal termination';
  x = [res.sol.itr.xx(:); to_vsdp_fmt(obj, res.sol.itr.barx)];
  y = res.sol.itr.y;
  z = vsdp.svec (obj, obj.c - obj.At*y, 1);
  f_objective = [obj.c'*x; obj.b'*y];
  obj.add_solution (sol_type, x, y, z, f_objective, solver_info);
end

end


function [subi, subj, subk, subl, val] = to_mosek_fmt (x)
% TO_MOSEK_FMT  Translate the VSDP quantities 'At' and 'c' to MOSEK format.
%
%   The output are four to five vectors of equal length representing data
%   for the non-zero values of the VSDP quantities
%
%      subi - constraint index (only relevant for 'At')
%      subj - cone       index
%      subk - row        index
%      subl - column     index
%      val  - non-zero  values
%
%   See also vsdp.

% Compute the lower triangular matrix in each cell.
x = cellfun(@(x) tril (vsdp.smat([], x, 1)), x, 'UniformOutput', false);

% Get the non-zero entries including the indices.
[subk, subl, val] = cellfun(@find, x, 'UniformOutput', false);

% Compute index vectors for the constraints (subi) and the cones (subj).
[ii, jj] = meshgrid (1:size(x,2), 1:size(x,1));
subi = cellfun (@(x, i) i * ones(size(x)), subk, num2cell(ii), ...
  'UniformOutput', false);
subj = cellfun (@(x, j) j * ones(size(x)), subk, num2cell(jj), ...
  'UniformOutput', false);

% Vectorize quantities
subi = vertcat (subi{:})'; % constraint index
subj = vertcat (subj{:})'; % cone       index
subk = vertcat (subk{:})'; % row        index
subl = vertcat (subl{:})'; % column     index
val  = vertcat (val{:})';  % non-zero  values
end


function x = to_vsdp_fmt (obj, x)
% TO_VSDP_FMT  Translate the MOSEK solution 'X' to VSDP format.
%
%   See also vsdp, vsdp.svec.

if (isempty (x))
  return;
end

% Offset to the semidefinite cones.
sdp_offset = obj.K.idx.s(1,1) - 1;
x = [zeros(sdp_offset, 1); x(:)];

% Index vector for sorting.
[~,~,~,vlidx] = vsdp.sindex (obj);

S = warning ('off', 'VSDP:svec:justScale');
x = vsdp.svec (obj, x([(1:sdp_offset)'; vlidx]), 2);
warning (S);

% Strip non-semidefinite portion.
x(1:sdp_offset) = [];
end
