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
[~, prob.barc.subj, ...
  prob.barc.subk, prob.barc.subl, prob.barc.val] = to_mosek_fmt(c);
[prob.bara.subi, prob.bara.subj, ...
  prob.bara.subk, prob.bara.subl, prob.bara.val] = to_mosek_fmt(A);

prob.bardim = obj.K.s;

prob.a = sparse(obj.m,0);
prob.c = [];
prob.blc = obj.b';
prob.buc = obj.b';

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
  x = to_vsdp_fmt(obj, res.sol.itr.barx);
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
%    The computations are similar to those in vsdp.svec.  The ordering of
%    the entries inside the MOSEK solution for each SDP cone are shown
%    below in Matrix 'J'.  Those have to be translated to the VSDP ordering
%    shown in Matrix 'I' below.  See the following example for N = 5:
%
%            [ 1            ]         [ 1            ]
%            [ 2  3         ]         [ 2  7         ]
%        I = [ 4  5  6      ]     J = [ 3  8 13      ]
%            [ 7  8  9 10   ]         [ 4  9 14 19   ]
%            [11 12 13 14 15]         [ 5 10 15 20 25]
%
%    The difference to 'lidx' from vsdp.svec is that 'lidx' addresses the
%    lower triangular part for a quadratic NxN matrix, whereas 'idx' below
%    just addresses the vectorized
%
%   See also vsdp, vsdp.svec.

% Index vector for sorting.
idx = zeros(obj.n, 1);

% For each SDP cone
for j = 1:length(obj.K.s)
  N = obj.K.s(j);
  sidx = obj.K.idx.s(j,:);
  I = zeros (N);
  I(triu (true (N))) = 1:(N*(N+1)/2);
  I = I';
  I(I == 0) = [];
  % Insert indices shifted by the cone offset.
  idx(sidx(1):sidx(2),:) = I + sidx(1) - 1;
end
S = warning ('off', 'VSDP:svec:justScale');
x = vsdp.svec (obj, x(idx), 2);
warning (S);
end
