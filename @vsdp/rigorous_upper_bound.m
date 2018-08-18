function obj = rigorous_upper_bound (obj, ybnd)
% RIGOROUS_UPPER_BOUND  Rigorous upper bound for conic programming.
%
%   obj.rigorous_upper_bound()  Compute a rigorous upper bound of the dual
%      optimal value  b'*y  and a rigorous enclosure of a primal strict
%      feasible (near optimal) solution of a VSDP object 'obj'.
%
%   obj.rigorous_upper_bound(ybnd)  Optionally a priori known finite upper
%      bounds of the dual optimal solution 'y' can be provided to speed up
%      the computation time, but no rigorous enclosure of a primal strict
%      feasible (near optimal) solution is computed.  It is recommend to prefer
%      infinite or no bounds at all instead of unreasonable large bounds.  This
%      improves the quality of the lower bound in many cases, but may increase
%      the computational time.
%
%      The a priori bounds must fulfill the following dual boundedness
%      assumption:
%
%        -ybnd(i) <= y(i) <= ybnd(i),  i = 1:length(b).
%
%
%   See also vsdp, vsdp.rigorous_lower_bound.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% Validate dual upper bounds.
if (nargin == 2)
  ybnd = ybnd(:);
  if (length (ybnd) ~= length (obj.b))
    error ('VSDP:rigorous_upper_bound:badYBND', ...
      ['rigorous_upper_bound: The length of the upper bound vector ', ...
      '''ybnd'' must be %d, but is %d.'], length (obj.b), length (ybnd));
  end
else  % If dual upper bounds are not given, use infinite bounds.
  ybnd = inf(length (obj.b), 1);
end

% If the problem was not approximately solved before, do it now.
if (isempty (obj.solutions.approximate))
  warning ('VSDP:rigorous_upper_bound:noApproximateSolution', ...
    ['rigorous_upper_bound: The conic problem has no approximate ', ...
    'solution yet, which is now computed using ''%s''.'], obj.options.SOLVER);
  obj.solve (obj.options.SOLVER, 'Approximate');
end

% If any upper bound 'ybnd' is not finite, proceed with other algorithm.
if (~all (isfinite (ybnd)))
  if (nargin == 2)
    warning ('VSDP:rigorous_upper_bound:infiniteBounds', ...
      ['rigorous_upper_bound: At least one element in the bounds ', ...
      '''ybnd'' was not finite.  Using algorithm for infinite bounds.']);
  end
  obj = rigorous_upper_bound_infinite_bounds (obj);
  return;
end
rub = tic;

% Step 1: Project approximate solution into the respective cones, e.g.
%         compute x^+.

x = vsdp_indexable (intval (obj.solutions.approximate.x), obj);

% LP cones
x.l = max (x.l, 0);

% Second-order cones
for j = 1:length(K.q)
  xq = x.q(j);
  % Very simple projection: 'x(1) >= ||x(2:end)||' holds, if 'x(1)' is set to
  % the maximum of both sites of the inequality.
  xq(1) = max (xq(1), sup (norm (xq(2:end))));
  x.q(j) = xq;
end

% SDP cones
for j = 1:length(K.s)
  % Find minimal eigenvalue.  For smat, remember that 'x' is scaled by mu = 2.
  E_min = min (inf_ (vsdp.verify_eigsym (vsdp.smat ([], x.s(j), 1/2))));
  % If the matrix in not positive semidefinite, perform a simple cone
  % projection, by shifting the diagonal by the minimal eigenvalue.
  if (E_min < 0)
    xs = x.s(j);
    idx = cumsum (1:K.s(j));  % Index vector for diagonal entries.
    xs(idx) = xs(idx) - E_min;
    x.s(j) = xs;
  end
end

% Step 2: Compute rigorous enclosure on '|A * x^+ - b|'.
ru = abs (obj.At' * x - obj.b);
% Step 3: Compute rigorous upper bound of the dual optimal value.  As 'x' and
%         'ru' are interval quantities, everything operation is rigorous.
fU = sup (obj.c' * x + ru' * ybnd);
solver_info.iter = 0;
solver_info.termination = 'Normal termination';
solver_info.elapsed_time = toc(rub);
obj.add_solution ('Rigorous upper bound', nan, [], nan, [nan, fU], solver_info);
end

function obj = rigorous_upper_bound_infinite_bounds (obj)
% Create helper structures.
[~, num_of_bounds, vidx, sdp_matrix] = obj.rigorous_lower_cone_bound();
% Bound and perturbation parameter.  alpha = obj.options.ALPHA
%
%  epsilon = epsilon + (alpha.^k) .* dl,  where dl < 0.
%
k         = zeros (num_of_bounds, 1);  % Counter for perturbation.
epsilon   = zeros (num_of_bounds, 1);  % Factor  for perturbation.
x_epsilon = zeros (obj.n, 1);          % 'epsilon' translated to 'x'.

x = obj.solutions.approximate.x;
rub = tic;
iter = 0;
while (iter <= obj.options.ITER_MAX)
  % Step 1: Compute rigorous enclosure for 'A * x = b'.
  x = vsdp.verify_uls (obj, obj.At', obj.b, x);
  if ((~isintval (x) || any (isnan (x))))
    error ('VSDP:rigorous_upper_bound:noUlsEnclosure', ...
      ['rigorous_upper_bound: Could not find a rigorous solution for the ', ...
      'linear system of constraints.']);
  end
  
  % Step 2: Verified lower bounds on 'x' for each cone.
  lb = obj.rigorous_lower_cone_bound (x, 1/2, false);
    
  % Step 3: Cone feasibility check and upper bound computation:
  %
  % If all lower bounds 'lb' on 'x' are non-negative, there are no cone
  % violations and 'x' is an enclosure of a primal strict feasible (near
  % optimal) solution.
  idx = (lb < 0);
  
  if (~any (idx))  % If no violations.
    fU = sup (obj.c' * x);
    solver_info.termination = 'Normal termination';
    break;  % SUCCESS
  end
  
  % Step 4: Perturb midpoint problem, such that
  %
  %    P(epsilon) = (mid([A]), mid([b]), mid([c]) + c_epsilon)
  %
  k      (idx) = k      (idx) + 1;
  epsilon(idx) = epsilon(idx) - (obj.options.ALPHA .^ k (idx)) .* lb(idx);
  if (~all (isfinite (epsilon)))
    error ('VSDP:rigorous_upper_bound:infinitePerturbation', ...
      'rigorous_upper_bound: Perturbation exceeds finite range.');
  end
  x_epsilon(vidx) = [epsilon(1:(obj.K.f + obj.K.l + length (obj.K.q)),1); ...
    sdp_matrix * epsilon((end - length (obj.K.s) + 1):end,1)];
  
  % Display short perturbation statistic.
  iter = iter + 1;
  if (obj.options.VERBOSE_OUTPUT)
    fprintf ('\n\n');
    fprintf ('--------------------------------------------------\n');
    fprintf ('  VSDP.RIGOROUS_UPPER_BOUND  (iteration %d)\n', iter);
    fprintf ('--------------------------------------------------\n');
    fprintf ('  Violated cones    (lb < 0): %d\n',      sum (idx));
    fprintf ('  Max. violation     min(lb): %+.2e\n',   min (lb));
    fprintf ('  Perturbation  max(epsilon): %+.2e\n\n', max (epsilon));
    fprintf ('  Solve perturbed problem using ''%s''.\n', obj.options.SOLVER);
    fprintf ('--------------------------------------------------\n\n');
  end
  
  % Step 5: Solve perturbed problem.
  %
  % Set perturbation parameters.
  obj.perturbation.b = (x_epsilon' * obj.At)';  % b := b - A * x(epsilon)
  obj.perturbation.c = [];
  
  obj.solve (obj.options.SOLVER, 'Rigorous upper bound');
  sol = obj.solutions.rigorous_upper_bound;
  if ((~strcmp (sol.solver_info.termination, 'Normal termination')) ...
      || isempty (sol.x) || any (isnan (sol.x)) || any (isinf (sol.x)))
    error ('VSDP:rigorous_upper_bound:unsolveablePerturbation', ...
      ['rigorous_upper_bound: Conic solver could not find a solution ', ...
      'for perturbed problem']);
  end
  % Store last successful solver info and new dual solution.
  solver_info = sol.solver_info;
  x = sol.x + x_epsilon;  % Undo perturbation.
end

%TODO
if (iter == obj.options.ITER_MAX)
  error ('maximum number of iterations reached');
end

% Update solver info and store solution.
solver_info.elapsed_time = toc(rub);
solver_info.iter = iter;
obj.add_solution ('Rigorous upper bound', x, [], lb, [nan, fU], solver_info);
end
