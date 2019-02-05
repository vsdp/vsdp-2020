function obj = check_dual_infeasible (obj, x)
% CHECK_PRIMAL_INFEASIBLE  Dual infeasibility check for VSDP conic program.
%
%   obj.check_dual_infeasible()  Check VSDP object 'obj' to be dual infeasible
%                                by using the approximate solution.
%
%   obj.check_dual_infeasible(x)  Check VSDP object 'obj' to be dual infeasible
%                                 by using 'x'.
%
%       Use a theorem of alternatives to claim a conic program to be primal
%       infeasible (see https://vsdp.github.io/references.html#Jansson2007
%       for proofs and details):
%
%       Let 'x' satisfy
%
%         (1) A*x = 0,  x in K,  and  c'*x < 0
%
%       then the system
%
%         (2) z := c - A'*y,  z in K^*
%
%       has no solution.  Either (1) or (2) can be satisfied, but not both.
%
%
%   See also vsdp, vsdp.rigorous_lower_bound, vsdp.rigorous_upper_bound,
%            vsdp.check_primal_infeasible.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (1, 2);

if (nargin == 1)
  % If the problem was not approximately solved yet, do it now.
  if (isempty (obj.solutions.approximate))
    if (~isempty (obj.options.SOLVER))
      warning ('VSDP:check_dual_infeasible:noApproximateSolution', ...
        ['check_dual_infeasible: The conic problem has no approximate ', ...
        'solution yet, which is now computed using ''%s''.'], ...
        obj.options.SOLVER);
      obj.solve (obj.options.SOLVER, 'Approximate');
    else
      obj.solve ();  % Interactive mode.
    end
  end
  x = obj.solutions.approximate.x;
else  % Use given solution.
  x = x(:);
  if (length (x) ~= obj.n)
    error ('VSDP:check_primal_infeasible:badApproximateSolution', ...
      ['check_primal_infeasible: The approximate vectorized primal ', ...
      'solution ''x'' must be a vector of length %d, but is %d.'], obj.n, ...
      length (x));
  end
end
if (isempty (x) || ~all (isfinite (x)))
  error ('VSDP:check_dual_infeasible:badApproximateSolution', ...
    ['check_dual_infeasible: The approximate primal solution ''x'' is ', ...
    'empty or contains infinite entries.']);
end

cdi = tic;

% Step 1: Compute beta.
beta = obj.c' * x;
% Step 2: If 'beta' is not negative, we cannot claim anything about the
%         dual problem feasibility.
if (sup (beta) < 0)
  % Step 2: Find a rigorous inclusion for the system
  %
  %         A  * x = 0
  %         c' * x = beta
  %
  rhs = intval ([zeros(obj.m, 1); inf_(beta)]);
  x = vsdp.verify_uls (obj, [obj.At'; obj.c'], rhs, x);
  if (~isintval (x) || ~all (isfinite (x)))
    error ('VSDP:check_dual_infeasible:ulsNotSolvable', ...
      'check_dual_infeasible: Cannot solve the underdetermined linear system.');
  end
  lb = obj.rigorous_lower_cone_bound (x, 1/2, false);
  if (all (lb >= 0))
    is_infeasible = true;
    solver_info.elapsed_time = toc(cdi);
    solver_info.termination  = 'Normal termination';
    obj.add_solution ('Certificate dual infeasibility', x, [], lb, ...
      [nan, is_infeasible], solver_info);
    return;
  end
end

is_infeasible = false;
solver_info.elapsed_time = toc(cdi);
solver_info.termination  = 'Normal termination';
obj.add_solution ('Certificate dual infeasibility', [], [], [], ...
  [nan, is_infeasible], solver_info);

end
