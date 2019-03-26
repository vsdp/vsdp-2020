function obj = check_primal_infeasible (obj, y)
% CHECK_PRIMAL_INFEASIBLE  Primal infeasibility check for VSDP conic program.
%
%   obj.check_primal_infeasible()  Check VSDP object 'obj' to be primal
%                                  infeasible by using the approximate
%                                  solution.
%
%   obj.check_primal_infeasible(y)  Check VSDP object 'obj' to be primal
%                                   infeasible by using 'y'.
%
%       Use a theorem of alternatives to claim a conic program to be primal
%       infeasible (see https://vsdp.github.io/references.html#Jansson2007
%       for proofs and details):
%
%       Let 'y' satisfy
%
%         (1) -A'*y in K^*  and  y'*b > 0
%
%       then the system
%
%         (2) A*x = b, x in K
%
%       has no solution.  Either (1) or (2) can be satisfied, but not both.
%
%   Example:
%
%       EPSILON = -0.01;
%       DELTA = 0.1;
%       blk(1,:) = {'s'; 2};
%       C{1,1} = [0 0; 0 0];
%       A{1,1} = [1 0; 0 0];
%       A{2,1} = [0 1; 1 DELTA];
%       b = [EPSILON; 1];
%
%       obj = vsdp(blk,A,C,b).solve ('sdpt3') ...
%                            .rigorous_lower_bound () ...
%                            .rigorous_upper_bound () ...
%                            .check_primal_infeasible () ...
%                            .check_dual_infeasible ()
%
%   See also vsdp, vsdp.rigorous_lower_bound, vsdp.rigorous_upper_bound,
%            vsdp.check_dual_infeasible.
%

% Copyright 2004-2019 Christian Jansson (jansson@tuhh.de)

narginchk (1, 2);

if (nargin == 1)
  % If the problem was not approximately solved yet, do it now.
  if (isempty (obj.solutions.approximate))
    if (~isempty (obj.options.SOLVER))
      warning ('VSDP:check_primal_infeasible:noApproximateSolution', ...
        ['check_primal_infeasible: The conic problem has no approximate ', ...
        'solution yet, which is now computed using ''%s''.'], ...
        obj.options.SOLVER);
      obj.solve (obj.options.SOLVER, 'Approximate');
    else
      obj.solve ();  % Interactive mode.
    end
  end
  y = obj.solutions.approximate.y;
else  % Use given solution.
  y = y(:);
  if (length (y) ~= obj.m)
    error ('VSDP:check_primal_infeasible:badApproximateSolution', ...
      ['check_primal_infeasible: The approximate dual solution ''y'' ', ...
      'must be a vector of length %d, but is %d.'], obj.m, length (y));
  end
end
if (isempty (y) || ~all (isfinite (y)))
  error ('VSDP:check_primal_infeasible:badApproximateSolution', ...
    ['check_primal_infeasible: The approximate dual solution ''y'' is ', ...
    'empty or contains infinite entries.']);
end

cpi = tic;

% Ensure, that the approximate dual solution 'y' solves the free variable part,
% see
%
%   https://vsdp.github.io/references.html#Anjos2007
%
% for details.
if (obj.K.f > 0)
  % Compute rigorous enclosure for underdetermined linear interval system of
  % dimension (K.f x m):
  %
  %    At.f * y = 0
  %
  y = vsdp.verify_uls (obj, obj.At(1:obj.K.f,:), zeros(obj.K.f, 1), y);
  if ((~isintval (y) || ~all (isfinite (y))))
    error ('VSDP:check_primal_infeasible:noBoundsForFreeVariables', ...
      ['check_primal_infeasible: Could not find a rigorous solution of ', ...
      'the linear system of free variables.']);
  end
else
  y = intval (y);
end

% Step 1: If  b'*y  is not positive, we cannot claim anything about the
%         primal problem feasibility.
if (inf (obj.b' * y) > 0)
  % Step 2: Rigorous lower cone bound for '-At*y' with free variables.
  lb = obj.rigorous_lower_cone_bound (-obj.At * y, 1, false);
  % Step 3: If '-At*y' is in the cone, 'y' is a rigorous certificate of primal
  %         infeasibility.
  if (all (lb >= 0))
    is_infeasible = true;
    solver_info.elapsed_time = toc(cpi);
    solver_info.termination  = 'Normal termination';
    obj.add_solution ('Certificate primal infeasibility', [], y, lb, ...
      [is_infeasible, nan], solver_info);
    return;
  end
end

is_infeasible = false;
solver_info.elapsed_time = toc(cpi);
solver_info.termination  = 'Normal termination';
obj.add_solution ('Certificate primal infeasibility', [], [], [], ...
  [is_infeasible, nan], solver_info);

end
