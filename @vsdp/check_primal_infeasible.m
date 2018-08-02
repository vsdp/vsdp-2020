function obj = check_primal_infeasible (obj)
% CHECK_PRIMAL_INFEASIBLE  Infeasibility check for VSDP conic program.
%
%   obj.check_primal_infeasible()  Check VSDP object 'obj' to be infeasible.
%
%       Using a theorem of alternatives to claim a conic program to be primal
%       infeasible (see https://vsdp.github.io/references.html#Jansson2007 for
%       proofs and details):
%
%       Let y satisfy
%       
%         (1) A'*y in K^*  and  y'*b < 0
%
%       then the system
%
%         (2) A*x = b, x in K
%
%       has no solution.  Either (1) or (2) can be satisfied, but not both.
%
%   Example:
%
%       %  min <[0 0; 0 0], X>
%       % s.t. <[1 0; 0 0], X> = [e];
%       %      <[0 1; 1 d], X> = [1];
%       %                   X in K
%
%       e = -0.01;  % Infeasible, because X(1,1) <= 0!
%       d = 0.1;
%       A1 = [1 0; 0 0];
%       A2 = [0 1; 1 d];
%       b = [e; 1];
%       c = zeros(4,1);
%       K.s = 2;
%       A = [A1(:), A2(:)];  % vectorize
%       obj = vsdp(A, b, c, K).check_primal_infeasible()
%
%   See also vsdp, vsdp.rigorous_lower_bound, vsdp.rigorous_upper_bound.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (1, 1);
% If the problem was not approximately solved before, do it now.
if (isempty (obj.solutions('Approximate')))
  warning ('VSDP:check_primal_feasible:noApproximateSolution', ...
    ['check_primal_feasible: The conic problem has no approximate ', ...
    'solution yet, which is now computed using ''%s''.'], obj.options.SOLVER);
  obj.solve (obj.options.SOLVER, 'Approximate');
end
y = obj.solutions('Approximate').y;
if (isempty (y) || ~all (isfinite (y)))
    error ('VSDP:check_primal_feasible:badApproximateSolution', ...
    ['check_primal_feasible: The approximate dual solution ''y'' is ', ...
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
    error ('VSDP:check_primal_feasible:noBoundsForFreeVariables', ...
      ['check_primal_feasible: Could not find a rigorous solution of ', ...
      'the linear system of free variables.']);
  end
else
  y = intval (y);
end

% Step 1: If  y'*b  is not negative, we cannot claim anything about the
%         primal problem feasibility.
if (sup (obj.b' * y) < 0)
  % Step 2: Rigorous lower cone bound for At*y with free variables.
  lb = obj.rigorous_lower_cone_bound (obj.At * y, 1, false);
  % Step 3: If 'At*y' is in the cone, 'y' is a rigorous certificate of primal
  %         infeasibility.
  if (all (lb >= 0))
    is_infeasible = true;
    solver_info.elapsed_time = toc(cpi);
    solver_info.termination  = 'Normal termination';
    obj.add_solution ('Certificate primal infeasibility', [], y, lb, ...
      [is_infeasible, nan], solver_info);
  end
end

is_infeasible = false;
solver_info.elapsed_time = toc(cpi);
solver_info.termination  = 'Normal termination';
obj.add_solution ('Certificate primal infeasibility', [], [], [], ...
  [is_infeasible, nan], solver_info);

end
