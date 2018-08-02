function obj = check_dual_infeasible (obj)
% CHECK_PRIMAL_INFEASIBLE  Dual infeasibility check for VSDP conic program.
%
%   obj.check_dual_infeasible()  Check VSDP object 'obj' to be dual infeasible.
%
%       Using a theorem of alternatives to claim a conic program to be primal
%       infeasible (see https://vsdp.github.io/references.html#Jansson2007 for
%       proofs and details):
%
%       Let x satisfy
%       
%         (1) A*x = 0,  x in K,  and  c'*x < 0
%
%       then the system
%
%         (2) z := c - A'*y,  z in K^*
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
%       obj = vsdp(A, b, c, K).check_dual_infeasible()
%
%       K.s = 2;
%       A1 = [-1 0; 0  0];
%       A2 = [ 0 0; 0 -1];
%        c = [0 1 1 0]';
%        b = [-1; 0];
%       A = [A1(:), A2(:)];  % vectorize
%       obj = vsdp(A, b, c, K).check_dual_infeasible()
%
%       K.s = 3;
%       A1 = [1 0 0; 0 0 0; 0 0 0];
%       A2 = [0 0 1; 0 1 0; 1 0 0];
%        c = [1 0 0; 0 1 0; 0 0 0];
%        b = [0; 1];
%       A = [A1(:), A2(:)];  % vectorize
%       obj = vsdp(A, b, c(:), K).check_dual_infeasible()
%
%   See also vsdp, vsdp.rigorous_lower_bound, vsdp.rigorous_upper_bound,
%            vsdp.check_primal_infeasible.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (1, 1);
% If the problem was not approximately solved before, do it now.
if (isempty (obj.solutions('Approximate')))
  warning ('VSDP:check_dual_feasible:noApproximateSolution', ...
    ['check_dual_feasible: The conic problem has no approximate ', ...
    'solution yet, which is now computed using ''%s''.'], obj.options.SOLVER);
  obj.solve (obj.options.SOLVER, 'Approximate');
end
x = obj.solutions('Approximate').x;
if (isempty (x) || ~all (isfinite (x)))
    error ('VSDP:check_dual_feasible:badApproximateSolution', ...
    ['check_dual_feasible: The approximate primal solution ''x'' is ', ...
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
    error ('VSDP:check_dual_feasible:ulsNotSolvable', ...
    'check_dual_feasible: Cannot solve the underdetermined linear system.');
  end
  lb = obj.rigorous_lower_cone_bound (x, 1/2, false);
  if (all (lb >= 0))
    is_infeasible = true;
    solver_info.elapsed_time = toc(cdi);
    solver_info.termination  = 'Normal termination';
    obj.add_solution ('Certificate dual infeasibility', x, [], lb, ...
      [nan, is_infeasible], solver_info);
  end
end

is_infeasible = false;
solver_info.elapsed_time = toc(cdi);
solver_info.termination  = 'Normal termination';
obj.add_solution ('Certificate dual infeasibility', [], [], [], ...
  [nan, is_infeasible], solver_info);

end
