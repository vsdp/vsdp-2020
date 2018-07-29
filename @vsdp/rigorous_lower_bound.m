function [fL,y,dl,info] = rigorous_lower_bound (obj, xbnd)
% VSDPLOW  Verified lower bound for conic programming.
%
%   [fL,y,dl,info] = VSDPLOW(At,b,c,K,[],y0) Computes a verified lower bound of
%      the primal optimal value and a rigorous enclosure of dual strict feasible
%      (near optimal) solutions of a conic problem in the standard primal-dual
%      form.  This form and the block-diagonal format (A,b,c,K) is explained in
%      'mysdps.m'.
%
%         'y0'     A dual feasible (eps-optimal) solution of the same dimension
%                  as input b.  This solution can be computed using 'mysdps'.
%
%      The output is:
%
%         'fL'     Verified lower bound of the primal optimal value.
%
%         'y'      Rigorous enclosure of dual strict feasible solutions.
%
%         'dl'     Verified lower bounds of eigenvalues or spectral values of
%                  `z = c - A' * y`.
%
%         'info'   Struct containing further information.
%           - iter  The number of iterations.
%
%   VSDPLOW(A,b,c,K,x0,y0,z0) optionally provide the other approximate
%      solutions of 'mysdps' (x0 and z0).
%
%   VSDPLOW(A,b,c,K,[],y0,[],xu) optionally provide known finite upper bounds
%      of the eigenvalues or spectral values of the primal optimal solution x.
%      We recommend to use infinite bounds xu(j) instead of unreasonable large
%      bounds xu(j).  This improves the quality of the lower bound in many
%      cases, but may increase the computational time.
%
%   VSDPLOW(A,b,c,K,[],y0,[],[],opts) optionally provide a structure for
%      additional parameter settings, explained in vsdpinit.
%
%   See also mysdps, vsdpinit, vsdpup, vsdpinfeas.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% Check, if primal upper bounds are given, otherwise use default upper bounds.
num_of_bounds = obj.K.f + obj.K.l + length(obj.K.q) + length(obj.K.s);
if ((nargin < 2) || isempty (xbnd))
  xbnd = inf (num_of_bounds, 1);
else
  xbnd = xbnd(:);
  if (length (xbnd) ~= num_of_bounds)
    error ('VSDP:rigorous_lower_bound:badXBND', ...
      ['rigorous_lower_bound: The length of the upper bound vector ', ...
      '''xbnd'' must be %d, but is %d.'], num_of_bounds, length (xbnd));
  end
end

% Bound and pertubation parameter.  alpha = obj.options.ALPHA
%
%  epsilon = epsilon + (alpha.^k) .* dl_,  where dl_ < 0.
%
dl        = zeros (num_of_bounds, 1);
k         = zeros (num_of_bounds, 1);
epsilon   = zeros (num_of_bounds, 1);  % factor for perturbation
c_epsilon = zeros (obj.n, 1);          % 'epsilon' translated to 'c'.

% If the problem was not approximately solved before, do it now.
if (isempty (obj.solution))
  warning ('VSDP:rigorous_lower_bound:noApproximateSolution', ...
    ['rigorous_lower_bound: The conic problem has no approximate ', ...
    'solution yet, which is now computed using ''%s''.'], obj.options.SOLVER);
  obj.solve (obj.options.SOLVER);
end
old_solution = obj.solution;

% Algorithm for both finite and infinite upper bounds `xu`.
iter = 0;
while (iter <= obj.options.ITER_MAX)
  y = obj.solution.y;
  iter = iter + 1;
  
  idx = obj.K.idx.f;
  % If infinite upper bounds for free variables are given.  Ensure, that
  % the approximate dual solution 'y' solves the free variable part, see
  %
  %   https://vsdp.github.io/references.html#Anjos2007
  %
  % for details.
  if (any (isinf (xbnd(idx(1):idx(end)))))
    % Compute rigorous enclosure for underdetermined linear interval system of
    % dimension (K.f x m):
    %
    %    At.f * y = c.f
    %
    y = vsdp.verify_uls (obj, obj.At(idx(1):idx(end),:), ...
      obj.c(idx(1):idx(end)), y);
    if (~isintval (y) || any (isnan (y)))
      tidy_up (obj, old_solution);
      error ('VSDP:rigorous_lower_bound:noBoundsForFreeVariables', ...
        ['rigorous_lower_bound: Could not find a verified solution of the ', ...
        'linear  system of free variables.']);
    end
  end
  
  % Step 1: Compute rigorous enclosure [d] for  c - At*y.
  d = obj.c - obj.At * intval (y);
  
  % Step 2: Verified lower bounds on cone eigenvalues
  %
  % Free variables: if infinite upper bounds for free variables are given, we
  %                 ensured above, that the solution is contained, thus we
  %                 assume no defect!
  if (any (isinf (xbnd(idx(1):idx(end)))))
    dl(idx(1):idx(end)) = zeros (obj.K.f, 1);
  else
    dl(idx(1):idx(end)) = mag (d(idx(1):idx(end)));
  end
  
  % LP cones
  idx = obj.K.idx.l;
  dl(idx(1):idx(end)) = inf_ (d(idx(1):idx(end)));
  
  % Second-order cones
  if (~isempty (obj.K.q))
    idx = obj.K.idx.q;
    % dl = inf (d(1) - ||d(2:end)||)
    dl(idx(1,1) + (1:length(obj.K.q))) = ...
      inf_ (d(idx(:,1)) - norm (d((idx(:,1)+1):idx(:,end))));
  end
  
  % SDP cones
  offset = obj.K.f + obj.K.l + length(obj.K.q);
  for j = 1:length(obj.K.s)
    idx = obj.K.idx.s(j,:);
    E_ = inf_ (vsdp.verify_eigsym (vsdp.smat ([], d(idx(1):idx(end)), 1)));
    % Weight the smallest eigenvalue lower bound E_ by the number of negative
    % eigenvalues.
    dl(j + offset) = min (E_) * sum (E_ < 0);
  end
  
  % Step 3: Cone feasibility check and lower bound computation:
  %
  % a) If the defect  d = c - A' * y  lies in each cone, then all lower bounds
  %    'dl' on 'd' are non-negative, 'y' is dual feasible, and  inf_ (b' * y)
  %    is the rigorous lower bound 'fL' on the primal objective value.
  %
  % b) If there is a cone violation 'dl(dl < 0)' and there exist a finite
  %    upper bound `xu(dl < 0)` on 'x' for each violated constraint, then
  %    a correction term has to be added to 'fL', which in case a) is zero.
  %
  % The correction term `defect` refers to the defect of the free variables.
  %
  idx = (dl < 0);
  if (~any (idx) || ~any (isinf (xbnd(idx))))
    fL = inf_ (obj.b' * y + dl(dl < 0)' * xbnd(dl < 0));
    tidy_up (obj, old_solution);
    return;  % SUCCESS
  end
  
  % Step 4: Perturb midpoint problem, such that
  %
  %    P(epsilon) = (mid([A]), mid([b]), mid([c]) + c_epsilon)
  %
  k      (idx) = k      (idx) + 1;
  epsilon(idx) = epsilon(idx) - (obj.options.ALPHA .^ k (idx)) .* dl(idx);
  if (~all (isfinite (epsilon)))
    tidy_up (obj, old_solution);
    error ('VSDP:rigorous_lower_bound:infinitePertubation', ...
      'rigorous_lower_bound: Perturbation exceeds finite range.');
  end
  c_epsilon(vidx) = epsilon((obj.K.f+1):end);
  
  % Index vector for perturbation.  In case of semidefinite programs, only
  % the diagonal elements have to be perturbed.  Those are easily obtained:
  vidx = vsdp.sindex (obj.K);
  vidx = vidx(:,1); % Get only diagonal entries of SDP cones.
  if (obj.K.l > 0)
    % In case of linear cones, each variable has to be pertubed.
    idx = obj.K.idx.l;
    vidx(idx(1):idx(end)) = true;
  end
  if (~isempty (obj.K.idx.q))
    % In case of second-order cones, only the first element is pertubed.
    vidx(obj.K.idx.q(:,1)) = true;
  end
  
  
  % Display short pertubation statistic.
  if (obj.options.VERBOSE_OUTPUT)
    fprintf ('\n\n');
    fprintf ('--------------------------------------------------\n');
    fprintf ('  VSDP.RIGOROUS_LOWER_BOUND  (iteration %d)\n', iter);
    fprintf ('--------------------------------------------------\n');
    fprintf ('  Violated cones (dl < 0): %d\n',          sum (idx));
    fprintf ('  Max. cone violation:     %+.2e\n',       min (dl_));
    fprintf ('  Max. pertubation:        %+.2e\n\n',     max (epsilon));
    fprintf ('  Solve pertubed problem using ''%s''.\n', obj.options.SOLVER);
    fprintf ('--------------------------------------------------\n\n');
  end
  
  % Step 5: Solve perturbed problem
  obj.pertubation.c = c_epsilon;
  obj.solve (obj.options.SOLVER);
  if ((obj.solution.info ~= 0) || isempty (y) || any (isnan (y)) ...
      || any (isinf (y)))
    tidy_up (obj, old_solution);
    error ('VSDP:rigorous_lower_bound:unsolveablePertubation', ...
      ['rigorous_lower_bound: Conic solver could not find a solution ', ...
      'for perturbed problem']);
  end
end

if (iter == VSDP_OPTIONS.ITER_MAX)
  error ('VSDP:VSDPLOW', 'VSDPLOW: maximum number of iterations reached');
end
obj.solution = old_solution;
fL = -Inf;
y = NaN;
dl = NaN;

end

function tidy_up (obj, old_solution)
% TIDY_UP  Revert changes to the VSDP object.
obj.pertubation.c = [];
obj.solution = old_solution;
end