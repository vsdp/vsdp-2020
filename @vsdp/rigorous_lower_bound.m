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
nbnd = obj.K.l + length(obj.K.q) + length(obj.K.s);  % Number of bounds.
if ((nargin < 2) || isempty (xbnd))
  xbnd = inf (nbnd, 1);
else
  xbnd = xbnd(:);
  if (length (xbnd) ~= nbnd)
    error ('VSDP:rigorous_lower_bound:badXBND', ...
      ['rigorous_lower_bound: The length of the upper bound vector ', ...
      '''xbnd'' must be %d, but is %d.'], size (obj.K.blk, 1), length (xbnd));
  end
end

epsj = ones(nbnd,1);  % factor for perturbation
ce   = sparse(obj.n,1);        % perturbation for c
pertS = cell(length(obj.K.s),1);  % diagonal perturbations of semidefinite blocks

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

% Save rounding mode.
rnd = getround();
setround(0);

% Algorithm for both finite and infinite upper bounds `xu`.
info.iter = 0;
while (info.iter <= obj.options.ITER_MAX)
  info.iter = info.iter + 1;
  
  % If infinite upper bounds for free variables are given.  Ensure, that
  % the approximate dual solution 'y' solves the free variable part, see
  %
  %   https://vsdp.github.io/references.html#Anjos2007
  %
  % for details.
  if ((obj.K.f > 0) && (isinf (max (xbnd(1:obj.K.f)))))
    % Compute rigorous enclosure for underdetermined linear interval system of
    % dimension (K.f x m):
    %
    %    At.f * y = c.f
    %
    y = verifylss (obj.At(1:obj.K.f,:), obj.c(1:obj.K.f));
    if (~isintval (y) || any (isnan (y)))
      error ('VSDP:rigorous_lower_bound:noBoundsForFreeVariables', ...
        ['rigorous_lower_bound: Could not find a verified solution of the ', ...
        'linear  system of free variables.']);
    end
  end
  
  % Step 1: Compute rigorous enclosure [d] for  c - At*y.
  d = obj.c - obj.At * intval (y);
  
  % Step 2: Verified lower bounds on cone eigenvalues
  %
  % Free variables
  dl.f = mag (d(1:obj.K.f));
  
  % LP cones
  idx = obj.K.idx.l;
  dl.l = inf_ (d(idx(1):idx(end)));
  
  % Second-order cones
  idx = obj.K.idx.q;
  % dl = inf (d(1) - ||d(2:end)||)
  dl.q = inf_ (d(idx(:,1)) - norm (d((idx(:,1)+1):idx(:,end))));
  
  % SDP cones
  dl.s = cell (obj.K.s, 1);
  for j = 1:length(dl.s)
    idx = obj.K.idx.s(j,:);
    dl.s{j} = veigsym (d(idx(1):idx(end)));
    if (min (dl.s{j}) < 0)
      pertS{j} = epsj(ofs) * pertS{j}(:);
    end
  end
  
  % Concatenate all lower bounds to a single vector.
  dl_ = [dl.f; dl.l; dl.q; dl.s];
  
  % Step 3: Cone feasibility check and lower bound computation:
  %
  % a) If the defect  d = c - A' * y  lies in each cone, then all lower bounds
  %    'dl' on 'd' are non-negative, 'y' is dual feasible, and  inf_ (b' * y)
  %    is the rigorous lower bound 'fL' on the primal objective value.
  %
  % b) If there is a cone violation `dl(dl < 0)` and there exist a finite
  %    upper bound `xu(dl < 0)` on `x` for each violated constraint, then
  %    a correction term has to be added to `fL`, which in case a) is zero.
  %
  % The correction term `defect` refers to the defect of the free variables.
  %
  if (all (dl_ >= 0) || ~any (isinf (xbnd(dl_ < 0))))
    fL = inf_ (b' * y + dl_(dl_ < 0)' * xbnd(dl_ < 0));
    return;  % SUCCESS
  end
  
  % Step 4: Perturb midpoint problem.
  %
  % The perturbed midpoint problem has the form
  %
  %    P(e) = (mid([A]), mid([b]), mid([c]) + ce)
  idx = 1:K.l+length(K.q);
  
  epsj(dl_ < 0) = epsj(dl_ < 0) * (1 + obj.options.ALPHA);
  
  if isempty(idx)
    ceps = ceps + sparse(pertI,1,cat(1,pertS{:}),obj.n,1);
  else
    ceps = ceps + sparse(pertI,1,cat(1,epsj(idx).*min(dl(idx),0),pertS{:}), ...
      obj.n,1);
  end
  if (~all(isfinite(ceps)))
    warning ('VSDP:VSDPLOW', 'VSDPLOW: perturbation extended range');
    break;
  end
  % Update perturbation factors.
  
  
  % Step 5: Solve perturbed problem
  [~,x0,y,z0,info_solver] = mysdps(At,b,c+ceps,K,x0,y,z0);
  if ((info_solver ~= 0) || isempty(y) || any(isnan(y)) || any(isinf(y)))
    warning ('VSDP:VSDPLOW', ...
      'VSDPLOW: conic solver could not find solution for perturbed problem');
    break;
  end
end

setround(rnd); % reset rounding mode

if (info.iter == VSDP_OPTIONS.ITER_MAX)
  error ('VSDP:VSDPLOW', 'VSDPLOW: maximum number of iterations reached');
end

fL = -Inf;
y = NaN;
dl = NaN;

end
