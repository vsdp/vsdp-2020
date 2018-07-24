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
if isempty(xbnd)
  xbnd = inf (size (obj.K.blk, 1), 1);
else
  xbnd = xbnd(:);
  if (length (xbnd) ~= size (obj.K.blk, 1))
    error ('VSDP:rigorous_lower_bound:badXBND', ...
      ['rigorous_lower_bound: The length of the upper bound vector ', ...
      '''xbnd'' must be %d, but is %d.'], size (obj.K.blk, 1), length (xbnd));
  end
end


nc = K.l + length(K.q) + length(K.s);  % number of cone constraints
dl = -inf(nc,1);    % dual lower bounds / trace bounds
epsj = ones(nc,1);  % factor for perturbation
ceps = sparse(obj.n,1);        % perturbation for c
pertS = cell(length(K.s),1);  % diagonal perturbations of semidefinite blocks

% Index vector for perturbation entries
vidx = vsdp.sindex (obj.K);
pertI = find (vidx(:,1));

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
  %   a) LP cones
  %
  idx = obj.K.idx.l;
  dl.l = inf_ (d(idx));
  % If any LP lower bound is negative: perturb all constraints close to
  % zero inside the LP cone using a heuristic.
  if (any (dl.l) < 0)
    pval = [min(dl.l), max(dl.l) * (-1e-13)];
    pval = abs (min (pval));
    % Take the minimum of absolute and relative pertubation.
    dl.l(dl.l < pval) = min (dl.l(dl.l < pval) + pval, 1.05 * pval);
  end
  
  %
  %   b) Second-order cones
  %
  dl.q = cell (obj.K.q, 1);
  for j = 1:length(dl.q)
    idx = obj.K.idx.q(j,:);
    % dl = inf (z(1) - ||z(2:end)||)
    dl.q{j} = inf_ (d(idx(1)) - norm (intval (d((idx(1)+1):idx(end)))));
  end
  
  %
  %   c) SDP cones
  %
  dl.s = cell (obj.K.s, 1);
  for j = 1:length(dl.s)
    idx = obj.K.idx.s(j,:);
    dl.s{j} = veigsym (d(idx(1):idx(end)));
    if (min (dl.s{j}) < 0)
      pertS{j} = epsj(ofs) * pertS{j}(:);
    end
  end
  
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
  if (all (dl >= 0) || ~any (isinf (xbnd(dl < 0))))
    fL = inf_ (b' * y + dl(dl < 0)' * xbnd(dl < 0) - xbnd(1:K.f)' * (mag (d(1:K.f))));
    return;         % SUCCESS
  end
  
  % Step 4: Perturb problem
  idx = 1:K.l+length(K.q);
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
  epsj(dl < 0) = epsj(dl < 0) * (1 + obj.options.ALPHA);
  
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
