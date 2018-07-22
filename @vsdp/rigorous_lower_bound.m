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
I = []; % index vector inside loop
info.iter = 0;
while (info.iter <= obj.options.ITER_MAX)
  info.iter = info.iter + 1;
  setround(1);  % default for rigorous computation in steps 1-3
  
  % Step 1: Defect computation, free variables handling
  
  % If infinite upper bounds for free variables are given.
  if ((K.f > 0) && (isinf (max (xbnd(1:K.f)))))
    % Solve dual linear constraints rigorously
    [y0,I] = vuls([], [], At(1:K.f,:), c(1:K.f), [], [], y0, I);
    if (~isintval (y0))
      warning ('VSDP:VSDPLOW', ...
        'VSDPLOW: could not find solution of dual equations');
      break;
    end
  end
  
  % Compute rigorous enclosure for `z = c - A' * y`  (with free variables).
  z = c - At * y0;
  
  % Compute defect by free variables, if finite upper bounds for free variables
  % are given.
  defect = 0;
  if ((K.f > 0) && (isfinite (max (xbnd(1:K.f)))))
    defect = xbnd(1:K.f)' * (mag (z(1:K.f)));
  end
  
  % Step 2: Verified lower bounds on cone eigenvalues
  
  % LP cones
  if (K.l > 0)
    idx = K.f+1:K.f+K.l;
    zl = inf_ (z(idx));
    % If any element of the LP cone is negative all constraints close to zero
    % will be perturbed.
    if (any (zl) < 0)
      pert = min ((-1e-13) * max (zl), min (zl));  % TODO: reference?
      idx = zl > -min (zl);
      zl(idx) = min (zl(idx) + pert, 1.05 * pert);
    end
    dl(1:K.l) = zl;
  end
  
  % Second-order cones
  for j = 1:length(obj.K.q)
    idx = obj.K.idx.q(j,:);
    % dl = inf (z(1) - ||z(2:end)||)
    dl(K.l+j) = inf_ (z(idx(1)) - norm (intval (z((idx(1)+1):idx(end)))));
  end
  
  % SDP cones
    for j = 1:length(obj.K.s)
    idx = obj.K.idx.s(j,:);
    [lmin,dl(ofs),pertS{j}] = veigsym (z(idx(1):idx(end)); % TODO
    if (lmin > 0)
      dl(ofs) = lmin;
    end
    pertS{j} = epsj(ofs) * pertS{j}(:);
  end
  
  % Step 3: Cone feasibility check and lower bound computation:
  %
  % a) If the defect  d = c - A' * y  lies in each cone, then all lower bounds
  %    'dl' on 'd' are non-negative, `y` is dual feasible, and `inf_ (b' * y)`
  %    is the rigorous lower bound `fL` on the primal objective value.
  %
  % b) If there is a cone violation `dl(dl < 0)` and there exist a finite
  %    upper bound `xu(dl < 0)` on `x` for each violated constraint, then
  %    a correction term has to be added to `fL`, which in case a) is zero.
  %
  % The correction term `defect` refers to the defect of the free variables.
  %
  if (all (dl >= 0) || ~any (isinf (xbnd(dl < 0))))
    fL = inf_ (b' * y0 + dl(dl < 0)' * xbnd(dl < 0) - defect);
    y = y0;
    setround(rnd);  % reset rounding mode
    return;         % SUCCESS
  end
  
  % Step 4: Perturb problem
  setround(0);  % no code for rigorous computations
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
  [~,x0,y0,z0,info_solver] = mysdps(At,b,c+ceps,K,x0,y0,z0);
  if ((info_solver ~= 0) || isempty(y0) || any(isnan(y0)) || any(isinf(y0)))
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
