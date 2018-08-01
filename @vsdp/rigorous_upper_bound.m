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
  if (length (ybnd) ~= length (b))
    error ('VSDP:rigorous_upper_bound:badYBND', ...
      ['rigorous_upper_bound: The length of the upper bound vector ', ...
      '''ybnd'' must be %d, but is %d.'], length (b), length (ybnd));
  end
else  % If dual upper bounds are not given, use infinite bounds.
  ybnd = inf(length (b), 1);
end

% If the problem was not approximately solved before, do it now.
if (isempty (obj.solutions('Approximate solution')))
  warning ('VSDP:rigorous_upper_bound:noApproximateSolution', ...
    ['rigorous_upper_bound: The conic problem has no approximate ', ...
    'solution yet, which is now computed using ''%s''.'], obj.options.SOLVER);
  obj.solve (obj.options.SOLVER, 'Approximate solution');
end
x = vsdp_indexable (intval (obj.solutions('Approximate solution').x), obj);

% If any upper bound 'ybnd' is not finite, proceed with other algorithm.
if (~all (isfinite (ybnd)))
  if (nargin == 2)
    warning ('VSDP:rigorous_upper_bound:infiniteBounds', ...
    ['rigorous_upper_bound: At least one element in the bounds ', ...
    '''ybnd'' was not finite.  Using algorithm for infinite bounds.']);
  end
  obj = rigorous_upper_bound_infinite_bounds (obj, x);
  return;
end
rub = tic;

% Step 1: Project approximate solution into the respective cones, e.g.
%         compute x^+.

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

function obj = rigorous_upper_bound_infinite_bounds (obj, x)
lb = -inf(nc,1);    % dual lower bounds
epsj = ones(nc,1);  % factor for perturbation
xeps = sparse(length(x),1);   % perturbation for x
pertS = cell(length(K.s),1);  % diagonal perturbations of semidefinite blocks

% Index vector for perturbation entries
pertI = ones(sum(K.s),1);
pertI(cumsum(K.s(1:end-1))+1) = 1-K.s(1:end-1);
pertI = [ones(K.l+(~isempty(K.q)),1); K.q(1:end-1); cumsum(pertI)];
pertI(1) = K.f + 1;
pertI = cumsum(pertI);

% Save rounding mode
rnd = getround();
setround(0);

info.iter = 0;
err_msg = '';
while (info.iter <= VSDP_OPTIONS.ITER_MAX)
  info.iter = info.iter + 1;
  setround(1);  % default rounding for verification part
  
  % Step 1: Compute rigorous enclosure for `A * x = b`.
  [x, I] = vuls([], [], A, b, [], [], x, I);
  if (~isintval(x))
    err_msg = 'could not find solution of primal equations';
    break;
  end
  
  % Step 2: Verified lower bounds for each cone.
  
  % LP cones
  if (K.l > 0)
    ind = (K.f + 1):(K.f + K.l);
    xl = -inf(x(ind));
    % Treat all linear variables as one cone: if any element is negative,
    % all constraints close to zero will be perturbated.
    xl_max = max(xl);
    if (xl_max > 0)  % there is a negative bound
      ind = (xl > (-xl_max / 2));
      xl(ind) = xl(ind) + xl_max;
    end
    lb(1:K.l) = -xl;
  end
  
  % Second-order cones
  ind = K.f + K.l;
  for j = 1:length(K.q)
    ind = (ind(end) + 2):(ind(end) + K.q(j));
    % lb = inf (z(1) - ||z(2:end)||)
    lb(K.l+j) = inf_ (z(ind(1)-1) - norm (intval (z(ind))));
  end
  
  % SDP cones
  blke = K.f + K.l + sum(K.q) + sum(K.s.*(K.s+1))/2;
  for j = length(K.s):-1:1
    ofs = K.l + length(K.q) + j; % FIXME: K.f?
    dind = cumsum(1:K.s(j));  % index for diagonal entries
    blk3 = dind(end);
    vx = x(blke-blk3+1:blke) / 2;
    vx(dind) = 2 * vx(dind);  % regard mu=2 for x
    [lmin,lb(ofs),pertS{j}] = bnd4sd(vx, VSDP_OPTIONS.FULL_EIGS_ENCLOSURE);
    if (lmin > 0)
      lb(ofs) = lmin;
    end
    pertS{j} = epsj(ofs) * pertS{j}(:);
    blke = blke - blk3;
  end
  
  % Step 3: Cone feasibility check and upper bound computation:
  %
  % If all lower bounds `lb` are non-negative, there are no cone violations
  % and the rigorous upper bound `fU` on the dual objective value can be
  % computed.
  %
  if (all (lb >= 0))
    fU = sup (c' * x);
    x = export_vsdp(imported_fmt,K,x);
    setround(rnd);  % reset rounding mode
    return; % SUCCESS
  end
  
  % Step 4: Perturb problem
  setround(0);  % no code for rigorous computations
  ind = 1:K.l+length(K.q);
  if isempty(ind)
    xeps = xeps + sparse(pertI,1,cat(1,pertS{:}),length(c),1);
  else
    xeps = xeps + sparse(pertI,1,cat(1,epsj(ind).*min(lb(ind),0),pertS{:}), ...
      length(c),1);
  end
  if (~all(isfinite(xeps)))
    err_msg = 'perturbation extended range';
    break;
  end
  epsj(lb < 0) = epsj(lb < 0) * (1 + VSDP_OPTIONS.ALPHA); % update perturbation factor
  
  % Step 5: Solve perturbed problem
  [~,x0,y0,z0,info_solver] = mysdps(A,b+(xeps'*A)',c,K,x0,y0,z0);
  if ((info_solver ~= 0) || isempty(x0) || any(isnan(x0)) || any(isinf(x0)))
    err_msg = 'conic solver could not find solution for perturbed problem';
    break;
  end
  x = x0 - xeps;  % undo perturbation
end

setround(rnd); % reset rounding mode

if (info.iter == VSDP_OPTIONS.ITER_MAX)
  err_msg = 'maximum number of iterations reached';
end

if (VSDP_OPTIONS.VERBOSE_OUTPUT)
  disp(['VSDPUP: ', err_msg]);
end

fU = Inf;
x = NaN;
lb = NaN;
end
