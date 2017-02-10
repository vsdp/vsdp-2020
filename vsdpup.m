function [fU,x,lb,info] = vsdpup(A,b,c,K,x,y0,z0,yu,opts)
% VSDPUP  Verified upper bound for semidefinite-quadratic-linear programming.
%
%   [fU,x,lb,info] = VSDPUP(A,b,c,K,x0) Computes a verified upper bound of the
%      primal optimal value and a rigorous enclosure of dual strict feasible
%      (near optimal) solutions of a conic problem in the standard primal-dual
%      form.  This form and the block-diagonal format (A,b,c,K) is explained in
%      'mysdps.m'.
%
%         'x0'     A primal feasible (eps-optimal) solution of the same
%                  dimension as input c.  This solution can be computed using
%                  'mysdps'.
%
%      The output is:
%
%         'fU'     Verified upper bound of the dual optimal value.
%
%         'x'      Rigorous enclosure of primal strict feasible solutions.
%
%         'lb'     Verified lower bounds of the eigenvalues or spectral values
%                  of x with respect to K.
%
%         'info'   Struct containing further information.
%           - iter  The number of iterations.
%
%   VSDPUP(A,b,c,K,x0,y0,z0) optionally provide the other approximate
%      solutions of 'mysdps' (y0 and z0).
%
%   VSDPUP(A,b,c,K,x0,[],[],yu) optionally provide known finite upper bounds
%      of the dual optimal solution y.  The following dual boundedness
%      assumption is assumed: an optimal dual solution y satisfies
%
%         -yu(i) <= y(i) <= yu(i),  i = 1:m.
%
%      We recommend to use infinite bounds yu(i) instead of unreasonable large
%      bounds yu(i).  This improves the quality of the lower bound in many 
%      cases, but may increase the computational time.
%
%   VSDPUP(A,b,c,K,x0,[],[],[],opts) optionally provide a structure for
%      additional parameter settings, explained in vsdpinit.
%
%   See also mysdps, vsdpinit, vsdplow, vsdpinfeas.

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

% check input
narginchk(5,9);
if isempty(A) || isempty(b) || isempty(c) || isempty(K) || isempty(x)
  error('VSDP:VSDPUP','non-empty input arguments are required');
end
if ((nargin < 6) || any(isnan(y0)))
  y0 = [];
end
if ((nargin < 7) || any(isnan(z0)))
  z0 = [];
end
if (nargin < 8)
  yu = [];
end
if (nargin < 9)
  opts = [];
end

VSDP_OPTIONS = vsdpinit(opts);

[A,Arad,b,brad,c,crad,K,x,y0,z0,IF] = import_vsdp(A,b,c,K,x,y0,z0);

% Check if approximation are applicable
if (any(isnan(x)))
  error('VSDP:VSDPUP', 'approximate solution x contains NaN');
end

% Get problem dimensions
dim3 = length(c);
nc = K.l + length(K.q) + length(K.s);  % number of cone constraints

% Check dual upper bound
yu = yu(:);
if (isempty(yu))
  yu = inf(length(b),1);
elseif (length(yu) ~= length(b))
  error('VSDP:VSDPUP', 'upper bound vector has wrong dimension');
end

% Save rounding mode
rnd = getround();
setround(0);

% Algorithm with finite dual bounds yu, projection into cone
if (isfinite(max(yu)))
  setround(1);  % default rounding mode for verification code
  x = full(x);  % faster for verification code

  % Projection of linear part
  x(K.f+1:K.f+K.l) = max(x(K.f+1:K.f+K.l),0);

  % Projection of second-order cone part
  blke = K.f + K.l;
  for j = 1:length(K.q)
    xj1 = x(blke+1);  % xq(1)
    xjn = x(blke+2:blke+K.q(j));  % xq(2:end)
    xjn = sqrtsup(xjn'*xjn);      % sup(||xq(2:end)||)
    x(blke+1) = max(xj1,xjn);     % very simple projection
    blke = blke + K.q(j);
  end

  % Projection of semidefinite part, force to psd cone
  blke = K.f + K.l + sum(K.q);
  for j = 1:length(K.s)
    dind = cumsum(1:K.s(j));  % index for diagonal entries
    vx = x(blke+1:blke+dind(end));
    vx(dind) = 2*vx(dind);    % regard mu=2 for x
    lambdaj = bnd4sd(vx);
    if (lambdaj < 0)  % simple projection into feasible space
      x(dind+blke) = x(dind+blke) + (-lambdaj)/2;
    end
    blke = blke + dind(end);
  end

  % x <in> K:  now regard defect = |A*x-b|
  defect = resmag(x',A,b',1,0,Arad,brad',0);
  % fU = sup(x'*c + defect*yu)
  fU = prodsup(x',c,0,crad) + defect*yu;
  x = NaN;
  lb = NaN;
  info.iter = 1;
  
  setround(rnd); % reset rounding mode
  return;
end

% Algorithm with infinite dual bounds yu
x0 = x;  % starting point
I = [];  % index vector inside loop
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

info.iter = 0;
err_msg = '';
while (info.iter <= VSDP_OPTIONS.ITER_MAX)
  info.iter = info.iter + 1;
  setround(1);  % default rounding for verification part

  % Step 1: Compute rigorous enclosure for A*x = b
  [x, I] = vuls([],[],struct('mid',A,'rad',Arad),...
    struct('mid',b,'rad',brad),[],[],x,I);
  if ~isstruct(x)
    err_msg = 'could not find solution of primal equations';
    break;
  else
    xrad = sparse(x.rad);
    x = full(x.mid);
  end

  % Step 2: Verified lower bounds on cone eigenvalues
  
  % Bound for linear variables
  if (K.l > 0)
    ind = (K.f + 1):(K.f + K.l);
    xl = xrad(ind) - x(ind);  % -inf(xl)
    % Treat all linear variables as one cone: if any element is negative,
    % all constraints close to zero will be perturbated.
    xl_max = max(xl);
    if (xl_max > 0)  % there is a negative bound
      ind = (xl > (-xl_max / 2));
      xl(ind) = xl(ind) + xl_max;
    end
    lb(1:K.l) = -xl;
  end
  % Eigenvalue bound for second-order cone variables
  ind = K.f + K.l;
  for j = 1:length(K.q)
    ind = (ind(end) + 2):(ind(end) + K.q(j));
    xj1 = xrad(ind(1)-1) - x(ind(1)-1); % -inf(xq(1))
    xjn = abs(x(ind)) + xrad(ind);      % sup(abs(xq(2:end)))
    xjn = sqrtsup(xjn'*xjn);  % sup(||xq(2:end)||)
    lb(K.l+j) = -(xj1+xjn);   % inf(xq(1) - ||xq(2:end)||)
  end
  % Eigenvalue bound for semidefinite cone
  blke = K.f + K.l + sum(K.q) + sum(K.s.*(K.s+1))/2;
  for j = length(K.s):-1:1
    ofs = K.l + length(K.q) + j; % FIXME: K.f?
    dind = cumsum(1:K.s(j));  % index for diagonal entries
    blk3 = dind(end);
    vx = struct('mid',x(blke-blk3+1:blke) / 2,'rad',xrad(blke-blk3+1:blke) / 2);
    vx.mid(dind) = 2 * vx.mid(dind);  % regard mu=2 for x
    vx.rad = vx.rad + vx.rad.*sparse(dind,1,1,blk3,1);
    [lmin,lb(ofs),pertS{j}] = bnd4sd(vx,1,VSDP_OPTIONS.FULL_EIGS_ENCLOSURE);
    if (lmin > 0)
      lb(ofs) = lmin;
    end
    pertS{j} = epsj(ofs) * pertS{j}(:);
    blke = blke - blk3;
  end

  % Step 3: Cone feasibility check, computing upper bound
  lbi = (lb < 0);
  if (~any(lbi))
    fU = prodsup(x',c,xrad',crad);  % sup(x'*c)
    x = export_vsdp(IF,K,midrad(full(x),full(xrad)));
    setround(rnd);  % reset rounding mode
    return; % SUCCESS
  end

  % Step 4: Perturb problem
  setround(0);  % no code for rigorous computations
  ind = 1:K.l+length(K.q);
  if isempty(ind)
    xeps = xeps + sparse(pertI,1,cat(1,pertS{:}),dim3,1);
  else
    xeps = xeps + sparse(pertI,1,cat(1,epsj(ind).*min(lb(ind),0),pertS{:}), ...
      dim3,1);
  end
  if (~all(isfinite(xeps)))
    err_msg = 'perturbation extended range';
    break;
  end
  epsj(lbi) = epsj(lbi) * (1 + VSDP_OPTIONS.ALPHA); % update perturbation factor

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
