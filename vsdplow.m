function [fL,y,dl,info] = vsdplow(A,b,c,K,x0,y,z0,xu,opts)
% VSDPLOW  Verified lower bound for semidefinite-quadratic-linear programming.
%
%   [fL,y,dl,info] = VSDPLOW(A,b,c,K,[],y0) Computes a verified lower bound of
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
%                  z = c-A'*y.
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

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

% check input
narginchk(6,9);
if isempty(A) || isempty(b) || isempty(c) || isempty(K) || isempty(y)
  error('VSDP:VSDPLOW', 'non-empty input arguments are required');
end
if ((nargin < 7) || any(isnan(z0)))
  z0 = [];
end
if (nargin < 8)
  xu = [];
end
if (nargin < 9)
  opts = [];
end

VSDP_OPTIONS = vsdpinit(opts);

[A,b,c,K,x0,y,z0] = import_vsdp(A,b,c,K,x0,y,z0);

% Check if approximations are applicable
if (any(isnan(y)))
  error('VSDP:VSDPLOW', 'approximate solution y contains NaN');
end
if (any(isnan(x0)))
  x0 = [];
end

% Get problem dimensions
dim3 = length(c);
nc = K.l + length(K.q) + length(K.s);  % number of cone constraints

% Check primal upper bound
xu = xu(:);
if isempty(xu)
  xu = inf(K.f + nc,1);
elseif (length(xu) ~= (K.f + nc))
  error('VSDP:VSDPLOW', 'upper bound vector has wrong dimension');
end

dl = -inf(nc,1);    % dual lower bounds / trace bounds
epsj = ones(nc,1);  % factor for perturbation
ceps = sparse(dim3,1);        % perturbation for c
pertS = cell(length(K.s),1);  % diagonal perturbations of semidefinite blocks

% Extract free part
xuf = xu(1:K.f);
xu(1:K.f) = [];

% Index vector for perturbation entries
pertI = ones(sum(K.s),1);
pertI(cumsum(K.s(1:end-1))+1) = 1 - K.s(1:end-1);
pertI = [ones(K.l+(~isempty(K.q)),1); K.q(1:end-1); cumsum(pertI)];
pertI(1) = K.f + 1;
pertI = cumsum(pertI);

% save rounding mode
rnd = getround();
setround(0);

% Algorithm with finite or infinite primal bounds xu
I = []; % index vector inside loop
info.iter = 0;
err_msg = '';
while (info.iter <= VSDP_OPTIONS.ITER_MAX)
  info.iter = info.iter + 1;
  setround(1);  % default for rigorous computation in steps 1-3
  
  % Step 1: Defect computation, free variables handling
  if ((K.f > 0) && (isinf(max(xuf))))
    % Solve dual linear constraints rigorously
    [y,I] = vuls([], [], A(1:K.f,:), c(1:K.f), [], [], y, I);
    if (~isa(y, 'intval'))
      err_msg = 'could not find solution of dual equations';
      break;
    end
  end
  
  % Compute rigorous enclosure for z = c - A*y  (with free variables)
  [z,zrad] = spdotK(mid(c),1,A,-mid(y),2);  % for point matrices
  zrad = zrad + rad(c);  % regard radii of other parameters
  if any(rad(y))
    zrad = zrad + abs(A)*rad(y);
  end
  if ~isempty(find(rad(A),1))
    zrad = zrad + rad(A)*(abs(mid(y))+rad(y));
  end
  
  defect = 0;  % defect by free variables
  if ((K.f > 0) && (isfinite(max(xuf))))  % given finite upper bounds
    defect = xuf' * (abs(z(1:K.f)) + zrad(1:K.f));
  end
  
  % Step 2: Verified lower bounds on cone eigenvalues
  
  % Bound for linear variables
  if (K.l > 0)
    ind = K.f+1:K.f+K.l;
    zjn = zrad(ind) - z(ind);  % -inf(zl)
    % interpret all lp vars as one cone:
    %  if any element is negative all constraints close to zero will be
    %  perturbated
    zj1 = max(zjn);
    if (zj1 > 0)  % there is a negative bound
      ind = zjn > -zj1;
      zj1 = max((-1e-13)*min(zjn),zj1);
      zjn(ind) = min(zjn(ind)+zj1,1.05*zj1);
    end
    dl(1:K.l) = -zjn;
  end
  % Eigenvalue bound for second-order cone variables
  ind = K.f + K.l;
  for j = 1:length(K.q)
    ind = ind(end)+2:ind(end)+K.q(j);
    zj1 = zrad(ind(1)-1) - z(ind(1)-1); % -inf(zq(1))
    zjn = abs(z(ind)) + zrad(ind);      % sup(abs(zq(2:end)))
    zjn = sqrtsup(zjn'*zjn);  % sup(||zq(2:end)||)
    dl(K.l+j) = -(zj1+zjn);   % inf(zq(1)-||zq(2:end)||)
  end
  % Eigenvalue bound for semidefinite cone
  blke = K.f + K.l + sum(K.q) + sum(K.s.*(K.s+1))/2;
  for j = length(K.s):-1:1
    ofs = K.l + length(K.q) + j;
    blks = blke - K.s(j)*(K.s(j)+1)/2 + 1;
    [lmin,dl(ofs),pertS{j}] = bnd4sd(midrad(z(blks:blke), zrad(blks:blke)), ...
      VSDP_OPTIONS.FULL_EIGS_ENCLOSURE);
    if (lmin > 0)
      dl(ofs) = lmin;
    end
    pertS{j} = epsj(ofs) * pertS{j}(:);
    blke = blks - 1;
  end
  
  % Step 3: Cone feasibility check, computing lower bound
  dli = find(dl < 0);
  if ~any(isinf(xu(dli)))
    % inf(min(dl,0)*xu + b'*y - defect)
    fL = -(sum(dl(dli)'*(-xu(dli))) + prodsup(-mid(b)',mid(y),rad(b)',rad(y)) + defect);
    setround(rnd);  % reset rounding mode
    return; % SUCCESS
  end
  
  % Step 4: Perturb problem
  setround(0);  % no code for rigorous computations
  ind = 1:K.l+length(K.q);
  if isempty(ind)
    ceps = ceps + sparse(pertI,1,cat(1,pertS{:}),dim3,1);
  else
    ceps = ceps + sparse(pertI,1,cat(1,epsj(ind).*min(dl(ind),0),pertS{:}), ...
      dim3,1);
  end
  if (~all(isfinite(ceps)))
    err_msg = 'perturbation extended range';
    break;
  end
  epsj(dli) = epsj(dli) * (1 + VSDP_OPTIONS.ALPHA); % update perturbation factor
  
  % Step 5: Solve perturbed problem
  [~,x0,y,z0,info_solver] = mysdps(A,b,c+ceps,K,x0,y,z0);
  if ((info_solver ~= 0) || isempty(y) || any(isnan(y)) || any(isinf(y)))
    err_msg = 'conic solver could not find solution for perturbed problem';
    break;
  end
end

setround(rnd); % reset rounding mode

if (info.iter == VSDP_OPTIONS.ITER_MAX)
  err_msg = 'maximum number of iterations reached';
end

if (VSDP_OPTIONS.VERBOSE_OUTPUT)
  disp(['VSDPLOW: ', err_msg]);
end

fL = -Inf;
y = NaN;
dl = NaN;

end
