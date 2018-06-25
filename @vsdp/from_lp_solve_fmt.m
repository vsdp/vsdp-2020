function [obj, pd] = from_lp_solve_fmt (A, b, c, e, lb, ub)
% FROM_LP_SOLVE_FMT  Import LP problem data from LP_SOLVE format.
%
%   [obj, pd] = vsdp.FROM_LP_SOLVE_FMT (A, b, c, e, lb, ub);
%
%   The LP_SOLVE problem format is:
%
%         max  c'*x
%         s.t. A*x <e> b,
%              lb <= x <= ub.
%
%   The problem data is:
%
%       'A'  double(m,n)  linear constraint matrix
%       'b'  double(m,1)  right-hand side vector
%       'c'  double(n,1)  primal objective function vector
%       'e'  double(m,1)  sense of the inequalities:
%            -1 - less or equal than
%             0 - equality
%             1 - greater or equal than
%       'lb' double(n,1)  vector of lower bounds
%       'ub' double(n,1)  vector of upper bounds
%
%   For more information on the LP_SOLVE format, see:
%
%     http://lpsolve.sourceforge.net/5.5/MATLAB.htm
%
%   Note, that if the problem gets dualized, e.g. if `pd == 'd'`, the optimal
%   value must be negated to get the right sign.
%
%   See also from_mps_file.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk(4, 6);
b = b(:);
c = c(:);

[m, n] = size(A);
if ((m ~= length(b)) || (n ~= length(c)))
  error ('VSDP:FROM_LP_SOLVE_FMT', ...
    'from_lp_solve_fmt: c or b not compatible with A');
elseif m~=length(e)
  error ('VSDP:FROM_LP_SOLVE_FMT', ...
    'from_lp_solve_fmt: dimension of flag vector "e" does not match');
end

if (nargin < 6)
  ub = inf(n,1);
end
if (nargin < 5)
  lb = -inf(n,1);
end

% resolve fixed bounds
idx = find(lb==ub);  % index for fixed bounds
nfx = length(idx);  % number of fixed bounds
if nfx>0
  A = [A; sparse(1:nfx,idx,1,nfx,n)];
  e = [e; zeros(nfx,1)];
  b = [b; lb(idx)];
  lb(idx) = -inf;
  ub(idx) = inf;
  m = m + nfx;
end

% indices
indl = find(e<0);  % inequalities less than  (lower)
indu = find(e>0);  % inequalities greater than  (upper)

Ai = [A(indl,:); -A(indu,:)];  % part of coefficient matrix for inequalities
A([indl; indu],:) = [];  % part of coefficient matrix for equality constraints

bi = [b(indl); -b(indu)];
b([indl; indu]) = [];
mi = length(bi);

indl = find(lb>-inf);  % index for (applied) lower bounds
nl = length(indl);  % number of lower bounds

if all(lb(indl)==0) && min(ub)==inf && n+mi<m+nl
  % primal form with free + slack variables
  idx = [find(lb==-inf); indl];
  A = [A(:,idx)' Ai(:,idx)'; ...
    sparse(mi,m-mi) speye(mi)];  % free vars first + add slack vars
  b = [b; bi];
  c = [-c(idx); sparse(mi,1)];
  K.f = n - nl;
  K.l = nl + mi;
  pd = 'd';
else
  % extend dual form
  indu = find(ub<inf);
  nu = length(indu);
  A = [A; Ai; sparse(1:nl,indl,-1,nl,n); sparse(1:nu,indu,1,nu,n)];
  Ai = c;  % Ai used as place-holder
  c = [b; bi; -lb(indl); ub(indu)];
  b = Ai;
  K.f = m - mi;
  K.l = mi + nl + nu;
  pd = 'p';
end

% Recursive call to default VSDP constructor.
obj = vsdp (A, b, c, K);

end
