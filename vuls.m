function [X, I, J] = vuls(A,a,B,b,xl,xu,x0,I,opts)
% VULS  Verification for underdetermined linear systems.
%
%         A * x <= a,
%         B * x == b,     or if B is transposed:  x' * B == b'
%         xl <= x <= xu.
%
%   The input A (m*n matrix), a (m-vector), B (p*n or n*p matrix), and b
%   (p-vector) can be real or interval quantities; the simple bounds xl and xu
%   must be real, but may be infinite; the approximate solution x0 must be real.
%   The optional input vector I must contain p indices out of {1,...,n} such
%   that the submatrix B(:,I) is nonsingular; if opts.VERIFY_FULL_LSS is true
%   the full non-symmetric lss enclosure algorithm will be applied.
%
%   The output is:
%
%      X   a box (n-interval vector), containing for every real
%          input (A,a,B,b)  within the interval input data a  solution
%          x of the above system, provided J is empty. Especially,
%          existence of solutions is verified, and moreover
%          X is computed close to x0 in a specified manner; for details see
%          [Jansson2004].
%
%      I   index vector such that the p*p submatrix B(:,I) (or in the
%          transposed case: B(I,:)) is nonsingular.
%
%          X := nan(n,1), I = [], J = [];
%          if existence of solutions cannot be proved, and verified finite
%          bounds cannot be computed. This is the case, if
%          (i) B has no full rank, or (ii) the linear interval solver
%          VERIFYLSS cannot compute rigorous bounds, or (iii) the box
%          of simple bounds [xl,xu] has no appropriate interior. Otherwise,
%          J returns the vector of indices of violated inequalities
%          for X:
%             J.ineqlin:  violated row indices of A * X <= b,
%             J.lower: violated indices of  xl <= X ,
%             J.upper: violated indices of  X <= xu .
%
%   See also vsdpinfeas, vsdplow, vsdpup.

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

% check parameter
if nargin<7 || isempty(x0) || ~isnumeric(x0)
  error('VULS: not enough input parameter or no initial solution set');
end
x0 = x0(:);
if ~isnumeric(a)
  error('VULS: inappropriate a');
end
a = a(:);
isIntval = 1;
if isstruct(A) || isstruct(B) || isstruct(b)
  isIntval = 0;
end
% read interval b
if isnumeric(b)
  b = b(:);  brad = sparse(size(b,1),1);
elseif isa(b,'intval') || all(isfield(b,{'mid','rad'}))
  brad = b.rad(:);  b = b.mid(:);
else
  error('VULS: unexpected data format of b');
end
% dimensions
m = length(a);
n = length(x0);
p = length(b);
% read interval A
if isnumeric(A) && ~isempty(A)
  Arad = sparse(size(A,1),size(A,2));
elseif isa(A,'intval') || all(isfield(A,{'mid','rad'}))
  Arad = A.rad;  A = A.mid;
elseif ~isempty(A)
  error('VULS: unexpected data format of A');
end
if all(size(A)==[m n])
  A = A';  Arad = Arad';
elseif ~isempty(A) && (any(size(A)~=[n m]) || any(size(Arad)~=[n m]))
  error('VULS: dimension of A does not match to a and x0');
end
% read interval B
if isnumeric(B) && ~isempty(B)
  Brad = sparse(max(size(B)),min(size(B)));
elseif isa(B,'intval') || all(isfield(B,{'mid','rad'}))
  Brad = B.rad;  B = B.mid;
else
  error('VULS: cannot import data format of B');
end
if all(size(B)==[p n])
  B = B';
end
if all(size(Brad)==[p n])
  Brad = Brad';
elseif ~isempty(B) && (any(size(B)~=[n p]) || any(size(Brad)~=[n p]))
  error('VULS: dimension of B does not match to b and x0');
end
% read bounds for x0
xl = xl(:);
if length(xl)~=n
  xl = -inf(n,1);
end
xu = xu(:);
if length(xu)~=n
  xu = inf(n,1);
end
% check index and opts parameter
if nargin<8 || length(I)~=p
  I = [];
end
if nargin<9
  opts = [];
end


VSDP_OPTIONS = vsdpinit(opts);

% preparation

% initial output
X = nan(n,1);
J.ineqlin = [];
J.lower = [];
J.upper = [];

% projection of x0 into the interior of [xl,xu]
if any(xu<xl)  % appropriate bounds ?
  disp('VULS: simple bounds are not feasible');
  return;
end
x0 = max(xl,x0);
x0 = min(xu,x0);

% enclosure X
if ~isempty(B)
  % determine the basis with lu decomposition
  if isempty(I)
    % bias towards greater x
    [~,I] = sort(abs(x0),'descend');
    if issparse(B)  % use UMFPACK
      [~,~,rp,~] = lu(B(I,:),[0.95 0.75],'vector');
    else
      [~,~,rp] = lu(B(I,:),'vector');
    end
    I = sort(I(rp(1:p)));  % sort for faster access
    clear dum rp idx isort xt;
  end
  % solving underdetermined interval system
  if n~=p
    x0(I) = 0;  % used for x0(N)
    % bI = b - B(N,:)' * x0(N)
    [bI,bIrad] = spdotK(b',1,B,-x0,3);
    rnd = getround();
    setround(1);
    bIrad = bIrad + brad + (abs(x0)'*Brad)';
    % basis of B
    B = B(I,:)';  Brad = Brad(I,:)';
    % reset rounding
    setround(rnd);
  else
    B = B';  Brad = Brad';
  end
  % verify lss
  if (VSDP_OPTIONS.VERIFY_FULL_LSS)
    bI = midrad(full(bI),full(bIrad));
    B = midrad(full(B),full(Brad));
    XI = verifylss(B,bI);
  else
    bI = midrad(sparse(bI),sparse(bIrad));
    B = midrad(sparse(B),sparse(Brad));
    bI = B'*bI;
    B = B'*B;
    B = {B.inf, B.sup};  % {inf, sup}
    B = infsup(max(B{1},B{1}'),min(B{2},B{2}'));  % use symmetry
    %B = midrad(tril(B.mid)+tril(B.mid,-1)',min(B.rad,B.rad'));
    XI = verifylss(B,bI);
  end
  % found an inclusion?
  if ~any(isnan(XI)) && ~isempty(XI)
    X = intval(x0);
    X(I) = XI;
  else
    return;
  end
end

% test of inequalities
if ~isempty(a)
  A = midrad(full(A),full(Arad));
  J.ineqlin = find(~(A*X <= a));
end
J.lower = find(~(xl<=X));
J.upper = find(~(X<=xu));

if ~isIntval
  X = struct('mid',X.mid,'rad',X.rad);
end

end
