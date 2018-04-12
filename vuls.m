function [X, I, J] = vuls(A,a,B,b,xl,xu,x0,I,opts)
% VULS  Verification for underdetermined linear systems.
%
%         A * x <= a,
%         B * x == b,
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

narginchk (7, 9);

a = a(:);
if (~(isfloat (a) || isintval (a)) || ~isreal (a))
  error('VSDP:VULS', 'VULS: bad vector `a`.');
end
b = b(:);
if (~(isfloat (b) || isintval (b)) || ~isreal (b))
  error('VSDP:VULS', 'VULS: bad vector `b`.');
end
xl = xl(:);
if (~isfloat (xl) || ~isreal (xl))
  error ('VSDP:VULS', 'VULS: bad lower bound `xl`.');
end
xu = xu(:);
if (~isfloat (xu) || ~isreal (xu))
  error ('VSDP:VULS', 'VULS: bad upper bound `xu`.');
end
x0 = x0(:);
if (~isfloat (x0) || ~isreal (x0))
  error ('VSDP:VULS', 'VULS: bad initial solution set `x0`.');
end

n = length (x0);
m = length (a);
p = length (b);

if (~(isfloat (A) || isintval (A)) || ~isreal (A) || (any (size(A) ~= [m, n])))
  error('VSDP:VULS', 'VULS: bad matrix `A`.');
end

if (~(isfloat (B) || isintval (B)) || ~isreal (B) || (any (size(B) ~= [p, n])))
  error('VSDP:VULS', 'VULS: bad matrix `B`.');
end


% Check index and opts parameter
if ((nargin < 8) || (length(I) ~= p))
  I = [];
end
if (nargin < 9)
  opts = [];
end

VSDP_OPTIONS = vsdpinit(opts);

% Default output.
X = nan(n,1);
J.ineqlin = [];
J.lower = [];
J.upper = [];


% Check bounds `xl`, `xu` for sanity.
if (length (xl) ~= n)
  warning ('VSDP:VULS', 'VULS: using infinite bounds `xl`.');
  xl = -inf(n,1);
end
if (length (xu) ~= n)
  warning ('VSDP:VULS', 'using infinite bounds `xu`.');
  xu = inf(n,1);
end
if (any (xu < xl))
  error ('VSDP:VULS', 'simple bounds are infeasible: `xu >= xl` is violated.');
end
% Projection of `x0` into the interior of [xl, xu]
x0 = max (xl, x0);
x0 = min (xu, x0);

% Compute a verified enclosure `X` for `B * x = b`.
if (~isempty (B))
  % Determine basis by using an LU-decomposition.
  if (isempty (I))
    % Bias towards greater `x0`.
    [~, I] = sort (abs (x0), 'descend');
    if (issparse (B))
      % Use UMFPACK with threshold.
      [~,~,row_perm,~] = lu (mid (B(I,:)), [0.95 0.75], 'vector');
    else
      [~,~,row_perm] = lu (mid (B(I,:)), 'vector');
    end
    % Sort indices for faster access.
    I = sort (I(row_perm(1:p)));
  end
  % Solve underdetermined interval system 
  if (n ~= p)
    % Subtract from right-hand side `b` the non-basis part.
    b = b - B(~I,:)' * x0(~I);
    % Reduce `B` to it's basis.
    B = B(I,:)';
  else
    B = B';
  end

  if (VSDP_OPTIONS.VERIFY_FULL_LSS)
    b = full (b);
    B = full (B);
  else
    % Use normal equations, e.g. solve `(B'*B) * x == (B'*b)`.
    b = B'*b;
    B = B'*B;
  end
  
  % Found an inclusion using INTLAB.
  XI = verifylss (B, b);
  if ((~any (isnan (XI))) && (~isempty (XI)))
    X = intval(x0);
    X(I) = XI;
  else
    return;
  end
end

% Compute inequality violations.
if (~isempty (a))
  J.ineqlin = find (~(full (A) * X <= a));
end
J.lower = find (~(xl <= X));
J.upper = find (~(X <= xu));


end
