function x = verify_uls (obj, A, b, x0)
% VERIFY_ULS  Verification for underdetermined linear interval systems.
%
%   Let a linear systems of equations with inexact input data
%
%      [A] * x = [b]
%
%   with an interval matrix A (n x n) and an interval right-hand side b
%   (n x 1).  The aim of this function is to compute an interval vector [x]
%   (n x 1) containing the solution set
%
%      SolSet([A],[b]) := { x in R^n : Ax = b  for some  A in [A], b in [b] }.
%
%   If all A in [A] are nonsingular, then the solution set is bounded and
%   satisfies, by definition, the property
%
%   (i)  for all A in [A], for all b in [b] exists some x in [x] : Ax = b.
%
%   For the computation of enclosures in the case of large linear systems, we
%   refer to <https://vsdp.github.io/references.html#Rump2013a>.
%
%   The computation of rigorous lower and upper bounds for the optimal value
%   requires considering a modified problem.  There, a nonsquare interval
%   matrix [A] (m x n) with m < n and a right-hand side [b] (m x 1) are given,
%   and the goal is to compute an interval vector [x] such that property (i) is
%   fulfilled.
%
%   Since m < n (in most cases m is much smaller than n), the solution set
%   SolSet([A],[b]) is in general unbounded, whereas property (i) requires
%   finding only an enclosure of a part of the solution set.
%
%   Obviously, there are many possibilities for computing such a part of the
%   solution set.  We need to compute such an enclosure [x] with respect to a
%   given vector x0 and an index set I as subset of {1:n}, and proceed as
%   follows:
%
%     1. Find a basis index set I (m x 1), such that A(:,I) is regular.
%     2. Set [b] := [b] - [A(:,~I)] * x0(~I), where ~I := {1:n} \ I.
%     3. set [A] := [A(:,I)].
%     4. Compute an enclosure [x] of the solution set with square interval
%        matrix and right-hand side such that  SolSet([A],[b]) in [x]  by using
%        an algorithm for square linear interval systems.
%
%   For more theoretical background, see:
%
%     [1] https://vsdp.github.io/references.html#Jansson2004
%     [2] https://vsdp.github.io/references.html#Jansson2007a
%
%   Example:
%
%       A  = [0 1 0 infsup(0.9,1.1)];
%       b  = 2;
%       x0 = [0 1 0 1]';
%       x  = vsdp.verify_uls([], A, b, x0);
%
%   See also vsdpinfeas, vsdplow, vsdpup.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (4, 4);

if (~isempty (obj) && ~isa (obj, 'vsdp'))
  error ('VSDP:vuls:badObj', ...
    'vuls: The first argument must be empty or a VSDP object.');
end

b = b(:);
if (~(isfloat (b) || isintval (b)) || ~isreal (b))
  error ('VSDP:vuls:badB', ...
    'vuls: data type of right-hand side vector ''b'' unsupported.');
end
m = length (b);

x0 = x0(:);
if (~isfloat (x0) || ~isreal (x0))
  error ('VSDP:vuls:badX0', ...
    'vuls: data type of initial solution ''x0'' unsupported.');
end
n = length (x0);

if (m >= n)
  warning ('VSDP:vuls:notUnderdetermined', ...
    'vuls: The system is not underdetermined m >= n.');
end

if (~(isfloat (A) || isintval (A)) || ~isreal (A))
  error ('VSDP:vuls:badA', ...
    'vuls: data type of matrix ''A'' unsupported.');
end
if (any (size (A) ~= [m, n]))
  error ('VSDP:vuls:badA', ...
    'vuls: ''A'' must be a %d x %d matrix, but is a %d x %d matrix.', ...
    m, n, size (A,1), size (A,2));
end


% Step 1: Determine basis index set 'I'.
if (isempty (obj))
  I = [];
else
  I = obj.cache('I');
end
if (isempty (I))
  % Bias towards greater entries in 'x0'.
  [~,I] = sort (abs (x0), 'descend');
  % Reorder the columns in 'A' and transpose.
  At = mid (A(:,I))';
  % Now the LU-decomposition will perform a "row"-pivoting 'p', which
  % essentially becomes a "column"-pivoting that is wanted.
  if (issparse (At))
    % Use UMFPACK with threshold.
    [~,~,p,~] = lu (At, [0.95 0.75], 'vector');
  else
    [~,~,p] = lu (At, 'vector');
  end
  % Sort indices for faster access.
  I = sort (I(p(1:m)));
  if (~isempty (obj))
    obj.cache(I);
  end
end

% In the following the entries x0(I) shall be replaced by a verified
% enclosure.  Therefore, set them explicitly to zero, such that the
% reduction of the right-hand side 'b' can happen in a straight forward
% way.
x0(I) = 0;

% Steps 2+3: Prepare square nonsingular interval system.
%
% Subtract from right-hand side the non-basis part.
b = b - A * x0;
% Reduce matrix to nonsingular basis.
A = A(:,I);

% Step 4: Compute an enclosure of the solution set.
xI = verifylss (A, b);
if ((~any (isnan (xI))) && (~isempty (xI)))
  % The idea is not to replace the entire approximate solution x0 by the
  % verified enclosure xI.  The approximate solution is "extended" to contain
  % the solution set.
  x = intval(x0);
  x(I) = xI;
else
  x = nan(n,1);
end

end
