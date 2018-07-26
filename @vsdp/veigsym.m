function lambda = veigsym (A)
% VEIGSYM  Verified enclosure for all eigenvalues of a symmetric matrix.
%
%   lambda = VEIGSYM(A) a verified enclosure for all eigenvalues of matrix 'A'
%      in form of an interval vector 'lambda' is computed.  The matrix 'A' must
%      be a full or sparse symmetric matrix.  'A' is allowed to be a real or
%      interval quantity.
%
%   Example:
%
%       A = [1 2 3;
%            2 1 4;
%            3 4 5];
%       A = midrad(A, 0.01*ones(3));
%       lambda = veigsym(A);
%
%   See also vsdplow, vsdpup, vsdpinfeas.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

if ~min(min( mid(A) == mid(A)'))
  error('VSDP:VEIGSYM', 'matrix must be symmetric')
end

% Main routine using Weyl's Perturbation Theorem
% see Lecture Notes on Optimization with Result Verification
[V, D] = eig (full (mid (A)));
E = A - V * intval(D) * V';
r = abss(norm(E,inf));
lambda = midrad(diag(D),r);

end
