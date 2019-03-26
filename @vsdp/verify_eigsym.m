function E = verify_eigsym (A)
% VERIFY_EIGSYM  Verified enclosure for all eigenvalues of a symmetric matrix.
%
%   E = VERIFY_EIGSYM(A) a verified enclosure for all eigenvalues of matrix 'A'
%      in form of an interval vector 'E' is computed.  The matrix 'A' must
%      be a full or sparse symmetric matrix.  'A' is allowed to be a real or
%      interval quantity.
%
%   For more theoretical background, see:
%
%     https://vsdp.github.io/references.html#Jansson2007 (p. 19)
%
%   Example:
%
%       % Matrix A has eigenvalues {0, 3, 6}.
%       A = [ 3 -2  2;
%            -2  4  0;
%             2  0  2];
%       A = midrad (A, 0.01*ones(3));
%       E = vsdp.verify_eigsym (A);
%
%   See also vsdp, rigorous_lower_bound, rigorous_upper_bound,
%            check_primal_infeasible, check_dual_infeasible.

% Copyright 2004-2019 Christian Jansson (jansson@tuhh.de)

if (~issymmetric (mid (A)))
  error ('VSDP:verify_eigsym:notSymmetric', ...
    'verify_eigsym: Matrix must be symmetric.')
end

% Make use of Weyl's Perturbation Theorem.
[V, D] = eig (full (mid (A)));
E = A - V * intval (D) * V';
r = mag (norm (E, inf));
E = midrad (diag (D), r);

end
