function A = svec (obj, A, mu, param)
% SVEC  Symmetric vectorization operator.
%
%   a = vsdp.svec( [], A)
%   a = vsdp.svec(  K, A)
%   a = vsdp.svec(obj, A)
%   a = vsdp.svec(  ... , mu, param)
%
%   with required arguments:
%
%      'K' or 'obj'  Cone structure or VSDP object.  For details, see
%                    'help vsdp.validate_cone'
%      'A'           Quadratic matrix or vectorized quadratic matrices.
%
%   and optional arguments:
%
%      'mu'     (default = sqrt(2)) Scaling factor for off-diagonal elements.
%      'param'  Additional parameters for 'A'.
%         - 'unsymmetric' (default)  Assume matrices have no symmetry and both
%                         triangular offdiagonal parts  have to be taken into
%                         account.
%         - 'symmetric'   Assume matrices are symmetric and only regard the
%                         upper triangular part.
%
%   For the trace inner product of symmetric matrices X and Y holds:
%
%        <X,Y> := trace(X*Y) = svec(K,X,sqrt(2))' * svec(K,Y,sqrt(2))
%                            = svec(K,X,1)'       * svec(K,Y,2).
%
%   Don't use the default 'mu = sqrt(2)' for verified computations.  For a
%   better comprehension, assume a symmetric 3x3 matrix 'A':
%
%                              [a]
%                              [B]
%                              [D]                                 [a   ]
%         [a b d]              [b]                                 [b*mu]
%     A = [B c e]  ==>  A(:) = [c]  -->  vsdp.SVEC(A(:), mu)  -->  [c*mu] = a
%         [D E f]              [E]  <--  vsdp.smat(a,  1/mu)  <--  [d   ]
%                              [d]                                 [e*mu]
%                              [e]                                 [f   ]
%                              [f]
%
%   Only non-redundant coefficients (non-capital letters) are stored and the
%   off-diagonal elements are scaled by factor 'mu'.
%
%   See also vsdp.sindex, vsdp.smat.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)


narginchk(2, 4);

if (~isfloat (A) && ~isintval (A))
  error ('VSDP:svec:badTypeA', ...
    'svec: ''A'' must be (vectorized) floating-point or interval matrix.');
end

% Check optional parameters.
if (nargin < 3)
  mu = sqrt(2);
else
  if (~isscalar (mu) || ~isnumeric (mu))
    error ('VSDP:svec:badMu', ...
      'svec: ''mu'' must be a positive scalar.');
  end
  mu = double (mu);
end
if (nargin < 4)
  isSymmetric = false;
else
  param = validatestring (param, {'symmetric', 'unsymmetric'});
  isSymmetric = strcmp (param, 'symmetric');
end

% Determine how to vectorize the input:
% a) square double/sparse/intval matrix A
if (isempty (obj))
  A = A(:);
  K.s = length (A);
  if (K.s ~= round (K.s))
    error ('VSDP:svec:badA', ...
      'svec: With empty first parameter ''A'' must be a square matrix.');
  end
elseif (isa (obj, 'vsdp'))
  [K, N, n] = deal (obj.K, obj.N, obj.n);
else
  [K, N, n] = validate_cone (obj)
end

% There is nothing to do, if the matrix has already the length of the condensed
% semidefinite cones 'n' or if there are no semidefinite cones at all.
if (~isfield (K, 's') || (size(A, 1) == n))
  return;
end

if (~isSymmetric)
  [~, midx, lidx] = vsdp.sindex (K.s);
  % Compute average of lower and upper off diagonal elements and store them in
  % the upper part.
  A(midx(:,2)) = (A(midx(:,2)) + A(lidx)) / 2;
else
  [~, midx] = vsdp.sindex (K.s);
end
% Scale off diagonal elements.
A(midx(:,2)) = A(midx(:,2)) * mu;
% Drop lower triangular elements.
A = A(midx(:,1) | midx(:,2));
return;

end
