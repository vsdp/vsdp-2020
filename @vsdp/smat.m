function A = smat (obj, a, mu)
% SMAT  Inverse operator of symmetric vectorization (vsdp.smat).
%
%   A = vsdp.smat( [], a)
%   A = vsdp.smat(  K, a)
%   A = vsdp.smat(obj, a)
%   A = vsdp.smat(  ... , mu)
%
%   with required arguments:
%
%      'K' or 'obj'  Cone structure or VSDP object.  For details, see
%                    'help vsdp.validate_cone'
%      'a'           smat-vectorized matrices.
%
%   and optional arguments:
%
%      'mu'     (default = 1/sqrt(2)) Scaling factor for off-diagonal elements.
%
%   For a better comprehension, assume a symmetric 3x3 matrix 'A' and see the
%   documentation of 'help vsdp.smat':
%
%                              [a]
%                              [B]
%                              [D]                                   [a   ]
%         [a b d]              [b]                                   [b*mu]
%     A = [B c e]  ==>  A(:) = [c]  -->  vsdp.smat([],A(:), mu) -->  [c*mu] = a
%         [D E f]              [E]  <--  vsdp.smat([],a,  1/mu) <--  [d   ]
%                              [d]                                   [e*mu]
%                              [e]                                   [f   ]
%                              [f]
%
%   See also vsdp.sindex, vsdp.smat, vsdp.sscale.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk(2, 3);
if (~isfloat (a) && ~isintval (a))
  error ('VSDP:smat:badTypeA', ...
    'smat: ''A'' must be (vectorized) floating-point or interval matrix.');
end

% Check optional parameters.
if (nargin < 3)
  mu = 1/sqrt(2);
else
  if (~isscalar (mu) || ~isnumeric (mu))
    error ('VSDP:smat:badMu', ...
      'smat: ''mu'' must be a positive scalar.');
  end
  mu = double (mu);
end

% Determine how to vectorize the input:
% a) square double/sparse/intval matrix A
if (isempty (obj))
  n = length (a);
  K.s = sqrt(0.25 + 2*n) - 0.5;  % because n = K.s*(K.s+1)/2
  N = K.s^2;
  if (~isvector (a) || (K.s ~= round (K.s)))
    error ('VSDP:smat:badA', ...
      ['smat: With empty first parameter ''a'' must be ', ...
      'a svec-vectorized quantity.']);
  end
elseif (isa (obj, 'vsdp'))
  [K, N, n] = deal (obj.K, obj.N, obj.n);
else  % Otherwise 'obj == K' should be a valid cone structure.
  [K, N, n] = vsdp.validate_cone (obj);
end

% There is nothing to do, if there are no semidefinite cones at all or if the
% the matrix has already the length of the condensed semidefinite cones 'n'.
if (~isfield (K, 's') || (size(a, 1) == N))
  return;
end

if (size (a, 1) ~= n)
  error ('VSDP:smat:badA', ...
    ['smat: ''size(A,1)'' must be %d (or %d in condensed format) ', ...
    'but is %d.'], N, n, size (a, 1));
end

if (issparse (a))
  A = sparse (N, size (a, 2));
else
  A = zeros (N, size (a, 2));
end
if (isintval (a))
  A = intval(A);
end

[vidx, midx] = vsdp.sindex (K.s);
% Scale off diagonal elements.
A(midx(:,1),:) = a(vidx(:,1),:);
A(midx(:,2),:) = a(vidx(:,2),:) * mu;
A(midx(:,2),:) = a(vidx(:,2),:) * mu;

end
