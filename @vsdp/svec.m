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
%         - 'symmetric'   (default)  Assume matrices are symmetric and only
%                         regard the upper triangular part.
%         - 'unsymmetric' Assume matrices have no symmetry and both triangular
%                         off-diagonal parts have to be taken into account.
%
%   For the trace inner product of symmetric matrices X and Y holds:
%
%        <X,Y> := trace(X*Y) = vsdp.svec(K,X  )' * vsdp.svec(K,Y  )
%                            = vsdp.svec(K,X,1)' * vsdp.svec(K,Y,2).
%
%   Don't use the default 'mu = sqrt(2)' for verified computations.  For a
%   better comprehension, assume a symmetric 3x3 matrix 'A':
%
%                              [a]
%                              [B]
%                              [D]                                   [a   ]
%         [a b d]              [b]                                   [b*mu]
%     A = [B c e]  ==>  A(:) = [c]  -->  vsdp.svec([],A(:), mu) -->  [c*mu] = a
%         [D E f]              [E]  <--  vsdp.smat([],a,  1/mu) <--  [d   ]
%                              [d]                                   [e*mu]
%                              [e]                                   [f   ]
%                              [f]
%
%   Only non-redundant coefficients (non-capital letters) are stored and the
%   off-diagonal elements are scaled by factor 'mu'.
%
%   See also vsdp, vsdp.sindex, vsdp.smat.
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
  isSymmetric = true;
else
  param = validatestring (param, {'symmetric', 'unsymmetric'});
  isSymmetric = strcmp (param, 'symmetric');
end

% Determine how to vectorize the input:
% a) square double/sparse/intval matrix A
if (isempty (obj))
  if (isvector (sup (A))) % "sup (A)" due to Octave bug with classes.
    % Assume A was already svec-vectorized.
    K.s = sqrt(0.25 + 2*length (A)) - 0.5;  % because n = K.s*(K.s+1)/2
    % If this is not the case, assume A = A(:).
    if (K.s ~= round (K.s))
      K.s = sqrt (length (A));
    end
  else  % Assume Matrix input.
    A = A(:);
    K.s = sqrt (length (A));
  end
  try
    [K, N, n] = vsdp.validate_cone (K);
  catch ME
    switch (ME.identifier)
      case 'VSDP:validate_cone:needPositiveIntegers'
        error ('VSDP:svec:badA', ...
          'svec: With empty first parameter ''A'' must be a square matrix.');
      otherwise
        rethrow (ME);
    end
  end
  obj = K;  % Validated cone to forward.
elseif (isa (obj, 'vsdp'))
  [K, N, n] = deal (obj.K, obj.N, obj.n);
else  % Otherwise 'obj == K' should be a valid cone structure.
  [K, N, n] = vsdp.validate_cone (obj);
  obj = K;  % Validated cone to forward.
end

% There is nothing to do, if there are no semidefinite cones at all.
if (isempty (K.s))
  return;
end

% If the matrix has already the condensed cone length 'n', just scale the off
% diagonal elements by 'mu'.
if (size (A, 1) == n)
  if (mu ~= 1)
    warning ('VSDP:svec:justScale', ...
      ['svec: Input ''A'' is already condensed vectorized, just scale ', ...
      'off-diagonal elements.']);
    idx = ~vsdp.sindex (obj);  % Inverted for off-diagonal elements.
    A(idx,:) = A(idx,:) * mu;
  end
  return;
end

% Otherwise scale the upper triangular matrix under consideration of the
% symmetry and drop the lower triangular elements.
if (size (A, 1) ~= N)
  error ('VSDP:svec:badA', ...
    ['svec: ''size(A,1)'' must be %d (or %d in condensed format) ', ...
    'but is %d.'], N, n, size (A, 1));
end

if (~isSymmetric)
  [~, midx, lidx] = vsdp.sindex (obj);
  % Compute average of lower and upper off diagonal elements and store them in
  % the upper part.
  A(midx(:,2),:) = (A(midx(:,2),:) + A(lidx,:)) / 2;
else
  [~, midx] = vsdp.sindex (obj);
end
% Scale off diagonal elements.
if (mu ~= 1)
  A(midx(:,2),:) = A(midx(:,2),:) * mu;
end
% Preserve non-semidefinite entries.
midx(1:(K.f + K.l + sum (K.q))) = true;
% Drop lower triangular elements.
A = A((midx(:,1) | midx(:,2)),:);

end
