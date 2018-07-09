function A = svec (obj, A, mu, param)
% SVEC  Symmetric vectorization operator.
%
%   A = vsdp.svec( [], A)
%   A = vsdp.svec(  K, A)
%   A = vsdp.svec(obj, A)
%   A = vsdp.svec(  ... , mu, opt)
%
%      'K' or 'obj'  Cone structure or VSDP object.
%      'A'           Quadratic matrix to vectorize.
%
%      'mu'     (default = 1) Scaling factor for off-diagonal elements.
%      'param'  Additional parameters for 'A'.
%         - 'unsymmetric' (default)  Assume matrices have no symmetry and both
%                         triangular parts  have to be taken into account.
%         - 'symmetric'   Assume matrices are symmetric and only regard the
%                         upper triangular part.
%
%      For the trace inner product of symmetric matrices X and Y holds:
%
%        <X,Y> := trace(X*Y) = svec(X,K,sqrt(2))' * svec(Y,K,sqrt(2))
%                            = svec(X,K,1)'       * svec(Y,K,2).
%
%      Don't use 'mu = sqrt(2)' for verified computations.  For better
%      comprehension, assume a symmetric 3x3 matrix 'A':
%
%                                 [a]
%                                 [B]
%                                 [C]                            [a   ]
%            [a b c]              [b]                            [b*mu]
%        A = [B d e]  ==>  A(:) = [d]  ==>  vsdp.SVEC(A(:),mu) = [c*mu]
%            [C E f]              [E]                            [d   ]
%                                 [c]                            [e*mu]
%                                 [e]                            [f   ]
%                                 [f]
%
%      Only non-redundant coefficients are stored and the off-diagonal elements
%      are scaled by factor 'mu'.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)


narginchk(2, 4);

% Check optional parameter
if (nargin < 3)
  mu = 1;
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
% a) square double matrix A
if (isempty (obj) && ismatrix (A) && isfloat (A) && all (size (A) == size (A')))
  [~, midx] = vsdp.index(size (A, 1));
  A = A(midx(:,1) | midx(:,3));
  return;
end

% There is nothing to do, if the stored matrix has already the length of the
% condensed semidefinite cones.
if (size(obj.At, 1) == obj.n)
  return;
end

% Compute the index for the upper triangular part of each semidefinite cone
% block.
num_of_sdp_blocks = length(K.s);
if (isempty (obj.svec_idx.upper))
  obj.svec_idx.upper = cell (num_of_sdp_blocks + 1, 1);
  obj.svec_idx.upper{1} = true(nos,1);
  for k = 1:num_of_sdp_blocks
    obj.svec_idx.upper{k+1} = reshape(triu(true(K.s(k))),[],1);
  end
  obj.svec_idx.upper = vertcat (obj.svec_idx.upper{:});
end


% Compute the index for the lower triangular part of each semidefinite cone
% block.
if (~isSymmetric && isempty (obj.svec_idx.lower))
  obj.svec_idx.lower = cell (num_of_sdp_blocks + 1, 1);
  obj.svec_idx.lower{1} = (1:nos)';
  blks = nos + 1;
  for k = 2:num_of_sdp_blocks+1
    nk = K.s(k-1);
    obj.svec_idx.lower{k} = nk * ones(nk*(nk+1)/2,1);
    obj.svec_idx.lower{k}(cumsum([1 1:nk-1])) = [blks 1:-nk:1-nk*(nk-2)];
    obj.svec_idx.lower{k} = cumsum(obj.svec_idx.lower{k});
    blks = blks + nk*nk;
  end
  obj.svec_idx.lower = vertcat (obj.svec_idx.lower{:});
end


% Remove all rows that are not indexed, regard `mu`.
if (isSymmetric)
  A = 0.5 * (A(obj.svec_idx.upper,:) + A(obj.svec_idx.lower,:));
else
  A = A(obj.svec_idx.upper,:);
end

% Scale if necessary.
if (mu ~= 1)
  A = sscale (vA, K, mu);
end



for j = 1 : l
  %Multiplication of lower and upper part
  Aj = A{j};
  Ajl= mult * tril(Aj,-1);
  Aju = Ajl';
  Aj = Ajl + Aju + diag(diag(Aj));
  if isintval(Aj)
    vA = intval(vA);
  end
  blocksize = size(Aj,1);
  blocklength = blocksize*(blocksize+1)/2;
  vAende = vAende + blocklength;
  Index =  repmat((1:blocksize),blocksize,1);
  Jndex = Index';
  vA(vAstart:vAende,1) = Aj(Index<=Jndex);
  %     vA(vAstart:vAende,1) = Aj(find(tril(ones(blocksize))));
  vAstart = vAstart + blocklength;
end

end


