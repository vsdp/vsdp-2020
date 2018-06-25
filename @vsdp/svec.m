function A = svec (obj, A, mu, isSymmetric)
% SVEC  Symmetric vectorization operator for VSDP objects (internal function).
%
%   A = obj.svec(A, mu, isSymmetric)
%
%      'A'            Quadratic matrix to vectorize.
%      'mu'           Scaling factor for off-diagonal elements.
%      'isSymmetric'  Symmetry assumption for 'A'.
%         - false  Assume matrices have no symmetry and both triangular parts
%                  have to be taken into account.
%         - true   Assume matrices are symmetric and only regard the upper
%                  triangular part.
%
%      For the trace inner product of symmetric matrices X and Y holds:
%
%        <X,Y> := trace(X*Y) = svec(X,K,sqrt(2))' * svec(Y,K,sqrt(2))
%                            = svec(X,K,1)'       * svec(Y,K,2).
%
%      Don't use mu = sqrt(2) for verified computations.  For better
%      comprehension, assume a 3x3 double matrix `A`:
%
%                                 [a]
%                                 [b]
%                                 [c]                    [a   ]
%            [a b c]              [b]                    [b*mu]
%        A = [b d e]  ==>  A(:) = [d]  ==>  svec(A,mu) = [c*mu]
%            [c e f]              [e]                    [d   ]
%                                 [c]                    [e*mu]
%                                 [e]                    [f   ]
%                                 [f]
%
%      Only non-redundant coefficients are stored and the off-diagonal elements
%      are scaled by factor `mu`.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)


narginchk(2, 3);
if ((nargin < 2) || isempty(mu))
  mu = 1;
end
if ((nargin < 3) || isempty(isSymmetric))
  isSymmetric = true;
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
