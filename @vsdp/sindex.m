function [vidx, midx, lidx] = sindex (obj)
% SINDEX  Compute indices for symmetric matrices of a cone structure 'K'.
%
%   [vidx, midx, lidx] = vsdp.sindex (K);
%
%      'obj'  Cone structure, or VSDP object.  See 'help vsdp.validate_cone'
%             for details.
%
%       'vidx'  Indices of the diagonal entries for the symmetric vectorized
%               (vsdp.svec) matrices.
%      '~vidx'  Indices of the upper tridiagonal entries for the symmetric
%                   vectorized (vsdp.svec) matrices.
%      'midx(:,1)'  Indices of the diagonal entries for the full matrices.
%      'midx(:,2)'  Indices of the upper tridiagonal entries for the full
%                   matrices.
%      'lidx'       Indices of the lower tridiagonal entries for the full
%                   matrices in the order given by 'midx(:,2)'.
%
%   For a better comprehension of this function, consider the following
%   quadratic 4x4 (K.s = 4) matrix 'A' with Fortran indices 'idx' and a
%   vsdp.svec-vectorized Matrix 'Avec':
%
%             [a b d g]            [ 1  5  9 13]
%             [B c e h]            [ 2  6 10 14]
%         A = [D E f i]      idx = [ 3  7 11 15]
%             [G H I j]            [ 4  8 12 16]
%
%      Avec = [a b c d e f g h i j]'
%
%   The goal of symmetric vectorization (vsdp.svec) is to store only the upper
%   triangular part in 'Avec', because the entries B,D,G,... are redundant in
%   the symmetric case.
%
%   To extract the upper part, one has to determine the indices in 'idx', in
%   this example 'Avec == A([1, 5:6, 9:11, 13:16])'.  For reasons of memory
%   efficiency, this function uses logical indexing rather than index vectors,
%   as the order is given natrually by the Fortran indices.
%
%   By computing the index matrices and vector
%
%      [vidx, midx, lidx] = vsdp.sindex(K);   % with K.s = 4
%
%   the main diagonal is indexed by:
%
%      Avec( vidx) == A(midx(:,1)) == [a c f j]'
%
%   and the triangular upper matrix elements are indexed by:
%
%      Avec(~vidx) == A(midx(:,2)) == [b d e g h i]'
%
%   Ocasionally, it is important to extract the lower triangular part of 'A' as
%   well.  But the Fortran index order by logical indexing would yield
%   [B D G E H I].  To obtain the same order of the elements as in
%   'A(midx(:,2))', use the index vector 'lidx' for the lower triangular part:
%
%      A(lidx) == [B D E G H I]'
%
%   NOTE: in general one should avoid to compute 'midx' and 'lidx', as it is an
%   expensive operation.  Working on a vectorized representation of 'A' is
%   more favorable.
%
%   Example:
%
%      A = [11, 12, 13, 14; ...
%           21, 22, 23, 24; ...
%           31, 32, 33, 34; ...
%           41, 42, 43, 44];
%      K.s = 4;
%      [~,midx,lidx] = vsdp.sindex(K);
%      disp ( A(midx(:,1))' )  % Diagonal elements of A
%      disp ( A(midx(:,2))' )  % Upper triangular part of A
%      disp ( A(lidx)' )       % Lower triangular part of A
%
%   See also vsdp.svec, vsdp.smat.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

if (isa (obj, 'vsdp'))
  K = obj.K;
else
  K = vsdp.validate_cone (obj);
end

idxs = cell(length(K.s), max (1, nargout));
switch (nargout)
  case {0, 1}
    for i = 1:length(K.s)
      idxs{i} = idx (K.s(i));
    end
  case 2
    for i = 1:length(K.s)
      [idxs{i,1}, idxs{i,2}] = idx (K.s(i));
    end
  case 3
    for i = 1:length(K.s)
      [idxs{i,1}, idxs{i,2}, idxs{i,3}] = idx (K.s(i));
    end
end
offset = K.f + K.l + sum (K.q);
vidx = false (offset, 1);
if (sum (K.s) > 0)  % Preserve logical data type in case of empty arrays [].
  vidx = [vidx; vertcat(idxs{:,1})];
end
if (nargout > 1)
  midx = false (offset, 2);
  if (sum (K.s) > 0)  % Preserve logical data type in case of empty arrays [].
    midx = [midx; vertcat(idxs{:,2})];
  end
end
if (nargout == 3)
  offset = cumsum ([0; K.s(1:end-1).^2]) + offset;
  for i = 1:length(offset)
    idxs{i,3} = idxs{i,3} + offset(i);
  end
  lidx = vertcat (idxs{:,3});
end

end

function [vidx, midx, lidx] = idx (N)
% IDX  Performs the index computation for a single matrix.
%

vidx = false(N * (N + 1) / 2, 1);
vidx(cumsum(1:N),1) = true;

if (nargout >= 2)
  midx = false(N^2, 2);
  d = diag (true (N, 1));
  midx(:,1) = d(:);
  d = triu (true (N), 1);
  midx(:,2) = d(:);
end

if (nargout == 3)
  % The following few lines of code are very tricky and should be explained.
  % Functionally, the code is equivalent to:
  %
  %    lidx = reshape (1:N^2, N, N);
  %    lidx = tril (lidx, -1)';
  %    lidx(lidx == 0) = [];
  %    lidx = lidx(:);
  %
  % But is by factor 2 slower and requires much more memory.  Assume the
  % following NxN matrix of indices in logical order I and Fortran order J for
  % N = 5:
  %
  %          [ 1            ]         [ 1            ]
  %          [ 2  3         ]         [ 2  7         ]
  %      I = [ 4  5  6      ]     J = [ 3  8 13      ]
  %          [ 7  8  9 10   ]         [ 4  9 14 19   ]
  %          [11 12 13 14 15]         [ 5 10 15 20 25]
  %
  % The desired vector 'lidx' is:
  %
  %   Position:       1 2     4           7            == cumsum ([1, 1:(N-2)]))
  %   lidx =        [ 2 3 8   4   9 14    5    10 15 20]
  %
  %        = cumsum([ 2 1 N (1-N) N  N (1-N-N)  N  N  N])
  %
  
  lidx = N * ones (N * (N - 1) / 2, 1);
  if (~isempty (lidx))  % catch N == 1
    lidx(cumsum ([1, 1:(N-2)])) = [2, (1:-N:(1-N*(N-3)))];
    lidx = cumsum (lidx);
  end
end

end
