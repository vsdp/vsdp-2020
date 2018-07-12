function [vidx, midx, lidx] = sindex (d)
%INDEX  Compute indices for symmetric matrices of dimensions 'd'.
%
%   [vidx, midx, lidx] = sindex (d);
%
%      'd'  scalar or vector containing the matrix dimensions.
%
%      'vidx(:,1)'  Indices of the diagonal entries for the symmetric vectorized
%                   (svec) matrices.
%      'vidx(:,2)'  Indices of the upper tridiagonal entries for the symmetric
%                   vectorized (svec) matrices.
%      'midx(:,1)'  Indices of the diagonal entries for the full d*d matrices.
%      'midx(:,2)'  Indices of the upper tridiagonal entries for the full d*d
%                   matrices.
%      'lidx'       Indices of the lower tridiagonal entries for the full d*d
%                   matrices in the order of 'midx(:,2)'.
%
%   For a better comprehension of this function, consider the following
%   quadratic 4x4 (d = 4) matrix 'A' with Fortran indices 'idx' and a
%   svec-vectorized Matrix 'a':
%
%             [a b d g]            [ 1  5  9 13]
%             [B c e h]            [ 2  6 10 14]
%         A = [D E f i]      idx = [ 3  7 11 15]
%             [G H I j]            [ 4  8 12 16]
%
%      Avec = [a b c d e f g h i j]'
%
%   The goal of symmetric vectorization (svec) is to store only the upper
%   triangular part in 'Avec', because the entries B,D,G,... are redundant in
%   the symmetric case.
%
%   To extract the upper part, one has to determine the indices in 'idx', in
%   this example 'Avec == A([1, 5:6, 9:11, 13:16])'.  For reasons of memory
%   efficiency, this function uses logical indexing rather than index vectors,
%   as the order is given natrually by the Fortran indices.
%
%   By computing the logical index matrices
%
%      [vidx, midx, lidx] = vsdp.sindex(4);
%
%   the expression
%
%      Avec(vidx(:,1)) == A(midx(:,1)) == [a c f j]'
%
%   extracts main diagonal and
%
%      Avec(vidx(:,2)) == A(midx(:,2)) == [b d e g h i]'
%
%   the triangular upper matrix elements.  Ocasionally, it is important to
%   extract the lower triangular part of 'A' as well.  But the Fortran index
%   order by logical indexing would yield [B D G E H I].  To obtain the
%   same order of the elements as in 'A(midx(:,2))', use the index vector
%   'lidx' for the lower triangular part:
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
%      [~,midx,lidx] = vsdp.sindex(4);
%      disp ( A(midx(:,1))' )  % Diagonal elements of A
%      disp ( A(midx(:,2))' )  % Upper triangular part of A
%      disp ( A(lidx)' )       % Lower triangular part of A
%
%   See also vsdp.svec, vsdp.smat.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

if (isempty (d) || ~isnumeric (d) || ~isvector (d))
  error ('VSDP:index:badD', ...
    'index: ''d'' must be a numeric vector or scalar.');
end
d = d(:);

idxs = cell(length(d), max (1, nargout));
switch (nargout)
  case {0, 1}
    for i = 1:length(d)
      idxs{i} = idx (d(i));
    end
  case 2
    for i = 1:length(d)
      [idxs{i,1}, idxs{i,2}] = idx (d(i));
    end
  case 3
    for i = 1:length(d)
      [idxs{i,1}, idxs{i,2}, idxs{i,3}] = idx (d(i));
    end
end
vidx = vertcat (idxs{:,1});
if (nargout > 1)
  midx = vertcat (idxs{:,2});
end
if (nargout == 3)
  offset = cumsum ([0; d(1:end-1).^2]);
  for i = 1:length(offset)
    idxs{i,3} = idxs{i,3} + offset(i);
  end
  lidx = vertcat (idxs{:,3});
end

end

function [vidx, midx, lidx] = idx (dim)
% IDX  Performs the index computation for a single matrix.
%

vidx = false(dim * (dim + 1) / 2, 2);
vidx(cumsum(1:dim),1) = true;
vidx(:,2) = ~(vidx(:,1));

if (nargout >= 2)
  midx = false(dim^2, 2);
  d = diag (true (dim, 1));
  midx(:,1) = d(:);
  d = triu (true (dim), 1);
  midx(:,2) = d(:);
end

if (nargout == 3)
  lidx = reshape (1:dim^2, dim, dim);
  lidx = tril (lidx, -1)';
  lidx(lidx == 0) = [];
  lidx = lidx(:);
end

end
