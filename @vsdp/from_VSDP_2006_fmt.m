function obj = from_VSDP_2006_fmt (blk, A, C, b, X0, y0, Z0)
% FROM_VSDP_2006_FMT  Import SDP problem data from VSDP 2006 format.
%
%   obj = vsdp.from_VSDP_2006_fmt (blk, A, C, b, X0, y0, Z0)
%
%      The VSDP 2006 block-diagonal structure format is:
%
%         min  sum(j=1:n| <  C{j}, X{j}>)
%         s.t. sum(j=1:n| <A{i,j}, X{j}>) = b(i)  for i = 1:m
%              X{j} must be positive semidefinite for j = 1:n
%
%      The problem data of the block-diagonal structure:
%
%         'blk'  cell(n,2)
%         'A'    cell(m,n)
%         'C'    cell(n,1)
%         'b'  double(m,1)
%
%      The j-th block C{j} and the blocks A{i,j}, for i = 1:m, are real
%      symmetric matrices of common size s_j, and blk(j,:) = {'s'; s_j}.
%
%      The blocks C{j} and A{i,j} must be stored as individual matrices in
%      dense or sparse format.
%
%      The optional initial guess format is:
%
%         'X0'   cell(n,1)
%         'y0' double(m,1)
%         'Z0'   cell(n,1)
%
%   See also @vsdp/to_VSDP_2006_fmt.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (4, 7);

% Translate cone structure `blk` to `K`.
if (~iscell (blk) || isempty (blk{1,2}))
  error ('VSDP:FROM_VSDP_2006_FMT', ...
    'from_VSDP_2006_fmt: bad cone structure `blk`');
else
  K.s = cat (1, blk{:,2});
end

% index vector to speed up svec
Ivec = [];

% Translate constraint matrix cells `A` to vectorized transposed matrix `At`.
if (iscell (A) && ~isempty (A{1}))
  [m, n] = size(A);
  At = cell(n,1);
  for i = n:-1:1
    At{i} = reshape (cat (2, A{:,i}), [], m);
    A(:,i) = [];
  end
  [At, Ivec] = vsdp.vsvec (cat (1, At{:}), K, 1, 1);
end

% Translate primal objective matrix cell `C` to vector `c`.
if (iscell (C) && ~isempty (C{1}))
  C = cellfun (@(x) x(:), C, 'UniformOutput', false);
  [c, Ivec] = vsdp.vsvec (cat (1, C{:}), K, 1, 1, Ivec);  % mu = 1
end

% Translate primal solution matrix cell `X0` to vector `x0`.
if ((nargin >= 5) && iscell (X0) && ~isempty (X0{1}))
  X0 = cellfun (@(x) x(:), X0, 'UniformOutput', false);
  [x0, Ivec] = vsdp.vsvec (cat (1, X0{:}), K, 2, 1, Ivec);  % mu = 2
else
  x0 = [];
end

if (nargin < 6)
  y0 = [];
end

% Translate slack variable matrix cell `Z0` to vector `z0`.
if ((nargin == 7) && iscell (Z0) && ~isempty (Z0{1}))
  Z0 = cellfun (@(x) x(:), Z0, 'UniformOutput', false);
  z0 = vsdp.vsvec (cat (1, Z0{:}), K, 1, 1, Ivec);  % mu = 1
else
  z0 = [];
end

% Recursive call to default VSDP constructor.
obj = vsdp (At, b, c, K, x0, y0, z0);

end
