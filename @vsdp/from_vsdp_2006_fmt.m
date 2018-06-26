function obj = from_vsdp_2006_fmt (blk, A, C, b, X0, y0, Z0)
% FROM_VSDP_2006_FMT  Import SDP problem data from VSDP 2006 format.
%
%   obj = vsdp.FROM_VSDP_2006_FMT (blk, A, C, b, X0, y0, Z0)
%
%   The VSDP 2006 block-diagonal structure format is:
%
%      min  sum(j=1:n| <  C{j}, X{j}>)
%      s.t. sum(j=1:n| <A{i,j}, X{j}>) = b(i)  for i = 1:m
%           X{j} must be positive semidefinite for j = 1:n
%
%   The problem data of the block-diagonal structure:
%
%      'blk'  cell(n,2)
%      'A'    cell(m,n)
%      'C'    cell(n,1)
%      'b'  double(m,1)  Identical to current format.
%
%   The j-th block C{j} and the blocks A{i,j}, for i = 1:m, are real symmetric
%   matrices of common size s_j, and blk(j,:) = {'s'; s_j}.
%
%   The blocks C{j} and A{i,j} must be stored as individual matrices in dense
%   or sparse format.
%
%   The optional initial guess format is:
%
%      'X0'   cell(n,1)
%      'y0' double(m,1)  Identical to current format.
%      'Z0'   cell(n,1)
%
%   Example:
%
%       blk(1,:) = {'s'; 2};
%       A{1,1} = [0 1; 1 0];
%       A{2,1} = [1 1; 1 1];
%         C{1} = [1 0; 0 1];
%            b = [1; 2.0001];
%
%   See also to_vsdp_2006_fmt.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (4, 7);

if (~iscell (blk) || isempty (blk{1,2}))
  error ('VSDP:FROM_VSDP_2006_FMT:no_cone_structure', ...
    'from_VSDP_2006_fmt: bad cone structure `blk`');
end
% Translate cone structure `blk` to `K`.
K.s = cat (1, blk{:,2});

% Translate constraint matrix cells `A` to vectorized transposed matrix `At`.
if (iscell (A) && ~isempty (A{1}))
  % Vectorize all cells of transposed 'A'.
  A = cellfun (@(x) x(:), A', 'UniformOutput', false);
  for i = 1:size(A, 1)
    % Concatenate all vectorized objective and constraint matrices horizontally
    % into the first cell, e.g. for each block i
    %
    %                [ |        |  ]
    %       A{i,1} = [ A1, ..., Am ]
    %                [ |        |  ]
    %
    A{i,1} = horzcat (A{i,:});
  end
  A = vertcat (A{:,1});
else
  A = [];
end

% Translate primal objective matrix cell `C` to vector `c`.
if (iscell (C) && ~isempty (C{1}))
  C = cellfun (@(x) x(:), C, 'UniformOutput', false);
  c = vertcat (C{:});
else
  c = [];
end

% Translate primal solution matrix cell `X0` to vector `x0`.
if ((nargin >= 5) && iscell (X0) && ~isempty (X0{1}))
  X0 = cellfun (@(x) x(:), X0, 'UniformOutput', false);
  x0 = vertcat (X0{:});
else
  x0 = [];
end

% Translate slack variable matrix cell `Z0` to vector `z0`.
if ((nargin == 7) && iscell (Z0) && ~isempty (Z0{1}))
  Z0 = cellfun (@(x) x(:), Z0, 'UniformOutput', false);
  z0 = vertcat (Z0{:});
else
  z0 = [];
end

% Recursive call to default VSDP constructor.
obj = vsdp (A, b, c, K, x0, y0, z0);

end
