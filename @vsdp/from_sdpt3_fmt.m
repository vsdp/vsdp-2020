function obj = from_sdpt3_fmt (blk, At, b, c, x0, y0, z0)
% FROM_SDPT3_FMT  Import conic problem from a SDPT3 format.
%
%   obj = vsdp.FROM_SDPT3_FMT (blk, At, b, c, x0, y0, z0);
%
%   The problem data of the SDPT3-4.0 is:
%
%      'n'  Number of cones (semidefinite 's', quadratic 'q', linear 'l',
%                            unrestricted (free variables) 'u')
%      'blk'  cell(n,2)
%      'At'   cell(n,1)  Each cell  double(n_j,m)
%      'c'    cell(n,1)
%      'b'  double(m,1)  Identical to current VSDP format.
%
%   In case of a semidefinite cone j, the matrix At{j} has the condensed format,
%   but is scaled with factor 'sqrt(2)'.
%
%   The optional initial guess format is:
%
%      'x0'    cell(n,1)
%      'y0'  double(m,1)  Identical to current VSDP format.
%      'z0'    cell(n,1)
%
%   The following example with block diagonal structure was taken from [2] with
%   vectorization, described in [1].
%
%   Example:
%
%      C{1}    = [ 0 0 0 0;
%                  0 0 0 0;
%                  0 0 1 2;
%                  0 0 2 1];
%      blk{1,1} = 's'; blk{1,2} = [2 2];
%      At{1} = [ ...
%         svec(blk, [ 1 0 0 0;
%                     0 1 0 0;
%                     0 0 0 0;
%                     0 0 0 0]), ...
%         svec(blk, [ 1 0 0 0;
%                     0 0 0 0;
%                     0 0 1 0;
%                     0 0 0 0]), ...
%         svec(blk, [ 0 1 0 0;
%                     1 0 0 0;
%                     0 0 0 1;
%                     0 0 1 0]), ...
%         svec(blk, [ 0 0 0 0;
%                     0 1 0 0;
%                     0 0 0 0;
%                     0 0 0 1])];
%       b = [1; 1; 1; 1];
%       obj = vsdp.FROM_SDPT3_FMT (blk, At, b, C);
%
%
% For more information on the SDPT3-4.0 format, see:
%
%   [1] http://www.math.nus.edu.sg/~mattohkc/sdpt3/guide4-0-draft.pdf
%   [2] http://www.math.nus.edu.sg/~mattohkc/sdpt3/sdpexample.html
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk(4, 7);

% Convert cone structure.
if (isempty (blk))
  error('VSDP:from_sdpt3_fmt:noConeStructure', ...
    'from_sdpt3_fmt: Cannot convert cone structure.');
end
if (size (blk, 2) > 2)
  error('VSDP:from_sdpt3_fmt:lowRankStructureUnsupported', ...
    'from_sdpt3_fmt: VSDP does not support low rank structures.');
end

% Cell indices of respective cones.
idx_f = ([blk{:, 1}] == 'u');
idx_l = ([blk{:, 1}] == 'l');
idx_q = ([blk{:, 1}] == 'q');
idx_s = ([blk{:, 1}] == 's');

K.f = sum ([blk{idx_f, 2}]);
K.l = sum ([blk{idx_l, 2}]);
K.q = [blk{idx_q, 2}];
K.s = [blk{idx_s, 2}];

% Ensure cone order.
sort_idx = [find(idx_f), find(idx_l), find(idx_q), find(idx_s)];
if (length (sort_idx) ~= size (blk, 1))
  error('VSDP:from_sdpt3_fmt:unsupportedCones', ...
    'from_sdpt3_fmt: unsupported ''blk'' cell format or empty cells');
end
c  =  c(sort_idx);
At = At(sort_idx);

% Convert cells to matrices.
c  = vsdp.cell2mat (c(:));
At = cell2mat (At(:));  % No vectorization needed, use native 'cell2mat'.

% Treat optional parameter.
if (nargin >= 5)
  x0 = x0(sort_idx);
  x0 = vsdp.cell2mat (x0(:));
else
  x0 = [];
end
if (nargin < 6)
  y0 = [];
end
if (nargin == 7)
  z0 = z0(sort_idx);
  z0 = vsdp.cell2mat (z0(:));
else
  z0 = [];
end

obj = vsdp (At, b, c, K, x0, y0, z0);

end
