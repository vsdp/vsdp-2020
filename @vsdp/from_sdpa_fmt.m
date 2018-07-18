function obj = from_sdpa_fmt (bLOCKsTRUCT, c, F, x0, X0, Y0)
% FROM_SDPA_FMT  Import conic problem from SDPA-M format.
%
%   obj = vsdp.FROM_SDPA_FMT (bLOCKsTRUCT, c, F, x0, X0, Y0);
%
%   The primal SDPA-M format is:
%
%      min  c' * x
%      s.t. X = sum(i=1:m| F(i) * x(i)) - F(0)
%           X positive semidefinite
%
%   The problem data is:
%
%      'mDIM'   double(1,1)  Number of primal variables.
%      'nBLOCK' double(1,1)  Number of SDP blocks in 'F'.
%
%      'bLOCKsTRUCT'  double(nBLOCK,1)  Block structure of 'F'.
%          bLOCKsTRUCT(i) > 0: symmetric matrix
%          bLOCKsTRUCT(i) < 0: diagonal  matrix
%      'c'   double(mDIM,1)     Objective function vector.
%      'F'   cell(nBLOCK,mDIM)  Constraint matrices.
%      'x0'  double(mDIM,1)     Approximation of initial solution 'x'.
%      'X0'  cell(nBLOCK,mDIM)  Approximation of initial solution 'X'.
%      'Y0'  cell(nBLOCK,mDIM)  Approximation of initial solution 'Y'.
%
%   Example 1 from [1] with a single 2x2 SDP block and three constraints:
%
%      m = 3;
%      mDIM = m + 1;
%      nBLOCK = 1;
%      bLOCKsTRUCT = [2];
%      c = [48, -8, 20];
%      F = cell(nBLOCK, mDIM);
%      F{1,1} = [-11,  0 ;  0, 23];
%      F{1,2} = [ 10,  4 ;  4,  0];
%      F{1,3} = [  0,  0 ;  0, -8];
%      F{1,4} = [  0, -8 ; -8,  2];
%      obj = vsdp.FROM_SDPA_FMT (bLOCKsTRUCT, c, F);
%
%   For more information on the SDPA-M format, see:
%
%     [1] https://sourceforge.net/projects/sdpa/files/sdpa-m/sdpamManual.pdf
%
%   See also from_sdpa_file.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk(3, 6);
if (isempty (bLOCKsTRUCT))
  error ('VSDP:FROM_SDPAM_FMT:no_block_structure', ...
    'from_sdpam_fmt: ''bLOCKsTRUCT'' has to be set.');
end
if (size (F, 1) ~= length (bLOCKsTRUCT))
  error ('VSDP:FROM_SDPAM_FMT:bad_cone_dimension', ...
    ['from_sdpam_fmt: According to ''bLOCKsTRUCT'' the ''size(F,1)'' ', ...
    'must be ''%d'' but is ''%d''.'], length (bLOCKsTRUCT), size (F, 1));
end
if (size (F, 2) ~= (length (c) + 1))
  error ('VSDP:FROM_SDPAM_FMT:bad_cone_dimension', ...
    ['from_sdpam_fmt: According to ''c'' the ''size(F,2)'' ', ...
    'must be ''%d'' but is ''%d''.'], length (c) + 1, size (F, 2));
end

idx_lp  = (bLOCKsTRUCT < 2); % Indices of diagonal (LP) and 1x1 blocks.
idx_sdp = ~idx_lp;           % Indices of semidefinite blocks.

K.l = sum (abs (bLOCKsTRUCT(idx_lp)));
K.s = bLOCKsTRUCT(idx_sdp);

b = -c(:);  % This is true for the VSDP format.  Here 'c' is the input!

% Put diagonal (LP) blocks to top and vectorize cells to matrix.
F = [F(idx_lp,:); F(idx_sdp,:)];
F = -vsdp.cell2mat (F);  % F = [c, At]

% Split F into c and A.
At = F(:,2:end);
c = F(:,1);       % Creates new 'c'!

% Treat optional paramter.
if (nargin > 3)
  y = x0;
else
  y = [];
end

if (nargin > 4)
  % Put diagonal (LP) blocks to top and concatenate.
  z = vsdp.cell2mat ([X0(idx_lp); X0(idx_sdp)]);
else
  z = [];
end

if (nargin > 5)
  % Put diagonal (LP) blocks to top and concatenate.
  x = vertcat([Y0(idx_lp); Y0(idx_sdp)]);
else
  x = [];
end

% The VSDP constructor cares for the condensed semidefinite variables.
obj = vsdp (At, b, c, K, x, y, z);

end
