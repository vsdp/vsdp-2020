function obj = from_sdpam_fmt (bLOCKsTRUCT, c, F, x0, X0, Y0)
% FROM_SDPAM_FMT  Import conic problem from SDPA-M format.
%
%   obj = vsdp.FROM_SDPAM_FMT (bLOCKsTRUCT, c, F, x0, X0, Y0);
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
%      mDIM = 3;
%      nBLOCK = 1;
%      bLOCKsTRUCT = [2];
%      c = [48, -8, 20];
%      F = cell(nBLOCK, mDIM+1);
%      F{1,1} = [-11,  0 ;  0, 23];
%      F{1,2} = [ 10,  4 ;  4,  0];
%      F{1,3} = [  0,  0 ;  0, -8];
%      F{1,4} = [  0, -8 ; -8,  2];
%
%   For more information on the SDPA-M format, see:
%
%     [1] https://sourceforge.net/projects/sdpa/files/sdpa-m/sdpamManual.pdf
%
%   See also from_sdpa_file.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk(3, 6);
if (isempty (bLOCKsTRUCT))
  error('VSDP:SDPAM2VSDP','"bLOCKsTRUCT" has to be set');
end

% Default values.
At = [];
b = -c(:);  % This is true for the VSDP format.
c = [];
x = [];
y = [];
z = [];

idx_lp  = (bLOCKsTRUCT < 2); % Indices of diagonal (LP) and 1x1 blocks.
idx_sdp = ~idx_lp;           % Indices of semidefinite blocks.

K.l = sum (abs (bLOCKsTRUCT(idx_lp)));
K.s = bLOCKsTRUCT(idx_sdp);

if (~isempty (F))
  % Vectorize all cells.
  F = cellfun (@(x) x(:), F, 'UniformOutput', false);
  % Put diagonal (LP) blocks to top.
  F = [F(idx_lp,:); F(idx_sdp,:)];
  
  for i = 1:size(F, 1) % == nBLOCK
    % Concatenate all vectorized objective and constraint matrices
    % horizontally into the first cell, e.g. for each block i
    %
    %                [ |   |        |  ]
    %       F{i,1} = [ C , A1, ..., Am ]
    %                [ |   |        |  ]
    %
    F{i,1} = horzcat (F{i,:});
  end
  % Finally concatenate all first vectorized blocks vertically and drop the
  % rest of F.
  F = -vertcat (F{:,1});  % F = [c, At]
  
  % Split F into c and A.
  c = F(:,1);
  At = F(:,2:end);
end

% Convert x0 to y.
if (nargin > 3)
  y = x0;
end

% Convert X0 to z.
if ((nargin > 4) && (~isempty (X0)))
  % Vectorize all cells.
  X0 = cellfun (@(x) x(:), X0, 'UniformOutput', false);
  % Put diagonal (LP) blocks to top and concatenate.
  z = vertcat([X0(idx_lp); X0(idx_sdp)]);
end

% Convert Y0 to x.
if ((nargin > 5) && (~isempty (Y0)))
  % Vectorize all cells.
  Y0 = cellfun (@(x) x(:), Y0, 'UniformOutput', false);
  % Put diagonal (LP) blocks to top and concatenate.
  x = vertcat([Y0(idx_lp); Y0(idx_sdp)]);
end

obj = vsdp (At, b, c, K, x, y, z);

end
