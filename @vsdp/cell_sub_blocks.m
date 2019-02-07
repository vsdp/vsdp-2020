function x = cell_sub_blocks (X, blk)
% CELL_SUB_BLOCKS  Vectorize cell blocks with SDP sub blocks.
%
%   x = cell_sub_blocks (X, K)  If there are no SDP sub-blocks, this function
%                               is equivalent to vsdp.cell2mat.
%
%   A SDP sub-block is
%
%     blk{j,1} = 's';
%     blk{j,2} = [ s_j1, s_j2, ..., s_jp ];
%
%
%              { [       ]                }               [     |    ]
%              { [ X{j1} ]                }               [ X{j1}(:) ]
%              { [       ]                }               [     |    ]
%       X{j} = {            ...           }   ==>  x(j) = [          ]
%              {                [       ] }               [     |    ]
%              {                [ X{jp} ] }               [ X{jp}(:) ]
%              {                [       ] }               [     |    ]
%
%   See also vsdp, vsdp.cell2mat.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

if (nargin == 2 && iscell (X) && iscell (blk))
  idx = find (strcmp (blk(:,1), 's') & (cellfun (@length, blk(:,2)) > 1));
  for j = 1:length (idx)
    Xj = X{idx(j)};
    % Create sub-block indices: [{1:3}, {4:5}, ...].
    subidx = mat2cell (1:length(Xj), 1, blk{idx(j),2});
    % Extract the sub blocks to cell array.
    Xj = cellfun(@(x) Xj(x,x), subidx, 'UniformOutput', false);
    % Vectorize and write back.
    X{idx(j)} = cell2mat (cellfun (@(x) x(:), Xj(:), 'UniformOutput', false));
  end
  x = cell2mat (cellfun (@(x) x(:), X, 'UniformOutput', false));
else
  error ('VSDP:CELL_SUB_BLOCKS:noCell', ...
    'cell_sub_blocks: Inputs must be a cell arrays.');
end
end
