function x = cell2mat (X)
% CELL2MAT  Vectorizes all cells and concatenates them to matrix.
%
%    In the conversion methods of VSDP an often reoccurring pattern is to
%    convert some `A' = cell(n,m)` into a vectorized matrix `At = double(N,m)`,
%    where `N >> n` depends on the content of the cells:
%
%         { [        ]     [        ] }            [     |             |     ]
%         { [ X{1,1} ] ... [ X{1,m} ] }            [ X{1,1}(:) ... X{1,m}(:) ]
%         { [        ]     [        ] }            [     |             |     ]
%     X = {                           }   ==>  x = [                         ]
%         { [        ]     [        ] }            [     |             |     ]
%         { [ X{n,1} ] ... [ X{n,m} ] }            [ X{n,1}(:) ... X{n,m}(:) ]
%         { [        ]     [        ] }            [     |             |     ]
%
if (iscell (X))
  x = cell2mat (cellfun (@(x) x(:), X, 'UniformOutput', false));
else
  error ('VSDP:CELL2MAT:no_cell', ...
    'cell2mat: ''%s'' must be a cell array.', inputname(1));
end
end