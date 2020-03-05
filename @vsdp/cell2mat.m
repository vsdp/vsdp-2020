function x = cell2mat (X)
% CELL2MAT  Vectorize all cells and concatenates them to matrix.
%
%   x = cell2mat (X)  In the conversion methods of VSDP an often reoccurring
%                     pattern is to convert some `A' = cell(n,m)` into a
%                     vectorized matrix `At = double(N,m)`, where `N >> n`
%                     depends on the content of the cells:
%
%         { [        ]     [        ] }            [     |             |     ]
%         { [ X{1,1} ] ... [ X{1,m} ] }            [ X{1,1}(:) ... X{1,m}(:) ]
%         { [        ]     [        ] }            [     |             |     ]
%     X = {                           }   ==>  x = [                         ]
%         { [        ]     [        ] }            [     |             |     ]
%         { [ X{n,1} ] ... [ X{n,m} ] }            [ X{n,1}(:) ... X{n,m}(:) ]
%         { [        ]     [        ] }            [     |             |     ]
%
%   See also vsdp.
%

% Copyright 2004-2020 Christian Jansson (jansson@tuhh.de)

if (iscell (X))
  x = cell2mat (cellfun (@(x) x(:), X, 'UniformOutput', false));
else
  error ('VSDP:CELL2MAT:noCell', ...
    'cell2mat: ''%s'' must be a cell array.', inputname(1));
end
end
