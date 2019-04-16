function obj = columns (obj, names)
% COLUMNS  Use specified columns by header names or indices.

% Copyright 2018-2019 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

if (isnumeric (names) || islogical (names))
  obj.cdata_view = obj.cdata_view(:,names);
elseif (iscell (names))
  names = names(:)';
  idx = zeros (size (names));
  for i = 1:length (names)
    if (isempty (name2idx (obj, names{i})))
      obj.cdata_view(:,end + 1) = new_column (obj, names{i});
    end
    idx(i) = name2idx (obj, names{i});
  end
  obj.cdata_view = obj.cdata_view(:,idx);
else
  error ('VSDP_BENCHMARK:EXPORTER:columns', ...
    'columns: Input must be numeric indices or a cell array of strings.');
end
end


function new_col = new_column (obj, col)
% NEW_COLUMN  Compute a dependent column of existing table.

new_col = cell (size (obj.cdata_view, 1), 1);
new_col{1,1} = col;

if (strncmp (col, 'mu_', 3))
  % Compute accuracy mu <https://vsdp.github.io/references.html#Jansson2006>,
  % that is 'mu_<op1>_<op2>'.
  cols = strsplit (col(4:end), '_');
  fop = @(a, b) (a - b) ./ max (1, (abs (a) + abs (b)) ./ 2);
elseif (any (col == '/'))  % '<op1>/<op2>'
  cols = strsplit (col, '/');
  fop = @rdivide;
elseif (any (col == '-'))  % '<op1>-<op2>'
  cols = strsplit (col, '-');
  fop = @minus;
else
  warning ('VSDP_BENCHMARK:export:nonExistingColumn', ...
    'export: Ignore column ''%s''.', col);
  new_col = [];
  return;
end

% Resolve column number.
cols{1} = name2idx (obj, cols{1});
cols{2} = name2idx (obj, cols{2});

% Perform operation for non-empty values.
idx = cellfun (@(x) ~isempty(x), obj.cdata_view(2:end,[cols{1}, cols{2}]));
idx = [false; idx(:,1) & idx(:,2)];
new_col(idx) = num2cell(fop ([obj.cdata_view{idx,cols{1}}], ...
  [obj.cdata_view{idx,cols{2}}]));
end


function idx = name2idx (obj, name)
idx = find (strcmp (name, obj.cdata_view(1,:)));
end
