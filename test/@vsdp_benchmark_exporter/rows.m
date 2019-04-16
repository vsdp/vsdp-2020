function obj = rows (obj, op)
% ROWS  Select specified rows by 'key=val' string or indices.
%
%    In case of a string, 'key' must match a column header name and 'val' is a
%    regular expression to be matched.
%

% Copyright 2018-2019 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

if (isnumeric (op) || islogical (op))
  obj.cdata_view = obj.cdata_view(op,:);
elseif (ischar (op))
  op = strsplit (op, '=');
  idx = name2idx (obj, op{1});
  idx = cellfun (@(x) ~isempty (regexp (x, op{2}, 'once')), ...
    obj.cdata_view(2:end,idx));
  obj.cdata_view = obj.cdata_view([true; idx],:);
else
  error ('VSDP_BENCHMARK:EXPORTER:rows', ...
    'rows: Input must be numeric indices or a ''key=val'' string.');
end
end

function idx = name2idx (obj, name)
idx = find (strcmp (name, obj.cdata_view(1,:)));
end
