function str = to_markdown (obj)
% TO_MARKDOWN  Export 'obj.cdata_view' to Markdown format.

% Copyright 2018-2019 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

% Ensure string entries.
cdata = obj.to_cell_strings ().cdata_view;

% Add row below.
cdata{end+1,1} = [];
for i = 1:size (cdata, 2)
  % Determine max. column width.
  col_width = max (cellfun (@length, cdata(:,i)));
  % Stretch cells left aligned to column max. column width.
  cdata([1,3:end],i) = cellfun ( ...
    @(str) sprintf (['%-', num2str(col_width), 's'], ...
    str), cdata(1:end-1,i), 'UniformOutput', false);
  cdata{2,i} = repmat ('-', 1, col_width);
end
% Concatenate linewise.  Add dummy columns to get first and last bar.
for i = 1:size (cdata, 1)
  cdata{i,1} = strtrim (strjoin ([{''}, cdata(i,:), {''}], ' | '));
end
cdata = cdata(:,1);
str = strjoin (cdata, '\n');
end
