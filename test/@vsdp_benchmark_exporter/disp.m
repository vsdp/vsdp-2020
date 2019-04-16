function disp (obj)
% DISP  Display short export overview.

% Copyright 2018-2019 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

fprintf ('\n VSDP Benchmark Export\n\n');
fprintf ('  Data directory:\n    %s\n\n', obj.data_dir);

if (any (strcmp (obj.cdata_view(1,:), 'lib')))
  benchmarks = unique (obj.columns ({'lib'}).cdata_view(2:end));
  fprintf ('  Benchmarks: %s\n\n', strjoin (benchmarks, ', '));
end
if (any (strcmp (obj.cdata_view(1,:), 'sname')))
  solver = unique (obj.columns ({'sname'}).cdata_view(2:end));
  fprintf ('  Solver:     %s\n\n', strjoin (solver, ', '));
end

[m, n] = size (obj.cdata_view);

if (any ([m, n] > [10, 4]))
  truncation_msg = sprintf(' (%d x %d  truncated)', m, n);
else
  truncation_msg = '';
end

fprintf ('  Current view%s:\n\n', truncation_msg);
obj.cdata_view = obj.cdata_view( ...
  1:min(10,size(obj.cdata_view,1)),...
  1:min(4,size(obj.cdata_view,2)));
str = obj.to_markdown ();
str = regexprep (str, '\n', '\n    ');
fprintf('    %s\n\n', str);
end
