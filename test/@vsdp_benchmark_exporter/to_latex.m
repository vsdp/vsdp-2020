function str = to_latex (obj, body_only)
% TO_LATEX  Export 'obj.cdata_view' to LaTeX markup (longtable environment).

% Copyright 2018-2019 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

% Ensure string entries.
cdata = obj.to_cell_strings ().cdata_view;

% Save table head.
thead = strjoin (cdata(1,:), '\n& ');

% Process table body.
for i = 2:size (cdata, 1)
  cdata{i,1} = sprintf ('%s \\\\', strjoin (cdata(i,:), '\n& '));
end
cdata = cdata(2:end,1);
str = sprintf ('%s\n', strjoin (cdata, '\n'));

% Improve output for special numerical values.
str = strrep (str, '& Inf',  '& {$+\infty$}');
str = strrep (str, '& -Inf', '& {$-\infty$}');
str = strrep (str, '& NaN',  '& {NaN}');

if ((nargin == 2) && strcmp (body_only, 'body_only'))
  thead = [];
  header = [];
  footer = [];
else
  thead = sprintf ('{%s}\n\\toprule\n%s \\\\\n\\toprule\n', ...
    repmat('c', 1, size (cdata(1,:), 2)), thead);
  header = strjoin ({ ...
    '\documentclass{article}', ...
    '\usepackage{booktabs}', ...
    '\usepackage{longtable}', ...
    '\usepackage{mathtools}', ...
    '', ...
    '\begin{document}', ...
    '\begin{longtable}'}, '\n');
  footer = strjoin ({ ...
    '\bottomrule', ...
    '\end{longtable}', ...
    '\end{document}'}, '\n');
end
str = sprintf ('%s%s%s%s', header, thead, str, footer);

% Fix unescaped underscores in entire output.
str = regexprep (str, '(?<!\\)_',  '\\_');
end
