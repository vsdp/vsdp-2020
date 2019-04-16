function obj = to_cell_strings (obj, use_columns)
% TO_CELL_STRINGS  Export numerical values in 'obj.cdata_view' to strings.
%
%   The second input 'use_columns' determines the string conversion.
%

% Copyright 2018-2019 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

if (all (all (cellfun (@ischar, obj.cdata_view))))
  return;  % All done.
end

if (nargin == 2)
  if ischar (use_columns)
    % Default columns to export.
    u_columns = { ...
      'lib',      '%s',   'Library';
      'name',     '%s',   'Name';
      'sname',    '%s',   'Solver';
      'm',        '%d',   '$m$';
      'n',        '%d',   '$n$';
      'K_f',      '%d',   '$K_f$';
      'K_l',      '%d',   '$K_l$';
      'K_q',      '%d',   '$K_q$';
      'K_s',      '%d',   '$K_s$';
      'fp',       '%e',   '$f_p$';
      'fd',       '%e',   '$f_d$';
      'fL',       '%.8e', '$\underline{f_p}$';
      'fU',       '%.8e', '$\overline{f_d}$';
      'fU-fL',    '%.8e', '$\overline{f_d} - \underline{f_p}$';
      'mu_fU_fL', '%.2e', '$\mu(\overline{f_d}, \underline{f_p})$';
      'ts',       '%.2f', '$t_s$';
      'tL',       '%.2f', '$\underline{t}$';
      'tL/ts',    '%.2e', '$\underline{t}/t_s$';
      'tU',       '%.2f', '$\overline{t}$';
      'tU/ts',    '%.2e', '$\overline{t}/t_s$'};
    switch (use_columns)
      case 'default'
        use_columns = u_columns(:,1:2);
      case 'default_with_headers'
        use_columns = u_columns;
      otherwise
        error ('VSDP_BENCHMARK:EXPORTER:to_cell_strings', ...
          'to_cell_strings: Use ''default'' or ''default_with_headers''.');
    end
    
    % Match header.
    idx = zeros (size (obj.cdata_view(1,:)));
    for i = 1:length(obj.cdata_view(1,:))
      idx(i) = find (strcmp (use_columns(:,1), obj.cdata_view(1,i)));
    end
    use_columns = use_columns(idx, :);
  end
else
  use_columns = [];
end

% Convert data to strings.
if (size (use_columns, 2) > 1)
  for j = 1:size(obj.cdata_view, 2)
    obj.cdata_view(2:end,j) = cellfun (...
      @(x) sprintf (use_columns{j,2}, x), obj.cdata_view(2:end,j), ...
      'UniformOutput', false);
  end
else
  obj.cdata_view = cellfun (@num2str, obj.cdata_view, ...
    'UniformOutput', false);
end

% Replace header.
if (size (use_columns, 2) == 3)
  obj.cdata_view(1,:) = use_columns(:,3);
end
end
