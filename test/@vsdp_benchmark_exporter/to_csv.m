function str = to_csv (obj)
% TO_CSV  Export 'obj.cdata_view' to comma separated values (CSV).

% Copyright 2018-2019 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

% Ensure string entries.
cdata = obj.to_cell_strings ().cdata_view;

for i = 1:size (cdata, 1)
  cdata{i,1} = strjoin (cdata(i,:), ',');
end
cdata = cdata(:,1);
str = strjoin (cdata, '\n');
end
