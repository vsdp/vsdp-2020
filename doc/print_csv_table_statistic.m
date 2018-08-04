function str = print_csv_table_statistic (filename)
% PRINT_CSV_TABLE_STATISTIC  Prints short statistic of a benchmark CSV table.
%
%   See also benchmark.

% Copyright 2016-2018 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

str = fileread(filename);
cstr = strsplit(str,'\n');
cstr(end) = []; % ignore newline at EOF
for i = 1:length(cstr)
  cstr{i} = strsplit(cstr{i},',');
end
carray = cell(length(cstr),length(cstr{1}));
for i = 1:length(cstr)
  carray(i,:) = cstr{i};
end

% Remove not required columns
idxFrom = @(str) find (strcmp (str, carray(1,:)));
carray(:,[idxFrom('Problem'), idxFrom('Equations'), idxFrom('Variables'), ...
  idxFrom('fU'), idxFrom('fL')]) = [];
idxFrom = @(str) find (strcmp (str, carray(1,:)));

% Compute tu/ts and tl/ts
ts = cellfun (@str2double, carray(2:end, idxFrom('ts')));
tu = cellfun (@str2double, carray(2:end, idxFrom('tu')));
tl = cellfun (@str2double, carray(2:end, idxFrom('tl')));
tu = tu ./ ts;
tl = tl ./ ts;
carray(2:end, idxFrom('tu')) = num2cell(tu);
carray(2:end, idxFrom('tl')) = num2cell(tl);
carray{1,     idxFrom('tu')} = 'tuts';
carray{1,     idxFrom('tl')} = 'tlts';

% Convert mufUfL to double
mufUfL = cellfun(@str2double,carray(2:end,idxFrom('mufUfL')));
carray(2:end,idxFrom('mufUfL')) = num2cell(mufUfL);

% Remove not required columns
carray(:,idxFrom('ts')) = [];
idxFrom  = @(strfind) find(strcmp(strfind, carray(1,:)));
rowsFrom = @(strfind) find(strcmp(strfind, carray(:,idxFrom('Solver'))));

% Create result table for each solver
solver = unique(carray(2:end,idxFrom('Solver')));
resultTable = cell(6,length(solver));
resultTable{1,1} = '               ';
resultTable{2,1} = 'med (mu(fU,fL))';
resultTable{3,1} = 'max (mu(fU,fL))';
resultTable{4,1} = 'min (mu(fU,fL))';
resultTable{5,1} = 'med (tu/ts)    ';
resultTable{6,1} = 'med (tl/ts)    ';
for i = 1:length(solver)
  mufUfL = cell2mat(carray(rowsFrom(solver{i}),idxFrom('mufUfL')));
  mufUfL(isnan(mufUfL)) = [];
  tuts = cell2mat(carray(rowsFrom(solver{i}),idxFrom('tuts')));
  tlts = cell2mat(carray(rowsFrom(solver{i}),idxFrom('tlts')));
  resultTable{1,i+1} = sprintf('%-8s',solver{i});
  resultTable{2,i+1} = sprintf('%-8.2e',median(mufUfL));
  resultTable{3,i+1} = sprintf('%-8.2e',max(mufUfL));
  resultTable{4,i+1} = sprintf('%-8.2e',min(mufUfL));
  resultTable{5,i+1} = sprintf('%-8.2f',median(tuts));
  resultTable{6,i+1} = sprintf('%-8.2f',median(tlts));
end

% Join output table
str = cell(size(resultTable,1),1);
for i = 1:length(str)
  str{i} = strjoin(resultTable(i,:), '    ');
end
str = strjoin(str, '\n');
end
