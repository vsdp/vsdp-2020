function obj = add_benchmark (obj, name, dir_pattern, name_fun)
% ADD_BENCHMARK  Generates a list of jobs for a benchmark.
%
%   name         The name of the benchmark library.
%   dir_pattern  A pattern for the 'dir()' function to extract the test cases.
%   name_fun     Function to extract the test case name from the file name.
%
%   If obj.BENCHMARK_DIR is set to the 'benchmarks' directory of the repository
%   from  https://github.com/vsdp/vsdp.github.io, one only has to provide the
%   name of the benchmark, that is 'DIMACS', 'ESC', 'RDM', 'SDPLIB',
%   'SPARSE_SDP', or 'all'.
%
%   See also vsdp_benchmark.
%

% Copyright 2004-2019 Christian Jansson (jansson@tuhh.de)

% Recursive call with well-known benchmarks.
if (nargin == 2)
  if (isempty (obj.BENCHMARK_DIR))
    error ('VSDP_BENCHMARK:add_benchmark:benchmarkDirUnknown', ...
      'add_benchmark: Please set ''obj.BENCHMARK_DIR''.');
  end
  switch (name)
    case 'all'
      obj.add_benchmark ('DIMACS');
      obj.add_benchmark ('ESC');
      obj.add_benchmark ('RDM');
      obj.add_benchmark ('SDPLIB');
      obj.add_benchmark ('SPARSE_SDP');
    case 'DIMACS'
      obj.add_benchmark ('DIMACS', ...
        fullfile (obj.BENCHMARK_DIR, 'DIMACS', 'data', '**', '*.mat.gz'), ...
        @(str) str(1:end - length('.mat.gz')));
    case 'ESC'
      S = warning ('off', 'VSDP_BENCHMARK:add_benchmark:notUniqueNames');
      obj.add_benchmark ('ESC', ...
        fullfile (obj.BENCHMARK_DIR, 'ESC', 'data', '*.dat-s.gz'), ...
        @(str) strtok (str, "_"));
      % Fix not unique test cases 'CH2'.
      idx = find (strcmp (obj.BENCHMARK(:,2), 'CH2'));
      for i = 1:length (idx)
        [~,str] = fileparts (obj.BENCHMARK{idx(i),3});
        str = strtok (str, 'STO');
        obj.BENCHMARK{idx(i),2} = str(1:end-1);
      end
      warning (S);
    case 'RDM'
      obj.add_benchmark ('RDM', ...
        fullfile (obj.BENCHMARK_DIR, 'RDM', 'data', '*.dat-s.gz'), ...
        @(str) strtok (str, "."));
    case 'SDPLIB'
      obj.add_benchmark ('SDPLIB', ...
        fullfile (obj.BENCHMARK_DIR, 'SDPLIB', 'data', '*.dat-s'), ...
        @(str) str(1:end - length('.dat-s')));
    case 'SPARSE_SDP'
      obj.add_benchmark ('SPARSE_SDP', ...
        fullfile (obj.BENCHMARK_DIR, 'SPARSE_SDP', 'data', '*.dat-s.gz'), ...
        @(str) str(1:end - length('.dat-s.gz')));
  end
  return;
end

list = {};
f = dir (dir_pattern);
for i = length(f):-1:1
  list(i,:) = {name, name_fun(f(i).name), fullfile(f(i).folder, f(i).name)};
end
obj.BENCHMARK = [obj.BENCHMARK; list];

% Check for redundant names.
idx = strcat (obj.BENCHMARK(:,1), obj.BENCHMARK(:,2));
if (length (idx) > length (unique (idx)))
  warning ('VSDP_BENCHMARK:add_benchmark:notUniqueNames', ...
    'add_benchmark: Names in library ''%s'' are not unique.', name);
end
end
