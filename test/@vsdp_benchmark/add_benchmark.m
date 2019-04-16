function add_benchmark (obj, name, dir_pattern, name_fun)
% ADD_BENCHMARK  Generates a list of jobs for a benchmark library.
%
%   name         The name of the benchmark library.
%   dir_pattern  A pattern for the 'dir()' function to extract the test cases.
%   name_fun     Function to extract the test case name from the file name.
%
%   See also vsdp_benchmark.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

f = dir (dir_pattern);
for i = length(f):-1:1
  list(i).lib   = name;
  list(i).name  = name_fun (f(i).name);
  list(i).file  = fullfile (f(i).folder, f(i).name);
  list(i).m = [];
  list(i).n = [];
  list(i).K_f = [];
  list(i).K_l = [];
  list(i).K_q = [];
  list(i).K_s = [];
  list(i).values = [];
end

% Check for redundant names.
if (length ({list.name}) > length (unique ({list.name})))
  warning ('VSDP_BENCHMARK:add_benchmark:notUniqueNames', ...
    'add_benchmark: Names in library ''%s'' are not unique.', name);
end
obj.BENCHMARK = [obj.BENCHMARK, list];

end
