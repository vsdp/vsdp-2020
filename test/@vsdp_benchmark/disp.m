function disp (obj, item)
% DISP  Display short information about VSDP benchmark object.
%
%   obj.disp()  Display an overview about the benchmark data.
%
%   obj.disp(item)  With a second argument 'item' display only information
%                   about one of 'system', 'solver', or 'benchmark'.
%
%   See also vsdp_benchmark.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

if (nargin > 1)
  item = validatestring (item, {'system', 'solver', 'benchmark'});
else
  item = 'all';
end

switch (item)
  case 'all'
    fprintf ('  VSDP benchmark\n');
    fprintf ('  --------------\n\n');
    fprintf ('%16s: ''%s''\n', 'RESULT_DIR', obj.RESULT_DIR);
    fprintf ('%16s: ''%s''\n\n',  'TMP_DIR', obj.TMP_DIR);
    disp_system ();
    fprintf ('\n');
    disp_solver (obj);
    fprintf ('\n');
    disp_benchmark (obj);
    fprintf ('\n\n');
  otherwise
    eval (sprintf ('disp_%s(obj)', item));
end
end


function disp_system (~)
% Information about the system in use.

% Very generic system information.
fprintf ('  System information:\n\n');
if (exist ('OCTAVE_VERSION', 'builtin'))
  fprintf ('    GNU Octave Version %s\n', version ());
else
  fprintf ('    MATLAB(TM) Version %s\n', version ());
end
fprintf ('%16s: ''%s''\n', 'Architecture', computer ());

% On a Unix (Linux) system, one is able to find out more about the used system.
if (isunix ())
  % Linux kernel information.
  [~, sys] = system ('uname -a');
  fprintf ('%16s: ''%s\n', '"uname -a"', sys(1:60));
  fprintf ('%19s%s''\n', '', sys(61:end-1));
  
  % RAM size.
  [~, ram] = system ('cat /proc/meminfo');
  % Last part of first line is total memory.
  [~, ram] = strtok (strtok (ram, sprintf ('\n')));
  fprintf ('%16s: %s\t%s\n', 'RAM', strtrim (ram), '(/proc/meminfo)');
  
  % CPU model and number of CPU cores.
  [~, cpu] = system ('cat /proc/cpuinfo');
  cpu_model = regexp (cpu, 'model name\s*: *([^\n\r]+)', 'tokens');
  cpu_model = [cpu_model{:}];
  cpu_number = regexp (cpu, 'processor\s*: *(\d+)', 'tokens');
  cpu_number = cellfun (@str2num, [cpu_number{:}]);
  if (isequal (cpu_model{:}))
    cpu_model = sprintf ('%s', cpu_model{1});
  else
    cpu_model = sprintf ('%s [unsure]', cpu_model{1});
  end
  if (cpu_number == ((1:length(cpu_number)) - 1))
    cpu_number = sprintf ('%d cores', length(cpu_number));
  else
    cpu_number = sprintf ('%d cores [unsure]', length(cpu_number));
  end
  fprintf ('%16s: ''%s''\n', 'CPU', cpu_model);
  fprintf ('%19s%s\t%s\n', '', cpu_number, '(/proc/cpuinfo)');
end
end


function disp_solver (obj)
% Information about the solvers.

fprintf ('  Solver information:\n\n');
for i = 1:length (obj.SOLVER)
  fprintf ('%16s: ''%s''\n', obj.SOLVER(i).name, obj.SOLVER(i).setup_dir);
end
end


function disp_benchmark (obj)
% Information about the benchmarks.

fprintf ('  Benchmark information (%3d test cases):\n\n', ...
  length (obj.BENCHMARK));
if (~isempty (obj.BENCHMARK))
  bms = {obj.BENCHMARK.lib};
  bm_unique = unique(bms);
  for i = 1:length (bm_unique)
    fprintf ('%16s: %3d test cases\n', bm_unique{i}, ...
      sum (cellfun (@(x) strcmp (x, bm_unique{i}), bms)));
  end
end

end
