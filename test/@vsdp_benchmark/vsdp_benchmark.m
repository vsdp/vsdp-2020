classdef vsdp_benchmark < handle
  % VSDP_BENCHMARK Benchmark class for VSDP.
  %
  %   See also vsdp.
  %
  
  % Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)
  
  properties
    % A cell array of strings with VSDP benchmarks to run.  The columns are:
    %
    %   obj.BENCHMARK(:,1)  Benchmark library name.
    %   obj.BENCHMARK(:,2)  Short name for the test case.
    %   obj.BENCHMARK(:,3)  System dependend full path of the original test
    %                       case data.
    %
    %  Default: cell (0, 3).
    BENCHMARK = cell (0, 3);
    
    % A cell array of strings with solvers to use for this benchmark.
    %
    % Default: {}.
    SOLVER = {}
    
    % Absolute path to the benchmark directory of
    %
    %   https://github.com/vsdp/vsdp.github.io.
    %
    % Default: [].
    BENCHMARK_DIR
    
    % Absolute directory path for persistent data storage.
    %
    % Default: [].
    RESULT_DIR
    
    % Temporary absolute directory path for data storage.
    %
    % Default: Output of 'tempname()'.
    TMP_DIR
  end
  
  methods
    function obj = vsdp_benchmark ()
      % Create temprorary directory.
      tmp_dir = tempname ();
      mkdir (tmp_dir);
      obj.TMP_DIR = tmp_dir;
    end
    
    
    function set.BENCHMARK_DIR (obj, p)
      obj.BENCHMARK_DIR = obj.check_dir (p);
    end
    
    
    function set.RESULT_DIR (obj, p)
      obj.RESULT_DIR = obj.check_dir (p);
      [~,~] = mkdir (fullfile (p, 'data'));
      [~,~] = mkdir (fullfile (p, 'log'));
    end
    
    
    function set.TMP_DIR (obj, p)
      obj.TMP_DIR = obj.check_dir (p);
    end
    
    
    function p = check_dir (~, p)
      % CHECK_DIR  Check 'p' to be an existing directory.
      %
      %   Return the absolute path 'p' or an emtpy char array.
      
      if (exist (p, "dir") == 7)
        p = what (p);  % Get absolute path.
        p = p.path;
      else
        warning ('VSDP_BENCHMARK:check_dir:noDir', ...
          'check_dir: The directory ''%s'' does not exist', p);
        p = '';
      end
    end
  end
end
