classdef vsdp_benchmark
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
  end
  
  
  properties (GetAccess = public, SetAccess = protected)
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
    function obj = vsdp_benchmark (b_dir, r_dir, sol, t_dir)
      % VSDP_BENCHMARK Constructor for VSDP benchmarks.
      %
      %   obj = vsdp_benchmark (b_dir, r_dir)
      %
      %      Where 'b_dir' is the benchmark directory (obj.BENCHMARK_DIR) of
      %
      %         https://github.com/vsdp/vsdp.github.io.
      %
      %      and 'r_dir' the result directory (obj.RESULT_DIR).
      %
      %   obj = vsdp_benchmark (..., sol)  The thrid optional parameter sets
      %                                    the obj.SOLVER list.
      %
      %   obj = vsdp_benchmark (..., t_dir)  The fourth optional parameter sets
      %                                      the obj.TMP_DIR directory.
      %
      %   See also vsdp, vsdp_benchmark_exporter.
      %
      
      narginchk (2, 4);
      obj.BENCHMARK_DIR = obj.check_dir (b_dir);
      obj.RESULT_DIR    = obj.check_dir (r_dir);
      [~,~] = mkdir (fullfile (obj.RESULT_DIR, 'data'));
      [~,~] = mkdir (fullfile (obj.RESULT_DIR, 'log'));
      
      if (nargin > 2)
        obj.SOLVER = sol;
      end
      if (nargin > 3)
        obj.TMP_DIR = obj.check_dir (t_dir);
      else
        tmp_dir = tempname ();
        mkdir (tmp_dir);
        obj.TMP_DIR = tmp_dir;
      end
    end
    
    
    function p = check_dir (~, p)
      % CHECK_DIR  Check 'p' to be an existing directory.
      %
      %   Return the absolute path 'p' or an emtpy char array.
      
      if (exist (p, "dir") == 7)
        p = what (p);  % Get absolute path.
        p = p.path;
      else
        error ('VSDP_BENCHMARK:check_dir:noDir', ...
          'check_dir: The directory ''%s'' does not exist', p);
      end
    end
  end
end
