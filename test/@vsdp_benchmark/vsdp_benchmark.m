classdef vsdp_benchmark < handle
  % VSDP_BENCHMARK Benchmark class for VSDP.
  %
  %   See also vsdp.
  %
  
  % Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)
  
  properties
    % A structure array with cached values about the VSDP benchmark:
    %
    %   lib     Benchmark library name.
    %   name    Short name for the test case.
    %   file    System dependend full path of the original test case data.
    %   m       Number of constraints.
    %   n       Number of cone variables.
    %   K_f     Has free variables.
    %   K_l     Has LP   variables.
    %   K_q     Has SOCP variables.
    %   K_s     Has SDP  variables.
    %   values  A structure array with solutions for each solver:
    %     name  Name of the approximate solver.
    %     fp    Primal approximate solution.
    %     fd    Dual   approximate solution.
    %     ts    Time for approximate solution.
    %     fL    Rigorous lower bound.
    %     tL    Time for rigorous lower bound.
    %     fU    Rigorous upper bound.
    %     tU    Time for rigorous upper bound.
    %
    % The data corresponding to the cached solution is for test case "j" and
    % solver "i":
    %
    %   % Test case data to create a VSDP object in various formats.
    %   obj.BENCHMARK(j).file
    %
    % The remaining data files are located inside subdirectory "data" of the
    % directory "RESULT_DIR", identified by a common ID prefix for "i" and "j":
    %
    %   ID = fprintf('%s_%s_%s', obj.BENCHMARK(j).lib, ...
    %     obj.BENCHMARK(j).name, obj.BENCHMARK(j).values(i).name);
    %
    %   % Data for "vsdp.add_solution('approximate', ...);".
    %   fprintf('%s_approximate_solution.mat', ID);
    %
    %   % Data for "vsdp.add_solution('Rigorous lower bound', ...);".
    %   fprintf('%s_%s_%s_rigorous_lower_bound.mat', ID);
    %
    %   % Data for "vsdp.add_solution('Rigorous upper bound', ...);".
    %   fprintf('%s_%s_%s_rigorous_upper_bound.mat', ID);
    %
    % Default: [].
    BENCHMARK
    
    % A structure array of initialized solver strings.
    %
    % The benchmark programm assumes, that all solver
    %
    % Default: structure with fields
    %
    %   name    = {'intlab', 'vsdp', 'csdp', 'mosek', 'sdpa', 'sdpt3', 'sedumi'}
    %   check_fun = { for each 'name' a function string to check functionality }
    %   setup_dir = { for each 'name' a setup directory       }
    %   setup_fun = { for each 'name' a setup function string }
    %
    % Default: [].
    SOLVER
    
    % Absolute directory path for persistent data storage.
    %
    % Default: Subdirectory 'result' of current directory.
    RESULT_DIR
    
    % Temporary absolute directory path for data storage.
    %
    % Default: Output of 'tempname()'.
    TMP_DIR
  end
  
  methods
    function obj = vsdp_benchmark (dir)
      if (nargin > 0)
        obj.RESULT_DIR = dir;
      else
        obj.RESULT_DIR = 'result';
      end
      
      % Create temprorary directory.
      obj.TMP_DIR = tempname ();
      
      % Add directory of this class to load path.
      addpath (fileparts (fileparts (mfilename ('fullpath'))));
    end
    
    
    function set.RESULT_DIR (obj, p)
      % Create result directory with "data" subdirectory.
      if (~isdir (fullfile (p, 'data')))
        mkdir (fullfile (p, 'data'));
      end
      obj.RESULT_DIR = obj.check_dir (p);
      obj.load_state ();
    end
    
    
    function set.TMP_DIR (obj, p)
      if (~isdir (p))
        mkdir (p);
      end
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
