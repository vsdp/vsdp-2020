classdef vsdp_benchmark_exporter
  % VSDP_BENCHMARK_EXPORTER
  %
  %   Helper class to export benchmarks.
  %
  
  % Copyright 2018-2019 Kai T. Ohlhus (kai.ohlhus@tuhh.de)
  
  properties
    % Current view of the data to be exported.
    %
    % Default: [].
    cdata_view
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
    
    % Imported test cases for formatted export.  This variable cannot be
    % modified.  Use 'obj.cdata_view' instead.
    %
    % Default: [].
    cdata
  end
  
  
  methods
    function obj = vsdp_benchmark_exporter (b_dir, r_dir, t_dir)
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
      %   obj = vsdp_benchmark (..., t_dir)  The thrid optional parameter sets
      %                                      the obj.TMP_DIR directory.
      %
      %
      %   See also vsdp, vsdp_benchmark.
      %
      
      narginchk (2, 3);
      obj.BENCHMARK_DIR = obj.check_dir (b_dir);
      obj.RESULT_DIR    = obj.check_dir (r_dir);
      
      if (nargin > 2)
        obj.TMP_DIR = obj.check_dir (t_dir);
      else
        tmp_dir = tempname ();
        mkdir (tmp_dir);
        obj.TMP_DIR = tmp_dir;
      end
    end
    
    function str = cache_file_name (obj)
      % CACHE_FILE_NAME  Returns the name of the cache file.
      %
      
      str = fullfile (obj.RESULT_DIR, 'benchmark_cache.mat');
    end
    
    function obj = create_cache (obj)
      % CREATE_CACHE  Cache the data of obj.RESULT_DIR.
      %
      
      obj.cdata = obj.get_test_cases ();
      obj = obj.get_data ();
      obj = obj.reset_view ();
      obj.save_cache ();
    end
    
    function save_cache (obj)
      cdata_cache = obj.cdata;
      save (obj.cache_file_name, '-v7', 'cdata_cache')
    end
    
    function obj = load_cache (obj)
      if (exist (obj.cache_file_name, 'file') == 2)
        load (obj.cache_file_name, 'cdata_cache');
        obj.cdata = cdata_cache;
        obj = obj.reset_view ();
      else
        obj = obj.create_cache ();
      end
    end
    
    function obj = reset_view (obj)
      obj.cdata_view = obj.cdata;
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
