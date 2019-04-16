classdef vsdp_benchmark_exporter
  % VSDP_BENCHMARK_EXPORTER
  %
  %   Helper class to export benchmarks.
  %
  
  % Copyright 2018-2019 Kai T. Ohlhus (kai.ohlhus@tuhh.de)
  
  properties
    cdata_view
  end
  
  properties (GetAccess = public, SetAccess = protected)
    data_dir
    cdata
  end
  
  properties (Access = protected)
    tmp_dir
  end
  
  methods
    function obj = vsdp_benchmark_exporter (data_dir, use_cache)
      % EXPORTER Constructor.
      %
      %    obj = vsdp_benchmark_exporter (data_dir)  Reads all the benchmarks
      %                                in data_dir and creates a cache in the
      %                                directory above.
      %
      %    obj = vsdp_benchmark_exporter (data_dir, false)  Same as above, but
      %                                does not use any cache.
      
      if (~exist (data_dir, 'dir'))
        error ('VSDP:vsdp_benchmark_exporter', ...
          'vsdp_benchmark_exporter: Directory ''%s'' does not exist.', ...
          data_dir);
      end
      obj.data_dir = data_dir;
      obj.tmp_dir = tempname ();
      mkdir (obj.tmp_dir);
      
      % Use cache if applicable.
      if (((nargin > 1) && use_cache) ...
          || (exist (obj.cache_file_name, 'file') == 2))
        obj = obj.load_cache ();
      else
        obj = obj.create_cache ();
      end
    end
    
    function str = cache_file_name (obj)
      str = fullfile (obj.data_dir, '..', 'benchmark_cache.mat');
    end
    
    function obj = create_cache (obj)
      % CREATE_CACHE  Cache the data of obj.data_dir.
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
      load (obj.cache_file_name, 'cdata_cache');
      obj.cdata = cdata_cache;
      obj = obj.reset_view ();
    end
    
    function obj = reset_view (obj)
      obj.cdata_view = obj.cdata;
    end
  end
end
