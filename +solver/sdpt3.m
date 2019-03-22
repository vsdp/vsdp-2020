classdef sdpt3 < handle
  % SDPT3  Solver proxy class (not the acutal solver!).
  %
  %   For information about SDPT3, see:
  %
  %      http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
  %      https://github.com/sqlp/sdpt3
  %
  %   See also vsdp.solve.
  %
  
  % Copyright 2004-2019 Christian Jansson (jansson@tuhh.de)
  
  methods (Static)
    function obj = solve (obj, sol_type)
      % SOLVE  Approximately solve conic problem instance with SDPT3.
      %
      %   See also vsdp.solve.
      %
      
      narginchk (1, 2);
      
      solver.sdpt3.install (true);  % Show errors
      
      if (nargin == 1)
        sol_type = 'Approximate';
      end
      [A, b, c] = obj.get_midpoint_problem_data (sol_type);
      
      % Should initial solution guess be taken into account?
      if ((obj.options.USE_INITIAL_GUESS) ...
          && (~isempty (obj.solution('Initial'))))
        isol = obj.solution('Initial');
        [x0, y0, z0] = deal (isol.x, isol.y, isol.z);
      else
        [x0, y0, z0] = deal ([], [], []);
      end
      
      % Should special solver options be taken into account?
      if (~isempty (obj.options.SOLVER_OPTIONS))
        OPTIONS = obj.options.SOLVER_OPTIONS;
      else
        OPTIONS = [];
      end
      
      % Adapt output verbosity.
      if (~obj.options.VERBOSE_OUTPUT)
        if (exist ('sqlpmain.m', 'file') == 2) % if SDPT3-4.0
          OPTIONS.printlevel = 0; % default: 3
        else
          OPTIONS.printyes = 0;   % default: 1
        end
      end
      
      % Prepare data for solver.
      A = vsdp.smat (obj, A, 1);
      c = vsdp.smat (obj, c, 1);
      [blk, A, c, b] = read_sedumi (A, b, c, obj.K);
      
      % Call solver.
      tic;
      [~, x, y, z, info] = sqlp (blk, A, c, b, OPTIONS, x0, y0, z0);
      solver_info.elapsed_time = toc;
      
      % Store solution.
      x = vsdp.svec (obj, vsdp.cell_sub_blocks (x(:), blk), 2);
      z = vsdp.svec (obj, vsdp.cell_sub_blocks (z(:), blk), 1);
      f_objective = [obj.c'*x; obj.b'*y];
      solver_info.name = 'sdpt3';
      if (isstruct (info))
        termcode = info.termcode;  % SDPT3-4.0 output
      else
        termcode = info(1);  % SDPT3-3.x output
      end
      switch (termcode)
        case 0
          solver_info.termination = 'Normal termination';
        case 1
          solver_info.termination = 'Primal infeasible';
        case 2
          solver_info.termination = 'Dual infeasible';
        otherwise
          solver_info.termination = 'Unknown';
      end
      
      obj.add_solution (sol_type, x, y, z, f_objective, solver_info);
    end
    
    function [f,l,q,s] = supported_cones ()
      f = true;  % free   variables.
      l = true;  % linear variables.
      q = true;  % second-order cones.
      s = true;  % semidefinite cones.
    end
    
    function spath = install (varargin)
      % Returns the path to the installed and usable solver.  Otherwise return
      % an empty array.  No error messages are thrown.
      %
      % By passing one or more arguments interactive installation actions
      % happen and, in case of failures, error messages are thrown.
      %
      
      sname          = 'sdpt3';
      is_available   = @() exist ('sqlp', 'file') == 2;
      get_path       = @() fileparts (which ('install_sdpt3'));
      installer_file = 'install_sdpt3.m';
      do_error       = false;
      spath = solver.registry.generic_install (sname, is_available, ...
        get_path, installer_file, do_error);
      
      % Return on sucess or non-interactive silent mode.
      if ((nargin == 0) || (~isempty (spath)))
        return;
      end
      
      % Interactive mode: Install solver from the internet.
      fprintf ('\n  Unable to find the SDPT3 solver (GPLv2 license).\n\n');
      url = 'https://github.com/sqlp/sdpt3/archive/master.zip';
      prompt = sprintf ('  Install from ''%s'' [y/n]? ', url);
      str = input (prompt, 's');
      if (str(1) ~= 'y' && str(1) ~= 'Y')
        error ('VSDP:SOLVER:sdpt3:installError', ...
          'solver.sdpt3.install: sdpt3 is not available.');
      end
      download_path = fullfile (vsdp.settings ('vsdp', 'path'), 'download');
      mkdir (download_path);
      zip_file = fullfile (download_path, 'sdpt3.zip');
      urlwrite (url, zip_file);
      if (exist (zip_file, 'file') ~= 2)
        error ('VSDP:SOLVER:sdpt3:installError', ...
          'solver.sdpt3.install: Could not download the sdpt3 installer');
      end
      unzip (zip_file, download_path);
      delete (zip_file);
      spath = fullfile (download_path, sname);
      movefile (fullfile (download_path, 'sdpt3-master'), spath);
      run (fullfile (spath, installer_file));
      % Final check.
      if (~is_available ())
        error ('VSDP:SOLVER:sdpt3:installError', ...
          'solver.sdpt3.install: SDPT3 is not available.');
      end
      spath = vsdp.settings (sname, 'path', spath);
    end
  end
end
