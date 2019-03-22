classdef sedumi < handle
  % SEDUMI  Solver proxy class (not the acutal solver!).
  %
  %   For information about SeDuMi, see:
  %
  %      http://sedumi.ie.lehigh.edu/
  %      https://github.com/sqlp/sedumi
  %
  %   See also vsdp.solve.
  %
  
  % Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)
  
  methods (Static)
    function obj = solve (obj, sol_type)
      % SOLVE  Approximately solve conic problem instance with SeDuMi.
      %
      %   See also vsdp.solve.
      %
      
      narginchk (1, 2);
      
      % Check solver availability.
      if (exist ('sedumi', 'file') ~= 2)
        error ('VSDP:solve_sedumi:notAvailable', ...
          ['solve_sedumi: SeDuMi does not seem to be ready.  ', ...
          'Did you run ''install_sedumi()'' inside the solver directory?\n\n', ...
          'To select another solver, run:  %s.solve()'], inputname(1));
      end
      
      if (nargin == 1)
        sol_type = 'Approximate';
      end
      [A, b, c] = obj.get_midpoint_problem_data (sol_type);
      
      % Should initial solution guess be taken into account?
      if (obj.options.USE_INITIAL_GUESS)
        warning ('VSDP:solve_sedumi:ignoreInitialGuess', ...
          ['solve_sedumi: SeDuMi does not support initial guesses (x0,y0,z0) ' ...
          'and proceeds without them.']);
      end
      
      % Should special solver options be taken into account?
      if (~isempty (obj.options.SOLVER_OPTIONS))
        pars = obj.options.SOLVER_OPTIONS;
      else
        pars = [];
      end
      
      % Adapt output verbosity.
      if (~obj.options.VERBOSE_OUTPUT)
        pars.fid = 0;
      end
      
      % Prepare data for solver.
      A = vsdp.smat (obj, A, 1);
      c = vsdp.smat (obj, c, 1);
      
      % Call solver.
      tic;
      [x, y, info] = sedumi (A, b, c, obj.K, pars);
      solver_info.elapsed_time = toc;
      
      % Store solution.
      x = vsdp.svec (obj, x, 2);
      z = vsdp.svec (obj, c - A*y, 1);
      f_objective = [obj.c'*x; obj.b'*y];
      solver_info.name = 'sedumi';
      switch (info.pinf + 2*info.dinf)
        case 0
          solver_info.termination = 'Normal termination';
        case 1
          solver_info.termination = 'Primal infeasible';
        case 2
          solver_info.termination = 'Dual infeasible';
        case 3
          solver_info.termination = 'Primal and dual infeasible';
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
      
      sname          = 'sedumi';
      is_available   = @() exist (sname, 'file') == 2;
      get_path       = @() fileparts (which ('install_sedumi'));
      installer_file = 'install_sedumi.m';
      do_error       = false;
      spath = solver.registry.generic_install (sname, is_available, ...
        get_path, installer_file, do_error);
      
      % Return on sucess or non-interactive silent mode.
      if ((nargin == 0) || (~isempty (spath)))
        return;
      end
      
      % Interactive mode: Install solver from the internet.
      fprintf ('\n  Unable to find the SeDuMi solver (GPLv2 license).\n\n');
      url = 'https://github.com/sqlp/sedumi/archive/master.zip';
      prompt = sprintf ('  Install from ''%s'' [y/n]? ', url);
      str = input (prompt, 's');
      if (str(1) ~= 'y' && str(1) ~= 'Y')
        error ('VSDP:SOLVER:sedumi:installError', ...
          'solver.sedumi.install: SeDuMi is not available.');
      end
      download_path = fullfile (vsdp.settings ('vsdp', 'path'), 'download');
      [~,~] = mkdir (download_path);
      zip_file = fullfile (download_path, 'sedumi.zip');
      urlwrite (url, zip_file);
      if (exist (zip_file, 'file') ~= 2)
        error ('VSDP:SOLVER:sedumi:installError', ...
          'solver.sedumi.install: SeDuMi is not available.');
      end
      unzip (zip_file, download_path);
      delete (zip_file);
      spath = fullfile (download_path, sname);
      movefile (fullfile (download_path, 'sedumi-master'), spath);
      run (fullfile (spath, installer_file));
      % Final check.
      if (~is_available ())
        error ('VSDP:SOLVER:sedumi:installError', ...
          'solver.sedumi.install: SeDuMi is not available.');
      end
      spath = vsdp.settings (sname, 'path', spath);
    end
  end
end
