classdef csdp < handle
  % CSDP  Solver proxy class (not the acutal solver!).
  %
  %   For more information about CSDP, see:
  %
  %      https://projects.coin-or.org/Csdp/
  %      https://github.com/coin-or/Csdp
  %
  %   See also vsdp.solve.
  %
  
  % Copyright 2004-2019 Christian Jansson (jansson@tuhh.de)
  
  methods (Static)
    function obj = solve (obj, sol_type)
      % SOLVE  Approximately solve conic problem instance with CSDP.
      %
      %   See also vsdp.solve.
      %
      
      narginchk (1, 2);
      
      solver.csdp.install (1);  % Show errors
      
      % Check cones.
      if (sum (obj.K.q) > 0)
        error ('VSDP:solve_csdp:unsupportedCone', ...
          'solve_csdp: Second-order cones (K.q) are not supported by CSDP.');
      end
      
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
        pars = obj.options.SOLVER_OPTIONS;
      else
        pars = [];
      end
      
      % Adapt output verbosity.
      if (~obj.options.VERBOSE_OUTPUT)
        pars.printlevel = 0;
      end
      
      % Prepare data for solver.
      [b, c, x0, y0, z0] = deal (full (b), full (c), full (x0), full (y0), ...
        full (z0));
      
      % Convert to SeDuMi-Format (same as CSDP format).
      A = vsdp.smat (obj, A, 1);
      c = vsdp.smat (obj, c, 1);
      K = obj.K;
      if (K.f > 0)
        warning('VSDP:solve_csdp:unsupportedCone', ...
          ['solve_csdp: CSDP supports free variables (K.f) by converting ' ...
          'them to the difference of positive variables.  The resulting ', ...
          'problem is ill-posed.']);
        [A, b, c, K] = convertf (A, b, c, K);  % CSDP function.
      end
      
      args = {A, b, c, K, pars};
      % CSDP requires all initial solutions to be not empty!
      if (~isempty (x0) && ~isempty(y0) && ~isempty (z0))
        args = [args, {x0, y0, z0}];
      end
      
      % If using Windows, switch to the 'csdp.exe' directory.
      if (ispc ())
        p = fileparts (which ('csdp.exe'));
        p = cd (p);
      end
      
      % Call solver.
      tic;
      [x, y, z, INFO] = csdp (args{:});
      solver_info.elapsed_time = toc;
      
      if (ispc ())
        cd (p);
      end
      
      % Store solution.
      x = vsdp.svec (obj, x, 2);
      z = vsdp.svec (obj, z, 1);
      f_objective = [obj.c'*x; obj.b'*y];
      solver_info.name = 'csdp';
      switch (INFO)
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
      q = false; % second-order cones.
      s = true;  % semidefinite cones.
    end
    
    function spath = install (varargin)
      % Returns the path to the installed and usable solver.  Otherwise return
      % an empty array.  No error messages are thrown.
      %
      % By passing one or more arguments interactive installation actions
      % happen and, in case of failures, error messages are thrown.
      %
      
      sname          = 'csdp';
      is_available   = @() exist (sname, 'file') == 2;
      get_path       = @() fileparts (which (sname));
      installer_file = [];
      do_error       = (nargin > 0);
      spath = solver.registry.generic_install (sname, is_available, ...
        get_path, installer_file, do_error);
    end
  end
end
