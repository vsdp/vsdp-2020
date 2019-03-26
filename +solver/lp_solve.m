classdef lp_solve < handle
  % LP_SOLVE  Solver proxy class (not the acutal solver!).
  %
  %   For more information about LP_SOLVE, see:
  %
  %      http://lpsolve.sourceforge.net/5.5/index.htm
  %
  %   See also vsdp.solve.
  %
  
  % Copyright 2004-2019 Christian Jansson (jansson@tuhh.de)
  
  methods (Static)
    function obj = solve (obj, sol_type)
      % SOLVE  Approximately solve conic problem instance with LP_SOLVE.
      %
      %   See also vsdp.solve.
      %
      
      narginchk (1, 2);
      
      solver.lp_solve.install (true);                   % Show errors
      solver.registry.check_cones (obj, 'lp_solve', 1); % Show errors
      
      if (nargin == 1)
        sol_type = 'Approximate';
      end
      [A, b, c] = obj.get_midpoint_problem_data (sol_type);
      
      % Should initial solution guess be taken into account?
      if (obj.options.USE_INITIAL_GUESS)
        warning ('VSDP:solve_lp_solve:ignoreInitialGuess', ...
          ['solve_lp_solve: LP_SOLVE does not support initial guesses ' ...
          '(x0,y0,z0) and proceeds without them.']);
      end
      
      % Prepare data for solver.
      [A, b, c] = deal (A', full (b), full (c));
      lbound = [ ...
        -inf(obj.K.f, 1); ...
        zeros(obj.K.l, 1)];           % lower bounds
      vtypes = zeros (length (b), 1); % variable types: 0 == equality
      
      % Call solver.
      tic;
      [~, x, y, stat] = lp_solve (-c, A, b, vtypes, lbound);
      solver_info.elapsed_time = toc;
      
      % Store solution.
      y = -y;
      z = c - A'*y;
      f_objective = [obj.c'*x; obj.b'*y];
      solver_info.name = 'lp_solve';
      switch (stat)
        case 0
          solver_info.termination = 'Normal termination';
        case 2
          solver_info.termination = 'Primal infeasible';
        case 3
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
      s = false; % semidefinite cones.
    end
    
    function spath = install (varargin)
      % Returns the path to the installed and usable solver.  Otherwise return
      % an empty array.  No error messages are thrown.
      %
      % By passing one or more arguments interactive installation actions
      % happen and, in case of failures, error messages are thrown.
      %
      
      sname          = 'lp_solve';
      is_available   = @() exist (sname, 'file') == 2;
      get_path       = @() fileparts (which (sname));
      installer_file = [];
      do_error       = (nargin > 0);
      spath = solver.registry.generic_install (sname, is_available, ...
        get_path, installer_file, do_error);
    end
  end
end
