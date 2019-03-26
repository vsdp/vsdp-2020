classdef glpk < handle
  % GLPK  Solver proxy class (not the acutal solver!).
  %
  %   For more information about GLPK, see:
  %
  %      https://www.gnu.org/software/glpk/
  %
  %   See also vsdp.solve.
  %
  
  % Copyright 2004-2019 Christian Jansson (jansson@tuhh.de)
  
  methods (Static)
    function obj = solve (obj, sol_type)
      % SOLVE  Approximately solve conic problem instance with GLPK.
      %
      %   See also vsdp.solve.
      %
      
      narginchk (1, 2);
      
      solver.glpk.install (true);                   % Show errors
      solver.registry.check_cones (obj, 'glpk', 1); % Show errors
      
      if (nargin == 1)
        sol_type = 'Approximate';
      end
      [A, b, c] = obj.get_midpoint_problem_data (sol_type);
      
      % Should initial solution guess be taken into account?
      if (obj.options.USE_INITIAL_GUESS)
        warning ('VSDP:solve_glpk:ignoreInitialGuess', ...
          ['solve_glpk: GLPK does not support initial guesses (x0,y0,z0) ' ...
          'and proceeds without them.']);
      end
      
      % Should special solver options be taken into account?
      if (~isempty (obj.options.SOLVER_OPTIONS))
        param = obj.options.SOLVER_OPTIONS;
      else
        param = [];
      end
      
      % Adapt output verbosity.
      if (~obj.options.VERBOSE_OUTPUT)
        param.msglev = 0;
      end
      
      % Prepare data for solver.
      [A, b, c] = deal (A', full (b), full (c));
      
      sense = 1;  % Minimization
      lbound = [ ...
        -inf( obj.K.f, 1); ...
        zeros(obj.K.l, 1)];                 % lower bound
      ubound = inf (length (c), 1);         % upper bound
      vtype = repmat ('C', length (c), 1);  % variable   types: continuous
      ctype = repmat ('S', length (b), 1);  % constraint types: equality
      
      % Call solver.
      tic;
      [x, ~, errnum, extra] = glpk ...
        (c, A, b, lbound, ubound, ctype, vtype, sense, param);
      solver_info.elapsed_time = toc;
      
      % Store solution.
      if (isfield (extra, 'lambda'))
        y = extra.lambda;
        z = obj.c - obj.At*y;
      end
      f_objective = [obj.c'*x; obj.b'*y];
      solver_info.name = 'glpk';
      switch (errnum)
        case 0
          solver_info.termination = 'Normal termination';
        case 10
          solver_info.termination = 'Primal infeasible';
        case 11
          solver_info.termination = 'Dual infeasible';
        case 15
          solver_info.termination = 'Primal and dual infeasible';
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
      
      sname          = 'glpk';
      is_available   = @() exist (sname, 'file') == 2;
      get_path       = @() fileparts (which (sname));
      installer_file = [];
      do_error       = (nargin > 0);
      spath = solver.registry.generic_install (sname, is_available, ...
        get_path, installer_file, do_error);
    end
  end
end
