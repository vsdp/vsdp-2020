classdef linprog < handle
  % LINPROG  Solver proxy class (not the acutal solver!).
  %
  %   For more information about LINPROG, see:
  %
  %      https://www.mathworks.com/help/optim/ug/linprog.html
  %
  %   See also vsdp.solve.
  %
  
  % Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)
  
  methods (Static)
    function obj = solve (obj, sol_type)
      % SOLVE  Approximately solve conic problem instance with LINPROG.
      %
      %   See also vsdp.solve.
      %
      
      narginchk (1, 2);
      
      solver.linprog.install (true);  % Show errors
      
      % Check cones.
      if ((sum (obj.K.q) > 0) || (sum (obj.K.s) > 0))
        error ('VSDP:solve:unsupportedCones', ...
          ['solve_linprog: Second-order cones (K.q) and semidefinite ', ...
          'cones (K.s) are not supported by LINPROG.']);
      end
      
      if (nargin == 1)
        sol_type = 'Approximate';
      end
      [A, b, c] = obj.get_midpoint_problem_data (sol_type);
      
      % Should initial solution guess be taken into account?
      if ((obj.options.USE_INITIAL_GUESS) ...
          && (~isempty (obj.solution('Initial'))))
        x0 = full (obj.solution('Initial').x);
      else
        x0 = [];
      end
      
      % Should special solver options be taken into account?
      options = optimoptions ('linprog', 'Algorithm', 'interior-point-legacy');
      if (~isempty (obj.options.SOLVER_OPTIONS))
        options = optimoptions (options, obj.options.SOLVER_OPTIONS);
      end
      
      % Adapt output verbosity.
      if (~obj.options.VERBOSE_OUTPUT)
        options = optimoptions (options, 'Display', 'off');
      end
      
      % Prepare data for solver.
      [A, b, c] = deal (A', full (b), full (c));
      lbound = [ ...
        -inf(obj.K.f, 1); ...
        zeros(obj.K.l, 1)];       % lower bound
      ubound = inf(length(c),1);  % upper bound
      
      % Call solver.
      tic;
      [x, ~, flag, ~, lambda] = linprog ...
        (c, [], [], A, b, lbound, ubound, x0, options);
      solver_info.elapsed_time = toc;
      
      % Store solution.
      if (isfield (lambda, 'eqlin'))
        y = -lambda.eqlin;
        z = c - A'*y;
      end
      f_objective = [obj.c'*x; obj.b'*y];
      solver_info.name = 'linprog';
      switch (flag)
        case 1
          solver_info.termination = 'Normal termination';
        case -2
          solver_info.termination = 'Primal infeasible';
        case -3
          solver_info.termination = 'Dual infeasible';
        case -5
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
      
      sname          = 'linprog';
      is_available   = @() exist (sname, 'file') == 2;
      get_path       = @() fileparts (which (sname));
      installer_file = [];
      do_error       = (nargin > 0);
      spath = solver.registry.generic_install (sname, is_available, ...
        get_path, installer_file, do_error);
    end
  end
end
