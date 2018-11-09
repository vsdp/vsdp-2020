classdef vsdp_solution < handle
  
  properties (SetAccess = protected, GetAccess = public)
    % Solution type.
    %
    % One of:
    %   - 'Initial',
    %   - 'Approximate',
    %   - 'Rigorous lower bound',
    %   - 'Rigorous upper bound',
    %   - 'Certificate primal infeasibility',
    %   - 'Certificate dual infeasibility'
    %
    sol_type
    x   % Primal solution quantity.
    y   % Dual   solution quantity.
    z   % Slack  or bound quantity.
    f_objective  % Primal and dual objective values or rigorous bounds.
    % Additional solver information.
    %
    % Structure with fields:
    %
    %   - 'name'          Name of the used conic solver.
    %   - 'elapsed_time'  Solver runtime in seconds.
    %   - 'termination'
    %        'Unknown'
    %        'Normal termination'
    %        'Primal infeasible'
    %        'Dual infeasible'
    %        'Primal and dual infeasible'
    %
    solver_info
  end
  
  methods
    function obj = vsdp_solution (sol_type, x, y, z, f_objective, solver_info)
      % VSDP_SOLUTION  Container class for VSDP solutions.
      %
      %   sol = vsdp_solution(sol_type, x, y, z)  Create a solution object.
      %
      %   sol = vsdp_solution(..., f_objective, solver_info)
      %         Create a solution object.
      %
      %   Example:
      %
      %      sol_type = 'Initial';
      %      x = [1 2 3]';  % Primal solution.
      %      y = [4 5]';    % Dual solution.
      %      z = [];
      %      sol = vsdp_solution (sol_type, x, y, z)
      %
      %
      %   See also vsdp.
      
      % Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)
      
      narginchk (4, 6);
      
      obj.sol_type = validatestring (sol_type, {'Initial', ...
        'Approximate', 'Rigorous lower bound', ...
        'Rigorous upper bound', 'Certificate primal infeasibility', ...
        'Certificate dual infeasibility'});
      
      if (~isempty (x))
        validateattributes (x, {'double', 'intval'}, {'vector'});
        obj.x = x(:);
      else
        obj.x = [];
      end
      if (~isempty (y))
        validateattributes (y, {'double', 'intval'}, {'vector'});
        obj.y = y(:);
      else
        obj.y = [];
      end
      if (~isempty (z))
        validateattributes (z, {'double', 'intval'}, {'vector'});
        obj.z = z(:);
      else
        obj.z = [];
      end
      if ((isequal (obj.sol_type, 'Initial') ...
          || isequal (obj.sol_type, 'Approximate')) ...
          && ~isempty (obj.x) && ~isempty (obj.z) ...
          && (length (obj.x) ~= length (obj.z)))
        error ('VSDP:vsdp_solution:sizeMissmatch', ...
          'vsdp_solution: The lengths of ''x'' and ''z'' must agree.');
      end
      
      if (nargin > 4)
        validateattributes (f_objective(:), {'double'}, {'size', [2, 1]});
        obj.f_objective = full (f_objective(:));
      else
        obj.f_objective = nan (2, 1);
      end
      if (nargin > 5)
        validateattributes (solver_info, {'struct'}, {});
        % Just use all the given fields.
        obj.solver_info = solver_info;
        % But ensure certain fields to be set.
        if (~isfield (obj.solver_info, 'name'))
          obj.solver_info.name = 'None';
        end
        if (isfield (obj.solver_info, 'elapsed_time'))
          validateattributes (obj.solver_info.elapsed_time, {'double'}, ...
            {'scalar'});
        else
          obj.solver_info.elapsed_time = nan;
        end
        if (isfield (obj.solver_info, 'termination'))
          obj.solver_info.termination = validatestring ( ...
            obj.solver_info.termination,  {'Unknown', 'Normal termination', ...
            'Primal infeasible', 'Dual infeasible', ...
            'Primal and dual infeasible'});
        else
          obj.solver_info.termination = 'Unknown';
        end
      else
        obj.solver_info = [];
      end
    end
    
    function varargout = size (obj, dim)
      % SIZE  Size of the solution instance.
      m = length (obj.y);
      n = length (obj.x);
      if (nargin == 1)
        if (nargout < 2)
          varargout = {[m, n]};
        elseif (nargout == 2)
          varargout = {m, n};
        else
          varargout(1:2) = {m, n};
          varargout(3:nargout) = {1};
        end
      else
        if (dim == 1)
          varargout = {m};
        elseif (dim == 2)
          varargout = {n};
        else
          varargout = {1};
        end
      end
    end
    
    function disp (obj)
      % DISP  Visualize the solution.
      if (~isempty (obj.solver_info))
        if (~strcmp (obj.solver_info.name, 'None'))
          solver_str = sprintf ('Solver ''%s'': ', obj.solver_info.name);
        else
          solver_str = '';
        end
        if (isfield (obj.solver_info, 'iter'))
          iter_str = sprintf (', %d iterations', obj.solver_info.iter);
        else
          iter_str = '';
        end
        fprintf ('      %s%s, %.1f seconds%s.\n\n', ...
          solver_str, obj.solver_info.termination, ...
          obj.solver_info.elapsed_time, iter_str);
      end
      switch (obj.sol_type)
        case 'Initial'
          obj.mem_info ();
        case 'Approximate'
          if (~all (isnan (obj.f_objective)))
            fprintf ('        c''*x = %.15e\n', obj.f_objective(1));
            fprintf ('        b''*y = %.15e\n', obj.f_objective(2));
          else
            obj.mem_info ();
          end
          fprintf ('\n');
        case 'Rigorous lower bound'
          fprintf ('          fL = %.15e\n',   obj.f_objective(1));
          fprintf ('\n');
        case 'Rigorous upper bound'
          fprintf ('          fU = %.15e\n',   obj.f_objective(2));
          fprintf ('\n');
        case 'Certificate primal infeasibility'
          fprintf ('      ');
          if (obj.f_objective(1))
            fprintf (['A certificate of primal infeasibility ''y'' was ', ...
              'found.\n      The conic problem is primal infeasible.\n\n']);
          else
            fprintf (['NO certificate of primal infeasibility was found.\n\n']);
          end
        case 'Certificate dual infeasibility'
          fprintf ('      ');
          if (obj.f_objective(2))
            fprintf (['A certificate of dual infeasibility ''x'' was ', ...
              'found.\n      The conic problem is dual infeasible.\n\n']);
          else
            fprintf (['NO certificate of dual infeasibility was found.\n\n']);
          end
        otherwise
          error ('VSDP_SOLUTION:disp:unknownType', ...
            'disp: Unknown solution type.');
      end
    end
    
    function mem_info (obj)
      % MEM_INFO  Display memory information.
      [x, y, z] = deal (obj.x, obj.y, obj.z);
      S = whos ('x', 'y', 'z');
      names = sprintf('%2s:\n', S.name);
      names = strsplit (names, '\n');
      names = char (names(1:end-1)');
      sizes = [S.size];
      sizes = [sizes(1:2:end)', sizes(2:2:end)'];
      d = numel (num2str (max (sizes(:,1))));
      sizes = sprintf(['[ %', num2str(d), 'd x %d ]\n'], sizes');
      sizes = strsplit (sizes, '\n');
      sizes = char (sizes(1:end-1)');
      classes = char ({S.class}');
      sparsity = {' ', 'sparse'};
      sparsity = char (sparsity(1 + [S.sparse]')');
      byte_size = repmat ({'Bytes'}, length (S), 1);
      bytes = [S.bytes]';
      byte_size(bytes > 1024^2) = {'MB'};
      bytes(bytes > 1024^2) = bytes(bytes > 1024^2) ./ 1024^2;
      byte_size(bytes > 1024^1) = {'KB'};
      bytes(bytes > 1024^1) = bytes(bytes > 1024^1) ./ 1024^1;
      bytes = num2str (bytes, '%.1f');
      byte_size = char (byte_size);
      space = repmat (' ', length (S), 1);
      str = [space, space, space, space, space, names, ...
        space, sizes, space, space, sparsity, space, classes, ...
        space, space, bytes, space, byte_size];
      disp (str);
    end
    
    function [X0, y0, Z0] = to_2006_fmt (sol, obj)
      % TO_2006_FMT  Export conic problem solution to VSDP 2006 format.
      %
      %   [X0, y0, Z0] = sol.to_2006_fmt (obj);  Export the solution 'sol' of
      %      a VSDP 2018 conic problem 'obj' to VSDP 2006 format.
      %
      %   See also vsdp, vsdp_solution.
      %
      
      % check supported cones
      if (any (obj.K.f > 0) || any (obj.K.l > 0) || ~isempty (obj.K.q))
        error ('VSDP_SOLUTION:to_2006_fmt:unsupported_cones', ...
          'to_2006_fmt: The VSDP 2006 format supports only SDP cones.');
      end
      
      if (isempty (sol.x))
        X0 = {};
      else
        X0 = cell (length (obj.K.s), 1);
        x0 = vsdp_indexable (sol.x, obj);
        for j = 1:length(obj.K.s)
          X0{j} = vsdp.smat ([], x0.s(j), 1/2);
        end
      end
      y0 = sol.y;  % Identical in both versions.
      if (isempty (sol.z))
        X0 = {};
      else
        Z0 = cell (length (obj.K.s), 1);
        z0 = vsdp_indexable (sol.z, obj);
        for j = 1:length(obj.K.s)
          Z0{j} = vsdp.smat ([], z0.s(j), 1);
        end
      end
    end
    
  end
  
end
