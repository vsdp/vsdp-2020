classdef vsdp_solution < handle
  
  properties (SetAccess = protected, GetAccess = public)
    % Solution type.
    %
    % One of:
    %   - 'Initial guess',
    %   - 'Approximate solution',
    %   - 'Rigorous lower bound',
    %   - 'Rigorous upper bound',
    %   - 'Certificate primal infeasibility',
    %   - 'Certificate dual infeasibility'
    %
    type
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
    %        'Primal and dual infeasibile'
    %
    solver_info
  end
  
  methods
    function obj = vsdp_solution (type, x, y, z, f_objective, solver_info)
      % VSDP_SOLUTION  Container class for VSDP solutions.
      %
      %   sol = vsdp_solution(type, x, y, z)  Create a solution object.
      %
      %   sol = vsdp_solution(..., f_objective, solver_info)
      %         Create a solution object.
      %
      %   Example:
      %
      %      type = 'Initial guess'
      %      x = [1 2 3]';  % Primal solution.
      %      y = [4 5]';    % Dual solution.
      %      z = [];
      %      sol = vsdp_approx_solution (type, x, y, z)
      %
      %
      %   See also vsdp.vsdp.
      
      % Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)
      
      narginchk (4, 6);
      
      obj.type = validatestring (type, {'Initial guess', ...
        'Approximate solution', 'Rigorous lower bound', ...
        'Rigorous upper bound', 'Certificate primal infeasibility', ...
        'Certificate dual infeasibility'});
      
      if (~isempty (x))
        validateattributes (x, {'double'}, {'vector'});
        obj.x = x(:);
      else
        obj.x = [];
      end
      if (~isempty (y))
        validateattributes (y, {'double'}, {'vector'});
        obj.y = y(:);
      else
        obj.y = [];
      end
      if (~isempty (z))
        validateattributes (z, {'double'}, {'vector'});
        obj.z = z(:);
      else
        obj.z = [];
      end
      if (nargin > 4)
        validateattributes (f_objective, {'double'}, {'size', [2, 1]});
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
          obj.solver_info.name = 'Unknown';
        end
        if (isfield (obj.solver_info, 'elapsed_time'))
          obj.solver_info.elapsed_time = validateattributes (elapsed_time, ...
            {'double'}, {'scalar'});
        else
          obj.solver_info.elapsed_time = nan;
        end
        if (isfield (obj.solver_info, 'termination'))
          obj.solver_info.termination = validatestring ( ...
            obj.solver_info.termination,  {'Unknown', 'Normal termination', ...
            'Primal infeasible', 'Dual infeasible', ...
            'Primal and dual infeasibile'});
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
      n = max (length (obj.x), length (obj.z));
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
      fprintf ('  [x] %s:\n', obj.type);
      switch (obj.type)
        case 'Initial guess'
          obj.mem_info ();
        case 'Approximate solution'
          if (~isempty (obj.solver_info))
            fprintf ('      ');
            fprintf ('Solver ''%s'': %s termination after %.1f seconds.\n', ...
              obj.solver_info.name, obj.solver_info.termination, ...
              obj.solver_info.elapsed_time);
          end
          if (~all (isnan (obj.f_objective)))
            fprintf ('\n      ');
            fprintf ('c''*x = %.15e',     obj.f_objective(1));
            fprintf ('\n      ');
            fprintf ('b''*y = %.15e\n\n', obj.f_objective(2));
          else
            obj.mem_info ();
          end
        case 'Rigorous lower bound'
        case 'Rigorous upper bound'
        case 'Certificate primal infeasibility'
        case 'Certificate dual infeasibility'
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
    
  end
  
end
