classdef vsdp_approx_solution < handle
  
  properties (SetAccess = protected, GetAccess = public)
    x   % Approximate primal solution or initial guess 'x0'.
    y   % Approximate dual   solution or initial guess 'y0'.
    z   % Approximate primal solution or initial guess 'z0'.
    % (Optional) Primal and dual objective values for the provided solutions.
    f_objective
    % (Optional) Solver termination-code of this solution.
    %
    %  -1: Unknown
    %   0: Normal termination
    %   1: Primal infeasible
    %   2: Dual   infeasible
    %   3: Primal and dual infeasibile
    %
    info
    solver  % (Optional) Solver used to compute this solution.
  end
  
  methods
    function obj = vsdp_approx_solution (x, y, z, f_objective, solver, info)
      % VSDP_APPROX_SOLUTION  Container class for VSDP solutions.
      %
      %   sol = vsdp_approx_solution(x,y,z)  Create an solution object for an
      %       approximation or initial guess.
      %
      %   ___ = vsdp_approx_solution( ... , f_objective, solver, info)
      %       Optionally, add the primal and dual objective values
      %       'f_objective', the 'solver' used to compute this solution and the
      %       solver termination code 'info'.
      %
      %   Example:
      %
      %      c = [1 1 1]';  % Primal objective vector.
      %      b = [2 2]';    % Right-hand side objective vector.
      %      x = [1 2 3]';  % Primal solution.
      %      y = [4 5]';    % Dual solution.
      %      z = [];
      %      f_objective = [c'*x; b'*y];
      %      sol = vsdp_approx_solution (x, y, z, f_objective, 'sdpt3', 0)
      %
      %
      %   See also vsdp.vsdp.
      
      % Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)
      
      narginchk (3, 6);
      
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
      if (~isempty (obj.x) && ~isempty (obj.z))
        if (length (obj.x) ~= length (obj.z))
          error ('VSDP:vsdp_approx_solution:sizeMissmatch', ...
            ['vsdp_approx_solution: The lengths of ''x'' and ''z'' ', ...
            'must be the same.']);
        end
      end
      if (nargin > 4)
        validateattributes (f_objective, {'double'}, {'vector'});
        obj.f_objective = f_objective(:);
      else
        obj.f_objective = nan (2, 1);
      end
      if (nargin > 4)
        validateattributes (solver, {'char'}, {'vector'});
        obj.solver = solver;
      else
        obj.solver = 'none';
      end
      if (nargin > 5)
        validateattributes (info, {'double'}, {'scalar'});
        obj.info = info;
      else
        obj.info = -1;  % Unknown
      end
    end
    
    function varargout = size (obj, dim)
      % SIZE  Size of VSDP_APPROX_SOLUTION instance.
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
      if (~all (isnan (obj.f_objective)))
        if (~strcmp (obj.solver, 'none'))
          str = sprintf ('          %s by ''%s'' (code %d)\n', ...
            obj.info_text (obj.info), obj.solver, obj.info);
        else
          str = '';
        end
        fprintf ('     [x]  Approximate solution:\n%s\n', str);
        fprintf ('               c''*x = %.15e\n',   obj.f_objective(1));
        fprintf ('               b''*y = %.15e\n\n', obj.f_objective(2));
      else
        obj.mem_info ();
      end
    end
    
    function mem_info (obj)
      % MEM_INFO Display memory information.
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
      str = [space, space, space, names, ...
        space, sizes, space, space, sparsity, space, classes, ...
        space, space, bytes, space, byte_size];
      disp (str);
    end
    
  end
  
  methods (Static)
    function str = info_text (val)
      % INFO_TEXT  Translate info code to human readable text.
      switch (val)
        case 0
          str = 'Normal termination';
        case 1
          str = 'Primal infeasible';
        case 2
          str = 'Dual infeasible';
        case 3
          str = 'Primal and dual infeasibile';
        otherwise
          str = 'Unknown';
      end
    end
  end
  
end
