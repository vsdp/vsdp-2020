classdef vsdp_options < handle
  
  properties
    % Growing factor for problem perturbation.
    %
    % Default: 0.5.
    ALPHA = 0.5;
    
    % Use algorithm for stronger complete eigenvalue enclosure.
    %
    % Default: true.
    FULL_EIGS_ENCLOSURE = true;
    
    % Maximum number of iterations for the problem perturbation.
    %
    % Default: 10.
    ITER_MAX = 10;
    
    % Assume symmetry of input matrices for SDP cones.
    %
    % Some solvers like SeDuMi allow non-symmetric input.  If true VSDP is
    % increasing the speed of data conversion by processing only a triangular
    % part of the SDP block matrices.  Default: false.
    SDP_ASSUME_SYMMETRY = false;
    
    % Minimum size of an SDP block.
    %
    % If block size is smaller than 'MIN_SDPBLK_SIZE' the blocks will be
    % grouped when transforming into SDPT3 format.
    %
    % Default: 2500.
    SDP_MIN_BLK_SIZE = 2500;
    
    % Approximate solver to use.
    %
    % Select one of the supported solvers:
    %   'sdpt3' (default), 'sedumi', 'sdpa', 'csdp', 'sdplr', 'lp_solve',
    %   'linprog', or 'glpk'.
    SOLVER = 'sdpt3';
    
    % Option structure for approximate solver.
    %
    % Default: [].
    SOLVER_OPTIONS = [];
    
    % Use initial starting point (x0, y0, z0) if given.
    %
    % Default: false.
    USE_STARTING_POINT = false;
    
    % Do not display VSDP warnings and minimize solver output.
    %
    % Default: true.
    VERBOSE_OUTPUT = true;
    
    % Apply the full non-symmetric matrix LSS enclosure is applied.
    %
    % Default: false.
    VERIFY_FULL_LSS = false;
  end
  
  methods
    function obj = vsdp_options (varargin)
      % VSDP_OPTIONS  Default options for VSDP.
      %
      %   obj = VSDP_OPTIONS()  Returns an opitons object with the fields:
      %
      %       obj.ALPHA               = 0.5
      %       obj.FULL_EIGS_ENCLOSURE = false
      %       obj.ITER_MAX            = 10
      %       obj.SDP_ASSUME_SYMMETRY = false
      %       obj.SDP_MIN_BLK_SIZE    = 2500
      %       obj.SOLVER              = 'sdpt3'
      %       obj.SOLVER_OPTIONS      = []
      %       obj.USE_STARTING_POINT  = false
      %       obj.VERBOSE_OUTPUT      = true
      %       obj.VERIFY_FULL_LSS     = false
      %
      %   obj = VSDP_OPTIONS('solver')  Optionally, pass the name of the
      %      approximate solver to use.  Type 'help vsdp_options.SOLVER' for a
      %      list of supported solvers.
      %
      %   Example:
      %
      %       obj = vsdp_options('sdpt3')
      %
      %       % Same as above.
      %       obj = vsdp_options()
      %       obj.SOLVER = 'sdpt3';
      %
      %   See also vsdp.vsdp.
      
      % Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)
      
      narginchk (0, 1);
      
      % Check for default solver by explicitly setting it again.
      if (nargin == 1)
        obj.SOLVER = varargin{1};
      else
        obj.SOLVER = obj.SOLVER;
      end
    end
    
    function set.ALPHA (obj, val)
      try
        validateattributes (val, {'numeric'}, {'positive', 'scalar'});
      catch
        error ('VSDP:vsdp_options:ALPHA', ...
          'vsdp_options: ALPHA must be a numeric positive scalar.');
      end
      obj.ALPHA = val;
    end
    
    function set.FULL_EIGS_ENCLOSURE (obj, bool)
      try
        validateattributes (bool, {'logical'}, {'scalar'});
      catch
        error ('VSDP:vsdp_options:FULL_EIGS_ENCLOSURE', ...
          'vsdp_options: FULL_EIGS_ENCLOSURE must be a logical scalar');
      end
      obj.FULL_EIGS_ENCLOSURE = bool;
    end
    
    function set.ITER_MAX (obj, val)
      try
        validateattributes (val, {'numeric'}, {'positive', 'scalar'});
      catch
        error ('VSDP:vsdp_options:ITER_MAX', ...
          'vsdp_options: ITER_MAX must be a numeric positive scalar.');
      end
      obj.ITER_MAX = val;
    end
    
    function set.SDP_ASSUME_SYMMETRY (obj, bool)
      try
        validateattributes (bool, {'logical'}, {'scalar'});
      catch
        error ('VSDP:vsdp_options:SDP_ASSUME_SYMMETRY', ...
          'vsdp_options: SDP_ASSUME_SYMMETRY must be a logical scalar');
      end
      obj.SDP_ASSUME_SYMMETRY = bool;
    end
    
    function set.SDP_MIN_BLK_SIZE (obj, val)
      try
        validateattributes (val, {'numeric'}, {'positive','scalar'});
      catch
        error ('VSDP:vsdp_options:SDP_MIN_BLK_SIZE', ...
          'vsdp_options: SDP_MIN_BLK_SIZE must be a numeric positive scalar.');
      end
      obj.SDP_MIN_BLK_SIZE = val;
    end
    
    function set.SOLVER (obj, str)
      try
        validateattributes (str, {'char'}, {});
      catch
        error ('VSDP:vsdp_options:SOLVER', ...
          'vsdp_options: SOLVER must be a string (char array).');
      end
      defaultSolver = 'sdpt3';
      defaultSolverAvailable = (exist ('sqlp', 'file') == 2);
      str = lower (str);
      switch (str)
        case defaultSolver
          setSolver = defaultSolverAvailable;
        case 'sedumi'
          setSolver = (exist ('sedumi', 'file') == 2);
        case 'sdpa'
          setSolver = ((exist ('sdpam', 'file') == 2) ...
            && ((exist ('mexsdpa', 'file') == 3) ...
            || (exist ('callSDPA', 'file') == 2)));
        case 'csdp'
          setSolver = (exist ('csdp', 'file') == 2);
        case 'sdplr'
          setSolver = (exist ('sdplr', 'file') == 2);
        case 'lp_solve'
          setSolver = (exist ('lp_solve', 'file') == 2);
        case 'linprog'
          setSolver = (exist ('linprog', 'file') == 2);
        case 'glpk'
          setSolver = (exist ('glpk', 'file') == 2);
        otherwise
          error ('VSDP:vsdp_options:SOLVER', ...
            'vsdp_options: Solver ''%s'' is not supported.', str);
      end
      if (setSolver)
        obj.SOLVER = str;
      elseif (~strcmp(obj.SOLVER, defaultSolver))
        % If another solver was successully set in the past keep it.
        warning ('VSDP:vsdp_options:SOLVER', ...
          ['vsdp_options: Solver ''%s'' is not available.  ', ...
          'Keep solver: ''%s''.'], str, obj.SOLVER);
      elseif (defaultSolverAvailable)
        warning ('VSDP:vsdp_options:SOLVER', ...
          ['vsdp_options: Solver ''%s'' is not available.  ', ...
          'Use default solver: ''%s''.'], str, defaultSolver);
        obj.SOLVER = defaultSolver;
      else
        error ('VSDP:vsdp_options:SOLVER', ...
          ['vsdp_options: Solver ''%s'' and the default solver ''%s'' ', ...
          'are not available.'], str, defaultSolver);
      end
    end
    
    function set.SOLVER_OPTIONS (obj, str)
      try
        validateattributes (str, {'char'}, {});
      catch
        error ('VSDP:vsdp_options:SOLVER_OPTIONS', ...
          'vsdp_options: SOLVER_OPTIONS must be a string (char array).');
      end
      obj.SOLVER_OPTIONS = str;
    end
    
    function set.USE_STARTING_POINT (obj, bool)
      try
        validateattributes (bool, {'logical'}, {'scalar'});
      catch
        error ('VSDP:vsdp_options:SOLVER_OPTIONS', ...
          'vsdp_options: USE_STARTING_POINT must be a logical scalar.');
      end
      obj.USE_STARTING_POINT = bool;
    end
    
    function set.VERBOSE_OUTPUT (obj, bool)
      try
        validateattributes (bool, {'logical'}, {'scalar'});
      catch
        error ('VSDP:vsdp_options:VERBOSE_OUTPUT', ...
          'vsdp_options: VERBOSE_OUTPUT must be a logical scalar');
      end
      obj.VERBOSE_OUTPUT = bool;
      % Get a list of all warning IDs on Linux:
      %   grep -R "warning ('VSDP:*"
      state = {'off', 'on'};
      %TODO: At the end add all warning IDs
      warning (state{bool + 1}, 'VSDP:vsdp_options:SOLVER')
    end
    
    function set.VERIFY_FULL_LSS (obj, bool)
      try
        validateattributes (bool, {'logical'}, {'scalar'});
      catch
        error ('VSDP:vsdp_options:VERIFY_FULL_LSS', ...
          'vsdp_options: VERIFY_FULL_LSS must be a logical scalar');
      end
      obj.VERIFY_FULL_LSS = bool;
    end
    
  end
end
