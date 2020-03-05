classdef vsdp_options < handle
  
  properties
    % Growing factor for problem perturbation.
    %
    % Default: 1.5.
    ALPHA = 1.5;
    
    % Cache intermediate values where possible.
    %
    % Default: true.
    CACHE = true;
    
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
    % For a complete list of available solvers for your conic problem 'obj' run
    %
    %   obj.solve()
    %
    % Default: ''.
    SOLVER = '';
    
    % Option structure for approximate solver.
    %
    % Default: [].
    SOLVER_OPTIONS = [];
    
    % Use initial guess (x0, y0, z0) if given.
    %
    % Default: false.
    USE_INITIAL_GUESS = false;
    
    % Do not display VSDP warnings and minimize solver output.
    %
    % Default: true.
    VERBOSE_OUTPUT = true;
  end
  
  methods
    function obj = vsdp_options (varargin)
      % VSDP_OPTIONS  Default options for VSDP.
      %
      %   obj = VSDP_OPTIONS()  Returns an opitons object with the fields:
      %
      %       obj.ALPHA               = 0.5
      %       obj.CACHE               = true
      %       obj.FULL_EIGS_ENCLOSURE = false
      %       obj.ITER_MAX            = 10
      %       obj.SDP_ASSUME_SYMMETRY = false
      %       obj.SDP_MIN_BLK_SIZE    = 2500
      %       obj.SOLVER              = ''
      %       obj.SOLVER_OPTIONS      = []
      %       obj.USE_INITIAL_GUESS  = false
      %       obj.VERBOSE_OUTPUT      = true
      %
      %   obj = VSDP_OPTIONS('solver')  Optionally, pass the name of the
      %      approximate solver to use.  Type 'help vsdp_options.SOLVER' for
      %      more information.
      %
      %   Example:
      %
      %       obj = vsdp_options('sdpt3')
      %
      %       % Same as above.
      %       obj = vsdp_options()
      %       obj.SOLVER = 'sdpt3';
      %
      %   See also vsdp.
      
      
      % Copyright 2004-2020 Christian Jansson (jansson@tuhh.de)
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
      obj.SOLVER = str;
    end
    
    function set.SOLVER_OPTIONS (obj, opts)
      % No check, permit arbitrary objects.
      obj.SOLVER_OPTIONS = opts;
    end
    
    function set.USE_INITIAL_GUESS (obj, bool)
      try
        validateattributes (bool, {'logical'}, {'scalar'});
      catch
        error ('VSDP:vsdp_options:SOLVER_OPTIONS', ...
          'vsdp_options: USE_INITIAL_GUESS must be a logical scalar.');
      end
      obj.USE_INITIAL_GUESS = bool;
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
    
  end
end
