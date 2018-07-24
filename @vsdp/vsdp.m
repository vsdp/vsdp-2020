classdef vsdp < handle
  
  properties (GetAccess = public, SetAccess = protected)
    % Condensed cone dimension.
    %
    % 'n = K.f + K.l + sum(K.q) + (sum(K.s .* (K.s + 1)) / 2)'.
    n
    % Uncondensed cone dimension.
    %
    % 'N = K.f + K.l + sum(K.q) + sum(K.s .* K.s)'.
    N
    m % Number of constraints.
  end
  
  properties (GetAccess = public, SetAccess = public)
    At  % Transposed condensed constraint matrix 'n x m'.
    b   % Right-hand side vector, or dual objective vector.
    c   % Primal objective vector.
    K   % Cone structure, 'K.f', 'K.l', 'K.q', 'K.s'.
    x = [];   % Approximate primal solution or initial guess 'x0'.
    y = [];   % Approximate dual   solution or initial guess 'y0'.
    z = [];   % Approximate primal solution or initial guess 'z0'.
    %TODO     'info'   Termination-code with
    %                   0: indication of optimality (normal termination),
    %                   1: indication of primal infeasibility,
    %                   2: indication of dual infeasibility,
    %                   3: indication of both primal and dual infeasibilities,
    %                  -1: otherwise.
    solution = [];
    options = vsdp_options ();
  end
  
  % (Un-)Vectorization of `obj.At`, `obj.C`, `obj.X`, and `obj.Z`.
  properties (Access = protected)
    % Store index vectors to speed-up operations.
    svec_idx = struct('lower', [], 'upper', []);
    smat_idx = struct('lower', [], 'upper', []);
  end
  
  methods (Static)
    % Static constructor methods for VSDP objects.
    [obj, pd] = from_lp_solve_fmt (A, b, c, e, lb, ub);
    [obj, pd] = from_mps_file (fname);
    obj = from_sdpa_file (fname, blksize);
    obj = from_sdpa_fmt (bLOCKsTRUCT, c, F, x0, X0, Y0);
    obj = from_sdpt3_fmt (blk, At, b, c, x0, y0, z0);
    obj = from_vsdp_2006_fmt (blk, A, C, b, X0, y0, Z0);
    x = cell2mat (X);
    A = svec (obj, A, mu, param);
    A = smat (obj, a, mu);
    [vidx, midx, lidx] = sindex (K);
    [K, N, n] = validate_cone (K);
    [x, I] = vuls (A, b, x0, I);
  end
  
  methods (Access = public)
    info (obj);
    obj = validate (obj);
    obj = solve (obj, solver);
    [blk, A, C, b, X0, y0, Z0] = to_vsdp_2006_fmt (obj);
    
    % Default class methods
    varargout = size (obj, dim);
    disp (obj);
  end
  
  methods (Access = protected)
    obj = solve_csdp     (obj);
    obj = solve_glpk     (obj);
    obj = solve_linprog  (obj);
    obj = solve_lp_solve (obj);
    obj = solve_sdpa     (obj);
    obj = solve_sdpt3    (obj);
    obj = solve_sedumi   (obj);
  end
  
  methods
    function obj = vsdp (varargin)
      % VSDP  Conic problem data class.
      %
      %      A conic problem in primal and dual standard form is:
      %
      %         (P)  min   c'*x          (D)  max  b'*y
      %              s.t. At'*x = b           s.t. z := c - At*y
      %                       x in K               z in K^*
      %
      %      where K is a cartesian product of the cones of free variables
      %      (K.f), LP (K.l), SOCP (K.q), and SDP (K.s).  K^* is the dual cone.
      %      For a theoretical introduction into verified conic programming see
      %
      %         https://vsdp.github.io/references#Jansson2009
      %
      %      The problem data of the block-diagonal structure:
      %
      %         'At' double(n,m)  Store transposed A, usually m <<= n.
      %         'b'  double(m,1)
      %         'c'  double(n,1)
      %         'K'  struct {'f','l','q','s'}
      %
      %      The parameters A, b, and c may be stored in dense or sparse format
      %      and can be interval quantities.
      %
      %      The VSDP constructor can be called directly
      %
      %         obj = VSDP (obj)
      %         obj = VSDP (At,  b,  c, K)              % 2012 Format
      %         obj = VSDP (    ...      , x0, y0, z0)  % 2012 Format
      %         obj = VSDP (blk, At, C, b, X0, y0, Z0)  % 2006 Format
      %
      %      or the problem data can be imported by calling one of the methods:
      %
      %         obj = VSDP.from_lp_solve_fmt  (...)
      %         obj = VSDP.from_sdpam_fmt     (...)
      %         obj = VSDP.from_sdpt3_fmt     (...)
      %         obj = VSDP.from_vsdp_2006_fmt (...)
      %         obj = VSDP.from_mps_file      (...)
      %         obj = VSDP.from_sdpa_file     (...)
      %
      %   Example:
      %
      %
      %
      %   See also  VSDP.from_lp_solve_fmt, VSDP.from_sdpa_fmt,
      %             VSDP.from_sdpt3_fmt,    VSDP.from_vsdp_2006_fmt,
      %             VSDP.from_mps_file,     VSDP.from_sdpa_file.
      
      % Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)
      
      narginchk (1, 7);
      
      switch (nargin)
        case 1
          obj = varargin{1};
          if (~isa (obj, 'vsdp'))
            error ('VSDP:vsdp:noVsdpObj', ...
              ['vsdp: When called with a single argument, ', ...
              '''obj'' must be a VSDP object.']);
          end
        case {4, 5, 6, 7}
          % If first parameter is of type cell, we assume, that the VSDP 2006
          % input format was used  ==>  call static VSDP constructor.
          if (iscell (varargin{1}))
            obj = vsdp.fromVSDP2006Fmt (varargin{:});
            return;
          end
          [obj.At, obj.b, obj.c, obj.K] = varargin{1:4};
          % Treat optional parameter of solution guess.
          if (nargin >= 5)
            obj.x = varargin{5};
          end
          if (nargin >= 6)
            obj.y = varargin{6};
          end
          if (nargin == 7)
            obj.z = varargin{7};
          end
        otherwise
          error ('VSDP:vsdp:badNumberOfArguments', ...
            ['vsdp: Bad number of arguments, ', ...
            '''help vsdp.vsdp'' for more information.']);
      end
      
      obj = obj.validate ();
    end
  end
end
