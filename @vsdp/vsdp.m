classdef vsdp < handle
  
  properties
    % Condensed cone dimension, e.g.
    % `n = K.f + K.l + sum(K.q) + (sum(K.s .* (K.s + 1)) / 2)`.
    n
    % Uncondensed cone dimension, e.g.
    % `N = K.f + K.l + sum(K.q) + sum(K.s .* K.s)`.
    N
    m   % Number of constraints.
    At  % Transposed condensed constraint matrix `n x m`.
    b   % Right-hand side vector, or dual objective vector.
    c   % Primal objective vector.
    K   % Cone structure, `K.f`, `K.l`, `K.q`, `K.s`.
    x0  % Approximate primal solution.
    y0  % Approximate dual solution.
    z0  % Approximate primal solution.
  end
  
  properties (Dependent)
    A  % Condensed constraint matrix `m x n`.
  end
  methods
    function A = get.A (obj)
      A = obj.At';
    end
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
    obj = from_sdpam_fmt (bLOCKsTRUCT, c, F, x0, X0, Y0);
    obj = from_sdpt3_fmt (blk, At, b, c, x0, y0, z0)
    obj = from_vsdp_2006_fmt (blk, A, C, b, X0, y0, Z0);
    x = cell2mat (X)
    obj = validate (At, b, c, K, x0, y0, z0);
  end
  
  % Public methods.
  methods
    A = svec (obj, A, mu, isSymmetric);
    A = smat (obj, A, mu, isSymmetric);
    [blk, A, C, b, X0, y0, Z0] = to_vsdp_2006_fmt (obj);
    
    % Default class methods
    varargout = size (obj, dim);
  end
  
  methods
    function obj = vsdp (At, b, c, K, x0, y0, z0)
      % VSDP  Conic problem data class.
      %
      %      A conic problem in primal and dual standard form is:
      %
      %         (P)  min  c'*x          (D)  max  b'*y
      %              s.t. A*x = b            s.t. z := c - A'*y
      %                   x in K                  z in K^*
      %
      %      where K is a cartesian product of the cones of free variables
      %      (K.f), LP (K.l), SOCP (K.q), and SDP (K.s).  K^* is the dual cone.
      %      For a theoretical introduction into verified conic programming see
      %      [Jansson2009].
      %
      %      The problem data of the block-diagonal structure:
      %
      %         'At' double(n,m)  Store transposed A, because m <<= n.
      %         'b'  double(m,1)
      %         'c'  double(n,1)
      %         'K'  struct {'f','l','q','s'}
      %
      %      The parameters A, b, and c may be stored in dense or sparse format
      %      and can be interval quantities.
      %
      %   obj = VSDP (A, b, c, K)
      %   obj = VSDP (    ...   , x0, y0, z0)
      %
      %   The class constructor accepts conic problem data in
      %   VSDP 2006/2012 format.
      %
      % Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)
      
      narginchk (4, 7);
      
      if (isempty(At) || isempty(b) || isempty(c) || isempty(K))
        error ('VSDP:VSDP', 'VSDP: empty input parameter');
      end
      if (nargin < 5)
        x0 = [];
      end
      if (nargin < 6)
        y0 = [];
      end
      if (nargin < 7)
        z0 = [];
      end
     
      % If first parameter is of type cell, we assume, that the VSDP 2006 input
      % format was used  ==>  Perform a recursive call to import that format.
      if (iscell (At))
        obj = vsdp.fromVSDP2006Fmt (At, b, c, K, x0, y0, z0);
        return;
      end
      
      obj = validate (At, b, c, K, x0, y0, z0);
    end
  end
end
