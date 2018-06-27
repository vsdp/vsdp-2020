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
  
  % Static constructor methods for VSDP objects.
  methods (Static)
    [obj, pd] = from_lp_solve_fmt (A, b, c, e, lb, ub);
    [obj, pd] = from_mps_file (fname);
    obj = from_sdpa_file (fname, blksize);
    obj = from_sdpam_fmt (bLOCKsTRUCT, c, F, x0, X0, Y0);
    obj = from_vsdp_2006_fmt (blk, A, C, b, X0, y0, Z0);
    x = cell2mat (X)
  end
  
  % Public methods.
  methods
    A = svec (obj, A, mu, isSymmetric);
    A = smat (obj, A, mu, isSymmetric);
    [blk, A, C, b, X0, y0, Z0] = to_vsdp_2006_fmt (obj);
  end
  
  methods
    function obj = vsdp (At, b, c, K, x, y, z)
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
      %   obj = VSDP (    ...   , x, y, z)
      %
      %   The class constructor accepts conic problem data in SeDuMi and
      %   VSDP 2006/2012 format.
      %
      % Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)
      
      % check input
      narginchk (4, 7);
      
      if (isempty(At) || isempty(b) || isempty(c) || isempty(K))
        error ('VSDP:VSDP', 'VSDP: empty input parameter');
      end
      if (nargin < 5)
        x = [];
      end
      if (nargin < 6)
        y = [];
      end
      if (nargin < 7)
        z = [];
      end
      
      % If first parameter is of type cell, we assume, that the VSDP 2006 input
      % format was used  ==>  Perform a recursive call to import that format.
      if (iscell (At))
        obj = vsdp.fromVSDP2006Fmt (At, b, c, K, x, y, z);
        return;
      end
      
      % Prepare cone structure.
      obj.K.f = 0;
      obj.K.l = 0;
      obj.K.q = [];
      obj.K.s = [];
      if (isfield (K, 'f'))
        obj.K.f = sum(K.f);
      end
      if (isfield (K, 'l'))
        obj.K.l = sum(K.l);
      end
      if (isfield (K, 'q'))
        obj.K.q = K.q(K.q > 0);
      end
      if (isfield (K, 's'))
        obj.K.s = K.s(K.s > 0);
      end
      
      % Determine uncondensed cone dimension.
      obj.N = obj.K.f + obj.K.l + sum(obj.K.q) ...
        + sum(obj.K.s .* obj.K.s);
      
      % Determine condensed cone dimension.
      obj.n = obj.K.f + obj.K.l + sum(obj.K.q) ...
        + (sum(obj.K.s .* (obj.K.s + 1)) / 2);
      
      % Prepare vector `b`.
      if (~isfloat (b) && ~isa (b, 'intval'))
        error ('VSDP:VSDP:wrongTypeB', ...
          'VSDP: Vector `b` has wrong data type.');
      end
      obj.b = b(:);
      obj.m = length(obj.b);
      
      % Check type of matrix `At`.
      if (~isfloat (At) && ~isa (At, 'intval'))
        error ('VSDP:VSDP:wrongTypeAt', ...
          'VSDP: Matrix `At` has wrong data type.');
      end
      
      % Check if any dimension of matrix `At` matches the number of constraints,
      % that is the length of vector `b`.
      if (~any(size(At) == obj.m))
        error ('VSDP:VSDP:badDimesionAt', ...
          'VSDP: No dimension of `At` matches the length of vector `b`.');
      end
      obj.At = At;
      
      % Ensure transposed format for `At` (n x m).
      if (size (obj.At, 2) ~= obj.m)
        obj.At = obj.At';
      end
      
      % Ensure compact vectorized format.
      if (size (obj.At, 1) > obj.n)
        At = obj.svec(At,1,false);
      elseif (size (obj.At, 1) ~= obj.n)
        error ('VSDP:VSDP:badDimesionAt', ...
          'VSDP: `Bad cone dimension `n of a `At` matches the length of vector `b`.');
      end
      
      % check size of A
      if any(size(At) ~= [n, m])
        error('VSDP:IMPORT_VSDP','wrong dimension of coefficient matrix "A"');
      end
      
      % prepare interval input for c
      if (~isfloat(c) && ~isa(c, 'intval'))
        error('VSDP:IMPORT_VSDP','cannot import primal objective "c"');
      end
      c = c(:);
      % compact vectorized format
      if (length(c) ~= n)
        [c,Ivec] = vsvec(c,K,1,0,Ivec);
      end
      obj.c = c;
      
      % prepare x
      if (~isempty(x))
        if (~isfloat(x))
          error('VSDP:IMPORT_VSDP','primal solution vector "x" has to be numeric');
        end
        x = x(:);
        % compact vectorized format, mu=2
        if (length(x) ~= n)
          [x,Ivec] = vsvec(x,K,1,0,Ivec);  % Ivec can only be used with mu=1
          x = sscale(x,K,2);
        end
      end
      obj.x = x;
      
      % prepare y
      if (~isempty(y))
        y = y(:);
        if (~isfloat(y) || (length(y) ~= m))
          error('VSDP:IMPORT_VSDP','cannot import dual solution vector "y"');
        end
      end
      obj.y = y;
      
      % prepare z
      if (~isempty(z))
        if (~isfloat(z))
          error('VSDP:IMPORT_VSDP','cannot import dual solution "z"');
        end
        z = z(:);
        % compact vectorized format
        if (size(z,1) ~= n)
          z = vsvec(z,K,1,0,Ivec);
        end
      end
      obj.z = z;
      
    end
    
    function varargout = size (obj, dim)
      %SIZE Size of a table.
      if (nargin == 1)
        if (nargout < 2)
          varargout = {[obj.m, obj.n]};
        elseif (nargout == 2)
          varargout = {obj.m, obj.n};
        else
          varargout(1:2) = {obj.m, obj.n};
          varargout(3:nargout) = {1};
        end
      else
        if (dim == 1)
          varargout = {obj.m};
        elseif (dim == 2)
          varargout = {obj.n};
        else
          varargout = {1};
        end
      end
    end
  end
end
