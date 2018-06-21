classdef vsdp < handle
  properties
    n   % Condensed cone dimension `K.f + K.l + K.`
    m   % Number of constraints
    At  % Transposed constraint matrix `n x m`
    b
    c
    K
    x
    y
    z
  end
  
  methods(Static)
    obj = fromVSDP2006Fmt (blk, A, C, b, X0, y0, Z0)
    [blk, A, C, b, X0, y0, Z0] = toVSDP2006Fmt (obj)
    [vA, I] = vsvec (A, K, mu, sflag, I)
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
      
      % Prepare cone structure and determine condensed cone dimension.
      obj.n = 0;
      obj.K.f = 0;
      obj.K.l = 0;
      obj.K.q = [];
      obj.K.s = [];
      if (isfield (K, 'f'))
        obj.K.f = sum(K.f);
        obj.n = obj.n + obj.K.f;
      end
      if (isfield (K, 'l'))
        obj.K.l = sum(K.l);
        obj.n = obj.n + obj.K.l;
      end
      if (isfield (K, 'q'))
        obj.K.q = K.q(K.q > 0);
        obj.n = obj.n + sum (obj.K.q);
      end
      if (isfield (K, 's'))
        obj.K.s = K.s(K.s > 0);
        obj.n = obj.n + (sum(obj.K.s .* (obj.K.s + 1)) / 2);
      end
      
      % prepare interval input for b
      if (~isfloat (b) && ~isa (b, 'intval'))
        error ('VSDP:VSDP','cannot import dual objective "b"');
      end
      obj.b = b(:);
      m = length(b);
      
      % prepare interval input for A
      if (~isfloat (At) && ~isa(At, 'intval'))
        error('VSDP:IMPORT_VSDP','cannot import coefficient matrix "A"');
      end
      % index vector to speed-up svec
      Ivec = [];
      % compact vectorized format
      if (any(size(At) - [n, m]) && any(size(At) - [m, n]))
        imported_fmt = 'SEDUMI';
        [At,Ivec] = vsvec(At,K.s,1,0,Ivec);
      end
      % transposed format
      if (size(At,1) == m)
        At = At';
      end
      % check size of A
      if any(size(At) ~= [n, m])
        error('VSDP:IMPORT_VSDP','wrong dimension of coefficient matrix "A"');
      end
      obj.At = At;
      
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
  end
end
