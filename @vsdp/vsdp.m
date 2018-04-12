classdef vsdp
  properties
    At
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
    function obj = vsdp (A, b, c, K, x, y, z)
      % VSDP  Conic problem data class.
      %
      %      A conic problem in primal and dual standard form is:
      %
      %         (P)  min  c'*x          (D)  max  b'*y
      %              s.t. A*x = b            s.t. z := c - A'*y
      %                   x in K                  z in K*
      %
      %      where K is a cartesian product of the cones of free variables
      %      (K.f), LP (K.l), SOCP (K.q), and SDP (K.s).  K* is the dual cone.
      %      For a theoretical introduction into verified conic programming see
      %      [Jansson2009].
      %
      %      The problem data of the block-diagonal structure:
      %
      %         'A'  double(n,m)
      %         'b'  double(m,1)
      %         'c'  double(n,1)
      %         'K'  struct {'f','l','q','s'}
      %
      %      The A, b, and c may be stored in dense or sparse format and can be
      %      interval quantities as well.
      %
      %   obj = VSDP (A, b, c, K)
      %   obj = VSDP (    ...   , x, y, z)
      %
      %   The class constructor accepts conic problem data in SEDUMI,
      %   VSDP 2006/2012 format.
      %      
      % Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)
      
      % check input
      narginchk (4, 7);
      
      if (isempty(A) || isempty(b) || isempty(c) || isempty(K))
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
      
      % If first parameter is of type cell, we assume, that the old VSDP 2006
      % input format was used.  Therefore perform a recursive call to import
      % that format.
      if (iscell (A))
        obj = vsdp.fromVSDP2006Fmt(A, b, c, K, x, y, z);
        return;
      end
      
      % index vector to speed-up svec
      Ivec = [];
      
      
      % prepare cone structure
      f = 0;
      l = 0;
      q = [];
      s = [];
      fields = isfield(K, {'f','l','q','s'});
      if fields(1)
        f = sum(K.f);
      end
      if fields(2)
        l = sum(K.l);
      end
      if fields(3)
        q = K.q(K.q > 0);
      end
      if fields(4)
        s = K.s(K.s > 0);
      end
      obj.K = struct('f', f, 'l', l, 'q', q(:), 's', s(:));
      % appropriate dimensions
      dim3 = f + l + sum(q) + (sum(s .* (s + 1)) / 2);
      
      
      % prepare interval input for b
      if (~isnumeric(b) && ~isa(b, 'intval'))
        error('VSDP:IMPORT_VSDP','cannot import dual objective "b"');
      end
      obj.b = b(:);
      m = length(b);
      
      % prepare interval input for A
      if (~isnumeric(A) && ~isa(A, 'intval'))
        error('VSDP:IMPORT_VSDP','cannot import coefficient matrix "A"');
      end
      % compact vectorized format
      if (any(size(A) - [dim3, m]) && any(size(A) - [m, dim3]))
        imported_fmt = 'SEDUMI';
        [A,Ivec] = vsvec(A,K,1,0,Ivec);
      end
      % transposed format
      if (size(A,1) == m)
        A = A';
      end
      % check size of A
      if any(size(A) ~= [dim3, m])
        error('VSDP:IMPORT_VSDP','wrong dimension of coefficient matrix "A"');
      end
      obj.At = A;
      
      % prepare interval input for c
      if (~isnumeric(c) && ~isa(c, 'intval'))
        error('VSDP:IMPORT_VSDP','cannot import primal objective "c"');
      end
      c = c(:);
      % compact vectorized format
      if (length(c) ~= dim3)
        [c,Ivec] = vsvec(c,K,1,0,Ivec);
      end
      obj.c = c;
      
      % prepare x
      if (~isempty(x))
        if (~isnumeric(x))
          error('VSDP:IMPORT_VSDP','primal solution vector "x" has to be numeric');
        end
        x = x(:);
        % compact vectorized format, mu=2
        if (length(x) ~= dim3)
          [x,Ivec] = vsvec(x,K,1,0,Ivec);  % Ivec can only be used with mu=1
          x = sscale(x,K,2);
        end
      end
      obj.x = x;
      
      % prepare y
      if (~isempty(y))
        y = y(:);
        if (~isnumeric(y) || (length(y) ~= m))
          error('VSDP:IMPORT_VSDP','cannot import dual solution vector "y"');
        end
      end
      obj.y = y;
      
      % prepare z
      if (~isempty(z))
        if (~isnumeric(z))
          error('VSDP:IMPORT_VSDP','cannot import dual solution "z"');
        end
        z = z(:);
        % compact vectorized format
        if (size(z,1) ~= dim3)
          z = vsvec(z,K,1,0,Ivec);
        end
      end
      obj.z = z;
      
    end
  end
end
