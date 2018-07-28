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
    At  % Transposed condensed constraint matrix 'n x m'.
    b   % Right-hand side vector, or dual objective vector.
    c   % Primal objective vector.
    K   % Cone structure, 'K.f', 'K.l', 'K.q', 'K.s'.
    cache_memory % Memory for caching expensive computation results.
    solution     % Approximate solution or initial guess.
    options      % Options for this problem instance.
  end
  
  methods (Static)
    % Static constructor methods for VSDP objects.
    obj = from_mps_file     (fname);
    obj = from_sdpa_file    (fname, blksize);
    obj = from_sdpa_fmt     (bLOCKsTRUCT, c, F, x0, X0, Y0);
    obj = from_sdpt3_fmt    (blk, At, b, c, x0, y0, z0);
    obj = from_2006_fmt     (blk, A,  C, b, X0, y0, Z0);
    obj = from_lp_solve_fmt (     A,  b, c,  e, lb, ub);
    
    % Other static methods.
    x = cell2mat (X);
    A = svec (obj, A, mu, param);
    A = smat (obj, a, mu);
    E = verify_eigsym (A);
    x = verify_uls (obj, A, b, x0);
    [vidx, midx, lidx] = sindex (K);
    [K, N, n] = validate_cone (K);
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
      %         'K'  struct{'f','l','q','s'}
      %
      %      The parameters A, b, and c may be stored in dense or sparse format
      %      and can be interval quantities.
      %
      %      The VSDP constructor can be called directly
      %
      %         obj2 = VSDP (obj1)                       % Copy Problem data.
      %         obj  = VSDP (At,  b,  c, K)              % 2012 Format
      %         obj  = VSDP (    ...      , x0, y0, z0)  % 2012 Format
      %         obj  = VSDP (blk, At, C, b, X0, y0, Z0)  % 2006 Format
      %
      %      or the problem data can be imported by calling one of the methods:
      %
      %         obj = VSDP.from_mps_file     (...)
      %         obj = VSDP.from_sdpa_file    (...)
      %         obj = VSDP.from_sdpa_fmt     (...)
      %         obj = VSDP.from_sdpt3_fmt    (...)
      %         obj = VSDP.from_2006_fmt     (...)
      %         obj = VSDP.from_lp_solve_fmt (...)
      %
      %   To use an initial guess (x0,y0,z0) type:
      %
      %      obj.add_solution (x0, y0, z0);
      %      obj.options.USE_STARTING_POINT = true;
      %
      %   Example:
      %
      %       A1 = [0 1;
      %             1 0];
      %       A2 = [1 1;
      %             1 1];
      %       C =  [1 0;
      %             0 1];
      %       K.s = 2;
      %       At = [A1(:), A2(:)];  % Vectorize data
      %       c  = C(:);
      %       b  = [1; 2.0001];
      %       obj = vsdp(At, b, c, K).solve()
      %
      %   See also  VSDP.from_lp_solve_fmt, VSDP.from_sdpa_fmt,
      %             VSDP.from_sdpt3_fmt,    VSDP.from_sdpa_file,
      %             VSDP.from_mps_file,     VSDP.from_2006_fmt.
      
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
          % Make a "clean" copy, only with the underlying data.
          [K.f, K.l, K.q, K.s] = deal (obj.K.f, obj.K.l, obj.K.q, obj.K.s);
          data = {obj.At, obj.b, obj.c, K};
          obj = vsdp (data{:});
        case 4
          % If first parameter is of type cell, we assume, that the VSDP 2006
          % input format was used  ==>  call static VSDP constructor.
          if (iscell (varargin{1}))
            obj = vsdp.fromVSDP2006Fmt (varargin{:});
            return;
          end
          [obj.At, obj.b, obj.c, obj.K] = varargin{1:4};
          obj.options = vsdp_options ();
          obj.solution = [];
          obj = obj.validate ();
        case {5, 6, 7}
          % First create VSDP problem instance.
          obj = vsdp (varargin{1:4});
          % Then add initial solution guess.
          obj.add_solution (varargin{5:end});
        otherwise
          error ('VSDP:vsdp:badNumberOfArguments', ...
            ['vsdp: Bad number of arguments, type ', ...
            '''help vsdp.vsdp'' for more information.']);
      end
    end
  end
end
