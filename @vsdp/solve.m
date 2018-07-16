function [obj,x,y,z,info] = mysdps(P,varargin)
% MYSDPS  Interface to several approximate conic solvers for VSDP.
%
%   [obj,x,y,z,info] = MYSDPS(P)
%   [obj,x,y,z,info] = MYSDPS(A,b,c,K)
%      Uses a preselected solver to get an approximate solution of a conic
%      problem `P = (A,b,c,K)` in the standard primal-dual form:
%
%         (primal)  min  c'*x          (dual)  max  b'*y
%                   s.t. A*x = b               s.t. z := c - A'*y
%                        x in K                     z in K^*
%
%      where K is a cartesian product of the cones of free variables (K.f),
%      LP (K.l), SOCP (K.q), and SDP (K.s).  K^* is the dual cone.  For a
%      theoretical introduction into verified conic programming see
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
%      The function returns:
%
%         'obj'    Primal and dual objective value [c'*x, b'*y].
%
%         'x,y,z'  An approximate optimal solution or a primal or dual
%                  infeasibility certificate.
%
%         'info'   Termination-code with
%                   0: indication of optimality (normal termination),
%                   1: indication of primal infeasibility,
%                   2: indication of dual infeasibility,
%                   3: indication of both primal and dual infeasibilities,
%                  -1: otherwise.
%
%   MYSDPS(...,x0,y0,z0) optionally provide an initial guess (x0,y0,z0)
%      to the approximate solver.
%
%   MYSDPS(...,[],[],[],opts) optionally provide a structure for additional
%      parameter settings, explained in vsdpinit.
%
%   Example:
%
%       A1 = [0 1;
%             1 0];
%       b1 = 1;
%       A2 = [1 1;
%             1 1];
%       b2 = 2.0001;
%       C = [1 0;
%            0 1];
%       K.s = 2;
%
%       % Vectorize data
%       A = [A1(:),A2(:)];
%       b = [b1; b2];
%       c = C(:);
%       [obj,x,y,z,info] = mysdps(A,b,c,K);
%
%   See also vsdpinit, vsdpup, vsdplow, vsdpinfeas.

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

% Check input
narginchk(1,8);
if (nargin < 8)
  opts = [];
end

VSDP_OPTIONS = vsdpinit(opts);
OPTIONS = VSDP_OPTIONS.SOLVER_OPTIONS;

if (~isa (P, 'vsdp'))
  P = vsdp(P, varargin{:});
end

% In case of interval data solve midpoint problem
A = mid(P.A);
b = mid(P.b);
c = mid(P.c);

% initialization of default output
obj = [inf, -inf];
y = [];
z = [];
info = -1;

% call solver for problem
switch (VSDP_OPTIONS.SOLVER)
  case 'sedumi'
    % Setup
    if (~VSDP_OPTIONS.VERBOSE_OUTPUT)
      OPTIONS.fid = 0;
    end
    [A,c] = vsdp2sedumi(A,c,[],[],K,opts);

    % Call solver
    [x,y,INFO] = sedumi(A,b,c,K,OPTIONS);

    % Transform results to VSDP format
    if (~isempty(x) && ~isempty(y))
      z = c - A*y;
      obj = [c'*x, b'*y];
    end
    info = INFO.pinf + 2*INFO.dinf;

  case 'sdpt3'
    % Setup
    if (~VSDP_OPTIONS.USE_STARTING_POINT || isempty(x0) || isempty(y0) ...
        || isempty(z0))
      x0 = [];
      y0 = [];
      z0 = [];
    end
    if (~VSDP_OPTIONS.VERBOSE_OUTPUT)
      if (exist('sqlpmain.m','file') == 2) % if SDPT3-4.0
        OPTIONS.printlevel = 0; % default: 3
      else
        OPTIONS.printyes = 0;   % default: 1
      end
    end
    [blk,A,c,x0,z0] = vsdp2sdpt3(K,A,c,x0,z0,opts);

    % Call solver
    [obj,x,y,z,INFO] = sqlp(blk,A,c,b,OPTIONS,x0,y0,z0);

    % Transform results to VSDP format
    [~,~,~,x,z] = sdpt2vsdp(blk,[],[],x,z);
    % save info codes
    if isstruct(INFO)
      info = INFO.termcode;  % SDPT3-4.0 output
    else
      info = INFO(1);  % SDPT3-3.x output
    end

  case 'sdpa'
    % Check cones
    if ((K.f > 0) || (sum(K.q) > 0))
      error('VSDP:MYSDPS', ...
        'SOCP (K.q) and free variables (K.f) are not supported by SDPAM');
    end

    % Setup
    if (~VSDP_OPTIONS.USE_STARTING_POINT || isempty(x0) || isempty(y0) ...
        || isempty(z0))
      x0 = [];
      y0 = [];
      z0 = [];
    end
    if (~VSDP_OPTIONS.VERBOSE_OUTPUT)
      OPTIONS.print = 'no';
    end
    [mDIM,nBLOCK,bLOCKsTRUCT,ct,F,x0,X0,Y0] = vsdp2sdpam(A,b,c,K,x0,y0,z0);

    % Call solver
    if (exist('mexsdpa','file') == 3)
      [~,x,X,Y,~] = sdpam(mDIM,nBLOCK,bLOCKsTRUCT,ct,F,x0,X0,Y0,OPTIONS);
    elseif (exist('callSDPA','file') == 2)
      [x,X,Y] = callSDPA(mDIM,nBLOCK,bLOCKsTRUCT,ct,F,x0,X0,Y0,OPTIONS);
    else
      error('VSDP:MYSDPS', 'You need to compile the SDPA MEX-interface.');
    end

    % Transform results to VSDP format
    [~,~,~,~,x,y,z] = sdpam2vsdp(bLOCKsTRUCT,[],[],x,X,Y);
    if (~isempty(x) && ~isempty(y))
      obj = [c'*x, b'*y];
      info = 0;
    end

  case 'csdp'
    % Check cones
    if (sum(K.q) > 0)
      error('VSDP:MYSDPS', 'SOCP (K.q) is not supported by CSDP');
    end
    if (K.f > 0)
      warning('VSDP:MYSDPS', ['CSDP supports free variables (K.f) by ' ...
        'converting them to the difference of positive variables.  The ' ...
        'resulting problem is ill-posed.']);
    end

    % Setup
    if (~VSDP_OPTIONS.VERBOSE_OUTPUT)
      OPTIONS.printlevel = 0;
    end
    OPTIONS.check = 0; % check problem data for symmetry.
    [A,c,x0,z0] = vsdp2sedumi(A,c,x0,z0,K); % that format can be read by csdp

    % Call solver
    if (VSDP_OPTIONS.USE_STARTING_POINT && ~isempty(x0) && ~isempty(y0) ...
        && ~isempty(z0))
      [x,y,z,INFO] = csdp(A,full(b),full(c),K,OPTIONS,full(x0),full(y0), ...
        full(z0));
    else
      [x,y,z,INFO] = csdp(A,full(b),full(c),K,OPTIONS);
    end

    % Transform results to VSDP format
    if (~isempty(x) && ~isempty(y))
      obj = [c'*x, b'*y];
    end
    if any(INFO == [0, 1, 2])
      info = INFO;
    end

  case 'sdplr'
    % Check cones
    if ((K.f > 0) || (sum(K.q) > 0))
      error('VSDP:MYSDPS', ...
        'SOCP (K.q) and free variables (K.f) are not supported by SDPLR');
    end

    % Setup
    [A,c,x0] = vsdp2sedumi(A,c,x0,[],K); % that format can be read by sdplr

    % Call solver
    if (VSDP_OPTIONS.USE_STARTING_POINT && ~isempty(x0) && ~isempty(y0))
      [x,y] = sdplr(A,b,c,K,OPTIONS,[],x0,y0);
    else
      [x,y] = sdplr(A,b,c,K,OPTIONS);
    end

    % Transform results to VSDP format
    if (~isempty(x) && ~isempty(y))
      z = c - A*y;
      obj = [c'*x, b'*y];
      info = 0;
    end

  case 'lp_solve'
    % Check cones
    if ((sum(K.q) > 0) || (sum(K.s) > 0))
      error('VSDP:MYSDPS', ...
        'SOCP (K.q) and SDP (K.s) are not supported by LPSOLVE');
    end
    % Call solver
    [~,x,y,stat] = lp_solve(-full(c),A',full(b),zeros(length(b),1), ...
      [-inf(K.f,1); zeros(K.l,1)]);
    % Transform results to VSDP format
    if (~isempty(x) && ~isempty(y))
      y = -y;
      z = c - A*y;
      obj = [c'*x, b'*y];
    end
    switch (stat)
      case 0
        info = 0; % normal termination
      case 2
        info = 1; % primal infeasible
      case 3
        info = 2; % dual infeasible (primal unbounded)
      otherwise
        info = -1; % an error occured
    end

  case 'linprog'
    % Check cones
    if ((sum(K.q) > 0) || (sum(K.s) > 0))
      error('VSDP:MYSDPS', ...
        'SOCP (K.q) and SDP (K.s) are not supported by LINPROG');
    end

    % Setup
    if (~VSDP_OPTIONS.USE_STARTING_POINT || isempty(x0))
      x0 = [];
    end
    if (isempty(OPTIONS)) % Avoid warning
      OPTIONS = optimoptions('linprog','Algorithm','interior-point-legacy');
    end
    if (~VSDP_OPTIONS.VERBOSE_OUTPUT)
      OPTIONS = optimoptions(OPTIONS, 'Display', 'off');
    end
    lbound = [-inf(K.f,1); zeros(K.l,1)]; % lower bound
    ubound = inf(length(c),1);            % upper bound

    % Call solver
    [x,~,flag,~,lambda] = ...
      linprog(full(c),[],[],A',full(b),lbound,ubound,x0,OPTIONS);

    % Transform results to VSDP format
    if (~isempty(x))
      y = -lambda.eqlin;
      z = c - A*y;
      obj = [c'*x, b'*y];
    end
    switch (flag)
      case 1
        info = 0; % normal termination
      case -2
        info = 1; % primal infeasible
      case -3
        info = 2; % dual infeasible (primal unbounded)
      case -5
        info = 3; % primal and dual infeasible
      otherwise
        info = -1; % an error occured
    end

  case 'glpk'
    % Check cones
    if ((sum(K.q) > 0) || (sum(K.s) > 0))
      error('VSDP:MYSDPS', ...
        'SOCP (K.q) and SDP (K.s) are not supported by GLPK');
    end

    % Setup
    sense = 1;  % minimization
    lbound = [-inf(K.f,1); zeros(K.l,1)]; % lower bound
    ubound = inf(length(c),1);            % upper bound
    vtype = repmat('C',length(c),1); % continuous variable
    ctype = repmat('S',length(b),1); % equality constraint

    % Call solver
    [x,~,errnum,extra] = ...
      glpk(full(c),A',full(b),lbound,ubound,ctype,vtype,sense);

    % Transform results to VSDP format
    if (~isempty(x))
      y = extra.lambda;
      z = c - A*y;
      obj = [c'*x, b'*y];
    end
    switch (errnum)
      case 0
        info = 0; % normal termination
      case 10
        info = 1; % primal infeasible
      case 11
        info = 2; % dual infeasible
      case 15
        info = 3; % primal and dual infeasible
      otherwise
        info = -1; % an error occured
    end

  otherwise
    error('VSDP:MYSDPS','The solver %s could not be detected!!', ...
      VSDP_OPTIONS.SOLVER);
end

[x,z] = export_vsdp(imported_fmt,K,x,z);

end
