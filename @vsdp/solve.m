function obj = solve (obj, solver)
% SOLVE  Approximately solve the given conic problem instance.
%
%   See 'help vsdp.vsdp' for a description of the conic problem instance.
%
%   See 'help vsdp_options.SOLVER' for a list of supported solvers.
%
%   To use an initial guess (x0,y0,z0) type:
%   
%      obj.x = x0;
%      obj.y = y0;
%      obj.z = z0;
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
%       At = [A1(:),  A2(:)];  % Vectorize data
%       c  = C(:);
%       b  = [1; 2.0001];
%       obj = vsdp(At, b, c, K).solve()
%       
%
%   See also vsdpinit, vsdpup, vsdplow, vsdpinfeas.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (1, 2);
if (nargin == 2)
  obj.options.SOLVER = solver;
end

% Dispatch to solver.
switch (obj.options.SOLVER)
  case 'sdpt3'
    obj = obj.solve_sdpt3 ();
  case 'sedumi'
    obj = obj.solve_sedumi ();
  case 'sdpa'
    obj = obj.solve_sdpa ();
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
    if (~obj.options.VERBOSE_OUTPUT)
      OPTIONS.printlevel = 0;
    end
    OPTIONS.check = 0; % check problem data for symmetry.
    [A,c,x0,z0] = vsdp2sedumi(A,c,x0,z0,K); % that format can be read by csdp
    
    % Call solver
    if (obj.options.USE_STARTING_POINT && ~isempty(x0) && ~isempty(y0) ...
        && ~isempty(z0))
      [x,y,z,INFO] = csdp(A,full(b),full(c),K,OPTIONS,full(x0),full(y0), ...
        full(z0));
    else
      [x,y,z,INFO] = csdp(A,full(b),full(c),K,OPTIONS);
    end
    
    % Transform results to VSDP format
    if (~isempty(x) && ~isempty(y))
      objt = [c'*x, b'*y];
    end
    if any(INFO == [0, 1, 2])
      info = INFO;
    end
    
  case 'sdplr'
    % Check cones
    if ((obj.K.f > 0) || (sum(obj.K.q) > 0))
      error ('VSDP:solve:unsupportedCones', ...
        ['solve: Free variables (K.f) and second-order cones (K.q) ', ...
        'are not supported by ''%s'''], obj.options.SOLVER);
    end
    
    % Setup
    [A,c,x0] = vsdp2sedumi(A,c,x0,[],K); % that format can be read by sdplr
    
    % Call solver
    if (obj.options.USE_STARTING_POINT && ~isempty(x0) && ~isempty(y0))
      [x,y] = sdplr(A,b,c,K,OPTIONS,[],x0,y0);
    else
      [x,y] = sdplr(A,b,c,K,OPTIONS);
    end
    
    % Transform results to VSDP format
    if (~isempty(x) && ~isempty(y))
      z = c - A*y;
      objt = [c'*x, b'*y];
      info = 0;
    end
    
  case 'lp_solve'
    % Check cones
    if ((sum(obj.K.q) > 0) || (sum(obj.K.s) > 0))
      error ('VSDP:solve:unsupportedCones', ...
        ['solve: Second-order cones (K.q) and semidefinite cones (K.s) ', ...
        'are not supported by ''%s'''], obj.options.SOLVER);
    end
    lbound = [-inf(K.f,1); zeros(K.l,1)]; % lower bound
    
    % Call solver
    [~,x,y,stat] = lp_solve ...
      (-full(c), A', full(b), zeros(length(b),1), lbound);
    % Transform results to VSDP format
    if (~isempty(x) && ~isempty(y))
      y = -y;
      z = c - A*y;
      objt = [c'*x, b'*y];
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
    if ((sum(obj.K.q) > 0) || (sum(obj.K.s) > 0))
      error ('VSDP:solve:unsupportedCones', ...
        ['solve: Second-order cones (K.q) and semidefinite cones (K.s) ', ...
        'are not supported by ''%s'''], obj.options.SOLVER);
    end
    
    % Setup
    OPTIONS = optimoptions ('linprog', 'Algorithm', 'interior-point-legacy');
    if (~obj.options.VERBOSE_OUTPUT)
      OPTIONS = optimoptions(OPTIONS, 'Display', 'off');
    end
    lbound = [-inf(K.f,1); zeros(K.l,1)]; % lower bound
    ubound = inf(length(c),1);            % upper bound
    
    % Call solver
    [x,~,flag,~,lambda] = linprog ...
      (full(c), [], [], A', full(b), lbound, ubound, x0, OPTIONS);
    
    % Transform results to VSDP format
    if (~isempty(x))
      y = -lambda.eqlin;
      z = c - A*y;
      objt = [c'*x, b'*y];
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
    if ((sum(obj.K.q) > 0) || (sum(obj.K.s) > 0))
      error ('VSDP:solve:unsupportedCones', ...
        ['solve: Second-order cones (K.q) and semidefinite cones (K.s) ', ...
        'are not supported by ''%s'''], obj.options.SOLVER);
    end
    
    % Setup
    sense = 1;  % Minimization
    lbound = [-inf(K.f,1); zeros(K.l,1)]; % lower bound
    ubound = inf(length(c),1);            % upper bound
    vtype = repmat ('C', length(c), 1);   % continuous variable
    ctype = repmat ('S', length(b), 1);   % equality constraint
    
    % Call solver
    [x,~,errnum,extra] = glpk ...
      (full(c), A', full(b), lbound, ubound, ctype, vtype, sense);
    
    % Transform results to VSDP format
    if (~isempty(x))
      y = extra.lambda;
      z = c - A*y;
      objt = [c'*x, b'*y];
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
    error ('VSDP:solve:unsupportedCones', ...
      'solve: The solver %s could not be detected.', obj.options.SOLVER);
end

end
