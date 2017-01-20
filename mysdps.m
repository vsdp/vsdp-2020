function [obj,x,y,z,info] = mysdps(A,b,c,K,x0,y0,z0,opts)
% MYSDPS  Interface to several conic solvers for VSDP.
%
%   [obj,x,y,z,info] = mysdps(A,b,c,K,x0,y0,z0,opts)
%     Uses a preselected solver to get an approx. solution of a conic
%     problem in the standard primal-dual form:
%
%    (P)  min  c'*x          (D)  max  b'*y
%         s.t. A*x = b           s.t. z := c - A'*y
%              x in K                 z in K*
%
%     where K is a cartesian product of the cones R+, SOCP, PSD.
%
%     For a theoretical introduction into verified conic programming see
%     [Jansson2009].
%
% Input:
% A: nA x m coefficient matrix in SeDuMi or VSDP internal format
% b: a M x 1 vector
% c: a nA x 1 vector, primal objective
% K: a structure with following fields
%     - K.f stores the number of free variables
%     - K.l is the number of nonnegative components
%     - K.q lists the lengths of socp blocks
%     - K.s lists the dimensions of semidefinite blocks
% x0: a nA x 1 vector - a primal feasible (eps-optimal) solution
% y0: a M x 1 vector - a dual feasible (eps-optimal) solution
% z0: a nA x 1 vector - a dual feasible (eps-optimal) solution (slack vars)
% opts: structure for additional parameter settings:
%     regarded fields:
%     	'SOLVER'  to select one of the supported solvers:
%               'sedumi','sdpt3','sdpa','csdp','sdplr', 'lp_solve','linprog'
%    	'USE_STARTING_POINT'  decides whether initial starting point shall
%                             be used or not
%                               -> default: true
%
% Output:
% obj: obj(1) - primal approx. optimal value
%      obj(2) - dual approx. optimal value
% x: a nC x 1 vector - approx. primal optimal solution
% y: a M x 1 vector - approx. dual optimal solution
% z: a nC x 1 vector - approx. dual optimal solution (slack vars)
% info: exit code of conic solver:
%     info = 0: normal termination
%     info = 1: problem indicated to be primal infeasible
%     info = 2: problem indicated to be dual infeasible
%     info = 3: problem indicated to be primal and dual infeasible
%     info = -1: an error occured
%
%   See also vsdpinit.

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

% check input parameter
if nargin<4 || isempty(A) || isempty(b) || isempty(c) || isempty(K)
  error('VSDP:MYSDPS','not enpugh input parameter');
elseif nargin<5
  x0 = [];  y0 = [];  z0 = [];  opts = [];
elseif nargin<6
  y0 = [];  z0 = [];  opts = [];
elseif nargin<7
  z0 = [];  opts = [];
elseif nargin<8
  opts = [];
end

global VSDP_OPTIONS;  % global option settings

solver = [];  % default: allow every solver
if isfield(opts,'SOLVER')
  solver = opts.SOLVER;
elseif isfield(VSDP_OPTIONS,'SOLVER')
  solver = VSDP_OPTIONS.SOLVER;
end

useSTART = true;  % use starting point
if isfield(opts,'USE_STARTING_POINT')
  useSTART = logical(opts.USE_STARTING_POINT);
elseif isfield(VSDP_OPTIONS,'USE_STARTING_POINT')
  useSTART = logical(VSDP_OPTIONS.USE_STARTING_POINT);
end

% import data
[A,~,b,~,c,~,K,x0,y0,z0,IF] = import_vsdp(A,b,c,K,x0,y0,z0);


% get solver number
solverLIST = {'sedumi','sdpt3','sdpa','csdp','sdplr', ...
  'lp_solve','linprog'};
nsolvers = length(solverLIST);
if ischar(solver)
  solverID = find(strncmpi(solver,solverLIST,4));
elseif isnumeric(solver)
  solverID = solver;
else
  solverID = [];
end

% if no solver selected -> allow all solvers
if isempty(solverID)
  solverID = 1:nsolvers;
end

% initialization of default output
obj = [inf -inf];
x = [];  y = [];  z = [];
info = -1;


% call solver for problem

if any(solverID==1) && exist('sedumi','file')==2
  % call sedumi
  
  % options parameter for sedumi
  OPTIONS = [];
  if isfield(opts,'SEDUMI_OPTIONS')
    OPTIONS = opts.SEDUMI_OPTIONS;
  elseif isfield(VSDP_OPTIONS,'SEDUMI_OPTIONS')
    OPTIONS = VSDP_OPTIONS.SEDUMI_OPTIONS;
  end
  % transform data into sedumi format
  [A,c] = vsdp2sedumi(A,c,[],[],K,opts);
  % call sedumi
  [x,y,INFO] = sedumi(A,b,c,K,OPTIONS);
  % transform results to vsdp format
  if ~isempty(x)
    obj(1) = c'*x;
  end
  if ~isempty(y)
    obj(2) = b'*y;
    z = c - A*y;
  end
  [x,z] = export_vsdp(IF,K,x,z);
  info = INFO.pinf + 2*INFO.dinf;
  
elseif any(solverID==2) && exist('sqlp','file')==2
  % call sdpt3
  
  % options parameter for sdpt3
  OPTIONS = [];
  if isfield(opts,'SDPT3_OPTIONS')
    OPTIONS = opts.SDPT3_OPTIONS;
  elseif isfield(VSDP_OPTIONS,'SDPT3_OPTIONS')
    OPTIONS = VSDP_OPTIONS.SDPT3_OPTIONS;
  end
  % use starting point ?
  if ~useSTART || isempty(x0) || isempty(y0) || isempty(z0)
    x0 = [];  y0 = [];  z0 = [];
  end
  % transform data into sdpt3 format
  [blk,A,c,x0,z0] = vsdp2sdpt3(K,A,c,x0,z0,opts);
  % call sdpt3 solver [with intial point]
  [obj,x,y,z,INFO] = sqlp(blk,A,c,b,OPTIONS,x0,y0,z0);
  % transform results to vsdp format + export data
  [~,~,~,x,z] = sdpt2vsdp(blk,[],[],x,z);
  [x,z] = export_vsdp(IF,K,x,z);
  % save info codes
  if isstruct(INFO)
    info = INFO.termcode;  % SDPT3-4.0 output
  else
    info = INFO(1);  % SDPT3-3.x output
  end
  
elseif any(solverID==3) && exist('sdpam','file')==2
  % call sdpa
  
  % check cones
  if K.f>0 || sum(K.q)>0
    error('VSDP:MYSDPS','Second order cone and free variables are not supported by SDPAM');
  end
  % use starting point ?
  if ~useSTART || isempty(x0) || isempty(y0) || isempty(z0)
    x0 = [];  y0 = [];  z0 = [];
  end
  % transform data into sdpam format
  [mDIM,nBLOCK,bLOCKsTRUCT,ct,F,x0,X0,Y0] = vsdp2sdpam(A,b,c,K,x0,y0,z0);
  % options parameter for sdpam
  OPTIONS = [];
  if isfield(opts,'SDPAM_OPTIONS')
    OPTIONS = opts.SDPAM_OPTIONS;
  elseif isfield(VSDP_OPTIONS,'SDPAM_OPTIONS')
    OPTIONS = VSDP_OPTIONS.SDPAM_OPTIONS;
  end
  % call csdp solver [with intial]
  [~,x0,X0,Y0,~] = sdpam(mDIM,nBLOCK,bLOCKsTRUCT,ct,F,x0,X0,Y0,OPTIONS);
  % transform results to vsdp format + export
  [~,~,~,~,x,y,z] = sdpam2vsdp(bLOCKsTRUCT,[],[],x0,X0,Y0);
  if ~isempty(x) && ~any(isnan(x))
    obj(1) = c'*x;
  end
  if ~isempty(y) && ~any(isnan(y))
    obj(2) = b'*y;
  end
  [x,z] = export_vsdp(IF,K,x,z);
  info = 0;
  
  
elseif any(solverID==4) && exist('csdp','file')==2
  % call csdp
  
  % check cones
  if sum(K.q)>0
    error('VSDP:MYSDPS','Second order cone is not supported by CSDP');
  elseif K.f>0
    disp(['Attention: CSDP supports free variables by converting ' ...
      'them to the difference of positive variables.\n' ...
      'The resulting problem is ill-posed.']);
  end
  % transform data into sedumi format which can be read by csdp
  [A,c,x0,z0] = vsdp2sedumi(A,c,x0,z0,K,struct('ALLOW_TRIANGULAR',false));
  % options parameter for csdp
  OPTIONS.check = 0;
  if isfield(opts,'CSDP_OPTIONS')
    OPTIONS = opts.CSDP_OPTIONS;
  elseif isfield(VSDP_OPTIONS,'CSDP_OPTIONS')
    OPTIONS = VSDP_OPTIONS.CSDP_OPTIONS;
  end
  % call csdp solver [with intial]
  if useSTART && ~isempty(x0) && ~isempty(y0) && ~isempty(z0)
    [x,y,z,INFO] = csdp(A,full(b),full(c),K,OPTIONS,full(x0),full(y0),full(z0));
  else
    [x,y,z,INFO] = csdp(A,full(b),full(c),K,OPTIONS);
  end
  % transform results to vsdp format + export
  if ~isempty(x) && ~any(isnan(x))
    obj(1) = c'*x;
  end
  if ~isempty(y) && ~any(isnan(y))
    obj(2) = b'*y;
  end
  [x,z] = export_vsdp(IF,K,x,z);
  if any(INFO==[0 1 2])
    info = INFO;
  end
  
elseif any(solverID==5) && exist('sdplr','file')==2
  % call sdplr
  
  % check cones
  if K.f>0 || sum(K.q)>0
    error('VSDP:MYSDPS','Second order cone and free variables are not supported by SDPLR');
  end
  % options parameter for sdplr
  OPTIONS = [];
  if isfield(opts,'SDPLR_OPTIONS')
    OPTIONS = opts.SDPLR_OPTIONS;
  elseif isfield(VSDP_OPTIONS,'SDPLR_OPTIONS')
    OPTIONS = VSDP_OPTIONS.SDPLR_OPTIONS;
  end
  % transform data into sedumi format, which can be read by sdplr
  [A,c,x0] = vsdp2sedumi(A,c,x0,[],K,struct('ALLOW_TRIANGULAR',false));
  % call csdp solver [with intial]
  if useSTART && ~isempty(x0) && ~isempty(y0)
    [x,y] = sdplr(A,b,c,K,OPTIONS,[],x0,y0);
  else
    [x,y] = sdplr(A,b,c,K,OPTIONS);
  end
  % transform results to vsdp format
  if ~isempty(x)
    obj(1) = c'*x;
  end
  if ~isempty(y)
    obj(2) = b'*y;
    z = c - A*y;
  end
  [x,z] = export_vsdp(IF,K,x,z);
  info = 0;
  
elseif any(solverID==6) && exist('lp_solve','file')==2
  % call lp_solve
  
  % test of not supported cones,
  if sum(K.q)>0 || sum(K.s)>0
    error('VSDP:MYSDPS','Second order and semidefinite cones are not supported by LPSOLVE');
  end
  % call lp_solve solver
  [~,x,y,stat] = lp_solve(-full(c),A',full(b),zeros(length(b),1),[-inf(K.f,1);zeros(K.l,1)]);
  if ~isempty(x)
    obj(1) = c'*x;
  end
  if ~isempty(y)
    y = -y;
    obj(2) = b'*y;
    z = c - A*y;
  end
  info = (stat==2) + 2*(stat==3);
  
elseif any(solverID==7) && exist('linprog','file')==2
  % call linprog
  
  % options parameter for linprog
  OPTIONS = optimset('LargeScale','on');
  if isfield(opts,'LINPROG_OPTIONS')
    OPTIONS = opts.LINPROG_OPTIONS;
  elseif isfield(VSDP_OPTIONS,'LINPROG_OPTIONS')
    OPTIONS = VSDP_OPTIONS.LINPROG_OPTIONS;
  end
  
  % test of not supported cones,
  if sum(K.q)>0 || sum(K.s)>0
    error('VSDP:MYSDPS','Second order and semidefinite cones are not supported by LINPROG');
  end
  % use starting point ?
  if ~useSTART || isempty(x0)
    x0 = [];
  end
  % call linprog solver [with intial point]
  [x,~,flag,~,lambda] = ...
    linprog(full(c),[],[],A',full(b),[-inf(K.f,1); zeros(K.l,1)],inf(length(c),1),x0,OPTIONS);
  % transform results to vsdp format
  if ~isempty(x)
    y = -lambda.eqlin;
    z = c - A*y;
    obj(1) = c'*x;
    obj(2) = b'*y;
  end
  info = [1 2] * (abs(flag)==[2; 3]);
  
else
  % no solver called
  warning('VSDP:MYSDPS','Solver not found!');
  disp('The selected solver could not be detected!');
end

end
