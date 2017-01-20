function [fU x lb info] = vsdpup(A,b,c,K,x,y0,z0,yu,opts)
%% VSDPLOW - verified upper bound for semidefinite-quadratic-linear programming
%    [fU x lb info] = vsdpup(A,b,c,K,x0,y0,z0)
%
%% >> Description: 
%     computes verified upper bound of primal optimal value and rigorous    
%     enclosure of primal strict feasible (near optimal) solution of a conic 
%     problem in the standard primal-dual form:
%
%    (P)  min  c'*x          (D)  max  b'*y
%         s.t. A*x = b            s.t. z := c - A'*y
%              x in K                  z in K*
% 
%     where K is a cartesian product of the cones R+, SOCP, PSD.
%     For theoretical introduction into verified conic programming see: 
%     C. Jansson. On Verified Numerical Computations in Convex Programming.
%     Japan J. Indust. Appl. Math., 26:337–363, 2009
%
%% >> Input:
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
% yu: upperbounds of the absolute value of dual optimal solution y
% opts: structure for additional parameter settings:
%     fields:
%     	'ITER_MAX'   maximum number of iterations that can be used to
%                    perturbate the problem and find a feasible solution
%    	'ALPHA'   growing factor for problem perturbation
%                       -> default: 0.5
%     	'FULL_EIGS_ENCLOSURE'  if true code for stronger complete
%                              eigenvalue enclosure will be applied
%                              the code is a little bit slower
%                                   -> default: true
%     	'SOLVER'  to select one of the supported solvers:
%               'sedumi','sdpt3','sdpa','csdp','sdplr', 'lp_solve','linprog'
%    	'USE_STARTING_POINT'  decides whether initial starting point shall
%                             be used or not
%                               -> default: false
%
%% >> Output:
% fU: verified upper bound of the primal optimal value
%  x: a nA x 1 vector - rigorous enclosure of primal strict feasible solution
% lb: verified lower bounds of eigenvalues or spectral values of X with
%     respect to K
% info.iter: number of iterations
%

%% ********************************************************************* %%
%% This file is part of VSDP by V. Haerter, C. Jansson and M. Lange      %%
%% Copyright (c) 2012, C. Jansson                                        %%
%%                     Technical University of Hamburg (TUHH)            %%
%%                     Institute for Reliable Computing (IRC)            %%
%% VSDP can be freely used for private and academic purposes.            %%
%% Commercial use or use in conjunction with a commercial program which  %%
%% requires VSDP or part of it to function properly is prohibited.       %%
%% ********************************************************************* %%

%% Last modified:  
% 31/07/10    V. Haerter, comments added
% 08/18/12    M. Lange, complete rewrite
% 18/08/12    M. Lange, preconditioning to find better basis indices
%
%%
% TODO: rotated quadratic cones
%


%% input parameter

% check number of input arguments
if nargin<5 || isempty(A) || isempty(b) || isempty(c) || isempty(K)
    error('VSDP:VSDPUP','more input arguments are required');
elseif nargin<6
    y0 = [];  z0 = [];  yu = [];  opts = [];
elseif nargin<7
    z0 = [];  yu = [];  opts = [];
elseif nargin<8
    yu = [];  opts = [];
elseif nargin<9
    opts = [];
end

global VSDP_OPTIONS;  % global options structure

% use starting point - default: false
if ~isfield(VSDP_OPTIONS,'USE_STARTING_POINT')
    VSDP_OPTIONS.USE_STARTING_POINT = false;
end
if ~isfield(opts,'USE_STARTING_POINT')
    opts.USE_STARTING_POINT = VSDP_OPTIONS.USE_STARTING_POINT;
end

% parameter list
optLIST = {'ITER_MAX','ALPHA','FULL_EIGS_ENCLOSURE'};
optfs = isfield(opts,optLIST);
goptfs = isfield(VSDP_OPTIONS,optLIST);

% read parameters
[ITER_MAX, ALPHA, FULL_EIGS_ENCLOSURE] = deal(10, 0.5, true);
for i = 1:length(optLIST)
    if optfs(i)
        eval([optLIST{i},' = opts.',optLIST{i},';']);
    elseif goptfs(i)
        eval([optLIST{i},' = VSDP_OPTIONS.',optLIST{i},';']);
    end
end


%%  Preliminary steps / Prealocations

% default output
fU = Inf;
lb = NaN;
info.iter = 0;

% rounding mode
rnd = getround();
setround(0);

% import data
[A,Arad,b,brad,c,crad,K,x,y0,z0,IF] = import_vsdp(A,b,c,K,x,y0,z0);

% check if approximations are applicable
if isempty(x) || any(isnan(x))
    warning('VSDP:VSDPUP','not applicable approximations given (NaN)');
    x = NaN;
    return;
elseif any(isnan(y0)) || any(isnan(z0))
    y0 = [];  z0 = [];
end

% get problem data dimensions
m = length(b);  % number of linear system constraints
dim3 = length(x);
nc = K.l + length(K.q) + length(K.s);  % number of cone constraints

% check upper dual bound
yu = yu(:);
if isempty(yu)
    yu = inf(m,1);
elseif length(yu)~=m
    error('VSDP:VSDPUP','upper bound vector has wrong dimension');
end


%% Algorithm with finite dual bounds yu, projection into cone
if max(yu)<inf
    setround(1);  % default rounding mode for verification code
    x = full(x);  % faster for verification code
    
    % projection of linear part
    x(K.f+1:K.f+K.l) = max(x(K.f+1:K.f+K.l),0);
    
    % projection of socp part
    blke = K.f + K.l;
    for j = 1:length(K.q)
        xj1 = x(blke+1);  % xq(1)
        xjn = x(blke+2:blke+K.q(j));  % xq(2:end)
        xjn = sqrtsup(xjn'*xjn);  % sup(||xq(2:end)||)
        x(blke+1) = max(xj1,xjn);  % very simple projection
        blke = blke + K.q(j);
    end
    
    % projection of sdp part, force to psd cone
    blke = K.f + K.l + sum(K.q);
    for j = 1:length(K.s)
        dind = cumsum(1:K.s(j));  % index for diagonal entries
        vx = x(blke+1:blke+dind(end));
        vx(dind) = 2*vx(dind);  % regard mu=2 for x
        lambdaj = bnd4sd(vx);
        if lambdaj<0  % simple projection into feasible space
            x(dind+blke) = x(dind+blke) + (-lambdaj)/2;
        end
        blke = blke + dind(end);
    end
    
    % x <in> K:  now regard defect = |A*x-b|
    defect = resmag(x',A,b',1,0,Arad,brad',0);
    % fU = sup(x'*c + defect*yu)
    fU = prodsup(x',c,0,crad) + defect*yu;
    x = NaN;  lb = NaN;
    % reset rounding mode and return
    setround(rnd);
    return;
end
    
%% Algorithm with infinite dual bounds yu
% variable declarations and allocations
x0 = x;  % starting point
lb = -inf(nc,1);  % dual lower bounds
epsj = ones(nc,1);  % factor for perturbation amount
xeps = sparse(length(x),1);  % wanted perturbation for x
pertS = cell(length(K.s),1);  % for diagonal perturbations of sdp blocks
I = [];  % basic indices

% create index vector for perturbation entries
pertI = ones(sum(K.s),1);
pertI(cumsum(K.s(1:end-1))+1) = 1-K.s(1:end-1);
pertI = [ones(K.l+(~isempty(K.q)),1); K.q(1:end-1); cumsum(pertI)];
pertI(1) = K.f + 1;
pertI = cumsum(pertI);


%% **** main loop ****
while info.iter<=ITER_MAX
    info.iter = info.iter + 1;
    setround(1);  % default rounding for verification part
    
    %% 1.step: compute rigorous enclosure for x
    [x, I] = vuls([],[],struct('mid',A,'rad',Arad),...
        struct('mid',b,'rad',brad),[],[],x,I,opts);
    if ~isstruct(x)
        disp('VSDPUP: could not find solution of primal equations');
        break;
    else
        xrad = sparse(x.rad);
        x = full(x.mid);
    end
    
    %% 2.step: verified lower bounds on cone eigenvalues
    if K.l>0  % bound for linear variables
        ind = K.f+1:K.f+K.l;
        xjn = xrad(ind) - x(ind);  % -inf(xl)
        % interprete all lp vars as one cone:
        %  if any element is negative all constraints close to zero will be
        %  perturbated
        xj1 = max(xjn);
        if xj1>0  % there is a negative bound
            ind = xjn > -0.5*xj1;
            xjn(ind) = xjn(ind) + xj1;
        end
        lb(1:K.l) = - xjn;
    end
    % eigenvalue bound for second-order cone variables
    ind = K.f + K.l;
    for j = 1:length(K.q)
        ind = ind(end)+2:ind(end)+K.q(j);
        xj1 = xrad(ind(1)-1) - x(ind(1)-1);  % -inf(xq(1))
        xjn = abs(x(ind)) + xrad(ind);  % sup(abs(xq(2:end)))
        xjn = sqrtsup(xjn'*xjn);  % sup(||xq(2:end)||)
        lb(K.l+j) = -(xj1+xjn);  % inf(xq(1)-||xq(2:end)||)
    end
    % eigenvalue bound for semidefinite cone
    blke = K.f + K.l + sum(K.q) + sum(K.s.*(K.s+1))/2;
    for j = length(K.s):-1:1
        ofs = K.l + length(K.q) + j;
        dind = cumsum(1:K.s(j));  % index for diagonal entries
        blk3 = dind(end);
        vx = struct('mid',0.5*x(blke-blk3+1:blke),...
            'rad',0.5*xrad(blke-blk3+1:blke));
        vx.mid(dind) = 2 * vx.mid(dind);  % regard mu=2 for x
        vx.rad = vx.rad + vx.rad.*sparse(dind,1,1,blk3,1);
        [lmin,lb(ofs),pertS{j}] = bnd4sd(vx,1,FULL_EIGS_ENCLOSURE);
        if lmin>0
            lb(ofs) = lmin;
        end
        pertS{j} = epsj(ofs) * pertS{j}(:);
        blke = blke - blk3;
    end
    
    %% 3.step: cone feasibility check, computing upper bound
    lbi = lb<0;
    if ~any(lbi)
        fU = prodsup(x',c,xrad',crad);  % sup(x'*c)
        x = export_vsdp(IF,K,midrad(full(x),full(xrad)));
        setround(rnd);
        return;
    end
    
    %% 4.step: create some perturbed problem
    setround(0);  % no code for rigorous computations
    ind = 1:K.l+length(K.q);
    if isempty(ind)
        xeps = xeps + sparse(pertI,1,cat(1,pertS{:}),dim3,1);
    else
        xeps = xeps + sparse(pertI,1,cat(1,epsj(ind).*min(lb(ind),0),...
            pertS{:}),dim3,1);
    end
    if any(isinf(xeps)) || any(isnan(xeps))
        disp('VSDPUP: perturbation extended range');
        break;
    end
    epsj(lbi) = epsj(lbi) * (1+ALPHA);  % update perturbation factor
        
    %% 5.step: solve the perturbed problem
    clear lbi ind vx x xrad;  % free some memory before calling solver
    [objt,x0,y0,z0,INFO] = mysdps(A,b+(xeps'*A)',c,K,x0,y0,z0,opts);
    % if could not find solution or primal infeasible: break
    if isempty(x0) || any(isnan(x0) | isinf(x0)) || any(INFO(1)==[1 (2) 3])
        disp('VSDPUP: could not find solution for perturbed problem');
        break;
    end
    x = x0 - xeps;  % undo perturbation 
end
% **** end of main loop ****
setround(rnd);

% write output
if info.iter==ITER_MAX
    disp('VSDPUP: maximum number of iterations reached');
end
fU = Inf; x = NaN; lb = NaN;

%________________________________End VSDPUP_______________________________