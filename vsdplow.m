function [fL, y, dl, info] = vsdplow(A,b,c,K,x0,y,z0,xu,opts)
%% VSDPLOW - verified lower bound for an SQLP problem
%    [fL y dl info] = vsdplow(A,b,c,K,[],y0)
%    [fL y dl info] = vsdplow(A,b,c,K,x0,y0,z0,xu,opts)
%
%% >> Description: 
%     computes verified lower bound of primal optimal value and rigorous    
%     enclosure of dual strict feasible (near optimal) solution of a conic 
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
% xu: upper bounds of eigenvalues or spectral values of primal optimal
%     solution x
% opts: structure for additional parameter settings:
%     regarded fields:
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
% fL: verified lower bound of the primal optimal value
%  y: an M x 1 vector - rigorous enclosure of dual strict feasible solution
% dl: verified lower bounds of eigenvalues or spectral values of z=c-A'*y
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
% 16/06/12    M. Lange, adapted for new funtions and data format
% 06/08/12    M. Lange, new helper functions for faster computation
%
%%
% TODO: rotated quadratic cones
%

%% input parameter

% check number of input arguments
if nargin<6 || isempty(A) || isempty(b) || isempty(c) || isempty(K)
    error('VSDP:VSDPLOW','more input arguments are required');
elseif nargin<8
    xu = [];  opts = [];
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

% initial output
fL = -Inf;
dl = NaN;
info.iter = 0;

% rounding mode
rnd = getround();
setround(0);

% import data
[A,Arad,b,brad,c,crad,K,x0,y,z0] = import_vsdp(A,b,c,K,x0,y,z0);

% check if approximations are applicable
if isempty(y) || any(isnan(y))
    warning('VSDP:VSDPLOW','not applicable approximations given (NaN)');
    y = NaN;
    return;
elseif nargin<7 || any(isnan(x0)) || any(isnan(z0))
    x0 = [];  z0 = [];
end

% get problem data dimensions
dim3 = length(c);  % dimension
nc = K.l + length(K.q) + length(K.s);  % number of cone constraints

% variable declarations and allocations
I = [];  % basic indices
xu = xu(:);
if isempty(xu)
    xu = inf(K.f+nc,1);
elseif size(xu,1)~=K.f+nc
    error('VSDP:VSDPLOW','upper bound vector has wrong dimension');
end
xuf = xu(1:K.f);  xu(1:K.f) = [];
yrad = sparse(size(y,1),1);  % interval radius y

dl = -inf(nc,1);  % dual lower bounds / trace bounds
epsj = ones(nc,1);  % factor for perturbation amount
ceps = sparse(dim3,1);  % perturbation for c
pertS = cell(length(K.s),1);  % for diagonal perturbations of sdp blocks

% extract free part
Af = A(1:K.f,:);  Afrad = Arad(1:K.f,:);
cf = c(1:K.f);  cfrad = crad(1:K.f);

% create index vector for perturbation entries
pertI = ones(sum(K.s),1);
pertI(cumsum(K.s(1:end-1))+1) = 1 - K.s(1:end-1);
pertI = [ones(K.l+(~isempty(K.q)),1); K.q(1:end-1); cumsum(pertI)];
pertI(1) = K.f + 1;
pertI = cumsum(pertI);


%% Algorithm with finite/infinite primal bounds xu
% **** main loop ****
while info.iter<=ITER_MAX
    info.iter = info.iter + 1;
    setround(1);  % default for rigorous computation in steps 1-3
    
    %% 1.step: defect computation, free variables handling
    if K.f>0 && max(xuf)==inf
        % solve dual linear constraints rigorously
        [y,I] = vuls([],[],struct('mid',Af,'rad',Afrad),...
            struct('mid',cf,'rad',cfrad),[],[],y,I,opts);
        if ~isstruct(y)
            disp('VSDPLOW: could not find solution of dual equations');
            break;
        else
            yrad = y.rad;
            y = y.mid;
        end
    end
    
    % compute rigorous enclosure for z = c - A*y  (with free variables)
    [z,zrad] = spdotK(c,1,A,-y,2);  % for point matrices
    zrad = zrad + crad;  % regard radii of other parameters
    if any(yrad)
        zrad = zrad + abs(A)*yrad;
    end
    if ~isempty(find(Arad,1))
        zrad = zrad + Arad*(abs(y)+yrad);
    end

    defect = 0;  % defect by free variables
    if K.f>0 && max(xuf)<inf  % upper bounds a-priory known
        defect = xuf' * (abs(z(1:K.f))+zrad(1:K.f));
    end
    
    %% 2.step: verified lower bounds on cone eigenvalues
    if K.l>0  % bound for linear variables
        ind = K.f+1:K.f+K.l;
        zjn = zrad(ind) - z(ind);  % -inf(zl)
        % interprete all lp vars as one cone:
        %  if any element is negative all constraints close to zero will be
        %  perturbated
        zj1 = max(zjn);
        if zj1>0  % there is a negative bound
            ind = zjn > -zj1;
            zj1 = max((-1e-13)*min(zjn),zj1); 
            zjn(ind) = min(zjn(ind)+zj1,1.05*zj1);
        end
        dl(1:K.l) = - zjn;
    end
    % eigenvalue bound for second-order cone variables
    ind = K.f + K.l;
    for j = 1:length(K.q)
        ind = ind(end)+2:ind(end)+K.q(j);
        zj1 = zrad(ind(1)-1) - z(ind(1)-1);  % -inf(zq(1))
        zjn = abs(z(ind)) + zrad(ind);  % sup(abs(zq(2:end)))
        zjn = sqrtsup(zjn'*zjn);  % sup(||zq(2:end)||)
        dl(K.l+j) = -(zj1+zjn);  % inf(zq(1)-||zq(2:end)||)
    end
    % eigenvalue bound for semidefinite cone
    blke = K.f + K.l + sum(K.q) + sum(K.s.*(K.s+1))/2;
    for j = length(K.s):-1:1
        ofs = K.l + length(K.q) + j;
        blks = blke - K.s(j)*(K.s(j)+1)/2 + 1;
        [lmin,dl(ofs),pertS{j}] = bnd4sd(struct('mid',z(blks:blke),...
            'rad',zrad(blks:blke)),1,FULL_EIGS_ENCLOSURE);
        if lmin>0
            dl(ofs) = lmin;
        end
        pertS{j} = epsj(ofs) * pertS{j}(:);
        blke = blks - 1;
    end
    
    %% 3.step: cone feasibility check, computing lower bound
    dli = find(dl<0);
    if ~any(isinf(xu(dli)))
        % inf(min(dl,0)*xu + b'*y - defect)
        fL = -(sum(dl(dli)'*(-xu(dli))) + prodsup(-b',y,brad',yrad) + defect);
        y = midrad(full(y),full(yrad));
        setround(rnd);  % reset rounding mode
        return;
    end
    
    %% 4.step: create some perturbed problem
    setround(0);  % no code for rigorous computations
    ind = 1:K.l+length(K.q);
    if isempty(ind)
        ceps = ceps + sparse(pertI,1,cat(1,pertS{:}),dim3,1);
    else
        ceps = ceps + sparse(pertI,1,cat(1,epsj(ind).*min(dl(ind),0),...
            pertS{:}),dim3,1);
    end
    if any(isinf(ceps)) || any(isnan(ceps))
        disp('VSDLOW: perturbation extended range');
        break;
    end
    epsj(dli) = epsj(dli) * (1+ALPHA);  % update perturbation factor
    
    %% 5.step: solve the perturbed problem
    clear dli ind z zrad;  % free some memory before calling solver
    [obj,x0,y,z0,INFO] = mysdps(A,b,c+ceps,K,x0,y,z0,opts);
    % if could not found solution or dual infeasible, break
    if isempty(y) || any(isnan(y)) || any(isinf(y)) || any(INFO(1)==[(1) 2 3])
        disp('VSDLOW: conic solver could not find solution for perturbed problem');
        break;
    end
end   
% **** end of main loop ****

% reset rounding mode
setround(rnd);

% write output
if info.iter==ITER_MAX
    disp('VSDPLOW: maximum number of iterations reached');
end
y = NaN; fL = -Inf; dl = NaN;

%________________________________End VSDPLOW______________________________