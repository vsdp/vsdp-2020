function [infeas x y] = vsdpinfeas(A,b,c,K,choose,x0,y0,z0,opts)
%% VSDPINFEAS - infeasibility check for semidefinite-quadratic-linear programming
%    [infeas x y] = vsdpinfeas(A,b,C,K,choose,x0,y0,z0)
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
% choose: 'p' or 'd' charachter for primal reps. dual infeasibility check
% x0: a nA x 1 vector - a primal feasible (eps-optimal) solution
% y0: a M x 1 vector - a dual feasible (eps-optimal) solution
% z0: a nA x 1 vector - a dual feasible (eps-optimal) solution (slack vars)
% opts: structure for additional parameter settings:
%     fields:
%     	'FULL_EIGS_ENCLOSURE'  if true code for stronger complete
%                              eigenvalue enclosure will be applied
%                              the code is a little bit slower
%                                   -> default: true
%
%% >> Output:
% infeas: 1 if primal infeasibility is proved
%         0 if couldn't prove infeasibility
%        -1 if dual infeasibility is proved
% x: a nC x 1 vector - a rigorous Farkas solution if dual  
%    infeasibility could be proved (primal unbounded ray)
% y: a M x 1 vector - a rigorous Farkas solution if primal
%    infeasiblity could be proved (dual unbounded ray)
%
%% EXAMPLE:
% EPS = -0.01; DELTA = 0.1;
%  min <[0 0; 0 0],X>
%  s.t. <[1 0; 0 0],X> = [EPS];
%       <[0 1; 1 DELTA],X> = [1];
%       X is PSD
%
% b = [EPS; 1];
% c = zeros(4,1);
% A = [1 0 0; 0 1 DELTA]';
% choose = 'p'; 
%
% [objt,x0,y0,z0,info] = mysdps(A,b,c,K);            
% [isinfeas x y] = vsdpinfeas(blk,A,c,b,choose,x0,y0,z0);
% isinfeas =
%      1
% y =
%  -100.5998
%    -0.0060
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
% 08/09/12    M. Lange, rewrite for new format 
%
%% ToDo
% - rotated quadratic cones
%


%% check input
if nargin<7
    error('VSDP:VSDPINFEAS','to few input arguments');
elseif nargin<9
    opts = [];
end

% default output
infeas = 0;
x = NaN;  y = NaN;

% check if approximations are applicable
if choose=='p' && (isempty(y0) || any(isnan(y0)))
    disp('VSDPINFEAS: not applicable approximations given (NaN)');
    return;
elseif choose=='d' && (isempty(x0) || any(isnan(x0)))
    disp('VSDPINFEAS: not applicable approximations given (NaN)');
    return;
end

% global options structure
global VSDP_OPTIONS;

% read option parameter
FULL_EIGS_ENCLOSURE = true;
if isfield(opts,'FULL_EIGS_ENCLOSURE')
    FULL_EIGS_ENCLOSURE = opts.FULL_EIGS_ENCLOSURE;
elseif isfield(VSDP_OPTIONS,'FULL_EIGS_ENCLOSURE')
    FULL_EIGS_ENCLOSURE = VSDP_OPTIONS.FULL_EIGS_ENCLOSURE;
end

% rounding mode
rnd = getround();
setround(0);

% import data
[A,Arad,b,brad,c,crad,K,x0,y0,z0,IF] = import_vsdp(A,b,c,K,x0,y0,[]);

        
%% check primal infeasibility
if choose=='p'
    
    %% 1.step: rigorous enclosure for z = -A*y
    if K.f>0
        % solve dual linear constraints rigorously
        y0 = vuls([],[],struct('mid',A(1:K.f,:),'rad',Arad(1:K.f,:)),...
            zeros(K.f,1),[],[],y0,[],opts);
        if ~isstruct(y0)
            disp('VSDPINFEAS: could not verify primal infeasibility');
            setround(rnd);
            return;
        else
            y0rad = y0.rad;
            y0 = y0.mid;
        end
    else
        y0rad = sparse(size(y0,1),1);
    end
    % default rounding mode for verification
    setround(1);
    % z = -A*y  (with free variables)
    [z,zrad] = resmidrad(0,0,A,y0,0,0,Arad,y0rad);
    
    %% 2.step: check positive cost
    alpha = prodsup(-b',y0,brad',y0rad);
    if alpha>=0
        disp('VSDPINFEAS: could not verify primal infeasibility');
        setround(rnd);
        return;
    end
    
    %% 3.step: verified lower bounds on cone eigenvalues
    dl = inf;
    % bound for linear variables
    if K.l>0
        ind = K.f+1:K.f+K.l;
        dl = - max(zrad(ind) - z(ind));
    end
    % eigenvalue bound for second-order cone variables
    ind = K.f + K.l;
    for j = 1:length(K.q)
        ind = ind(end)+2:ind(end)+K.q(j);
        zj1 = zrad(ind(1)-1) - z(ind(1)-1);  % -inf(zq(1))
        zjn = abs(z(ind)) + zrad(ind);  % sup(abs(zq(2:end)))
        zjn = sqrtsup(zjn'*zjn);  % sup(||zq(2:end)||)
        dl = min(dl,-(zj1+zjn));  % inf(zq(1)-||zq(2:end)||)
    end
    % finish if dl<0
    if dl<0
        disp('VSDPINFEAS: could not verify primal infeasibility');
        setround(rnd);
        return;
    end
    % eigenvalue bound for semidefinite cone
    blke = K.f + K.l + sum(K.q) + sum(K.s.*(K.s+1))/2;
    for j = length(K.s):-1:1
        blks = blke - K.s(j)*(K.s(j)+1)/2 + 1;
        lambdaMin = bnd4sd(struct('mid',z(blks:blke),...
            'rad',zrad(blks:blke)),1,FULL_EIGS_ENCLOSURE);
        if lambdaMin<0
            disp('VSDPINFEAS: could not verify primal infeasibility');
            setround(rnd);
            return;
        end
        dl = min(dl,lambdaMin);
        blke = blks - 1;
    end
    
    %% 4.step: final feasibility check, compute result
    if dl>=0
        infeas = 1;
        if any(y0rad)
            y = midrad(full(y0),full(y0rad));
        else
            y = y0;
        end
    end
    setround(rnd);
    return;


%% check dual infeasibility
elseif choose=='d'
    
    %% 1.step: check c'*x > 0
    [alpha alpharad] = spdotK(c,x0,3);
    setround(1);  % default rounding mode for verification
    alpharad = alpharad + crad'*x0;
    if alpharad>=-alpha
    	disp('VSDPINFEAS: could not verify dual infeasibility');
    	setround(rnd);
        return;
    end
    
    %% 2.step: inclusion for x such that A*x = 0
    vx = vuls([],[],struct('mid',A,'rad',Arad),0*b,[],[],x0,[],opts);
    if ~isstruct(vx)
    	disp('VSDPINFEAS: could not verify dual infeasibility');
    	setround(rnd);
        return;
    elseif prodsup(vx.mid',c,full(vx.rad)',crad)>=0  % c'*x < 0 ?
        % find inclusion such that c'x = alpha
        A = cat(2,A,c);
        Arad = cat(2,Arad,crad);
        b(end+1) = alpha;
        brad(end+1) = 0;
        vx = vuls([],[],struct('mid',A,'rad',Arad),...
            struct('mid',b,'rad',brad),[],[],x0,[],opts);
        if ~isstruct(vx)
            disp('VSDPINFEAS: could not verify dual infeasibility');
            setround(rnd);
            return;
        end
    end
    x0rad = sparse(vx.rad);
    x0 = full(vx.mid);
    
    %% 3.step: verified lower bounds on cone eigenvalues
    lb = inf;
    % bound for linear variables
    if K.l>0
        ind = K.f+1:K.f+K.l;
        lb = - max(x0rad(ind) - x0(ind));
    end
    % eigenvalue bound for second-order cone variables
    ind = K.f + K.l;
    for j = 1:length(K.q)
        ind = ind(end)+2:ind(end)+K.q(j);
        xj1 = x0rad(ind(1)-1) - x0(ind(1)-1);  % -inf(xq(1))
        xjn = abs(x0(ind)) + x0rad(ind);  % sup(abs(xq(2:end)))
        xjn = sqrtsup(xjn'*xjn);  % sup(||xq(2:end)||)
        lb = min(lb,-(xj1+xjn));  % inf(xq(1)-||xq(2:end)||)
    end
    % finish if lb<0
    if lb<0
        disp('VSDPINFEAS: could not verify dual infeasibility');
        setround(rnd);
        return;
    end
    % eigenvalue bound for semidefinite cone
    blke = K.f + K.l + sum(K.q) + sum(K.s.*(K.s+1))/2;
    for j = length(K.s):-1:1
        dind = cumsum(1:K.s(j));  % index for diagonal entries
        blk3 = dind(end);
        vx = struct('mid',0.5*x0(blke-blk3+1:blke),...
            'rad',0.5*x0rad(blke-blk3+1:blke));
        vx.mid(dind) = 2 * vx.mid(dind);  % regard mu=2 for x
        vx.rad = vx.rad + vx.rad.*sparse(dind,1,1,blk3,1);
        lambdaMin = bnd4sd(vx,1,FULL_EIGS_ENCLOSURE);
        if lambdaMin<0
            disp('VSDPINFEAS: could not verify dual infeasibility');
            setround(rnd);
            return;
        end
        lb = min(lb,lambdaMin);
        blke = blke - blk3;
    end
    
    %% 4.step: final feasibility check, compute result
    if lb>=0
        infeas = -1;
        x = export_vsdp(IF,K,midrad(full(x0),full(x0rad)));
    end
    setround(rnd);
    return;


%% wrong "choose" parameter
else
    
    error('VSDP:VSDPINFEAS','"choose" must be p or d');
    
end


%______________________________End VSDPINFEAS_____________________________