function [A,b,c,K,pd] = lp2vsdp(A,b,c,e,lb,ub)
%% LP2VSDP:  transforms problem data from LP_SOLVE format to VSDP3 format 
%    [At,b,c,K] = lp2vsdp(A,b,c,e,lb,ub)
%
%% >> Description:
% The LP_SOLVE solver solves MILP problems of the form:
%    max c'*x
%    s.t. A*x <> b
%         lb <= x <= ub
% Function converts this form into the standard primal/dual form of conic
% problems treated by VSDP3. Note, that if the problem gets dualized the
% optimal value must be negated to get the right sign.
% 
%% >> Input:
% A: m x n matrix representing linear constraints
% b: m x 1 vector of right sides for the inequality constraints
% c: n x 1 vector of coefficients for primal objective function
% e: m x 1 flag vector that determines the sense of the inequalities
%    -1 - less than
%     0 - equality
%     1 - greater than
%
%% >> Output: 
% A:  Matrix of linear equations
% b, c: - coefficients of primal/dual objective functions
% K: K.f, K.l, K.q, K.s converted to column vectors
% pd: a scalar with one of 'p' or 'd'. If 'p' problem was not dualized.
%     If 'd' problem was dualized and optimal value must be negated. 
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
% 31/08/10    V. Haerter, comments added
% 09/08/12    M. Lange, speed improvement and code reduction
% 10/08/12    M. Lange, extended primal form conversion
%
%%
% TODO: improvements:  - computations in primal form with interval
%                        inclusion
%                      - mirroring of variables <in> [-inf 0] to R+
%

%% check input
b = b(:);
c = c(:);
[m n] = size(A);
if m~=length(b) || n~=length(c)
    error('VSDP:LP2VSDP','c or b not compatible with A');
elseif m~=length(e)
    error('VSDP:LP2VSDP','dimension of flag vector "e" does not match');
end
if nargin<5
    lb = -inf(n,1);
    ub = inf(n,1);
elseif nargin<6 
    ub = inf(n,1);
end

% resolve fixed bounds
ind = find(lb==ub);  % index for fixed bounds
nfx = length(ind);  % number of fixed bounds
if nfx>0
    A = [A; sparse(1:nfx,ind,1,nfx,n)];
    e = [e; zeros(nfx,1)];
    b = [b; lb(ind)];
    lb(ind) = -inf;
    ub(ind) = inf;
    m = m + nfx;
end

% indices
indl = find(e<0);  % inequalities less than  (lower)
indu = find(e>0);  % inequalities greater than  (upper)

Ai = [A(indl,:); -A(indu,:)];  % part of coefficient matrix for inequalities
A([indl; indu],:) = [];  % part of coefficient matrix for equality constraints

bi = [b(indl); -b(indu)];
b([indl; indu]) = [];
mi = length(bi);

indl = find(lb>-inf);  % index for (applied) lower bounds
nl = length(indl);  % number of lower bounds

if all(lb(indl)==0) && min(ub)==inf && n+mi<m+nl
    % primal form with free + slack variables
    ind = [find(lb==-inf); indl];
    A = [A(:,ind)' Ai(:,ind)'; ...
        sparse(mi,m-mi) speye(mi)];  % free vars first + add slack vars
    b = [b; bi];
    c = [-c(ind); sparse(mi,1)];
    K.f = n - nl;
    K.l = nl + mi;
    pd = 'd';
else
    % extend dual form
    indu = find(ub<inf);
    nu = length(indu);
    A = [A; Ai; sparse(1:nl,indl,-1,nl,n); sparse(1:nu,indu,1,nu,n)];
    Ai = c;  % Ai used as place-holder
    c = [b; bi; -lb(indl); ub(indu)];
    b = Ai;
    K.f = m - mi;
    K.l = mi + nl + nu;
    pd = 'p';
end

%___________________________End of LP2VSDP____________________________