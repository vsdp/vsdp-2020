function [mDIM,nBLOCK,bLOCKsTRUCT,c,F,x0,X0,Y0] = vsdp2sdpam(A,b,c,K,x,y,z)
%% VSDP2SDPT3:  transforms problem data from VSDP3 to SDPAM format
%    [mDIM,nBLOCK,bLOCKsTRUCT,c,F,x0,X0,Y0] = vsdp2sdpt3(At,b,c,K)
%
%% >> Input:
% A: a nA3 x M Matrix,
%     whereas nA = dimf+diml+dimq+dims3
%     dimf: number of free variables: dimf = sum(K.f>0)
%     diml: number of nonnegative variables: diml = sum(K.l>0)
%     dimq: sum of all socp variables: dimq = sum_i(K.q(i))
%     dims3: sum of all sdp variables: dims3 = sum_i(K.s(i)*(K.s(i)+1)/2)
% b: M x 1 vector - right hand side of linear constraints
% c: nA3 x 1 vector - primal objective function
% K: a structure with following fields
%     - K.f stores the number of free variables
%     - K.l is the number of nonnegative components
%     - K.q lists the lengths of socp blocks
%     - K.s lists the dimensions of semidefinite blocks
% x: a nA3 x 1 vector - approx. primal optimal solution  - svec(mu=2)
% y: a M x 1 vector - approx. dual optimal solution
% z: a nA3 x 1 vector - approx. dual optimal solution (slack vars)
%
%% >> Output:
% mDIM: scalar - number of primal variables
% nBLOCK: scalar - number of blocks of F
% bLOCKsTRUCT: scalar - represents the block structure of F
% c: vector - coefficient vector
% F: cell array of coefficient matrices
% x0: vector for initial solution
% X0,Y0: cell arrays of initial points
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
% 08/08/11    V. Haerter, corrected errors in transforming linear blocks
% 14/08/12    M. Lange, rewrite for improved performance and new functions
%


%% check input
if nargin<4 || ~isstruct(K) || isempty(A) || isempty(b) || isempty(c)
    error('VSDP:VSDP2SDPAM','the given problem is incomplete');
elseif nargin<5
    x = [];  y = [];  z = [];
elseif nargin<6
    y = [];  z = [];
elseif nargin<7
    z = [];
end

c = c(:);  x = x(:);  y = y(:);  z = z(:);


%% prepare data

% prepare K
fields = isfield(K,{'f','l','q','s'});
if fields(1) && sum(K.f)>0
    error('VSDP:VSDP2SDPAM','free variables are not supported by SDPA');
elseif fields(3) && sum(K.q)>0
    error('VSDP:VSDP2SDPAM','SOCP cone is not supported by SDPA');
elseif fields(2)
    K.l = sum(K.l);
else
    K.l = 0;
end
if fields(4)
    K.s = reshape(K.s(K.s>0),[],1);
else
    K.s = [];
end

% non-compact vectorized format
xdim = ~isempty(x);
zdim = ~isempty(z);
[c Imat] = vsmat(c,K,1,1);
if zdim
    c = [vsmat(z,K,1,1,Imat) c];
end
if xdim
    c = [vsmat(x,K,0.5,1,Imat) c];
end

% A = [x z c A]
if size(A,1)<size(A,2)
    A = -[c vsmat(A,K,1,1,Imat)'];
else
    A = -[c vsmat(A,K,1,1,Imat)];
end

% initialize outputs
mDIM = length(b);
if K.l
    bLOCKsTRUCT = [-K.l K.s'];
else
    bLOCKsTRUCT = K.s';
end
nBLOCK = length(bLOCKsTRUCT);
c = full(-b);  % full instead of sparse to catch a bug in sdpa
x0 = y;
F = cell(nBLOCK,mDIM+1);
X0 = cell(nBLOCK,1);
Y0 = cell(nBLOCK,1);

% free some memory
clear x z b y;


%% transform linear part
if K.l
    At = A(1:K.l,:);
    Y0{1} = At(:,1:xdim);
    X0{1} = At(:,xdim+1:xdim+zdim);
    for i = 1:mDIM+1
        F{1,i} = At(:,xdim+zdim+i);
    end
end


%% transform sdp part
blkj = K.l;
for j = (K.l>0)+1:nBLOCK
    nj = K.s(j-(K.l>0));
    At = reshape(A(blkj+1:blkj+nj*nj,:),nj,[]);
    blkj = blkj + nj*nj;
    Y0{j} = At(:,1:xdim*nj);
    X0{j} = At(:,xdim*nj+1:(xdim+zdim)*nj);
    blki = xdim + zdim;
    for i = 1:mDIM+1;
        F{j,i} = At(:,blki+1:blki+nj);
        blki = blki + nj;
    end
end

% empty unused cell arrays
if ~zdim
    X0 = [];
end
if ~xdim
    Y0 = [];
end

%__________________________End of VSDP2SDPAM___________________________