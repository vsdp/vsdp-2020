function [blk, At, ct, xt, zt] = vsdp2sdpt3(K,A,c,x0,z0,opts)
%% VSDP2SDPT3:  transforms problem data from VSDP to SDPT3 format
%    [blk,At,ct,xt,zt] = vsdp2sdpt3(K,A,c,x0,z0,opts)
%
% except for K all input parameter are optional
%
%% >> Input:
% K: a structure with following fields
%     - K.f stores the number of free variables
%     - K.l is the number of nonnegative components
%     - K.q lists the lengths of socp blocks
%     - K.s lists the dimensions of semidefinite blocks
% A:  a nA3 x M Matrix,
%     whereas nA = dimf+diml+dimq+dims
%     dimf: number of free variables: dimf = sum(K.f>0)
%     diml: number of nonnegative variables: diml = sum(K.l>0)
%     dimq: sum of all socp variables: dimq = sum_i(K.q(i))
%     dims3: sum of sdp variables: dims3 = sum_i(K.s(i)*(K.s(i)+1)/2)
% c: nC x 1 vector - primal objective function
% x0: a nC x 1 vector - approx. primal optimal solution
% z0: a nC x 1 vector - approx. dual optimal solution (slack vars)
% opts: structure for additional parameter settings:
%     regarded fields:
%       'SDPT3_VERSION'    Version number of SDPT3 solver.
%                              - default: 4.0
%       'MIN_SDPBLK_SIZE'  Minimum size of an sdp block. For efficiency 
%                          reasons SDPT3 allows to group smaller sdp
%                          blocks. Every sdp block with a dimension that is
%                          smaller MIN_SDPBLK_SIZE is considered as a small
%                          block that should be grouped.
%                               - default: 50
%
%% >> Output:
% blk: a cell array describing the block diagonal structure of problem data
% At: a cell array with At{1} = A(1:dimf,:), At{2} = A(dimf+1:dimf+diml,:)
%     and the same for socp cone, and At{i} = [svec(Ai1) ... svec(Aim)]
%     for positive semidefinite cone constraint matrices, where Aij is the
%     matrix of j-th block and i-th constraint
% ct: a cell array of matrices of dimensions given by K
% xt: a cell array of matrices of dimensions given by K
% zt: a cell array of matrices of dimensions given by K
% Is: index vector which describes the reorganization of the sdp blocks to
%     group small sdp blocks
%
%%
% Note that the right hand side of the linear constraints (b) and the
% dual optimal solution vector (y) have the same format in VSDP and SDPT3.
%
% This is a private conversion function. The dimension of the input
% parameter is not checked. This has to be done in the calling function.
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
% 30/07/12    M. Lange, rewrite for new parameter and lower memory usage
%
%% ToDo
% - code reduction
%


%% input parameter

% check input
if nargin<1 || ~isstruct(K)
    error('VSDP:VSDP2SDPT3','not enough input parameter');
elseif nargin<2
    A = [];  c = [];  x0 = [];  z0 = [];  opts = [];
elseif nargin<3
    c = [];  x0 = [];  z0 = [];  opts = [];
elseif nargin<4
    x0 = [];  z0 = [];  opts = [];
elseif nargin<5
    z0 = []; opts = [];
elseif nargin<6
    opts = [];
end

% global VSDP setting
global VSDP_OPTIONS;

% version number
version = 4.0;  % assumed version number
if isfield(opts,'SDPT3_VERSION')
    version = opts.SDPT3_VERSION;
elseif isfield(VSDP_OPTIONS,'SDPT3_VERSION')
    version = VSDP_OPTIONS.SDPT3_VERSION;
end

% max blocksize for 'small' sdp blocks
blksize = [];
if version<4.0  % ignore MIN_SDPBLK_SIZE if SOLVER_VERSION < 4.0
    blksize = 0;
elseif isfield(opts,'MIN_SDPBLK_SIZE')
    blksize = opts.MIN_SDPBLK_SIZE;
elseif isfield(VSDP_OPTIONS,'MIN_SDPBLK_SIZE')
    blksize = VSDP_OPTIONS.MIN_SDPBLK_SIZE;
end

% create sdp block group if  blk3s > blk3size = blksize*(blksize+1)/2
if isempty(blksize)
    blksize = 2500;
else
    blksize = blksize * (blksize+1) / 2;
end


%% preparation

% get problem data dimensions
fields = isfield(K,{'f','l','q','s'});
if fields(1)
    K.f = sum(K.f);
else
    K.f = 0;
end
if fields(2)
    K.l = sum(K.l);
else
    K.l = 0;
end
if fields(3)
    K.q = reshape(K.q(K.q>0),[],1);
else
    K.q = [];
end
if fields(4)
    K.s = reshape(K.s(K.s>0),[],1);
else
    K.s = 0;
end
dim = K.f + K.l + sum(K.q) + sum(K.s.*K.s);
dim3 = K.f + K.l + sum(K.q) + sum(K.s.*(K.s+1))/2;
n = (K.f>0) + (K.l>0) + length(K.q) + length(K.s);

% convert to appropriate format
if length(A)~=dim3
    A = vsvec(A,K,sqrt(2),0);
else
    A = sscale(A,K,sqrt(2));
end
Imat = [];
if length(c)~=dim
    [c Imat] = vsmat(c,K,2,0);
end
if all(length(x0)~=[0 dim])  % vx0 = svec(x0(:),K,2), !! mu=2 !!
    [x0 Imat] = vsmat(x0,K,1,0,Imat);
end
if all(length(z0)~=[0 dim])
    z0 = vsmat(z0,K,2,0,Imat);
end

% row- or column-wise
idim = 1 + (size(A,2)==dim3);

% error if free variables are used with version < 4.0
if version<4.0 && K.f>0
    error('VSDP:VSDP2SDPT3', ...
        'SDPT3 version < 4.0 does not support unconstrained variables');
end

% create cell arrays for output
blk = cell(n,2);
At = cell(n,1);
ct = cell(n,1);
xt = cell(n,1);
zt = cell(n,1);
        
% reducing the size of the input data after building the corresponding
% cells saves memory and increases the speed for the remaining data
% access
% the reduction, however, can have a significant effect on cpu-time
% the following is a tradeoff for cpu and memory efficiency
mem3tol = dim3 / min(n+1,3+sqrt(n));

% cell index
bli = 1;

% position indices
blk2e = 0;
blk3e = 0;


%% transform unconstrained variables
if K.f>0
    blk{bli,1} = 'u';
    blk{bli,2} = K.f;
    if ~isempty(A) && idim==1
        At{bli} = A(blk3e+1:blk3e+K.f,:);
        blk3e = blk3e+K.f;
    elseif ~isempty(A)
        At{bli} = A(:,blk3e+1:blk3e+K.f)';
        blk3e = blk3e+K.f;
    end
    if ~isempty(c)
        ct{bli} = c(blk2e+1:blk2e+K.f);
    end
    if ~isempty(x0)
        xt{bli} = x0(blk2e+1:blk2e+K.f);
    end
    if ~isempty(z0)
        zt{bli} = z0(blk2e+1:blk2e+K.f);
    end
    blk2e = blk2e + K.f;
    bli = bli + 1;
end


%% transform linear variables
if K.l>0 
    blk{bli,1} = 'l';
    blk{bli,2} = K.l;
    if ~isempty(A) && idim==1
        At{bli} = A(blk3e+1:blk3e+K.l,:);
        blk3e = blk3e+K.l;
    elseif ~isempty(A)
        At{bli} = A(:,blk3e+1:blk3e+K.l)';
        blk3e = blk3e+K.l;
    end
    if ~isempty(c)
        ct{bli} = c(blk2e+1:blk2e+K.l);
    end
    if ~isempty(x0)
        xt{bli} = x0(blk2e+1:blk2e+K.l);
    end
    if ~isempty(z0)
        zt{bli} = z0(blk2e+1:blk2e+K.l);
    end
    blk2e = blk2e + K.l;
    bli = bli + 1;
end


%% transform socp part
if ~isempty(K.q)
    dimq = sum(K.q);
    blk{bli,1} = 'q';
    blk{bli,2} = K.q(:)';
    if ~isempty(A) && idim
        At{bli} = A(blk3e+1:blk3e+dimq,:);
        blk3e = blk3e+dimq;
    elseif ~isempty(A)
        At{bli} = A(:,blk3e+1:blk3e+dimq)';
        blk3e = blk3e+dimq;
    end
    if ~isempty(c)
        ct{bli} = c(blk2e+1:blk2e+dimq);
    end
    if ~isempty(x0)
        xt{bli} = x0(blk2e+1:blk2e+dimq);
    end
    if ~isempty(z0)
        zt{bli} = z0(blk2e+1:blk2e+dimq);
    end
    blk2e = blk2e + dimq;
    bli = bli + 1;
end


%% transform sdp part, consider small blocks
if ~isempty(K.s)
    
    % save index positions before processing sdp cells
    bls = bli;
    blk2s = blk2e;
    
    % starting index for current group
    kstart = 1;
    
    %% create blk
    blke = 0;
    for k = 1:length(K.s)
        blke = blke + K.s(k)*(K.s(k)+1)/2;
        
        % group great enough or end of sdp data
        if blke>blksize || k==length(K.s)
            blk{bli,1} = 's';
            blk{bli,2} = reshape(K.s(kstart:k),1,[]);
            blke = 0;
            kstart = k + 1;
            bli = bli + 1;
        end
    end
    
    
    %% create At cells for sdp part
    if ~isempty(A)
        % keep sdp part only
        if blk3e>0 && idim
            A(1:blk3e,:) = [];
            blk3e = 0;
        elseif blk3e>0
            A(:,1:blk3e) = [];
            blk3e = 0;
        end
        % create cell blocks
        if idim==1
            for k = bls:bli-1
                nk3 = sum(blk{k,2}.*(blk{k,2}+1))/2;
                At{k} = A(blk3e+1:blk3e+nk3,:);
                blk3e = blk3e + nk3;
                if blk3e>mem3tol
                    A = A(blk3e+1:end,:);
                    blk3e = 0;
                end
            end
        else
        	for k = bls:bli-1
                nk3 = sum(blk{k,2}.*(blk{k,2}+1))/2;
                At{k} = A(:,blk3e+1:blk3e+nk3)';
                blk3e = blk3e + nk3;
                if blk3e>mem3tol
                    A = A(:,blk3e+1:end);
                    blk3e = 0;
                end
            end
        end
    end  % transformation of A
    
    
    %% block - cell transformation for sdp part of c
    if ~isempty(c)
        blk2e = blk2s;  % reset index position
        for k = bls:bli-1
            nk = blk{k,2};
            lk = length(nk);
            if lk==1  % single sdp block
                ct{k} = reshape(c(blk2e+1:blk2e+nk*nk),nk,nk);
                ct{k} = 0.5 * (ct{k} + ct{k}');
                blk2e = blk2e + nk*nk;
            else  % block diagonal form of sdp group
                ct{k} = cell(1,lk);
                for i = 1:lk
                    nki = nk(i);
                    ct{k}{i} = reshape(c(blk2e+1:blk2e+nki*nki),nki,nki);
                    ct{k}{i} = 0.5 * sparse(ct{k}{i} + ct{k}{i}');
                    blk2e = blk2e + nki*nki;
                end
                nki = nk(end);
                ct{k} = blkdiag(ct{k}{:});
            end
        end
    end  % transformation of c
    
    
    %% block - cell transformation for sdp part of x0
    if ~isempty(x0)
        blk2e = blk2s;  % reset index position
        for k = bls:bli-1
            nk = blk{k,2};
            lk = length(nk);
            if lk==1  % single sdp block
                xt{k} = reshape(x0(blk2e+1:blk2e+nk*nk),nk,nk);
                xt{k} = 0.5 * (xt{k} + xt{k}');
                blk2e = blk2e + nk*nk;
            else  % block diagonal form of sdp group
                xt{k} = cell(1,lk);
                for i = 1:lk
                    nki = nk(i);
                    xt{k}{i} = reshape(x0(blk2e+1:blk2e+nki*nki),nki,nki);
                    xt{k}{i} = 0.5 * sparse(xt{k}{i} + xt{k}{i}');
                    blk2e = blk2e + nki*nki;
                end
                nki = nk(end);
                xt{k} = blkdiag(xt{k}{:});
                blk2e = blk2e + nki*nki;
            end
        end
    end  % transformation of x0

    
    %% block - cell transformation for sdp part of z0
    if ~isempty(z0)
        blk2e = blk2s;  % reset index position
        for k = bls:bli-1
            nk = blk{k,2};
            lk = length(nk);
            if lk==1  % single sdp block
                zt{k} = reshape(z0(blk2e+1:blk2e+nk*nk),nk,nk);
                zt{k} = 0.5 * (zt{k} + zt{k}');
                blk2e = blk2e + nk*nk;
            else  % block diagonal form of sdp group
                zt{k} = cell(1,lk);
                for i = 1:lk
                    nki = nk(i);
                    zt{k}{i} = reshape(z0(blk2e+1:blk2e+nki*nki),nki,nki);
                    zt{k}{i} = 0.5 * sparse(zt{k}{i} + zt{k}{i}');
                    blk2e = blk2e + nki*nki;
                end
                nki = nk(end);
                zt{k} = blkdiag(zt{k}{:});
                blk2e = blk2e + nki*nki;
            end
        end
    end  % transformation of z0
    
end  % transform sdp part
    

%% remove unused cells
blk(bli:end,:) = [];
At(bli:end) = [];
ct(bli:end) = [];
xt(bli:end) = [];
zt(bli:end) = [];

if isempty(At{1})
    At = [];
end
if isempty(c)
    ct = [];
end
if isempty(x0)
    xt = [];
end
if isempty(z0)
    zt = [];
end

%_____________________________End VSDP2SDPT3___________________________