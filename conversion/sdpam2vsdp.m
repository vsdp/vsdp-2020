function [A,b,c,K,x,y,z] = sdpam2vsdp(bLOCKsTRUCT,c,F,x0,X0,Y0)
%% VSDP2SDPT3:  transforms problem data from SDPAM to VSDP3 format
%    [A,b,c,K,x,y,z] = sdpam2vsdp(bLOCKsTRUCT,c,F,x0,X0,Y0)
%
%% >> Input:
% mDIM: scalar - number of primal variables
% nBLOCK: scalar - number of blocks of F
% bLOCKsTRUCT: scalar - represents the block structure of F
% c: vector - coefficient vector
% F: cell array of coefficient matrices
% x0: vector for initial solution
% X0,Y0: cell arrays of initial points
%
%% >> Output:
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
% 22/08/12    M. Lange, written
%


%% prepare data
if nargin<1 || isempty(bLOCKsTRUCT)
    error('VSDP:SDPAM2VSDP','"bLOCKsTRUCT" has to be set');
end

% initial output
A = [];  b = [];  x = [];  y = [];  z = [];

% indices of diagonal blocks - linear part
indl = find(bLOCKsTRUCT<2);

% indices of semidefinite blocks
inds = find(bLOCKsTRUCT>1);

% index vector to speed up svec
Ivec = [];


%% create K
K = struct('l',sum(abs(bLOCKsTRUCT(indl))),'s',reshape(bLOCKsTRUCT(inds),[],1));


%% convert c to b
if nargin>1
    b = -c;
end


%% convert F to c & A
if nargin>2 && ~isempty(F)
    % vectorize linear part, sort blocks
    F = [ cellfun(@(x) x(linspace(1,numel(x),length(x))),F(indl,:),...
        'UniformOutput',false); F(inds,:) ];
    
    % concat
    mdim = size(F,2);
    for i = 1:size(F,1)
        F{i,1} = reshape(cat(2,F{i,:}),[],mdim);
    end
    F = -cat(1,F{:,1});  % [c A] in SeDuMi format
    
    % compact VSDP format
    [F Ivec] = vsvec(F,K,1,1);
    
    % split into c and A
    c = F(:,1);
    A = F(:,2:end);
    clear F;
else
    c = [];
end


%% write y
if nargin>4
    y = x0;
end


%% convert X0 to z
if nargin>3 && ~isempty(X0) && nargout>6
    % vectorize + sort blocks
    X0 = [ cellfun(@(x) reshape(x(linspace(1,numel(x),length(x))),[],1),...
             X0(indl),'UniformOutput',false) ; ...
           cellfun(@(x) x(:),X0(inds),'UniformOutput',false) ];
    
    % concat + svec
    [z Ivec] = vsvec(cat(1,X0{:}),K,1,1,Ivec);
end


%% convert Y0 to x
if nargin>5 && ~isempty(Y0) && nargout>4
    % vectorize + sort blocks
    Y0 = [ cellfun(@(x) reshape(x(linspace(1,numel(x),length(x))),[],1),...
             Y0(indl),'UniformOutput',false) ; ...
           cellfun(@(x) x(:),Y0(inds),'UniformOutput',false) ];
    
    % concat + svec
    x = vsvec(cat(1,Y0{:}),K,2,1,Ivec);  % mu=2
end


%__________________________End of SDPAM2VSDP___________________________