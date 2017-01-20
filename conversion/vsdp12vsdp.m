function [K,A,c,x,z] = vsdp12vsdp(blk,At,C,Xt,Zt)
%% VSDP2SEDUMI:  transforms problem data from VSDP 0.1 to new VSDP format
%    [K,A,c,x,z] = vsdp12vsdp(blk,At,C,Xt,Zt)
%
%% >> Input:
% blk: a cell array describing the block diagonal structure of problem data
% At: a cell array with At{i,1:m} = [svec(Ai1) ... svec(Aim)]
%     for positive semidefinite cone constraint matrices, where Aij is the
%     matrix of j-th block and i-th constraint
% C: a cell array of matrices of dimensions given by blk
% Xt: a cell array of matrices, initial primal solution
% Zt: a cell array of matrices, initial dual solution
%
%% >> Output:
% A: a dims3 x M Matrix,
%     whereas dims3: sum of all sdp variables: dims3 = sum_i(K.s(i)*(K.s(i)+1)/2)
% c: dims3 x 1 vector for primal objective function
% x: dims3 x 1 vector, initial primal solution
% z: dims3 x 1 vector, initial dual solution (slack vars)
% K: a structure with following fields
%     - K.s lists the dimensions of semidefinite blocks
%

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

% create K structure
if nargin<1 || ~iscell(blk) || isempty(blk{1,2})
  error('VSDP:VSDP12VSDP','cannot read "blk"');
else
  K.s = cat(1,blk{:,2});
end

% initial output
A = [];
c = [];
x = [];
z = [];

% index vector to speed up svec
Ivec = [];

% coefficient matrix
if nargin>1 && iscell(At) && ~isempty(At{1})
  [mblk,nblk] = size(At);
  A = cell(nblk,1);
  for i = nblk:-1:1
    % not memory but runtime efficient
    A{i} = reshape(cat(2,At{:,i}),[],mblk);
    At(:,i) = [];
  end
  [A,Ivec] = vsvec(cat(1,A{:}),K,1,1);
end

% primal objective vector
if nargin>2 && iscell(C) && ~isempty(C{1})
  C = cellfun(@(x) x(:),C,'UniformOutput',false);
  [c,Ivec] = vsvec(cat(1,C{:}),K,1,1,Ivec);
  clear C;
end

% initial primal solution vector
if nargin>3 && iscell(Xt) && ~isempty(Xt{1})
  Xt = cellfun(@(x) x(:),Xt,'UniformOutput',false);
  [x,Ivec] = vsvec(cat(1,Xt{:}),K,2,1,Ivec);  % mu=2
  clear Xt;
end

% initial dual solution vector (slack variables)
if nargin>4 && iscell(Zt) && ~isempty(Zt{1})
  Zt = cellfun(@(x) x(:),Zt,'UniformOutput',false);
  z = vsvec(cat(1,Zt{:}),K,1,1,Ivec);
  clear Zt;
end

end
