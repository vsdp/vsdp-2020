function [K,At,c,x0,z0] = vsdp12vsdp(blk,A,C,X0,Z0)
% VSDP12VSDP  Convert problem data from VSDP 2006 to VSDP 2012 format.
%
%   [K,A,c,x0,z0] = VSDP12VSDP(blk,A,C,X0,Z0)
%      The VSDP 2006 block-diagonal structure format is:
%
%         min  sum(j=1:n| <  C{j}, X{j}>)
%         s.t. sum(j=1:n| <A{i,j}, X{j}>) = b(i)  for i = 1:m
%              X{j} must be positive semidefinite for j = 1:n
%
%      The problem data of the block-diagonal structure:
%
%         'blk'  cell(n,2)
%         'A'    cell(m,n)
%         'C'    cell(n,1)
%         'b'  double(m,1)
%
%      The j-th block C{j} and the blocks A{i,j}, for i = 1:m, are real
%      symmetric matrices of common size s_j, and blk(j,:) = {'s'; s_j}.
%
%      The blocks C{j} and A{i,j} must be stored as individual matrices in
%      dense or sparse format.
%
%      The optional initial guess format is:
%
%         'X0'   cell(n,1)
%         'Z0'   cell(n,1)
%
%   See also import_vsdp.

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

% create K structure
if nargin<1 || ~iscell(blk) || isempty(blk{1,2})
  error('VSDP:VSDP12VSDP','cannot read "blk"');
else
  K.s = cat(1,blk{:,2});
end

% initial output
At = [];
c = [];
x0 = [];
z0 = [];

% index vector to speed up svec
Ivec = [];

% coefficient matrix
if nargin>1 && iscell(A) && ~isempty(A{1})
  [mblk,nblk] = size(A);
  At = cell(nblk,1);
  for i = nblk:-1:1
    % not memory but runtime efficient
    At{i} = reshape(cat(2,A{:,i}),[],mblk);
    A(:,i) = [];
  end
  [At,Ivec] = vsvec(cat(1,At{:}),K,1,1);
end

% primal objective vector
if nargin>2 && iscell(C) && ~isempty(C{1})
  C = cellfun(@(x) x(:),C,'UniformOutput',false);
  [c,Ivec] = vsvec(cat(1,C{:}),K,1,1,Ivec);
  clear C;
end

% initial primal solution vector
if nargin>3 && iscell(X0) && ~isempty(X0{1})
  X0 = cellfun(@(x) x(:),X0,'UniformOutput',false);
  [x0,Ivec] = vsvec(cat(1,X0{:}),K,2,1,Ivec);  % mu=2
  clear Xt;
end

% initial dual solution vector (slack variables)
if nargin>4 && iscell(Z0) && ~isempty(Z0{1})
  Z0 = cellfun(@(x) x(:),Z0,'UniformOutput',false);
  z0 = vsvec(cat(1,Z0{:}),K,1,1,Ivec);
  clear Zt;
end

end
