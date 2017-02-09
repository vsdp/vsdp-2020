function [A,Arad,b,brad,c,crad,K,x,y,z,use_fmt] = import_vsdp(A,b,c,K,x,y,z)
% IMPORT_VSDP  Imports different supported formats.
%
%   [A,Arad,b,brad,c,crad,K,x,y,z,IF] = IMPORT_VSDP(A,b,c,K,x,y,z)
%   IMPORT_VSDP(...,x,y,z)
%
%   Input: problem data in SEDUMI, old VSDP, or new VSDP internal format.
%
% Output:
% A/Arad: an nA3 x M,
%     whereas nA3 = dimf+diml+dimq+dims3
%     dimf: number of free variables: dimf = sum(K.f>0)
%     diml: number of nonnegative variables: diml = sum(K.l>0)
%     dimq: sum of all socp variables: dimq = sum_i(K.q(i))
%     dims3: sum of all sdp variables: dims3 = K.s'*K.s/2
% b/brad: an M x 1 vector
% c/crad: an nA3 x 1 vector,
% K: a structure with following fields
%     - K.f stores the number of free variables
%     - K.l is the number of nonnegative components
%     - K.q lists the lengths of socp blocks
%     - K.s lists the dimensions of semidefinite blocks
% x: an nA3 x 1 vector - a primal feasible (eps-optimal) solution, mu=2
% y: an M x 1 vector - a dual feasible (eps-optimal) solution
% z: an nA3 x 1 vector - a dual feasible (eps-optimal) solution (slack vars)
%

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

% check input
narginchk(4,7);
if isempty(A) || isempty(b) || isempty(c) || isempty(K)
  error('VSDP:IMPORT_VSDP','not enough input parameter');
end
if (nargin < 5)
  x = [];
end
if (nargin < 6)
  y = [];
end
if (nargin < 7)
  z = [];
end

% default input format
use_fmt = [];  % if empty: default internal format

% import old vsdp format
if iscell(A)  % (A,b,c,K,x,y,z) -> (blk,A,C,b,Xt,yt,Zt)
  use_fmt = 'VSDP01';
  tmp = K;  % b
  [K,A,c,x,z] = vsdp12vsdp(A,b,c,x,z);  % (blk,A,C,Xt,Zt)
  b = tmp;
end

% index vector to speed-up svec
Ivec = [];


% prepare cone structure
f = 0;  l = 0;  q = [];  s = [];
fields = isfield(K,{'f','l','q','s'});
if fields(1)
  f = sum(K.f);
end
if fields(2)
  l = sum(K.l);
end
if fields(3)
  q = reshape(K.q(K.q>0),[],1);
end
if fields(4)
  s = reshape(K.s(K.s>0),[],1);
end
K = struct('f',f,'l',l,'q',q,'s',s);
% appropriate dimensions
dim3 = f + l + sum(q) + sum(s.*(s+1))/2;


% prepare interval input for b
if isnumeric(b)
  b = b(:);  brad = sparse(size(b,1),1);
elseif isa(b,'intval') || all(isfield(b,{'mid','rad'}))
  brad = b.rad(:);  b = b.mid(:);
else
  error('VSDP:IMPORT_VSDP','cannot import dual objective "b"');
end
% check size of b
m = size(b,1);
if size(brad,1)~=m
  error('VSDP:IMPORT_VSDP','rad(b) has not the same dimension as mid(b)');
end


% prepare interval input for A
if isnumeric(A)
  Arad = 0;
elseif isa(A,'intval') || all(isfield(A,{'mid','rad'}))
  Arad = A.rad;  A = A.mid;
else
  error('VSDP:IMPORT_VSDP','cannot import coefficient matrix "A"');
end
% compact vectorized format
if (any(size(A) - [dim3, m]) && any(size(A) - [m, dim3]))
  use_fmt = 'SEDUMI';
  [A,Ivec] = vsvec(A,K,1,0,Ivec);
end
if ~any(find(Arad,1))
  Arad = sparse(dim3,m);
elseif (any(size(Arad) - [dim3, m]) && any(size(Arad) - [m, dim3]))
  Arad = vsvec(Arad,K,1,0,Ivec);
end
% transposed format
if (size(A,1) == m)
  A = A';
end
if (size(Arad,1) == m)
  Arad = Arad';
end
% check size of A
if any(size(A)~=[dim3 m] | size(Arad)~=[dim3 m])
  error('VSDP:IMPORT_VSDP','wrong dimension of coefficient matrix "A"');
end


% prepare interval input for c
if isnumeric(c)
  c = c(:);  crad = 0;
elseif isa(c,'intval') || all(isfield(c,{'mid','rad'}))
  crad = c.rad(:);  c = c.mid(:);
else
  error('VSDP:IMPORT_VSDP','cannot import primal objective "c"');
end
% compact vectorized format
if size(c,1)~=dim3
  [c,Ivec] = vsvec(c,K,1,0,Ivec);
end
if size(crad,1)~=dim3 && any(find(crad,1))
  crad = vsvec(crad,K,1,0,Ivec);
elseif size(crad,1)~=dim3
  crad = sparse(dim3,1);
end


% prepare x
if ~isempty(x)
  if ~isnumeric(x)
    error('VSDP:IMPORT_VSDP','primal solution vector "x" has to be numeric');
  end
  x = x(:);
  % compact vectorized format, mu=2
  if size(x,1)~=dim3
    [x,Ivec] = vsvec(x,K,1,0,Ivec);  % Ivec can only be used with mu=1
    x = sscale(x,K,2);
  end
end


% prepare y
if ~isnumeric(y) || all(length(y)~=[0 m])
  error('VSDP:IMPORT_VSDP','cannot import dual solution vector "y"');
end
y = y(:);


% prepare z
if ~isempty(z)
  if ~isnumeric(z)
    error('VSDP:IMPORT_VSDP','cannot import dual solution "z"');
  end
  z = z(:);
  % compact vectorized format
  if size(z,1)~=dim3
    z = vsvec(z,K,1,0,Ivec);
  end
end

end
