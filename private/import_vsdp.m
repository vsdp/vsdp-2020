function [A,b,c,K,x,y,z,imported_fmt] = import_vsdp(A,b,c,K,x,y,z)
% IMPORT_VSDP  Import and verify problem data.
%
%   [A,b,c,K,x,y,z,imported_fmt] = IMPORT_VSDP(A,b,c,K)
%   IMPORT_VSDP(...,x,y,z)
%
%   Input: problem data in SEDUMI, old VSDP, or new VSDP internal format.
%
%   The function returns:
%
%      'A,b,c,K,x,y,z'  Validated arguments in VSDP's block-diagonal format,
%                       that is explained in 'mysdps.m'.
%
%      'imported_fmt'   The recognized imported format.  One of
%                              []: empty means default format.
%                        'VSDP01': vsdp-2006 format.
%                        'SEDUMI': SeDuMi format.
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
imported_fmt = [];  % if empty: default internal format

% import old vsdp format
if iscell(A)  % (A,b,c,K,x,y,z) -> (blk,A,C,b,Xt,yt,Zt)
  imported_fmt = 'VSDP01';
  tmp = K;  % b
  [K,A,c,x,z] = vsdp12vsdp(A,b,c,x,z);  % (blk,A,C,Xt,Zt)
  b = tmp;
end

% index vector to speed-up svec
Ivec = [];


% prepare cone structure
f = 0;
l = 0;
q = [];
s = [];
fields = isfield(K, {'f','l','q','s'});
if fields(1)
  f = sum(K.f);
end
if fields(2)
  l = sum(K.l);
end
if fields(3)
  q = K.q(K.q > 0);
end
if fields(4)
  s = K.s(K.s > 0);
end
K = struct('f', f, 'l', l, 'q', q(:), 's', s(:));
% appropriate dimensions
dim3 = f + l + sum(q) + (sum(s .* (s + 1)) / 2);


% prepare interval input for b
if (~isnumeric(b) && ~isa(b, 'intval'))
  error('VSDP:IMPORT_VSDP','cannot import dual objective "b"');
end
b = b(:);
m = length(b);

% prepare interval input for A
if (~isnumeric(A) && ~isa(A, 'intval'))
  error('VSDP:IMPORT_VSDP','cannot import coefficient matrix "A"');
end
% compact vectorized format
if (any(size(A) - [dim3, m]) && any(size(A) - [m, dim3]))
  imported_fmt = 'SEDUMI';
  [A,Ivec] = vsvec(A,K,1,0,Ivec);
end
% transposed format
if (size(A,1) == m)
  A = A';
end
% check size of A
if any(size(A) ~= [dim3, m])
  error('VSDP:IMPORT_VSDP','wrong dimension of coefficient matrix "A"');
end

% prepare interval input for c
if (~isnumeric(c) && ~isa(c, 'intval'))
  error('VSDP:IMPORT_VSDP','cannot import primal objective "c"');
end
c = c(:);
% compact vectorized format
if (length(c) ~= dim3)
  [c,Ivec] = vsvec(c,K,1,0,Ivec);
end

% prepare x
if (~isempty(x))
  if (~isnumeric(x))
    error('VSDP:IMPORT_VSDP','primal solution vector "x" has to be numeric');
  end
  x = x(:);
  % compact vectorized format, mu=2
  if (length(x) ~= dim3)
    [x,Ivec] = vsvec(x,K,1,0,Ivec);  % Ivec can only be used with mu=1
    x = sscale(x,K,2);
  end
end

% prepare y
if (~isempty(y))
  y = y(:);
  if (~isnumeric(y) || (length(y) ~= m))
    error('VSDP:IMPORT_VSDP','cannot import dual solution vector "y"');
  end
end


% prepare z
if (~isempty(z))
  if (~isnumeric(z))
    error('VSDP:IMPORT_VSDP','cannot import dual solution "z"');
  end
  z = z(:);
  % compact vectorized format
  if (size(z,1) ~= dim3)
    z = vsvec(z,K,1,0,Ivec);
  end
end

end
