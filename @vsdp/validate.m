function obj = validate (At, b, c, K, x0, y0, z0)
% VALIDATE  Validate the VSDP 2012/18 input format.
%
%
%

narginchk (4, 7);

% The cone structure constains all relevant metadata for the actual data,
% especially the vectorized cone sizes 'N' (uncondensed size) and 'n' condensed
% size).
[obj.K, obj.N, obj.n] = validate_cone (K);

% Prepare vector `b`.
if (~isfloat (b) && ~isintval (b))
  error ('VSDP:validate:wrongTypeB', ...
    'validate: Vector `b` has wrong data type.');
end
b = b(:);

obj.m = length(b);
if ((obj.m == obj.n) && (obj.n == obj.N))
  warning ('VSDP:validate:sizeAmbiguity', ...
    ['validate: The number of constraints ''m'' equals the number of ', ...
    'variables ''n''.  Make sure, that the constraint matrix ''At'' is ', ...
    'a transposed n-times-m constraint matrix!']);
end

% Check type of matrix `At`.
if (~isfloat (At) && ~isintval (At))
  error ('VSDP:validate:wrongTypeAt', ...
    'validate: Matrix ''At'' has wrong data type.');
end

% Check if any dimension of matrix ''At'' matches the number of constraints,
% that is the length of vector ''b''.
if (~any(size(At) == obj.m))
  error ('VSDP:validate:badConstraintDimesionAt', ...
    'validate: No dimension of ''At'' matches the length of vector ''b''.');
end
% Ensure transposed format for ''At'' (n x m).
if (size (At, 2) ~= obj.m)
  At = At';
end
% Ensure valid cone dimension.
if (size (At, 1) == obj.N)
  At = svec (K, At, 1, false);  % No symmetry assumed => false.
elseif (size (At, 1) ~= n)
  error ('VSDP:validate:badConeDimesionAt', ...
    'validate: `Bad cone dimension `n of a `At` matches the length of vector `b`.');
end

% prepare interval input for c
if (~isfloat(c) && ~isa(c, 'intval'))
  error('VSDP:IMPORT_VSDP','cannot import primal objective "c"');
end
c = c(:);
% compact vectorized format
if (length(c) ~= n)
  [c,Ivec] = vsvec(c,K,1,0,Ivec);
end
obj.c = c;

% prepare x
if (~isempty(x0))
  if (~isfloat(x0))
    error('VSDP:IMPORT_VSDP','primal solution vector "x" has to be numeric');
  end
  x0 = x0(:);
  % compact vectorized format, mu=2
  if (length(x0) ~= n)
    [x0,Ivec] = vsvec(x0,K,1,0,Ivec);  % Ivec can only be used with mu=1
    x0 = sscale(x0,K,2);
  end
end
obj.x = x0;

% prepare y
if (~isempty(y0))
  y0 = y0(:);
  if (~isfloat(y0) || (length(y0) ~= m))
    error('VSDP:IMPORT_VSDP','cannot import dual solution vector "y"');
  end
end
obj.y = y0;

% prepare z
if (~isempty(z0))
  if (~isfloat(z0))
    error('VSDP:IMPORT_VSDP','cannot import dual solution "z"');
  end
  z0 = z0(:);
  % compact vectorized format
  if (size(z0,1) ~= n)
    z0 = vsvec(z0,K,1,0,Ivec);
  end
end
obj.z = z0;
end

function [K, N, n] = validate_cone (K_in)
K.f = 0;
K.l = 0;
K.q = [];
K.s = [];
K.f_idx = [];
K.l_idx = [];
K.q_idx = [];
K.s_idx = [];

if (isfield (K_in, 'f'))
  K.f = sum(K_in.f);
  K.f_idx = 1:K.f;
end
if (isfield (K_in, 'l'))
  K.l = sum(K_in.l);
  K.l_idx = K.f + (1:K.l);
end
if (isfield (K_in, 'q'))
  K.q = K_in.q(K_in.q > 0);
  % Compute index starts and ends with offset.
  K.q_idx = K.f + K.l + [cumsum([1, K.q(1:end-1)]); cumsum(K.q)];
  % Create index ranges in cells.
  K.q_idx = cellfun (@(x) x(1):x(2), num2cell (K.q_idx, 1), ...
    'UniformOutput', false);
end
if (isfield (K_in, 's'))
  K.s = K_in.s(K_in.s > 0);
  % Compute index starts and ends with offset.
  K.s_idx = K.s .* (K.s + 1) / 2;
  K.s_idx = K.f + K.l + sum(K.q) + ...
    [cumsum([1, K.s_idx(1:end-1)]); cumsum(K.s_idx)];
  % Create index ranges in cells.
  K.s_idx = cellfun (@(x) x(1):x(2), num2cell (K.s_idx, 1), ...
    'UniformOutput', false);
end

% Determine (un-)condensed cone dimension.
N = K.f + K.l + sum(K.q) + sum (K.s .* K.s);
n = K.f + K.l + sum(K.q) + (sum (K.s .* (K.s + 1)) / 2);
end