function obj = validate (obj)
% VALIDATE  Validate the VSDP 2012/18 input format.
%
%
%   See also vsdp.vsdp.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (1, 1);

% The cone structure constains all relevant metadata for the actual data,
% especially the vectorized cone sizes 'N' (uncondensed size) and 'n' condensed
% size.
[obj.K, obj.N, obj.n] = vsdp.validate_cone (obj.K);

% Prepare vector `b`.
if (~isfloat (obj.b) && ~isintval (obj.b))
  error ('VSDP:validate:badTypeB', ...
    'validate: Vector ''b'' must be numeric or an interval.');
end
obj.b = obj.b(:);

obj.m = length (obj.b);
if ((obj.m == obj.n) && (obj.n == obj.N))
  warning ('VSDP:validate:sizeAmbiguity', ...
    ['validate: The number of constraints ''m'' equals the number of ', ...
    'variables ''n''.  Make sure, that the constraint matrix ''At'' is ', ...
    'a transposed n-times-m constraint matrix!']);
end

if (~isfloat (obj.At) && ~isintval (obj.At))
  error ('VSDP:validate:badTypeAt', ...
    'validate: Constraint matrix ''At'' must be numeric or an interval.');
end

% Check if any dimension of matrix ''At'' matches the number of constraints,
% that is the length of vector ''b''.
if (~any (size (obj.At) == obj.m))
  error ('VSDP:validate:badConstraintDimesionAt', ...
    'validate: No dimension of ''At'' matches the length of vector ''b''.');
end
% Ensure transposed format for ''At'' (n x m).
if (size (obj.At, 2) ~= obj.m)
  obj.At = obj.At';
end
% Ensure valid cone dimension.
if ((obj.n < obj.N) && (size (obj.At, 1) == obj.N))
  obj.At = vsdp.svec (K, At, 1, false);  % No symmetry assumed => false.
elseif (size (obj.At, 1) ~= obj.n)
  error ('VSDP:validate:badConeDimesionAt', ...
    ['validate: ''size(At,1)'' must be %d ', ...
    '(or %d in condensed format) but is %d.'], obj.N, obj.n, size (obj.At, 1));
end

% prepare interval input for c
if (~isfloat (obj.c) && ~isa (obj.c, 'intval'))
  error ('VSDP:validate:badTypeC', ...
    'validate: Primal objective vector ''c'' must be numeric or an interval.');
end
obj.c = obj.c(:);
% Ensure valid cone dimension.
if ((obj.n < obj.N) && (length (obj.c) == obj.N))
  [c,Ivec] = vsvec(c,K,1,0,Ivec);
elseif (length (obj.c) ~= obj.n)
  error ('VSDP:validate:badLengthC', ...
    ['validate: Length of primal objective vector ''c'' must be %d ', ...
    '(or %d in condensed format) but is %d.'], obj.N, obj.n, length (obj.c));
end

if (~isempty (obj.x))
  if (~isfloat (obj.x))
    error ('VSDP:validate:badTypeX', ...
      'validate: Primal approximate solution vector ''x'' must be numeric.');
  end
  obj.x = obj.x(:);
  % Condense SDP cones if necessary.
  if ((obj.n < obj.N) && (length (obj.x) == obj.N))
    [x0,Ivec] = vsvec(x0,K,1,0,Ivec);  % Ivec can only be used with mu=1
    x0 = sscale(x0,K,2);
  elseif (length (obj.x) ~= obj.n)
    error ('VSDP:validate:badLengthX', ...
      ['validate: Length of primal solution vector ''x'' must be %d ', ...
      '(or %d in condensed format) but is %d.'], obj.N, obj.n, length (obj.x));
  end
end

% prepare y
if (~isempty (obj.y))
  if (~isfloat (obj.y))
    error ('VSDP:validate:badTypeY', ...
      'validate: Dual approximate solution vector ''y'' must be numeric.');
  end
  obj.y = obj.y(:);
  if (length (obj.y) ~= obj.m)
    error ('VSDP:validate:badLengthY', ...
      ['validate: Length of dual solution vector ''y'' must be %d ', ...
      'but is %d.'], obj.m, length (obj.y));
  end
end

if (~isempty (obj.z))
  if (~isfloat (obj.z))
    error ('VSDP:validate:badTypeZ', ...
      'validate: Primal approximate solution vector ''z'' must be numeric.');
  end
  obj.z = obj.z(:);
  % Condense SDP cones if necessary.
  if ((obj.n < obj.N) && (length (obj.z) == obj.N))
    [z0,Ivec] = vsvec(z0,K,1,0,Ivec);  % Ivec can only be used with mu=1
    z0 = sscale(z0,K,2);
  elseif (length (obj.z) ~= obj.n)
    error ('VSDP:validate:badLengthZ', ...
      ['validate: Length of primal solution vector ''z'' must be %d ', ...
      '(or %d in condensed format) but is %d.'], obj.N, obj.n, length (obj.z));
  end
end

end
