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
% Ensure vectorzied condensed cone dimension.
try
  obj.At = vsdp.svec (obj, obj.At, 1, 'unsymmetric');
catch ME
  error ('VSDP:validate:badAt', ...
    'validate: Cannot vectorize ''At''.');
end

% Prepare 'c'.
if (~isfloat (obj.c) && ~isa (obj.c, 'intval'))
  error ('VSDP:validate:badTypeC', ...
    'validate: Primal objective vector ''c'' must be numeric or an interval.');
end
% Ensure vectorzied condensed cone dimension.
try
  obj.c = vsdp.svec (obj, obj.c(:), 1, 'unsymmetric');
catch ME
  error ('VSDP:validate:badC', ...
    'validate: Cannot vectorize ''c''.');
end

end
