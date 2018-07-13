function [K, N, n] = validate_cone (K_in)
% VALIDATE_CONE  Ensure a given VSDP cone strucuture to be valid.
%
%

narginchk(1, 1);

if (~isstruct (K_in))
  error ('VSDP:validate_cone:noStruct', ...
    'validate_cone: Input type must be a struct.');
end

K.f = 0;
K.l = 0;
K.q = [];
K.s = [];
K.idx.f = [];
K.idx.l = [];
K.idx.q = [];
K.idx.s = [];

if (isfield (K_in, 'f'))
  K.f = sum (K_in.f(K_in.f > 0));
  if (K.f > 0)
    K.idx.f = [1, K.f];
  end
end
if (isfield (K_in, 'l'))
  K.l = sum (K_in.l(K_in.l > 0));
  if (K.l > 0)
    K.idx.l = K.f + [1, K.l];
  end
end
if (isfield (K_in, 'q'))
  K.q = K_in.q(K_in.q > 0);
  K.q = K.q(:);
  if (~isempty (K.q))
    K.idx.q = K.f + K.l + [cumsum([1; K.q(1:end-1)]); cumsum(K.q)];
  end
end
if (isfield (K_in, 's'))
  K.s = K_in.s(K_in.s > 0);
  K.s = K.s(:);
  if (~isempty (K.s))
    K.idx.s = K.s .* (K.s + 1) / 2;
    K.idx.s = K.f + K.l + sum(K.q) + ...
      [cumsum([1; K.idx.s(1:end-1)]), cumsum(K.idx.s)];
  end
end

if (any ([K.f; K.l; K.q; K.s] ~= round([K.f; K.l; K.q; K.s])))
  error ('VSDP:validate_cone:needPositiveIntegers', ...
    'validate_cone: All cone dimensions must be positive intergers.');
end

% Determine (un-)condensed cone dimension.
N = K.f + K.l + sum(K.q) + sum (K.s .* K.s);
n = K.f + K.l + sum(K.q) + (sum (K.s .* (K.s + 1)) / 2);

end
