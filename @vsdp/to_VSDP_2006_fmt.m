function [blk, A, C, b, X0, y0, Z0] = to_vsdp_2006_fmt (obj)
% TO_VSDP_2006_FMT  Export conic problem data to VSDP 2006 format.
%
%   [blk, A, C, b, X0, y0, Z0] = obj.toVSDP2006Fmt();
%
%      `obj` is a vsdp class object.
%
%   See also from_vsdp_2006_fmt.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% check supported cones
if (any (obj.K.f > 0) || any (obj.K.l > 0) || any (obj.K.q > 0))
  error ('VSDP:TOVSDP2006FMT', ...
    'toVSDP2006Fmt: The VSDP 2006 format supports only SDP cones.');
end
% convert x
if ~isempty(x)
  x = vsvec(x,K,2,0);  % to verify that internal VSDP format is given
  blke = length(x);
  for j = length(K.s):-1:1
    nj = K.s(j);
    blks = blke - nj*(nj+1)/2 + 1;
    Xt{j}(triu(true(nj))) = 0.5 * x(blks:blke);  % mu==2
    Xt{j} = reshape(Xt{j},nj,nj);
    Xt{j} = Xt{j} + Xt{j}';
    blke = blks - 1;
  end
  X0 = Xt;
end

if ((nargout > 6) && ~isempty(z))
  blke = length(z);
  for j = length(K.s):-1:1
    nj = K.s(j);
    blks = blke - nj*(nj+1)/2 + 1;
    Zt{j}(triu(true(nj))) = z(blks:blke);
    Xt{j} = reshape(Xt{j},nj,nj);
    Zt{j} = Xt{j} + Xt{j}' - diag(sparse(diag(Xt{j})));
    blke = blks - 1;
  end
  Z0 = Zt;
else
  Z0 = [];
end

end
