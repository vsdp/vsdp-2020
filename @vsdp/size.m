function varargout = size (obj, dim)
% SIZE  Size of VSDP instance.
if (nargin == 1)
  if (nargout < 2)
    varargout = {[obj.m, obj.n]};
  elseif (nargout == 2)
    varargout = {obj.m, obj.n};
  else
    varargout(1:2) = {obj.m, obj.n};
    varargout(3:nargout) = {1};
  end
else
  if (dim == 1)
    varargout = {obj.m};
  elseif (dim == 2)
    varargout = {obj.n};
  else
    varargout = {1};
  end
end
end