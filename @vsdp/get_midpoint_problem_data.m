function [At, b, c] = get_midpoint_problem_data (obj, sol_type)
% GET_MIDPOINT_PROBLEM_DATA  Get midpoint data of the conic problem.
%
%   See also vsdp.
%

% Copyright 2004-2020 Christian Jansson (jansson@tuhh.de)

[At, b, c] = deal (mid (obj.At), mid (obj.b), mid (obj.c));
if ((nargin > 1) && ~strcmp (sol_type, "Approximate"))
  % Return the perturbed midpoint problem just in case no approximate
  % solution is requested.
  if (~isempty (obj.perturbation.b))
    b = b - obj.perturbation.b;
  end
  if (~isempty (obj.perturbation.c))
    c = c - obj.perturbation.c;
  end
end
end
