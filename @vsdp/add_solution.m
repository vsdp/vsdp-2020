function obj = add_solution (obj, varargin)
% ADD_SOLUTION  Add initial guess or approximate solution.
%
%   See also vsdp.vsdp.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

if (nargin == 2)
  if (~isa (varargin{1}, 'vsdp_approx_solution'))
    error ('VSDP:add_solution:badType', ...
      ['add_solution: A single argument must be a ''vsdp_approx_solution'' ', ...
      'object.']);
  else
    sol = varargin{1};
  end
else
  sol = vsdp_approx_solution (varargin{:});
end

[m, n] = sol.size ();
if ((m > 0) && (m ~= obj.m))
  error ('VSDP:add_solution:dimMissmatch', ...
    ['add_solution: Dimension missmatch.  The number of problem ', ...
    'constraints is m = %d but the solution has %d.'], obj.m, m);
end
if ((n > 0) && (n ~= obj.n))
  error ('VSDP:add_solution:dimMissmatch', ...
    ['add_solution: Dimension missmatch.  The number of problem ', ...
    'variables is n = %d but the solution has %d.\n\n', ...
    'Are the solutions properly symmetric vectorized?\n\n', ...
    '    %s\n    %s'], obj.n, n, ...
    'x = vsdp.svec (obj, x(:), 2, ''unsymmetric'');', ...
    'z = vsdp.svec (obj, z(:), 1, ''unsymmetric'');');
end
obj.solution = sol;

end