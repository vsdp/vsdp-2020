function obj = add_solution (obj, varargin)
% ADD_SOLUTION  Add solution object.
%
%   obj = add_solution (obj, sol) Add vsdp solution object 'sol'
%                                  to vsdp object 'obj'.
%
%   See also vsdp, vsdp_solution.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

if (nargin == 2)
  if (~isa (varargin{1}, 'vsdp_solution'))
    error ('VSDP:add_solution:badType', ...
      ['add_solution: A single argument must be a ''vsdp_solution'' ', ...
      'object.']);
  else
    sol = varargin{1};
  end
else
  sol = vsdp_solution (varargin{:});
end

[m, n] = sol.size ();
if ((m > 0) && (m ~= obj.m))
  error ('VSDP:add_solution:dimMismatch', ...
    ['add_solution: Dimension mismatch.  The number of problem ', ...
    'constraints is m = %d but the solution has %d.'], obj.m, m);
end
if ((n > 0) && (n ~= obj.n))
  error ('VSDP:add_solution:dimMismatch', ...
    ['add_solution: Dimension mismatch.  The number of problem ', ...
    'variables is n = %d but the solution has %d.\n\n', ...
    'Are the solutions properly symmetric vectorized?\n\n', ...
    '    %s\n    %s'], obj.n, n, ...
    'x = vsdp.svec (obj, x(:), 2, ''unsymmetric'');', ...
    'z = vsdp.svec (obj, z(:), 1, ''unsymmetric'');');
end

switch (sol.sol_type)
  case 'Initial'
    obj.solutions.initial = sol;
  case 'Approximate'
    obj.solutions.approximate = sol;
  case 'Rigorous lower bound'
    obj.solutions.rigorous_lower_bound = sol;
  case 'Rigorous upper bound'
    obj.solutions.rigorous_upper_bound = sol;
  case 'Certificate primal infeasibility'
    obj.solutions.certificate_primal_infeasibility = sol;
  case 'Certificate dual infeasibility'
    obj.solutions.certificate_dual_infeasibility = sol;
  otherwise
    error ('VSDP:add_solution:unknownType', ...
      'add_solution: Unknown solution type.');
end

end
