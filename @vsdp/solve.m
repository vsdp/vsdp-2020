function obj = solve (obj, solver)
% SOLVE  Approximately solve the given conic problem instance.
%
%   See 'help vsdp.vsdp' for a description of the conic problem instance.
%
%   See 'help vsdp_options.SOLVER' for a list of supported solvers.
%
%   To use an initial guess (x0,y0,z0) type:
%
%      obj.x = x0;
%      obj.y = y0;
%      obj.z = z0;
%      obj.options.USE_STARTING_POINT = true;
%
%   Example:
%
%       A1 = [0 1;
%             1 0];
%       A2 = [1 1;
%             1 1];
%       C =  [1 0;
%             0 1];
%       K.s = 2;
%       At = [A1(:),  A2(:)];  % Vectorize data
%       c  = C(:);
%       b  = [1; 2.0001];
%       obj = vsdp(At, b, c, K).solve()
%
%
%   See also vsdpinit, vsdpup, vsdplow, vsdpinfeas.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (1, 2);
if (nargin == 2)
  obj.options.SOLVER = solver;
end

% Dispatch to solver.
switch (obj.options.SOLVER)
  case 'sdpt3'
    obj = obj.solve_sdpt3 ();
  case 'sedumi'
    obj = obj.solve_sedumi ();
  case 'sdpa'
    obj = obj.solve_sdpa ();
  case 'csdp'
    obj = solve_csdp ();
  case 'lp_solve'
    obj = solve_lp_solve ();
  case 'linprog'
    obj = solve_linprog ();
  case 'glpk'
    obj = solve_glpk ();
  otherwise
    error ('VSDP:solve:unsupportedCones', ...
      'solve: The solver %s could not be detected.', obj.options.SOLVER);
end

end
