function obj = solve (obj, varargin)
% SOLVE  Guided approximate solution the given conic problem instance.
%
%   obj.SOLVE()  Let the user choose interactively the solver to be used.  A
%                valid decision is saved in 'obj.options.SOLVER'.
%   obj.SOLVE('solver')  Solve the problem instance with 'solver'.  The
%                        decision is saved in 'obj.options.SOLVER' if empty.
%
%   See 'help vsdp' for a description of the conic problem instance.
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
%       At = [A1(:), A2(:)];  % Vectorize data
%       c  = C(:);
%       b  = [1; 2.0001];
%       obj = vsdp(At, b, c, K).solve()  % Interactive solver choice
%       obj = vsdp(At, b, c, K).solve('sdpt3')
%       obj = vsdp(At, b, c, K).solve_sdpt3()
%
%   See also vsdp.

% Copyright 2004-2019 Christian Jansson (jansson@tuhh.de)

narginchk (1, 3);

% Get list of available solvers with appropriate cone support.
slist_cs = solver.registry.list_cone_support (obj);
slist = intersect (solver.registry.list_available (), slist_cs);

% Treat foreseeable error due to cone support.
if (isempty (slist) ...
    || ((nargin > 1) && (isempty (intersect (slist, varargin{1})))))
  if (isempty (slist_cs))
    err_msg = 'VSDP is not working for this problem.';
  elseif (nargin > 1)
    err_msg = sprintf ('''%s'' does not support this problem.  %s %s.', ...
      varargin{1}, ...
      'Please use one of:', strjoin (slist_cs, ', '));
  else
    err_msg = sprintf ('No available solver supports this problem.  %s %s.', ...
      'Please use one of:', strjoin (slist_cs, ', '));
  end
  error ('VSDP:solve:unsupportedSolver', 'solve: %s', err_msg);
end

if (nargin == 1)
  if (isempty (obj.options.SOLVER))
    sol = select_solver_by_menu (obj, slist);
  else
    sol = obj.options.SOLVER;
  end
else
  sol = varargin{1};
end

% Dispatch to solver
obj = eval (['solver.', sol, '.solve (obj, varargin{2:end});']);
if (isempty (obj.options.SOLVER))
  obj.options.SOLVER = sol;
end
end


function sol = select_solver_by_menu (obj, slist)
% SELECT_SOLVER_BY_MENU  Display a cmd user dialog to select a solver.

sidx = find (strcmp (obj.options.SOLVER, slist));
if (isempty (sidx))
  sidx = 1;
end
idx = -1;
% Prompt the user to choose and available and supported solver.
while ((idx < 1) || (idx > length (slist)))
  fprintf ('\n');
  for i = 1:length (slist)
    if (i == sidx)
      fprintf ('  [%2d] %s\n', i, slist{i});
    else
      fprintf ('   %2d  %s\n', i, slist{i});
    end
  end
  fprintf ('\n');
  str = input (sprintf ('Choose one of the above solvers [1-%d]: ', ...
    length (slist)), 's');
  if (isempty (str))
    idx = sidx;
  else
    idx = str2double (str);
  end
end
sol = slist{idx};
end
