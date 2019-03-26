function obj = solve (obj, varargin)
% SOLVE  Guided approximate solution the given conic problem instance.
%
%   obj.SOLVE()  Let the user choose interactively the solver to be used.  A
%                valid decision is saved in 'obj.options.SOLVER'.
%   obj.SOLVE('solver')  Solve the problem instance with 'solver'.  This call
%                        is identical to the call 'obj.solve_solver()'.  The
%                        decision is not saved in 'obj.options.SOLVER'.
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

% Generate a list of by VSDP supported solvers.
solver_list = solver.registry.list_available ();

supported_solvers_str = ['''', strjoin(solver_list, ''', '''), ''''];

if (~isempty (obj.options.SOLVER))
  use_solver = obj.options.SOLVER;
elseif (nargin == 1)  % Interactive mode.
  initial_value = find (strcmp (obj.options.SOLVER, solver_list));
  if (isempty (initial_value))
    initial_value = find (strcmp ('sdpt3', solver_list));
  end
  idx = -1;
  % Prompt the user to choose between the supported solvers.
  % If 'obj.options.SOLVER' was set properly, the option is surrounded by braces
  % and chosen when just hitting enter.
  while ((idx < 1) || (idx > length (solver_list)))
    fprintf ('\n');
    for i = 1:length (solver_list)
      if (i == initial_value)
        fprintf ('  [%2d] %s\n', i, solver_list{i});
      else
        fprintf ('   %2d  %s\n', i, solver_list{i});
      end
    end
    fprintf ('\n');
    str = input (sprintf ('Choose one of the above solvers [1-%d]: ', ...
      length (solver_list)), 's');
    if (isempty (str))
      idx = initial_value;
    else
      idx = str2double (str);
    end
  end
  use_solver = solver_list{idx};
  obj.options.SOLVER = use_solver;
else  % Non-interactive mode.
  try
    use_solver = validatestring (varargin{1}, solver_list);
  catch
    error ('VSDP:solve:unsupportedSolver', ...
      ['The solver ''%s'' is not supported by VSDP.  ', ...
      'Choose one of:\n\n    %s\n\nand set ''obj.options.SOLVER''.'], ...
      varargin{1}, supported_solvers_str);
  end
  if (isempty (obj.options.SOLVER))
    obj.options.SOLVER = use_solver;
  end
end

% Dispatch to solver.
obj = eval (['solver.', use_solver, '.solve (obj, varargin{2:end});']);

end
