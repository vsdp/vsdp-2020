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

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (1, 3);

% Generate a list of by VSDP supported solvers.
supported_solvers = methods (obj);
supported_solvers = supported_solvers(strncmp (supported_solvers, 'solve_', 6));
supported_solvers = cellfun (@(x) x(7:end), supported_solvers, ...
  'UniformOutput', false);  % Strip 'solve_' prefix.

supported_solvers_str = ['''', strjoin(supported_solvers, ''', '''), ''''];

if (nargin == 1)  % Interactive mode.
  initial_value = find (strcmp (obj.options.SOLVER, supported_solvers));
  idx = -1;
  % Prompt the user to choose between the supported solvers.
  % If 'obj.options.SOLVER' was set properly, the option is surrounded by braces
  % and chosen when just hitting enter.
  while ((idx < 1) || (idx > length (supported_solvers)))
    fprintf ('\n');
    for i = 1:length (supported_solvers)
      if (i == initial_value)
        fprintf ('  [%2d] %s\n', i, supported_solvers{i});
      else
        fprintf ('   %2d  %s\n', i, supported_solvers{i});
      end
    end
    fprintf ('\n');
    str = input (sprintf ('Choose one of the above solvers [1-%d]: ', ...
      length (supported_solvers)), 's');
    if (isempty (str))
      idx = initial_value;
    else
      idx = str2double (str);
    end
  end
  solver = supported_solvers{idx};
  obj.options.SOLVER = solver;
else  % Non-interactive mode.
  try
    solver = validatestring (varargin{1}, supported_solvers);
  catch
    error ('VSDP:solve:unsupportedSolver', ...
      ['The solver ''%s'' is not supported by VSDP.  ', ...
      'Choose one of:\n\n    %s'], varargin{1}, supported_solvers_str);
  end
end

% Dispatch to solver.
obj = eval (['obj.solve_', solver, '(varargin{2:end})']);

end
