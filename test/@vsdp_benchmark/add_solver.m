function bool = add_solver (obj, name, check_fun, setup_dir, setup_fun)
% ADD_SOLVER  Sets up a solver for a VSDP benchmark.
%
%   name      The name of the solver.  This name should match the definition in
%             VSDP.
%   check_fun Function to check if the solver is ready to use.
%
%   setup_dir (optional)  Directory where to setup the solver.
%   setup_fun (optional)  Function  call  to setup the solver.
%
%   Returns true on success.
%
%   See also vsdp_benchmark.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (3, 5);

% Check for solver to be ready.
if (check_fun ())
  bool = add_solver_data (obj, name, check_fun, setup_dir, setup_fun);
  return;
end

% Try to setup the solver.
bool = false;
if (nargin > 3)
  setup_dir = obj.check_dir (setup_dir);
  if (~isempty (setup_dir))
    OLD_DIR = cd (setup_dir);
    setup_fun ();
    cd (OLD_DIR);
  end
end

% Check for solver to be ready.
if (check_fun ())
  bool = add_solver_data (obj, name, check_fun, setup_dir, setup_fun);
else
  warning ('VSDP_BENCHMARK:add_solver:solverNotReady', ...
    'bm_setup: The solver ''%s'' is not available.', name);
end

end

function bool = add_solver_data (obj, name, check_fun, setup_dir, setup_fun)
solver.name = name;
% TODO: Octave bug #43215, cannot store function handles, thus store as string.
solver.check_fun = func2str (check_fun);
solver.setup_dir = setup_dir;
solver.setup_fun = func2str (setup_fun);
obj.SOLVER = [obj.SOLVER, solver];
bool = true;
end
