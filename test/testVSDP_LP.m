function testVSDP_LP (testCase)
LP_VSDP_2012_P09 (testCase);  % K.l = 5
LP_VSDP_2012_P11 (testCase);  % K.l = 2 and K.f = 1
end

function LP_VSDP_2012_P09 (testCase)
use_solvers = {'sedumi', 'sdpt3'};
if (exist ('OCTAVE_VERSION', 'builtin'))
  use_solvers{end+1} = 'glpk';
else
  use_solvers{end+1} = 'linprog';
end

A = [-1, 2,  0, 1, 1;
  0, 0, -1, 0, 2];
b = [2; 3];
c = [0, 2, 0, 3, 5];
K.l = 5;

x_sol = [0; 0.25; 0; 0; 1.5];
y_sol = [1; 2];

obj = vsdp (A, b, c, K);
obj.options.VERBOSE_OUTPUT = false;

for i = 1:length(use_solvers)
  obj.options.SOLVER = use_solvers{i};
  obj.solve (obj.options.SOLVER);
  sol = obj.solutions('Approximate solution');
  verifyEqual (testCase, sol.solver_info.termination, 'Normal termination');
  verifyEqual (testCase, full (sol.x), x_sol, ...
    'AbsTol', 1e-7, 'RelTol', eps ());
  verifyEqual (testCase, full (sol.y), y_sol, ...
    'AbsTol', 1e-7, 'RelTol', eps ());
  
  obj.rigorous_lower_bound ();
  sol = obj.solutions('Rigorous lower bound');
  verifyEqual (testCase, sol.solver_info.termination, 'Normal termination');
  verifyEqual (testCase, sol.f_objective(1) <= (obj.c' * x_sol), true);
  verifyEqual (testCase, isnan (sol.f_objective(2)), true);
  verifyEqual (testCase, isempty (sol.x), true);
  verifyEqual (testCase, isintval (sol.y), true);
  verifyEqual (testCase, isreal (sol.z), true);
  verifyEqual (testCase, all (sol.z) >= 0, true);
end
end

function LP_VSDP_2012_P11 (testCase)
use_solvers = {'sedumi', 'sdpt3'};
if (exist ('OCTAVE_VERSION', 'builtin'))
  use_solvers{end+1} = 'glpk';
else
  use_solvers{end+1} = 'linprog';
end

A = [2, 1, -1;
  -1, 1, 1];
b = [0.5; 1];
c = [-0.5, 1, 1];
K.f = 1;
K.l = 2;

x_sol = [-1/6; 5/6; 0];
y_sol = [ 1/6; 5/6];

obj = vsdp (A, b, c, K);
obj.options.VERBOSE_OUTPUT = false;

for i = 1:length(use_solvers)
  obj.options.SOLVER = use_solvers{i};
  obj.solve (obj.options.SOLVER);
  sol = obj.solutions('Approximate solution');
  verifyEqual (testCase, sol.solver_info.termination, 'Normal termination');
  verifyEqual (testCase, full (sol.x), x_sol, ...
    'AbsTol', 1e-7, 'RelTol', eps ());
  verifyEqual (testCase, full (sol.y), y_sol, ...
    'AbsTol', 1e-7, 'RelTol', eps ());
  
  obj.rigorous_lower_bound ();
  sol = obj.solutions('Rigorous lower bound');
  verifyEqual (testCase, sol.solver_info.termination, 'Normal termination');
  verifyEqual (testCase, sol.f_objective(1) <= (obj.c' * x_sol), true);
  verifyEqual (testCase, isnan (sol.f_objective(2)), true);
  verifyEqual (testCase, isempty (sol.x), true);
  verifyEqual (testCase, isintval (sol.y), true);
  verifyEqual (testCase, isreal (sol.z), true);
  verifyEqual (testCase, all (sol.z) >= 0, true);
end
end