function testVSDP_SDP (testCase)
% Original problem in https://vsdp.github.io/references#Borchers2017 (pages
% 8 -- 14).
SDP_VSDP_2012_P15 (testCase);  % K.s = [2; 3; 2]
end

function SDP_VSDP_2012_P15 (testCase)
use_solvers = {'sedumi', 'sdpt3', 'sdpa', 'csdp'};

A = [3, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
  0, 0, 0, 0, 3, 0, 1, 0, 4, 0, 1, 0, 5, 0, 0, 0, 1];
c = [-2; -1; -1; -2; -3; 0; -1; 0; -2; 0; -1; 0; -3; 0; 0; 0; 0];
b = [1; 2];
K.s = [2; 3; 2];

% In Borchers2017 "2/3" is rounded to 0.667 what is too inaccurate.  A short
% evaluation with pen&paper assures that 'tr(A2*x) = 2  ==>  x(4) = 2/3'.
x_sol = [0.125; 0.25; 0.125; 2/3; zeros(8,1)];
y_sol = [-0.75; -1];
z_sol = [0.25; -0.25; 0.25; 0; 0; 2; 0; 0; 2; 0.75; 0; 1];

obj = vsdp (A, b, c, K);
obj.options.VERBOSE_OUTPUT = false;

for i = 1:length(use_solvers)
  obj.options.SOLVER = use_solvers{i};
  obj.solve (obj.options.SOLVER);
  sol = obj.solutions.approximate;
  verifyEqual (testCase, sol.solver_info.termination, 'Normal termination');
  % May solvers have problems with x([4,7]), thus big tolerance 1e-4
  verifyEqual (testCase, full (sol.x), x_sol, ...
    'AbsTol', 1e-4, 'RelTol', eps ());
  verifyEqual (testCase, full (sol.y), y_sol, ...
    'AbsTol', 1e-7, 'RelTol', eps ());
  verifyEqual (testCase, full (sol.z), z_sol, ...
    'AbsTol', 1e-7, 'RelTol', eps ());
  
  obj.rigorous_lower_bound ();
  sol = obj.solutions.rigorous_lower_bound;
  verifyEqual (testCase, sol.solver_info.termination, 'Normal termination');
  verifyEqual (testCase, sol.f_objective(1) <= (obj.c' * x_sol), true);
  verifyEqual (testCase, isnan (sol.f_objective(2)), true);
  verifyEqual (testCase, isempty (sol.x), true);
  verifyEqual (testCase, isintval (sol.y), true);
  verifyEqual (testCase, isreal (sol.z), true);
  verifyEqual (testCase, all (sol.z) >= 0, true);
  
  obj.rigorous_upper_bound ();
  sol = obj.solutions.rigorous_upper_bound;
  verifyEqual (testCase, sol.solver_info.termination, 'Normal termination');
  verifyEqual (testCase, isnan (sol.f_objective(1)), true);
  verifyEqual (testCase, sol.f_objective(2) >= (obj.b' * y_sol), true);
  verifyEqual (testCase, isintval (sol.x), true);
  verifyEqual (testCase, isempty (sol.y), true);
  verifyEqual (testCase, isreal (sol.z), true);
  verifyEqual (testCase, all (sol.z) >= 0, true);
end
end