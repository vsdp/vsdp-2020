function testVSDP_SOCP (testCase)
% Original problem in https://vsdp.github.io/references#ElGhaoui1997 (page
% 1052).
SOCP_VSDP_2012_P13 (testCase);  % K.q = [5; 5]
end

function SOCP_VSDP_2012_P13 (testCase)
use_solvers = {'sedumi', 'sdpt3'};

A = [ ...
  -1, 0, 0,  0, 0,  0, 0,  0,  0,  0;
  +0, 0, 0,  0, 0, -1, 0,  0,  0,  0;
  +0, 3, 0, -2, 1,  0, 0, -1,  0,  0;
  +0, 1, 1,  5, 4,  0, 0,  0, -1,  0;
  +0, 4, 1,  3, 5,  0, 0,  0,  0, -1];
b = [-1, -1, 0, 0, 0];
c = [0, 0, 2, 1, 3, 0, 1, 0, 0, 0]';
K.q = [5; 5];

P = A(3:5,2:5)';
q = c(2:5);

obj = vsdp (A, b, c, K);
obj.options.VERBOSE_OUTPUT = false;

for i = 1:length(use_solvers)
  obj.options.SOLVER = use_solvers{i};
  obj.solve (obj.options.SOLVER);
  sol = obj.solutions('Approximate');
  verifyEqual (testCase, sol.solver_info.termination, 'Normal termination');
  verifyEqual (testCase, ...
    round (sol.y(1), 7) >= round (norm (q - P*sol.y(3:5)), 7), true);
  verifyEqual (testCase, ...
    round (sol.y(2), 7) >= round (norm ([1; sol.y(3:5)]), 7), true);
  
  obj.rigorous_lower_bound ();
  sol = obj.solutions('Rigorous lower bound');
  verifyEqual (testCase, sol.solver_info.termination, 'Normal termination');
  verifyEqual (testCase, isfinite (sol.f_objective(1)), true);
  verifyEqual (testCase, isnan (sol.f_objective(2)), true);
  verifyEqual (testCase, isempty (sol.x), true);
  verifyEqual (testCase, isintval (sol.y), true);
  verifyEqual (testCase, isreal (sol.z), true);
  verifyEqual (testCase, all (sol.z) >= 0, true);
  
  obj.rigorous_upper_bound ();
  sol = obj.solutions('Rigorous upper bound');
  verifyEqual (testCase, sol.solver_info.termination, 'Normal termination');
  verifyEqual (testCase, isnan (sol.f_objective(1)), true);
  verifyEqual (testCase, isfinite (sol.f_objective(2)), true);
  verifyEqual (testCase, isintval (sol.x), true);
  verifyEqual (testCase, isempty (sol.y), true);
  verifyEqual (testCase, isreal (sol.z), true);
  verifyEqual (testCase, all (sol.z) >= 0, true);
end
end