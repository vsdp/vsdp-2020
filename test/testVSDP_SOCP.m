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
  y = obj.solution.y;
  verifyEqual (testCase, obj.solution.info, 0);
  verifyEqual (testCase, ...
    round (y(1), 7) >= round (norm (q - P*y(3:5)), 7), true);
  verifyEqual (testCase, ...
    round (y(2), 7) >= round (norm ([1; y(3:5)]), 7), true);
end
end