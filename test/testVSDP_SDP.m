function testVSDP_SDP (testCase)
% Original problem in https://vsdp.github.io/references#Borchers2017 (pages
% 8 -- 14).
SDP_VSDP_2012_P15 (testCase);  % K.s = [2; 3; 2]
end

function SDP_VSDP_2012_P15 (testCase)
use_solvers = {'sedumi', 'sdpt3'};

A = [3, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
  0, 0, 0, 0, 3, 0, 1, 0, 4, 0, 1, 0, 5, 0, 0, 0, 1];
c = [-2; -1; -1; -2; -3; 0; -1; 0; -2; 0; -1; 0; -3; 0; 0; 0; 0];
b = [1; 2];
K.s = [2; 3; 2];

% Only 3 digits accuracy in Borchers2017, assume at least 4.
x_sol = [0.125; 0.25; 0.125; 0.66667; zeros(8,1)];
y_sol = [-0.75; -1];
z_sol = [0.25; -0.25; 0.25; 0; 0; 2; 0; 0; 2; 0.75; 0; 1];

obj = vsdp (A, b, c, K);
obj.options.VERBOSE_OUTPUT = false;

for i = 1:length(use_solvers)
  obj.options.SOLVER = use_solvers{i};
  obj.solve ();
  verifyEqual (testCase, full (obj.x), x_sol, ...
    'AbsTol', 1e-4, 'RelTol', eps ());
  verifyEqual (testCase, full (obj.y), y_sol, ...
    'AbsTol', 1e-7, 'RelTol', eps ());
  disp(use_solvers{i})
  verifyEqual (testCase, full (obj.z), z_sol, ...
    'AbsTol', 1e-7, 'RelTol', eps ());
end
end