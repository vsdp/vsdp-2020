function testVSDP_SDP (testCase)
% Original problem in https://vsdp.github.io/references#Borchers2017 (pages
% 8 -- 14).
SDP_VSDP_2012_P15 (testCase);      % K.s = [2; 3; 2]
% Original problem in https://vsdp.github.io/references#Jansson2007a (page 192).
SDP_Jansson2007a_P192 (testCase);  % K.s = 3
end

function SDP_VSDP_2012_P15 (testCase)
use_solvers = {'sedumi', 'sdpt3', 'sdpa', 'csdp', 'mosek'};
use_solvers = intersect (use_solvers, solver.registry.list_available ());
if (isempty (use_solvers))
  warning ('VSDP:testVSDP:noSolver', ...
    'testVSDP_LP: skipping test, no solver available.');
end

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

for i = 1:length(use_solvers)
  % Make a clean copy.
  obj = vsdp (obj);
  obj.options.VERBOSE_OUTPUT = false;
  obj.options.SOLVER = use_solvers{i};
  if (strcmp (use_solvers{i}, 'sdpa'))
    obj.options.ALPHA = 20;
  end
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


function SDP_Jansson2007a_P192 (testCase)
use_solvers = {'sedumi', 'sdpt3', 'csdp', 'mosek'};
use_solvers = intersect (use_solvers, solver.registry.list_available ());
if (isempty (use_solvers))
  warning ('VSDP:testVSDP:noSolver', ...
    'testVSDP_LP: skipping test, no solver available.');
end

DELTA   = 1e-4;
EPSILON = 2 * DELTA;

c = [  0;   1/2;    0;
      1/2; DELTA;   0;
       0;    0;   DELTA ];

At = {};
At{1} = [ 0; -1/2; 0;
        -1/2;  0;  0;
          0;   0;  0 ];
At{2} = [ 1; 0; 0;
          0; 0; 0;
          0; 0; 0 ];
At{3} = [ 0; 0; 1;
          0; 0; 0;
          1; 0; 0 ];
At{4} = [ 0; 0; 0;
          0; 0; 1;
          0; 1; 0 ];
At = [At{:}];

b = [1; EPSILON; 0; 0];

K.s = 3;

obj = vsdp (At, b, c, K);

% A priori error bounds
xu = 1e5;
yu = 1e5 * [1 1 1 1]';

% An evaluation of Jansson2007a (p. 193) yields for EPSILON = 2e-4
% the values X_22 = 5000, X_33 = 0 and y, repectivly.
x_sol = [ EPSILON;        -1; 0;
               -1; 1/EPSILON; 0;
                0;         0; 0 ];
x_sol = vsdp.svec (obj, x_sol, 2);
y_sol = [0; -1/(2*EPSILON); 0; 0];
z_sol = [ 1/(2*EPSILON);   1/2; 0;
                    1/2; DELTA; 0;
                      0;     0; DELTA ];
z_sol = vsdp.svec (obj, z_sol, 2);

for i = 1:length(use_solvers)
  % Make a clean copy.
  obj = vsdp (obj);
  obj.options.VERBOSE_OUTPUT = false;
  obj.options.SOLVER = use_solvers{i};
  obj.solve (obj.options.SOLVER);
  sol = obj.solutions.approximate;
  verifyEqual (testCase, sol.solver_info.termination, 'Normal termination');
  verifyEqual (testCase, full (sol.x(1:3)), x_sol(1:3), 'RelTol', 1e-5);
  verifyEqual (testCase, full (sol.x(4:6)), x_sol(4:6), 'AbsTol', 1e-1);
  verifyEqual (testCase, full (sol.y(2)), y_sol(2), 'RelTol', 1e-2);
  verifyEqual (testCase, full (sol.y([1,3,4])), y_sol([1,3,4]), 'AbsTol', 1e-2);
  verifyEqual (testCase, full (sol.z([1:3,6])), z_sol([1:3,6]), 'RelTol', 1);
  verifyEqual (testCase, full (sol.z(4:5)), z_sol(4:5), 'AbsTol', 1e-5);
  
  obj.rigorous_lower_bound ();
  sol = obj.solutions.rigorous_lower_bound;
  verifyEqual (testCase, sol.solver_info.termination, 'Normal termination');
  verifyEqual (testCase, sol.f_objective(1) <= (obj.c' * x_sol), true);
  verifyEqual (testCase, isnan (sol.f_objective(2)), true);
  verifyEqual (testCase, isempty (sol.x), true);
  verifyEqual (testCase, isintval (sol.y), true);
  verifyEqual (testCase, isreal (sol.z), true);
  verifyEqual (testCase, all (sol.z) >= 0, true);
  
  obj.rigorous_lower_bound (xu);
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
  
  obj.rigorous_upper_bound (yu);
  sol = obj.solutions.rigorous_upper_bound;
  verifyEqual (testCase, sol.solver_info.termination, 'Normal termination');
  verifyEqual (testCase, isnan (sol.f_objective(1)), true);
  verifyEqual (testCase, sol.f_objective(2) >= (obj.b' * y_sol), true);
  verifyEqual (testCase, isempty (sol.x), true);
  verifyEqual (testCase, isempty (sol.y), true);
  verifyEqual (testCase, isempty (sol.z), true);
  verifyEqual (testCase, all (sol.z) >= 0, true);
end
end
