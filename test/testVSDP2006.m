function tests = testVSDP2006()
% TESTVSDP2006  Runs a testsuite for VSDP (version 2006).
%
%   Example:
%
%       addpath ('legacy')
%       clc; table (runtests ('testVSDP2006')) % Matlab
%       clc; testVSDP2006;                     % Octave
%
%   See also demovsdp.
%

% Copyright 2016-2019 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

tests = functiontests(localfunctions);
end

function [blk,A,C,b] = example_feasible()
blk(1,:) = {'s'; 2};
A{1,1} = [0 1; 1 0];
A{2,1} = [1 1; 1 1];
C{1}   = [1 0; 0 1];
b = [1; 2.0001];
end

function [blk,A,C,b] = example_DELTA(DELTA)
% DELTA > 0: the optimal value is -0.5, and strong duality holds,
% DELTA = 0: ill-posed problem with nonzero duality gap,
% DELTA < 0: primal and dual infeasible.
C{1} = [ 0   1/2    0;
  1/2 DELTA   0;
  0    0   DELTA];
A{1,1} = [  0  -1/2 0;
  -1/2   0  0;
  0    0  0];
A{2,1} = [1 0 0;
  0 0 0;
  0 0 0];
A{3,1} = [0 0 1;
  0 0 0;
  1 0 0];
A{4,1} = [0 0 0;
  0 0 1;
  0 1 0];
b = [1; 2*DELTA; 0; 0];
blk(1,:) = {'s'; 3};
end

function [blk,A,C,b,DELTA] = example_DELTA_feasible()
DELTA = 1e-4;
[blk,A,C,b] = example_DELTA(DELTA);
end

function [blk,A,C,b,DELTA] = example_DELTA_infeasible()
DELTA = -1e-4;
[blk,A,C,b] = example_DELTA(DELTA);
end

function [blk,A,C,b] = example_primal_infeasible()
blk(1,:) = {'s'; 2};
A{1,1} = [1 0; 0 0];
A{2,1} = [0 1; 1 0.005];
C{1}   = [0 0; 0 0];
b = [-0.01; 1];
end

function testVSDPCHECK_example_feasible(testCase)
[blk,A,C,b] = example_feasible();
[m,n] = vsdpcheck(blk,A,C,b);
verifyEqual(testCase, m, 2)
verifyEqual(testCase, n, 1)
end

function testVSDPCHECK_example_DELTA_feasible(testCase)
[blk,A,C,b,~] = example_DELTA_feasible();
[m,n] = vsdpcheck(blk,A,C,b);
verifyEqual(testCase, m, 4)
verifyEqual(testCase, n, 1)
end

function testVSDPCHECK_example_DELTA_infeasible(testCase)
[blk,A,C,b,~] = example_DELTA_infeasible();
[m,n] = vsdpcheck(blk,A,C,b);
verifyEqual(testCase, m, 4)
verifyEqual(testCase, n, 1)
end

function testVSDPCHECK_example_primal_infeasible(testCase)
[blk,A,C,b] = example_primal_infeasible();
[m,n] = vsdpcheck(blk,A,C,b);
verifyEqual(testCase, m, 2)
verifyEqual(testCase, n, 1)
end

function testMYSDPS_example_feasible(testCase)
[blk,A,C,b] = example_feasible();
[objt,X,y,Z,info] = mysdps(blk,A,C,b);
verifyEqual(testCase, objt, [1, 1], 'RelTol', 1e-4)
verifyEqual(testCase, X{1}, ones(2)/2, 'RelTol', 1e-4)
verifyEqual(testCase, objt(1), trace(C{1} * X{1}), 'RelTol', 1e-4)
verifyEqual(testCase, objt(2), b'*y, 'RelTol', 1e-4)
verifyEqual(testCase, Z{1}, C{1} - A{1,1}*y(1) - A{2,1}*y(2), 'RelTol', 1e-4)
verifyEqual(testCase, info, 0)
end

function testMYSDPS_example_DELTA_feasible(testCase)
[blk,A,C,b,DELTA] = example_DELTA_feasible();
[objt,X,y,Z,info] = mysdps(blk,A,C,b);
verifyEqual(testCase, objt, [-0.5 -0.5], 'RelTol', DELTA)
testX = X{1};
verifyEqual(testCase, testX(1:2,1:2), [2*DELTA, -1; -1, 1/(2*DELTA)], ...
  'RelTol', DELTA)
verifyEqual(testCase, testX(3,3) >= 0, true)
verifyEqual(testCase, y(2), -1/(4*DELTA), 'RelTol', DELTA)
verifyEqual(testCase, y([1,3,4]), [0;0;0], 'AbsTol', DELTA)
verifyEqual(testCase, objt(1), trace(C{1} * X{1}), 'RelTol', DELTA)
verifyEqual(testCase, objt(2), b'*y, 'RelTol', DELTA)
D = C{1};
for j = 1:4
  D = D - y(j) * A{j,1};
end
verifyEqual(testCase, Z{1}, D, 'RelTol', DELTA)
verifyEqual(testCase, info, 0)
end

function testMYSDPS_example_DELTA_infeasible(testCase)
[blk,A,C,b,~] = example_DELTA_infeasible();
[~,~,~,~,info] = mysdps(blk,A,C,b);
% indication of primal infeasibility
verifyEqual(testCase, info, 1)
end

function testMYSDPS_example_primal_infeasible(testCase)
[blk,A,C,b] = example_primal_infeasible();
[~,~,~,~,info] = mysdps(blk,A,C,b);
% indication of primal infeasibility
verifyEqual(testCase, info, 1)
end

function testVSDPUP_example_feasible(testCase)
[blk,A,C,b] = example_feasible();
[~,Xt,yt,Zt,~] = mysdps(blk,A,C,b);
[fU,X,lb] = vsdpup(blk,A,C,b,Xt,yt,Zt);
verifyEqual(testCase, fU >= 1.0001, true)
verifyEqual(testCase, all (all (in (X{1}, midrad(ones(2)/2,5e-4)))), true)
verifyEqual(testCase, all (lb >= 0), true)
end

function testVSDPUP_example_DELTA_feasible(testCase)
[blk,A,C,b,DELTA] = example_DELTA_feasible();
[~,Xt,yt,Zt,~] = mysdps(blk,A,C,b);
[fU,X,lb] = vsdpup(blk,A,C,b,Xt,yt,Zt);
verifyEqual(testCase, fU >= -0.5, true)
midX = [2*DELTA, -1, 0; -1, 1/(2*DELTA), 0; 0, 0, DELTA];
verifyEqual(testCase, all (all (in (X{1}, midrad(midX,2*DELTA)))), true)
verifyEqual(testCase, all (lb >= 0), true)
end

function testVSDPUP_example_DELTA_feasible_finite_bnd(testCase)
[blk,A,C,b,~] = example_DELTA_feasible();
[~,Xt,yt,Zt,~] = mysdps(blk,A,C,b);
yu = 1e5 * [1 1 1 1]';
[fU,X,lb] = vsdpup(blk,A,C,b,Xt,yt,Zt,yu);
% In this case: VSDPUP computes rigorously the finite upper bound without
% verifying primal feasibility.
verifyEqual(testCase, fU >= -0.5, true)
verifyEqual(testCase, isnan(X), true)
verifyEqual(testCase, isnan(lb), true)
end

function testVSDPUP_example_DELTA_infeasible(testCase)
[blk,A,C,b,~] = example_DELTA_infeasible();
[~,Xt,yt,Zt,~] = mysdps(blk,A,C,b);
S = warning ('off', 'VSDP:rigorous_upper_bound:unsolveablePerturbation');
[fU,X,lb] = vsdpup(blk,A,C,b,Xt,yt,Zt);
warning (S);
verifyEqual(testCase, isinf(fU), true)
verifyEqual(testCase, isnan(X), true)
verifyEqual(testCase, isnan(lb), true)
end

function testVSDPUP_example_primal_infeasible(testCase)
[blk,A,C,b] = example_primal_infeasible();
[~,Xt,yt,Zt,~] = mysdps(blk,A,C,b);
S = warning ('off', 'VSDP:rigorous_upper_bound:unsolveablePerturbation');
[fU,X,lb] = vsdpup(blk,A,C,b,Xt,yt,Zt);
warning (S);
verifyEqual(testCase, isinf(fU), true)
verifyEqual(testCase, isnan(X), true)
verifyEqual(testCase, isnan(lb), true)
end

function testVSDPLOW_example_feasible(testCase)
[blk,A,C,b] = example_feasible();
[~,Xt,yt,Zt,~] = mysdps(blk,A,C,b);
[fL,Y,dl] = vsdplow(blk,A,C,b,Xt,yt,Zt);
verifyEqual(testCase, fL <= 1.0001, true)
verifyEqual(testCase, all (in (Y, midrad([-1; 1],1e-4))), true)
verifyEqual(testCase, all (dl >= 0), true)
end

function testVSDPLOW_example_DELTA_feasible(testCase)
[blk,A,C,b,DELTA] = example_DELTA_feasible();
[~,Xt,yt,Zt,~] = mysdps(blk,A,C,b);
[fL,Y,dl] = vsdplow(blk,A,C,b,Xt,yt,Zt);
verifyEqual(testCase, fL <= -0.5, true)
verifyEqual(testCase, all (in (Y(2), ...
  midrad(-1/(4*DELTA), mag(Y(2))*DELTA))), true)
verifyEqual(testCase, all (in (Y([1,3,4]), midrad([0;0;0], DELTA))), true)
verifyEqual(testCase, all (dl >= 0), true)
end

function testVSDPLOW_example_DELTA_feasible_finite_bnd(testCase)
[blk,A,C,b,DELTA] = example_DELTA_feasible();
[~,Xt,yt,Zt,~] = mysdps(blk,A,C,b);
xu = 1e5;
[fL,Y,dl] = vsdplow(blk,A,C,b,Xt,yt,Zt,xu);
verifyEqual(testCase, fL <= -0.5, true)
verifyEqual(testCase, all (in (Y(2), ...
  midrad(-1/(4*DELTA), mag(Y(2))*DELTA))), true)
verifyEqual(testCase, all (in (Y([1,3,4]), midrad([0;0;0], DELTA))), true)
verifyEqual(testCase, all (dl >= -DELTA), true) % Small cone violation.
end

function testVSDPLOW_example_DELTA_infeasible(testCase)
[blk,A,C,b,~] = example_DELTA_infeasible();
[~,Xt,yt,Zt,~] = mysdps(blk,A,C,b);
S = warning ('off', 'VSDP:rigorous_lower_bound:unsolveablePerturbation');
[fL,Y,dl] = vsdplow(blk,A,C,b,Xt,yt,Zt);
warning (S);
verifyEqual(testCase, isinf(fL), true)
verifyEqual(testCase, isnan(Y), true)
verifyEqual(testCase, isnan(dl), true)
end

function testVSDPINFEAS_example_feasible(testCase)
[blk,A,C,b] = example_feasible();
[isinfeas,X,Y] = vsdpinfeas(blk,A,C,b,'p');
verifyEqual(testCase, isinfeas, 0)
verifyEqual(testCase, isnan(X), true)
verifyEqual(testCase, isnan(Y), true)
[isinfeas,X,Y] = vsdpinfeas(blk,A,C,b,'d');
verifyEqual(testCase, isinfeas, 0)
verifyEqual(testCase, isnan(X), true)
verifyEqual(testCase, isnan(Y), true)
end

function testVSDPINFEAS_example_DELTA_feasible(testCase)
[blk,A,C,b,~] = example_DELTA_feasible();
[isinfeas,X,Y] = vsdpinfeas(blk,A,C,b,'p');
verifyEqual(testCase, isinfeas, 0)
verifyEqual(testCase, isnan(X), true)
verifyEqual(testCase, isnan(Y), true)
[isinfeas,X,Y] = vsdpinfeas(blk,A,C,b,'d');
verifyEqual(testCase, isinfeas, 0)
verifyEqual(testCase, isnan(X), true)
verifyEqual(testCase, isnan(Y), true)
end

function testVSDPINFEAS_example_DELTA_infeasible(testCase)
% Note: should result in 'isinfeas = 1' in any case, but this is out of scope
% of VSDP.
[blk,A,C,b,~] = example_DELTA_infeasible();
[isinfeas,X,Y] = vsdpinfeas(blk,A,C,b,'p');
verifyEqual(testCase, isinfeas, 0)
verifyEqual(testCase, isnan(X), true)
verifyEqual(testCase, isnan(Y), true)
[isinfeas,X,Y] = vsdpinfeas(blk,A,C,b,'d');
verifyEqual(testCase, isinfeas, 0)
verifyEqual(testCase, isnan(X), true)
verifyEqual(testCase, isnan(Y), true)
end

function testVSDPINFEAS_example_primal_infeasible(testCase)
[blk,A,C,b] = example_primal_infeasible();
[isinfeas,X,Y] = vsdpinfeas(blk,A,C,b,'p');
verifyEqual(testCase, isinfeas, 1)
verifyEqual(testCase, isnan(X), true)
% TODO: improving ray Y does not match with a) value from help string of
%       vsdpinfeas or b) value from PDF documentation
%verifyEqual(testCase, Y, [-100.5998; -0.0060], 'RelTol', 1e-2)
%verifyEqual(testCase, Y, [-101.0835012952675; -0.01083501295267454], ...
%  'RelTol', 1e-2)
[isinfeas,X,Y] = vsdpinfeas(blk,A,C,b,'d');
verifyEqual(testCase, isinfeas, 0)
verifyEqual(testCase, isnan(X), true)
verifyEqual(testCase, isnan(Y), true)
end
