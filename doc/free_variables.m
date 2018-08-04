%% Handling Free Variables
%
% Free variables occur often in practice.  Handling free variables in interior
% point algorithms is a pending issue (see for example
% </references#Andersen2002 [Andersen2002]>,
% </references#Anjos2007 [Anjos2007]>, and
% </references#Meszaros1998 [Meszaros1998]>).  Frequently users
% convert a problem with free variables into one with restricted variables by
% representing the free variables as a difference of two nonnegative variables.
% This approach increases the problem size and introduces ill-posedness, which
% may lead to numerical difficulties.
%
% For an example we consider the test problem _nb_L1_ from the DIMACS test
% library </references#Pataki2002 [Pataki2002]>.  The problem
% originates from side lobe minimization in antenna engineering.  This is a
% second order cone programming problem with 915 equality constraints, 793 SOCP
% blocks each of size 3, and 797 nonnegative variables.  Moreover, the problem
% has two free variables that are described as the difference of four
% nonnegative variables.  This problem can be loaded from the examples
% directory of VSDP.  As the computation is more expensive, only the results
% are reported here:
%
%
%   vsdpinit('sdpt3');
%   load(fullfile('examples','nb_L1.mat'));
%   [objt,xt,yt,zt,info] = mysdps(A,b,c,K);
%   objt
%
%   objt =
%    -13.012270628163670 -13.012270796164543
%
%
% SDPT3 solves the problem without warnings, although it is ill-posed according
% to Renegar's definition </references#Renegar1994 [Renegar1994]>.
%
% Now we try to get rigorous bounds using the approximation of SDPT3.
%
%
%   fL = vsdplow(A,b,c,K,xt,yt,zt)
%   fU = vsdpup (A,b,c,K,xt,yt,zt)
%
%   fL =
%     -Inf
%   fU =
%    -13.012270341861644
%
%
% These results reflect that the interior of the dual feasible solution set is
% empty.  An ill-posed problem has the property that the distance to primal or
% dual infeasibility is zero.  If as above the distance to dual infeasibility
% is zero, then there are sequences of dual infeasible problems with input data
% converging to the input data of the original problem. Each problem of the
% sequence is dual infeasible and thus has the dual optimal solution $-\infty$.
% Hence, the result $-\infty$ of |vsdplow| is exactly the limit of the optimal
% values of the dual infeasible problems and reflects the fact that the
% distance to dual infeasibility is zero.  This demonstrates that the infinite
% bound computed by VSDP is sharp, when viewed as the limit of a sequence of
% infeasible problems.  We have a similar situation if the distance to primal
% infeasibility is zero.
%
% If the free variables are not converted into restricted ones then the problem
% is well-posed and a rigorous finite lower bound can be computed.
%
%
%   load(fullfile('examples','nb_L1free.mat'));
%   [objt,xt,yt,zt,info] = mysdps(A,b,c,K);
%   objt
%
%   objt =
%    -13.012270619970032 -13.012270818869100
%
%
% By using the computed approximations we obtain the following rigorous bounds:
%
%
%   fL = vsdplow(A,b,c,K,xt,yt,zt)
%   fU = vsdpup (A,b,c,K,xt,yt,zt)
%
%   fL =
%    -13.012270819014953
%   fU =
%    -13.012270617556419
%
% Therefore, without splitting the free variables, we get rigorous finite lower
% and upper bounds of the exact optimal value with an accuracy of about eight
% decimal digits.  Moreover, verified interior solutions are computed for both
% the primal and the dual problem, proving strong duality.
%
% In Table <benchmark_dimacs_free_2012_12_12.html> we display rigorous bounds
% for the optimal value of eight problems contained in the DIMACS test library
% that have free variables (see </references#Anjos2007 [Anjos2007]>
% and </references#Kobayashi2007 [Kobayashi2007]>).  These
% problems have been modified by reversing the substitution of the free
% variables.  We have listed the results for the problems with free variables
% and for the same problems when representing the free variables as the
% difference of two nonnegative variables.  The table contains the rigorous
% upper bounds $fU$, the rigorous lower bounds $fL$, and the computing times
% measured in seconds for the approximate solution $t_s$, the lower bound
% $t_u$, and the upper bound $t_l$, respectively.  The table demonstrates the
% drastic improvement if free variables are not split.
%
% Independent of the transformation of the free variables the primal problems
% of the _nql_ instances are ill-posed.  The weak error bound of the optimal
% constraints.  A solution for the _qssp180_ instance is due to the large
% number of equality system with 130141 equality constraints and 261365
% variables has to be solved rigorously.  In the next version of VSDP the
% accuracy for such large problems will be improved.
%
