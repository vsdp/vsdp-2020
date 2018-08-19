%% Semidefinite Programming
%
% The primal standard form of a conic program with $n_{s}$ symmetric positive
% semidefinite cones
%
% $$\mathbb{S}^{s_{j}}_{+} := \left\{ X \in \mathbb{R}^{s_{j} \times s_{j}}
% \colon\; X = X^{T},\; v^{T} X v \geq 0,\; \forall v \in \mathbb{R}^{s_{j}}
% \right\},\quad j = 1,\ldots,n_{s}.$$
%
% is
%
% $$\begin{array}{lll}
% \text{minimize}
% & \sum_{j=1}^{n_{s}} \langle C_{j}, X_{j} \rangle & \\
% \text{subject to}
% & \sum_{j=1}^{n_{s}} \langle A_{ij}, X_{j} \rangle = b_{i},
% & i = 1,\ldots,m, \\
% & X_{j} \in \mathbb{S}^{s_{j}}_{+},
% & j = 1,\ldots,n_{s},
% \end{array}$$
%
% with symmetric $s_{j} \times s_{j}$ matrices $A_{ij}$ and $C_{j}$.
% The dual problem form is
%
% $$\begin{array}{ll}
% \text{maximize} & b^{T} y \\
% \text{subject to}
% & Z_{j} := C_{j} - \sum_{i=1}^{m} y_{i} A_{ij}
%   \in \mathbb{S}^{s_{j}}_{+},\quad j = 1, \ldots, n_{s}.
% \end{array}$$
%
%%


%% A feasible SDP
%
% We consider an example from the CSDP User's Guide
% <https://vsdp.github.io/references.html#Borchers2017 [Borchers2017]>:
%
% $$\begin{array}{lll}
% \text{minimize}
% & \sum_{j=1}^{3} \langle C_{j}, X_{j} \rangle & \\
% \text{subject to}
% & \sum_{j=1}^{3} \langle A_{ij}, X_{j} \rangle = b_{i},\quad
%      i = 1,2, \\
% & X_{1} \in \mathbb{S}^{2}_{+}, \\
% & X_{2} \in \mathbb{S}^{3}_{+}, \\
% & X_{3} \in \mathbb{S}^{2}_{+},
% \end{array}$$
%
% where $b = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$,
%
% $$\begin{array}{ccc}
%   C^{s_{1}}_{1} = \begin{pmatrix} -2 & -1 \\ -1 & -2 \end{pmatrix},
% & C^{s_{2}}_{2} =
%   \begin{pmatrix} -3 & 0 & -1 \\ 0 & -2 & 0 \\ -1 & 0 & -3 \end{pmatrix},
% & C^{s_{3}}_{3} = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}, \\
%   A^{s_{1}}_{1,1} = \begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix},
% & A^{s_{2}}_{1,2} =
%   \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix},
% & A^{s_{3}}_{1,3} = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, \\
%   A^{s_{1}}_{2,1} = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix},
% & A^{s_{2}}_{2,2} =
%   \begin{pmatrix} 3 & 0 & 1 \\ 0 & 4 & 0 \\ 1 & 0 & 5 \end{pmatrix},
% & A^{s_{3}}_{2,3} = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}.
% \end{array}$$
%
% In the vectorized format the corresponding coefficient matrix |At| and the
% primal objective vector |c| are
%

At{1} = [ 3; 1;
          1; 3;
          0; 0; 0;
          0; 0; 0;
          0; 0; 0;
          1; 0;
          0; 0 ];
At{2} = [ 0; 0;
          0; 0;
          3; 0; 1;
          0; 4; 0;
          1; 0; 5;
          0; 0;
          0; 1 ];
At = [At{:}];

b = [ 1;
      2 ];

c = [ -2; -1; 
      -1; -2;
      -3;  0; -1;
       0; -2;  0;
      -1;  0; -3;
       0;  0;
       0;  0];

%%
% And the cone structure |K| for this problem is
%

K.s = [2 3 2];

%%
% Now we compute approximate solutions by using |solve| and then verified
% error bounds by using |rigorous_lower_bound| and |rigorous_upper_bound|:
%

obj = vsdp (At, b, c, K);
obj.options.VERBOSE_OUTPUT = false;
obj.solve('sdpt3');
obj.rigorous_lower_bound();
obj.rigorous_upper_bound();

%%
% Finally, we get an overview about all the performed computations:
%

disp (obj)

%%
% To compare the approximate solution |X|, |y|, and |Z| with
% <https://vsdp.github.io/references.html#Borchers2017 [Borchers2017]> the
% vectorized solution quantities |x| and |z| have to be transformed back to
% matrices by using |vsdp.smat| and the appropriate scaling factor:
%

x = full(obj.solutions.approximate.x);
X = {vsdp.smat([], x(1:3),   1/2),
     vsdp.smat([], x(4:9),   1/2),
     vsdp.smat([], x(10:12), 1/2)}
y = obj.solutions.approximate.y
z = full(obj.solutions.approximate.z);
Z = {vsdp.smat([], z(1:3),   1),
     vsdp.smat([], z(4:9),   1),
     vsdp.smat([], z(10:12), 1)}

%%
% The compuation of the rigorous lower bounds involves the computation of the
% smallest eigenvalues |Zl(j)=| $\lambda_{\min}([Z_{j}])$ for $j = 1,2,3$.
%

Zl = obj.solutions.rigorous_lower_bound.z'
Y  = obj.solutions.rigorous_lower_bound.y

%%
% Since all |Zl >= 0| it is proven that all matrices $Z_{j}$ are in the
% interior of the cone
% $\mathbb{S}^{2}_{+} \times \mathbb{S}^{3}_{+} \times \mathbb{S}^{2}_{+}$
% and |Y| is a rigorous enclosure of a dual strict feasible (near optimal)
% solution.
%
% Analogous computations are performed for the rigorous upper bound.  Here
% lower bounds on the smallest eigenvalue of the primal solution are computed
% |Xl(j)=| $\lambda_{\min}([X_{j}])$ for $j = 1,2,3$.
%

Xl = obj.solutions.rigorous_upper_bound.z'

%%
% The matrix |X| is a rigorous enclosure of a primal strict feasible (near
% optimal) solution and can be restored from the vectorized quantity
% |obj.solutions.rigorous_upper_bound.x| as shown for the approximate solution.
% We omit the dispay of the interval matrix |X| for brevity.
%
% Since all |Xl| are positive, strict feasibility for the primal problem is
% proved.  Thus strong duality holds for this example.
%


%% An infeasible SDP
%
% Now we consider the following example
% (see <https://vsdp.github.io/references.html#Jansson2007a [Jansson2007a]>):
%
% $$\begin{array}{ll}
% \text{minimize} & \langle C(\delta), X \rangle \\
% \text{subject to}
% & \langle A_{1}, X \rangle = 1, \\
% & \langle A_{2}, X \rangle = \epsilon, \\
% & \langle A_{3}, X \rangle = 0, \\
% & \langle A_{4}, X \rangle = 0, \\
% & X \in \mathbb{S}^{3}_{+},
% \end{array}$$
%
% with Lagragian dual
%
% $$\begin{array}{ll}
% \text{maximize} & y_{1} + \epsilon y_{2} \\
% \text{subject to}
% & Z(\delta) := C(\delta) - \sum_{i = 1}^{4} A_{i} y_{i}
%   \in \mathbb{S}^{3}_{+}, \\
% & y \in \mathbb{R}^{4},
% \end{array}$$
%
% where
%

c = @(DELTA) ...
    [  0;   1/2;    0;
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

b = @(EPSILON) [1; EPSILON; 0; 0];

K.s = 3;

%%         
% The linear constraints of the primal problem form imply
%
% $$X(\epsilon) = \begin{pmatrix}
% \epsilon & -1 & 0 \\ -1 & X_{22} & 0 \\ 0 & 0 & X_{33}
% \end{pmatrix} \in \mathbb{S}^{3}_{+}$$
%
% iff $X_{22} \geq 0$, $X_{33} \geq 0$, and $\epsilon X_{22} - 1 \geq 0$.
% The conic constraint of the dual form is
%
% $$Z(\delta) = \begin{pmatrix}
% -y_{2} & \frac{1+y_{1}}{2} & -y_{3} \\
% \frac{1+y_{1}}{2} & \delta & -y_{4} \\
% -y_{3} & -y_{4} & \delta \end{pmatrix} \in \mathbb{S}^{3}_{+}.$$
%
% Hence, for
%
% * $\epsilon \leq 0$: the problem is primal infeasible $\hat{f_{p}} = +\infty$.
% * $\delta   \leq 0$: the problem is dual   infeasible $\hat{f_{d}} = -\infty$.
% * $\epsilon = \delta = 0$: the problem is ill-posed and there is a duality
%   gap with $\hat{f_{p}} = +\infty$ and $\hat{f_{d}} = -1$.
% * $\epsilon > 0$ and $\delta > 0$: the problem is feasible with
%   $\hat{f_{p}} = \hat{f_{d}} = -1 + \delta / \epsilon$.
%
% We start with the last feasible case and expect
% $\hat{f_{p}} = \hat{f_{d}} = -1 + 10$ with.
%

DELTA   = 10^(-3);
EPSILON = 10^(-4);

%%
%

obj = vsdp (At, b(EPSILON), c(DELTA), K);
obj.options.VERBOSE_OUTPUT = false;
obj.solve('sdpt3');
obj.rigorous_lower_bound();
obj.rigorous_upper_bound();

%%
%

disp (obj)

%%
% Nothing bad happened, as expected.
%
% Now we change the setting for primal infeasiblilty, what SDPT3 detects as
% well:
%

DELTA   = -10^(-3);
EPSILON = -10^(-4);
obj = vsdp (At, b(EPSILON), c(DELTA), K);
obj.options.VERBOSE_OUTPUT = false;
obj.solve('sdpt3');
obj.check_primal_infeasible();
disp (obj)

%%
% The value of the return parameter |info| confirms successful termination of
% the solver.  The first eight decimal digits of the primal and dual optimal
% values are correct, since $\hat{f}_{p} = \hat{f}_{d} = -0.5$, all components
% of the approximate solutions |xt| and |yt| have at least five correct decimal
% digits.  Nevertheless, successful termination reported by a solver gives no
% guarantee on the quality of the computed solution.
%
% For instance, if we apply SeDuMi to the same problem we obtain:
%

vsdpinit('sedumi');
[objt,xt,yt,zt,info] = mysdps(A,b,c,K);
objt, xt, yt, info  % zt:  hidden for brevity

%%
% SeDuMi terminates without any warning, but some results are poor.  Since the
% approximate primal optimal value is smaller than the dual one, weak duality
% is not satisfied. In other words, the algorithm is not backward stable for
% this example.  The CSDP-solver gives similar results:
%

vsdpinit('sdpt3'); %TODO: csdp
[objt,xt,yt,zt,info] = mysdps(A,b,c,K);
objt, xt, yt, info  % zt:  hidden for brevity

%%
% A good deal worse are the results that can be derived with older versions of
% these solvers, including SDPT3 and SDPA
% </references#Jansson2006 [Jansson2006]>.
%
% Reliable results can be obtained by the functions |vsdplow| and |vsdpup|.
% Firstly, we consider |vsdplow| and the approximate solver SDPT3.
%

vsdpinit('sdpt3');
[objt,xt,yt,zt,info] = mysdps(A,b,c,K);
[fL,y,dl] = vsdplow(A,b,c,K,xt,yt,zt)

%%
% the vector |y| is a rigorous interior dual $??$-optimal solution where we shall
% see that $?? \approx 2.27 \times 10^{-8}$.  The positivity of |dl| verifies
% that |y| contains a dual strictly feasible solution.  In particular, strong
% duality holds.  By using SeDuMi similar rigorous results are obtained.  But
% for the SDPA-solver we get
%

vsdpinit('sdpa');
[objt,xt,yt,zt,info] = mysdps(A,b,c,K);
[fL,y,dl] = vsdplow(A,b,c,K,xt,yt,zt)

%%
% Thus, an infinite lower bound for the primal optimal value is obtained and
% dual feasibility is not verified.  The rigorous lower bound strongly depends
% on the computed approximate solution and therefore on the used approximate
% conic solver.
%
% Similarly, a verified upper bound and a rigorous enclosure of a primal
% $??$-optimal solution can be computed by using the |vsdpup| function
% together with SDPT3:
%

vsdpinit('sdpt3');
[objt,xt,yt,zt,info] = mysdps(A,b,c,K);
[fU,x,lb] = vsdpup(A,b,c,K,xt,yt,zt)

%%
% The output |fU| is close to the dual optimal value $\hat{f}_{d} = -0.5$.
% The interval vector |x| contains a primal strictly feasible solution, see
% \eqref{OptSolSDPExp}, and the variable |lb| is a lower bound for the smallest
% eigenvalue of |x|.  Because |lb| is positive, Slater's condition is fulfilled
% and strong duality is verified once more.
%
% Summarizing, by using SDPT3 for the considered example with parameter
% $?? = 10^{-4}$, VSDP verified strong duality with rigorous bounds for the
% optimal value
% $$
% -0.500000007 \leq \hat{f}_{p} = \hat{f}_{d} \leq -0.499999994.
% $$
%
% The rigorous upper and lower error bounds of the optimal value show only
% modest overestimation.  Strictly primal and dual feasible solutions are
% obtained.  Strong duality is verified.  Moreover, we have seen that the
% quality of the rigorous results depends strongly on the quality of the
% computed approximations.
%
