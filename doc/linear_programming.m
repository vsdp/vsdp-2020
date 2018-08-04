%% Linear Programming
%
% In this section we describe how linear programming problems can be solved
% with VSDP.  In particular, two linear programming examples are considered
% in detail.  Each conic problem is fully described by the four variables
% $(A,b,c,K)$.  The first two quantities represent the affine constraints
% $Ax = b$.  The third is the primal objective vector |c|, and the last
% describes the underlying cone.  The cone |K| is a structure with four fields:
% |K.f|, |K.l|, |K.q|, and |K.s|.  The field |K.f| stores the number of free
% variables $n_{f}$, the field |K.l| stores the number of nonnegative variables
% $n_{l}$, the field |K.q| stores the dimensions $q_{1}, \ldots, q_{n_{q}}$ of
% the second order cones, and similarly |K.s| stores the dimensions $s_{1},
% \ldots, s_{n_{s}}$ of the semidefinite cones.  If a component of |K| is
% empty, then it is assumed that the corresponding cone do not occur.
%
% Consider the linear programming problem
% $$
% \label{LP1}
% \begin{array}{ll}
% \text{minimize}   & 2x_{2} + 3x_{4} + 5x_{5}, \\
% \text{subject to} &
% \begin{pmatrix}
% -1 & 2 &  0 & 1 & 1 \\
%  0 & 0 & -1 & 0 & 2
% \end{pmatrix} x = \begin{pmatrix} 2 \\ 3 \end{pmatrix}, \\
% & x \in \mathbb{R}^{5}_{+},
% \end{array}
% $$
% with its corresponding dual problem
% $$
% \label{LP1Dual}
% \begin{array}{ll}
% \text{maximize}   & 2 y_{1} + 3 y_{2}, \\
% \text{subject to} &
% z = \begin{pmatrix} 0 \\ 2 \\ 0 \\ 3 \\ 5 \end{pmatrix} -
% \begin{pmatrix}
% -1 &  0 \\
% 2 &  0 \\
% 0 & -1 \\
% 1 &  0 \\
% 1 &  2
% \end{pmatrix} y \in \mathbb{R}^{5}_{+},
% \end{array}
% $$
%
% The unique exact optimal solution is given by
% $x^{*} = \begin{pmatrix} 0 & 0.25 & 0 & 0 & 1.5 \end{pmatrix}^{T}$,
% $y^{*} = \begin{pmatrix} 1 &2 \end{pmatrix}^{T}$ with
% $\hat{f}_{p} = \hat{f}_{d} = 8$.
%
% The input data of the problem in VSDP are
%

A = [-1, 2,  0, 1, 1;
      0, 0, -1, 0, 2];
b = [2; 3];
c = [0; 2; 0; 3; 5];
K.l = 5;

%%
% By using |vsdpinit| the approximate conic solver can be set globally for all
% VSDP functions. Here we choose the solver SDPT3:
%

vsdpinit('sdpt3');

%%
% If we call the |mysdps| for problem \eqref{LP1} this problem will be solved
% approximately, yielding the output
%

format infsup long
[objt,xt,yt,zt,info] = mysdps(A,b,c,K)

%%
% The returned variables have the following meaning: |objt| stores the
% approximate primal and dual optimal value, |xt| represents the approximate
% primal solution in vectorized form, and |(yt,zt)| is the dual solution pair.
% The last returned variable |info| gives information on the success of the
% conic solver:
%
% * |info = 0|: successful termination,
% * |info = 1|: the problem might be primal infeasible,
% * |info = 2|: the problem might be dual infeasible,
% * |info = 3|: the problem might be primal and dual infeasible.
%
% With the approximate solution, a verified lower bound of the primal optimal
% value can be computed by the function |vsdplow|:
%

[fL,y,dl] = vsdplow(A,b,c,K,xt,yt,zt)

%%
% The output consists of
%
% # a verified lower bound for the optimal value stored in |fL|,
% # a rigorous interval enclosure of a dual feasible solution |y|, and
% # a componentwise lower bound of $z = c - A^{T} y$ stored in |dl|.
%
% Since |dl| is positive, the dual problem is strictly feasible, and the
% rigorous interval vector |y| contains a dual interior solution.  Here only
% some significant digits of this interval vector are displayed.  The upper
% and lower bounds of the interval |y| can be obtained by using the |sup| and
% |inf| routines of INTLAB.  For more information about the |intval| data type
% see </references#Rump1999 [Rump1999]>.
%
% Next we compute an upper bound for the optimal value by using |vsdpup|:
%

[fU,x,lb] = vsdpup(A,b,c,K,xt,yt,zt)

%%
% The returned variables have a similar meaning to those of |vsdplow|: |fU| is
% a verified upper bound for the optimal value, |x| is a rigorous interval
% enclosure of a primal feasible solution, and |lb| is a componentwise lower
% bound for |x|.  Since |lb| is a positive vector, hence contained in the
% positive orthant, the primal problem is strictly feasible.  As for the
% interval vector |y| also for the interval vector |x| only some significant
% digits are displayed.  The quantity |x| is proper interval vector.  This
% becomes clear when displaying the first component with |midrad|
%

midrad(x(1))

%%
% or |infsup|
%

infsup(x(1))

%%
% Summarizing, we have obtained a primal dual interval solution pair that
% contains a primal and dual strictly feasible $??$-optimal solution, where
% $?? = \frac{2 (fU-fL)}{|fU|+|fL|} \approx 4.83 \times 10^{-9}$.
%
% How a linear program with free variables can be solved with VSDP is
% demonstrated by the following example with one free variable $x_{3}$:
% $$
% \begin{array}{ll}
% \text{minimize}   & x_{1} + x_{2} - 0.5 x_{3}, \\
% \text{subject to}
% & x_{1} - x_{2} + 2 x_{3} = 0.5, \\
% & x_{1} + x_{2} -   x_{3} = 1, \\
% & x_{1} \geq 0, \\
% & x_{2} \geq 0.
% \end{array}
% $$
%
% The optimal solution pair of this problem is
% $x^{*} = \begin{pmatrix} \dfrac{5}{6} & 0 & -\dfrac{1}{6} \end{pmatrix}^{T}$,
% $y^{*} = \begin{pmatrix} \dfrac{1}{6} & \dfrac{5}{6}\end{pmatrix}^{T}$ with
% $\hat{f}_{p} = \hat{f}_{d} = \dfrac{11}{12}$.
%
% When entering a problem the order of the variables is important: Firstly free
% variables, secondly nonnegative variables, thirdly second order cone
% variables, and last semidefinite variables.  This order must be maintained
% in the matrix |A| as well as in the primal objective |c|.  In the given
% example, the free variable is $x_{3}$, the nonnegative variables are $x_{1}$,
% $x_{2}$.  Second order cone variables and semidefinite variables are not
% present.  Therefore, the problem data are
%

K.f = 1;  % number of free variables
K.l = 2;  % number of nonnegative variables
A = [ 2, 1, -1;   % first column corresponds to free variable x3
     -1, 1,  1];  % second and third to bounded x1, x2
c = [-0.5; 1; 1]; % the same applies to c
b = [0.5; 1];

%%
% Rigorous bounds for the optimal value can be optained with:
%

[objt,xt,yt,zt,info] = mysdps(A,b,c,K);
fL = vsdplow(A,b,c,K,xt,yt,zt)
fU = vsdpup (A,b,c,K,xt,yt,zt)
