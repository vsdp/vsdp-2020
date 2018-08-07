%% Linear Programming
%
% In this section we describe how linear programming problems can be solved
% with VSDP.  In particular, two linear programs are considered in detail.
%%

%% First example
%
% Consider the linear program in primal standard form
% $$\begin{array}{ll}
% \text{minimize}   & 2x_{2} + 3x_{4} + 5x_{5}, \\
% \text{subject to} &
% \begin{pmatrix}
% -1 & 2 &  0 & 1 & 1 \\
%  0 & 0 & -1 & 0 & 2
% \end{pmatrix} x = \begin{pmatrix} 2 \\ 3 \end{pmatrix}, \\
% & x \in \mathbb{R}^{5}_{+},
% \end{array}$$
% with its corresponding dual problem
% $$\begin{array}{ll}
% \text{maximize}   & 2 y_{1} + 3 y_{2}, \\
% \text{subject to} &
% z = \begin{pmatrix} 0 \\ 2 \\ 0 \\ 3 \\ 5 \end{pmatrix} -
% \begin{pmatrix}
% -1 &  0 \\
%  2 &  0 \\
%  0 & -1 \\
%  1 &  0 \\
%  1 &  2
% \end{pmatrix} y \in \mathbb{R}^{5}_{+}.
% \end{array}$$
%
% The unique exact optimal solution is given by
% $x^{*} = (0, 0.25, 0, 0, 1.5)^{T}$, $y^{*} = (1, 2)^{T}$ with
% $\hat{f_{p}} = \hat{f_{d}} = 8$.
%
% The input data of the problem in VSDP are
%

A = [-1, 2,  0, 1, 1;
      0, 0, -1, 0, 2];
b = [2; 3];
c = [0; 2; 0; 3; 5];
K.l = 5;

%%
% To create a VSDP object of the linear program data above, we call the VSDP
% class constructor and do not suppress the output.
%

obj = vsdp (A, b, c, K)

%%
% The output seems quite verbose in the beginning, but it contains all relevant
% information about the possibilities of VSDP, including the commands, that
% can be run to compute rigorous bounds or certificates for infeasibility.
% After these computations, the results are displayed as short summary in the
% same context.
%
% By calling the |solve| method on the VSDP object |obj|, we can compute an
% approximate solution |x|, |y|, |z|, for example by using SDPT3.  When calling
% |solve| without any arguments, the user is asked to choose one of the
% supported solvers.
%

obj.solve ('sdpt3');

%%
% The solver output is often quite verbose.  Especially for large problem
% instances it is recommended to display the solver progress.  To suppress
% solver messages, the following option can be set:
%

obj.options.VERBOSE_OUTPUT = false;

%%
% To permanently assign an approximate solver to a VSDP object, use the
% following option:
%

obj.options.SOLVER = 'sdpt3';

%%
% By simply typing the VSDP object's name, the user gets a short summary of
% the solution state.
%

obj

%%
% On success, one can obtain the approximate |x| and |y| solutions.
%

format short
x = obj.solutions.approximate.x
y = obj.solutions.approximate.y

%%
% The approximate solution is close to the optimal solution
% $x^{*} = (0, 0.25, 0, 0, 1.5)^{T}$, $y^{*} = (1, 2)^{T}$.
%
% With this approximate solution, a rigorous lower bound |fL| of the primal
% optimal value $\hat{f_{p}}$ can be computed by calling:
%

obj.rigorous_lower_bound ();

format long
fL = obj.solutions.rigorous_lower_bound.f_objective(1)

%%
% Similarly, a rigorous upper bound |fU| of the dual optimal value $\hat{f_{d}}$
% can be computed by calling:
%

obj.rigorous_upper_bound ();

fU = obj.solutions.rigorous_upper_bound.f_objective(2)

%%
% All this information is available in the summary of the VSDP object and must
% only be extracted if necessary.
%

obj

%%
% Despite the rigorous lower bound |fL|, the solution object
% |obj.solutions.rigorous_lower_bound| contain more information:
%

format short
format infsup
Y = obj.solutions.rigorous_lower_bound.y

%%
%

Dl = obj.solutions.rigorous_lower_bound.z

%%
% # |Y| is a rigorous interval enclosure of a dual feasible near optimal
%   solution and
% # |Dl| a lower bound of of each cone in $z = c - A^{*} y$.  For a linear
%   program this is a lower bound on each component of |z|.
%
% Since |Dl| is positive, the dual problem is strictly feasible, and the
% rigorous interval vector |Y| contains a dual interior solution.  Here only
% some significant digits of this interval vector are displayed.  The upper
% and lower bounds of the interval |Y| can be obtained by using the |sup| and
% |inf| routines of INTLAB.  For more information about the |intval| data type
% see <https://vsdp.github.io/references.html#Rump1999 [Rump1999]>.
%

%%
% The information returned by |rigorous_upper_bound()| is similar:
%

X = obj.solutions.rigorous_upper_bound.x

%%
%

Xl = obj.solutions.rigorous_upper_bound.z

%%
% # |X| is a rigorous interval enclosure of a primal feasible near optimal
%   solution and
% # |Xl| a lower bound of of each cone in |X|.  Again, for a linear program
%   this is a lower bound on each component of |X|.
%
% Since |Xl| is a positive vector, |X| is contained in the positive orthant and
% the primal problem is strictly feasible.
%
% Summarizing, we have obtained a primal dual interval solution pair with an
% accuracy measured by
% $$\mu(a, b) = \dfrac{a-b}{\max\{1.0, (|a| + |b|)/2\}},$$
% see <https://vsdp.github.io/references.html#Jansson2006 [Jansson2006]>.
%

format shorte
mu = (fU - fL) / max (1, (abs (fU) + abs(fL)) / 2)


%% Second example with free variables
%
% How a linear program with free variables can be solved with VSDP is
% demonstrated by the following example with one free variable $x_{3}$:
% $$\begin{array}{ll}
% \text{minimize}   & x_{1} + x_{2} - 0.5 x_{3}, \\
% \text{subject to}
% & x_{1} - x_{2} + 2 x_{3} = 0.5, \\
% & x_{1} + x_{2} -   x_{3} = 1, \\
% & x_{1} \geq 0, \\
% & x_{2} \geq 0.
% \end{array}$$
%
% The optimal solution pair of this problem is
% $x^{*} = (\frac{5}{6}, 0, -\frac{1}{6})^{T}$,
% $y^{*} = (\frac{1}{6}, \frac{5}{6})^{T}$ with
% $\hat{f_{p}} = \hat{f_{d}} = \frac{11}{12} \approx 9.166\ldots$.
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
% The whole VSDP computation can be done in a few lines of code:
%

obj = vsdp (A, b, c, K);
obj.options.VERBOSE_OUTPUT = false;
obj.solve ('sdpt3').rigorous_lower_bound ().rigorous_upper_bound ();

%%
% Yielding
%

obj
