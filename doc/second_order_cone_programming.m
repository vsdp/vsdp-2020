%% Second Order Cone Programming
%
% Let us consider a least squares problem taken from
% <https://vsdp.github.io/references.html#ElGhaoui1997 [ElGhaoui1997]>
%
% $$\|b - Ax\|_2 = \min_{y \in \mathbb{R}^{3}} \|b - Ay\|_2$$
%
% with singular matrix
%

A = [ 3 1 4 ; ...
      0 1 1 ; ...
     -2 5 3 ; ...
      1 4 5 ];

%%
% and right-hand side
%

b = [ 0 ; ...
      2 ; ...
      1 ; ...
      3 ];

%%
% This problem can be formulated as second-order cone program in dual standard
% form:
%
% $$\begin{array}{ll}
% \text{maximize}   & -y_{1} - y_{2}, \\
% \text{subject to}
% & y_{1} \geq \| (b - A y_{3:5} ) \|_{2}, \\
% & y_{2} \geq \begin{Vmatrix}\begin{pmatrix}
% 1 \\ y_{3} \\ y_{4} \\ y_{5} \end{pmatrix}\end{Vmatrix}_{2}, \\
% & y \in \mathbb{R}^{5}.
% \end{array}$$
%
% The two inequalities can be written in the form
% $$\begin{pmatrix} y_{1} \\ b - A y_{3:5} \end{pmatrix} \in \mathbb{L}^{5}
% \quad\text{and}\quad
% \begin{pmatrix} y_{2} \\ 1 \\ y_{3:5} \end{pmatrix} \in \mathbb{L}^{5}.$$
% Since
% $$\begin{pmatrix} y_{1} \\ q - P y_{3:5} \end{pmatrix} =
% \underbrace{\begin{pmatrix} 0 \\ 0 \\ 2 \\ 1 \\ 3 \end{pmatrix}}_{=c_{1}^{q}}
% - \underbrace{\begin{pmatrix}
% -1 & 0 &  0 & 0 & 0 \\
%  0 & 0 &  3 & 1 & 4 \\
%  0 & 0 &  0 & 1 & 1 \\
%  0 & 0 & -2 & 5 & 3 \\
%  0 & 0 &  1 & 4 & 5
% \end{pmatrix}}_{=(A_{1}^{q})^{T}}
% \begin{pmatrix} y_{1} \\ y_{2} \\ y_{3} \\ y_{4} \\ y_{5} \end{pmatrix},$$
% and
% $$\begin{pmatrix} y_{2} \\ 1 \\ y_{3:5} \end{pmatrix} =
% \underbrace{\begin{pmatrix} 0 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}}_{=c_{2}^{q}}
% - \underbrace{\begin{pmatrix}
% 0 & -1 &  0 &  0 &  0 \\
% 0 &  0 &  0 &  0 &  0 \\
% 0 &  0 & -1 &  0 &  0 \\
% 0 &  0 &  0 & -1 &  0 \\
% 0 &  0 &  0 &  0 & -1
% \end{pmatrix}}_{=(A_{2}^{q})^{T}}
% \begin{pmatrix} y_{1} \\ y_{2} \\ y_{3} \\ y_{4} \\ y_{5} \end{pmatrix},$$
% our dual problem \eqref{socpDual} takes the form
% $$\begin{array}{ll}
% \text{maximize}
% & \underbrace{\begin{pmatrix} -1 & -1 & 0 & 0 & 0 \end{pmatrix}}_{=b^{T}} y, \\
% \text{subject to}
% & z = \underbrace{\begin{pmatrix}
%                   0 \\ 0 \\ 2 \\ 1 \\ 3 \\ 0 \\ 1 \\ 0 \\ 0 \\ 0
%                   \end{pmatrix}}_{=c}
%     - \underbrace{\begin{pmatrix}
%                   -1 &  0 &  0 &  0 &  0 \\
%                    0 &  0 &  3 &  1 &  4 \\
%                    0 &  0 &  0 &  1 &  1 \\
%                    0 &  0 & -2 &  5 &  3 \\
%                    0 &  0 &  1 &  4 &  5 \\
%                    0 & -1 &  0 &  0 &  0 \\
%                    0 &  0 &  0 &  0 &  0 \\
%                    0 &  0 & -1 &  0 &  0 \\
%                    0 &  0 &  0 & -1 &  0 \\
%                    0 &  0 &  0 &  0 & -1
%                   \end{pmatrix}}_{=A^{T}} y \in K^{*},
% \end{array}$$
% where $K^{*} = \mathbb{L}^{5} \times \mathbb{L}^{5}$.
%
% We want to solve this problem with SeDuMi and enter the problem data of the
% primal problem.
%

c = [  0; 0; 2;  1; 3;  0; 1;  0;  0;  0];
A = [ -1, 0, 0,  0, 0,  0, 0,  0,  0,  0;
       0, 0, 0,  0, 0, -1, 0,  0,  0,  0;
       0, 3, 0, -2, 1,  0, 0, -1,  0,  0;
       0, 1, 1,  5, 4,  0, 0,  0, -1,  0;
       0, 4, 1,  3, 5,  0, 0,  0,  0, -1];
b = [-1; -1; 0; 0; 0];

%%
% Apart from the data |(A,b,c)|, the vector |q = [5;5]| of the second order
% cone block sizes must be forwarded to the structure |K|:
%

K.q = [5;5];

%%
% Now we compute approximate solutions by using |solve| and then verified
% error bounds by using |vsdplow| and |vsdpup|:
%

obj = vsdp (A, b, c, K);
obj.options.VERBOSE_OUTPUT = false;
obj.solve ('sdpt3');
obj.rigorous_lower_bound ();
obj.rigorous_upper_bound ();

%%
% The approximate primal and dual optimal values and the rigorous lower and
% upper bounds are
%

disp (obj)

%%
%

fL = obj.solutions.approximate.y

%%
% The quantities |x| and |y| are not displayed here. The two output vectors
% |lb| and |dl| provide rigorous lower bounds for the eigenvalues of these
% variables.  Since both vectors are positive
%

fL = obj.solutions.rigorous_lower_bound.f_objective(1)
fU = obj.solutions.rigorous_upper_bound.f_objective(2)

lb = obj.solutions.rigorous_lower_bound.z
dl = obj.solutions.rigorous_upper_bound.z

%%
% the primal and the dual set of feasible solutions have a non empty relative
% interior.  Thus Slater's constraint qualifications (*Strong Duality Theorem*)
% imply that the duality gap is zero.  The intervals $x$ and $y$ contain
% interior feasible $\epsilon$-optimal solutions, where
%

epsilon = 2 * (fU - fL) / (abs (fU) + abs (fL))

%%
% The conic programming allows to mix constraints of different types.
% For instance, one can add the linear inequality
% $\sum_{i=1}^{5} y_{i} \leq 3.5$ to the previous dual problem.  By using the
% standard form \eqref{cpPrim}, \eqref{cpDual}, and the condensed quantities
% \eqref{condensedX} it follows that we must add to the matrix |A| the column
% consisting of ones and to |c| the value 3.5.  We extend the input data as
% follows:
%

A = [[1; 1; 1; 1; 1], A];
c = [3.5; c];
K.l = 1;

%%
%

obj = vsdp (A, b, c, K).solve ('sdpt3');
obj.rigorous_lower_bound ();
obj.rigorous_upper_bound ();

%%
% Then we obtain
%

disp (obj)
