%% Second Order Cone Programming
%
% This section explains how to work with second order cone problems.  Due to
% \eqref{stdPrim} the second order cone problem in standard primal form is
% given by
% $$
% \label{socpPrim}
% \begin{array}{ll}
% \text{minimize}
% & \langle c^{f}, x^{f} \rangle +
%   \sum_{i=1}^{n_{q}} \langle c_{i}^{q}, x_{i}^{q} \rangle, \\
% \text{subject to}
% & A^{f} x^{f} + \sum_{i=1}^{n_{q}} A_{i}^{q} x_{i}^{q} = b, \\
% & x^{f} \in \mathbb{R}^{n_{f}}, \\
% & x_{i}^{q} \in \mathbb{L}^{q_{i}}, \quad i = 1, \ldots, n_{q}.
% \end{array}
% $$
% The corresponding dual problem is
% $$
% \label{socpDual}
% \begin{array}{ll}
% \text{maximize}   & b^{T} y, \\
% \text{subject to}
% & (A_{i}^{q})^{T} y + z_{i}^{q} = c_{i}^{q},
% \quad z_{i}^{q} \in \mathbb{L}^{q_{i}},
% \quad i = 1,\ldots,n_{q}, \\
% & (A^{f})^{T} y + z^{f} = c^{f},
% \quad z^{f} \in \mathbb{R}^{n_{f}}, \\
% & z^{f} = 0.
% \end{array}
% $$
%
% Let us consider the total least squares problem taken from
% </references#ElGhaoui1997 [ElGhaoui1997]>:
% $$
% \label{LSQexample}
% \begin{array}{ll}
% \text{maximize}   & -y_{1} - y_{2}, \\
% \text{subject to}
% & y_{1} \geq || (q - P y_{3:5} ) ||_{2}, \\
% & y_{2} \geq
%   \begin{Vmatrix}\begin{pmatrix} 1 \\ y_{3:5} \end{pmatrix}\end{Vmatrix}_{2},
% \end{array}
% $$
% where
% $$
% P = \begin{pmatrix}
%  3 & 1 & 4 \\
%  0 & 1 & 1 \\
% -2 & 5 & 3 \\
%  1 & 4 & 5
% \end{pmatrix},
% \quad q = \begin{pmatrix} 0 \\ 2 \\ 1 \\ 3 \end{pmatrix},
% \quad y \in \mathbb{R}^5.
% $$
% We want to transform this problem to the dual form \eqref{socpDual}.
% The two inequalities can be written in the form
% $$
% \begin{pmatrix} y_{1} \\ q - P y_{3:5} \end{pmatrix} \in \mathbb{L}^{5}
% \quad\text{and}\quad
% \begin{pmatrix} y_{2} \\ 1 \\ y_{3:5} \end{pmatrix} \in \mathbb{L}^{5}.
% $$
% Since
% $$
% \begin{pmatrix} y_{1} \\ q - P y_{3:5} \end{pmatrix} =
% \underbrace{\begin{pmatrix} 0 \\ 0 \\ 2 \\ 1 \\ 3 \end{pmatrix}}_{=c_{1}^{q}}
% - \underbrace{\begin{pmatrix}
% -1 & 0 &  0 & 0 & 0 \\
%  0 & 0 &  3 & 1 & 4 \\
%  0 & 0 &  0 & 1 & 1 \\
%  0 & 0 & -2 & 5 & 3 \\
%  0 & 0 &  1 & 4 & 5
% \end{pmatrix}}_{=(A_{1}^{q})^{T}}
% \begin{pmatrix} y_{1} \\ y_{2} \\ y_{3} \\ y_{4} \\ y_{5} \end{pmatrix},
% $$
% and
% $$
% \begin{pmatrix} y_{2} \\ 1 \\ y_{3:5} \end{pmatrix} =
% \underbrace{\begin{pmatrix} 0 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}}_{=c_{2}^{q}}
% - \underbrace{\begin{pmatrix}
% 0 & -1 &  0 &  0 &  0 \\
% 0 &  0 &  0 &  0 &  0 \\
% 0 &  0 & -1 &  0 &  0 \\
% 0 &  0 &  0 & -1 &  0 \\
% 0 &  0 &  0 &  0 & -1
% \end{pmatrix}}_{=(A_{2}^{q})^{T}}
% \begin{pmatrix} y_{1} \\ y_{2} \\ y_{3} \\ y_{4} \\ y_{5} \end{pmatrix},
% $$
% our dual problem \eqref{socpDual} takes the form
% $$
% \label{SOCPexample}
% \begin{array}{ll}
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
% \end{array}
% $$
% where $K^{*} = \mathbb{L}^{5} \times \mathbb{L}^{5}$.
%
% We want to solve this problem with SeDuMi and enter the problem data of the
% primal problem.
%

clear A b c K

vsdpinit('sedumi');
c = [0; 0; 2; 1; 3; 0; 1; 0; 0; 0];
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
% Now we compute approximate solutions by using |mysdps| and then verified
% error bounds by using |vsdplow| and |vsdpup|:
%

[objt,xt,yt,zt,info] = mysdps(A,b,c,K);
[fL,y,dl] = vsdplow(A,b,c,K,xt,yt,zt);
[fU,x,lb] = vsdpup (A,b,c,K,xt,yt,zt);

%%
% The approximate primal and dual optimal values and the rigorous lower and
% upper bounds are
%

objt, fL, fU

%%
% The quantities |x| and |y| are not diplayed here. The two output vectors
% |lb| and |dl| provide rigorous lower bounds for the eigenvalues of these
% variables.  Since both vectors are positive
%

lb, dl

%%
% the primal and the dual set of feasible solutions have a non empty relative
% interior.  Thus Slater's constraint qualifications (*Strong Duality Theorem*)
% imply that the duality gap is zero.  The intervals $x$ and $y$ contain
% interior feasible $??$-optimal solutions, where
% $?? = \dfrac{2 (fU-fL)}{|fU|+|fL|} \approx 1.85 \times 10^{-9}$.
%
% The conic program \eqref{cpPrim} allows to mix constraints of different types.
% Let us, for instance, add the linear inequality
% $\sum_{i=1}^{5} y_{i} \leq 3.5$ to the previous dual problem.  By using the
% standard form \eqref{cpPrim}, \eqref{cpDual}, and the condensed quantities
% \eqref{condensedX} it follows that we must add to the matrix |A| the column
% consisting of ones and to |c| the value 3.5.  We extend the input data as
% follows:
%

A = [[1; 1; 1; 1; 1], A];
c = [3.5; c];
K.l = 1;
[objt,xt,yt,zt,info] = mysdps(A,b,c,K);
fL = vsdplow(A,b,c,K,xt,yt,zt);
fU = vsdpup (A,b,c,K,xt,yt,zt);

%%
% Then we obtain
%

fL, fU
