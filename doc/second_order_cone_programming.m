%% Second-order Cone Programming
%
%%

%%
% Consider a least squares problem from
% <https://vsdp.github.io/references.html#ElGhaoui1997 [ElGhaoui1997]>
%
% $$\left\|b_{data} - A_{data}\,\hat{y}\right\|_2
% = \min_{y_{3:5} \in \mathbb{R}^{3}}
% \left\|b_{data} - A_{data}\,y_{3:5}\right\|_2$$
%
% with singular matrix
%

A_data = [ 3 1 4 ;
           0 1 1 ;
          -2 5 3 ;
           1 4 5 ];

%%
% and right-hand side
%

b_data = [ 0 ;
           2 ;
           1 ;
           3 ];

%%
% This problem can be formulated as second-order cone program in dual standard
% form:
%
% $$\begin{array}{ll}
% \text{maximize}   & -y_{1} - y_{2}, \\
% \text{subject to}
% & y_{1} \geq \| (b_{data} - A_{data}\,y_{3:5} ) \|_{2}, \\
% & y_{2} \geq
% \begin{Vmatrix}\begin{pmatrix} 1 \\ y_{3:5} \end{pmatrix}\end{Vmatrix}_{2}, \\
% & y \in \mathbb{R}^{5}.
% \end{array}$$
%
% The two inequality constraints can be written as second-order cone vectors
%
% $$\begin{pmatrix} y_{1} \\ b_{data} - A_{data}\,y_{3:5} \end{pmatrix}
% \in \mathbb{L}^{5} \quad\text{and}\quad
% \begin{pmatrix} y_{2} \\ 1 \\ y_{3:5} \end{pmatrix} \in \mathbb{L}^{5}.$$
%
% Both vectors can be expressed as matrix-vector product of $y$
%
% $$\underbrace{\begin{pmatrix} 0 \\ b_{data} \end{pmatrix}}_{=c_{1}^{q}}
% - \underbrace{\begin{pmatrix}
% -1 & 0 & 0 & 0 & 0 \\
%  0 & 0 & ( & A_{data} & )
% \end{pmatrix}}_{=(A_{1}^{q})^{T}}
% \begin{pmatrix} y_{1} \\ y_{2} \\ y_{3} \\ y_{4} \\ y_{5} \end{pmatrix}
% \in \mathbb{L}^{5}$$
% and
%
% $$\underbrace{\begin{pmatrix} 0 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}}_{=c_{2}^{q}}
% - \underbrace{\begin{pmatrix}
% 0 & -1 &  0 &  0 &  0 \\
% 0 &  0 &  0 &  0 &  0 \\
% 0 &  0 & -1 &  0 &  0 \\
% 0 &  0 &  0 & -1 &  0 \\
% 0 &  0 &  0 &  0 & -1
% \end{pmatrix}}_{=(A_{2}^{q})^{T}}
% \begin{pmatrix} y_{1} \\ y_{2} \\ y_{3} \\ y_{4} \\ y_{5} \end{pmatrix}
% \in \mathbb{L}^{5}.$$
%
% With these formulations, the dual problem takes the form
% $$\begin{array}{ll}
% \text{maximize}
% & \underbrace{\begin{pmatrix} -1 & -1 & 0 & 0 & 0 \end{pmatrix}}_{=b^{T}} y, \\
% \text{subject to}
% & z = \underbrace{\begin{pmatrix}
%                   0 \\ b_{data} \\ 0 \\ 1 \\ 0 \\ 0 \\ 0
%                   \end{pmatrix}}_{=c}
%     - \underbrace{\begin{pmatrix}
%                   -1 &  0 &  0 &  0 &  0 \\
%                    0 &  0 &  ( & A_{data} & ) \\
%                    0 & -1 &  0 &  0 &  0 \\
%                    0 &  0 &  0 &  0 &  0 \\
%                    0 &  0 & -1 &  0 &  0 \\
%                    0 &  0 &  0 & -1 &  0 \\
%                    0 &  0 &  0 &  0 & -1
%                   \end{pmatrix}}_{=A^{T}} y \in K^{*}, \\
% & y \in \mathbb{R}^{5}.
% \end{array}$$
% where $K^{*} = \mathbb{L}^{5} \times \mathbb{L}^{5}$.
%
% We want to solve this problem with SeDuMi and enter the problem data of the
% primal problem.
%

At = zeros (10, 5);
At(1,1) = -1;
At(2:5, 3:5)  = A_data;
At(6,2) = -1;
At(8:10, 3:5) = -eye(3);
b = [-1 -1 0 0 0]';
c = [ 0 b_data' 0 0 0 0 0]';

%%
% Apart from the data |(At,b,c)|, the vector |q = [5;5]| of the second-order
% cone block sizes must be forwarded to the structure |K|:
%

K.q = [5;5];

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
% Now we analyze the resulting regularized least squares solution |y_SOCP =|
% $y_{3:5}$

y_SOCP = obj.solutions.approximate.y(3:5)

%%
% and compare it to a naive least squares solution |y_LS|
%

y_LS = A_data \ b_data

%%
%

[                  norm(y_SOCP)                    norm(y_LS);
 norm(b_data - A_data * y_SOCP)  norm(b_data - A_data * y_LS)]


%%
% The conic programming allows to mix constraints of different types.
% For instance, one can add the linear inequality
% $\sum_{i=1}^{5} y_{i} \leq 3.5$ to the previous dual problem.  We extend the
% input data as follows:
%

At = [1 1 1 1 1; At];
c =  [3.5      ; c ];
K.l = 1;

%%
% Note that the order of the cone variables matters for |At| and |c|:
%
% # |K.f| Free              variables
% # |K.l| Linear            variables
% # |K.q| Second-order cone variables
% # |K.s| Semidefinite      variables
%

obj = vsdp (At, b, c, K);
obj.options.VERBOSE_OUTPUT = false;
obj.solve('sdpt3');
obj.rigorous_lower_bound();
obj.rigorous_upper_bound();

%%
% Then we obtain
%

disp (obj)
