%% A Priori Bounds
%
% In many practical applications the order of the magnitude of a primal or dual
% optimal solution is known a priori.  This is the case in many combinatorial
% optimization problems, or, for instance, in truss topology design where the
% design variables such as bar volumes can be roughly bounded.  If such bounds
% are available they can speed up the computation of guaranteed error bounds
% for the optimal value substantially, see
% <https://vsdp.github.io/references.html#Jansson2006 [Jansson2006]>.
%
%%

%%
% For linear programming problems the upper bound for the variable $x^{l}$
% is a vector $\bar{x}$ such that $x^{l} \leq \bar{x}$.  For second
% order cone programming the upper bounds for block variables $x_{i}^{q}$
% with $i = 1,\ldots,n_{q}$ can be entered as a vector of upper bounds
% $\overline{\lambda}_{i}$ of the largest eigenvalues
% $$\lambda_{\max}(x_{i}^{q}) = (x_{i}^{q})_{1} + ||(x_{i}^{q})_{:}||_{2}.$$
% Similarly, for semidefinite programs upper bounds for the primal variables
% $X_{j}^{s}$ can be entered as a vector of upper bounds of the largest
% eigenvalues $\lambda_{\max}(X_{j}^{s})$, $j = 1,\ldots,n_{s}$. An upper bound
% $\bar{y}$ for the dual optimal solution $y$ is a vector which is
% componentwise larger then $|y|$. Analogously, for conic programs with free
% variables the upper bound can be entered as a vector $\bar{x}$ such that
% $|x^{f}| \leq \bar{x}$.
%
% As an example, we consider the previous SDP problem with an upper bound
% $xu = 10^{5}$ for $\lambda_{\max}(X)$.
%

c = [  0;   1/2;      0;
      1/2; 10^(-3);   0;
       0;    0;     10^(-3) ];

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

b = [1; 10^(-4); 0; 0];

K.s = 3;

obj = vsdp (At, b, c, K);
obj.options.VERBOSE_OUTPUT = false;

%%
% Now we compute approximate solutions by using |solve| and then verified
% error bounds by using |rigorous_lower_bound| and |rigorous_upper_bound|:
%

xu = 1e5;
yu = 1e5 * [1 1 1 1]';

obj.solve('sedumi') ...
   .rigorous_lower_bound(xu) ...
   .rigorous_upper_bound(yu)

%%
% yielding also a reasonable bound.
%
