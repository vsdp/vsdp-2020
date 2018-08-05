%% Conic Programming
%
%%

%% Primal and dual form
%
% VSDP can handle three self dual convex cones $\mathcal{K}$, that often occur
% in practical applications.  These are:
%
% * The non-negative orthant:
%   $$\mathbb{R}^{n}_{+} := \{ x \in \mathbb{R}^{n} \colon\; x_{i} \geq 0,
%   \; i = 1, \ldots, n \}.$$
% * The Lorentz cone (see
%   <https://vsdp.github.io/references.html#Alizadeh2003 [Alizadeh2003]>):
%   $$\mathbb{L}^{n} := \{ x \in \mathbb{R}^{n} \colon x_{1}
%   \geq \|x_{2:n}\|_{2}\}.$$
% * The cone of symmetric positive semidefinite matrices:
%   $$\mathbb{S}^{n}_{+} := \left\{ X \in \mathbb{R}^{n \times n} \colon\;
%   X = X^{T},\; v^{T} X v \geq 0,\; \forall v \in \mathbb{R}^{n} \right\}.$$
%
% If a quantity is in the interior of one of the above cones, the definitions
% above must hold with strict inequalities.
%
% By $\langle c, x \rangle := c^{T} x$ the usual Euclidean inner product of
% vectors in $\mathbb{R}^{n}$ is denoted.  For symmetric matrices
% $X, Y \in \mathbb{S}^{n}$ the inner product is given by
% $\langle X,Y \rangle := \text{trace}(XY)$.
%
% Let $A^{f}$ and $A^{l}$ be a $m \times n_{f}$ and a $m \times n_{l}$ matrix,
% respectively, and let $A_{i}^{q}$ be $m \times q_{i}$ matrices for
% $i = 1,\ldots,n_{q}$.  Let $c^{f} \in \mathbb{R}^{n_{f}}$,
% $c^{l} \in \mathbb{R}^{n_{l}}$, $c_{i}^{q} \in \mathbb{R}^{q_i}$, and
% $b \in \mathbb{R}^{m}$.  Moreover, let $A_{1,j}^{s}, \ldots, A_{m,j}^{s}$,
% $C_{j}^{s}$, $X_{j}^{s}$ be symmetric $(s_{j} \times s_{j})$ matrices for
% $j = 1, \ldots, n_{s}$.
%
% Now we can define the conic semidefinite-quadratic-linear programming
% problem in primal standard form:
% $$\begin{array}{ll}
% \text{minimize} &
% \langle c^{f}, x^{f} \rangle + \langle c^{l}, x^{l} \rangle +
% \sum_{i=1}^{n_{q}} \langle c_{i}^{q}, x_{i}^{q} \rangle +
% \sum_{j=1}^{n_{s}} \langle C_{j}^{s}, X_{j}^{s} \rangle \\
% \text{subject to} &
% A^{f} x^{f} + A^{l} x^{l} + \sum_{i=1}^{n_{q}} A_{i}^{q} x_{i}^{q} +
% \sum_{j=1}^{n_{s}}\mathcal{A}_{j}^{s}(X_{j}^{s}) = b,
% \end{array}$$
% where $x^{f} \in \mathbb{R}^{n_{f}}$ are "free variables",
% $x^{l} \in \mathbb{R}^{n_{l}}_{+}$ are "non-negative variables",
% $x_{i}^{q} \in \mathbb{L}^{q_i}$, $i = 1, \ldots, n_{q}$, are "second-order
% cone (SOCP) variables", and finally $X_{j}^{s} \in \mathbb{S}^{s_{j}}_{+}$,
% $j = 1, \ldots, n_{s}$ "positive semidefinite (SDP) variables".  The linear
% operator
% $$\mathcal{A}_{j}^{s}(X_{j}^{s}) :=
% \begin{pmatrix}
% \langle A_{1j}^{s}, X_{j}^{s} \rangle \\
% \vdots \\
% \langle A_{mj}^{s}, X_{j}^{s} \rangle
% \end{pmatrix}$$
% maps the symmetric matrices $X_{j}^{s}$ to $\mathbb{R}^{m}$.  The adjoint
% linear operator is
% $$(\mathcal{A}_{j}^{s})^{*} y := \sum_{k=1}^{m} A_{kj}^{s} y_{k}.$$
%
% The dual problem associated with the primal standard form is
% $$\begin{array}{ll}
% \text{maximize} & b^{T} y \\
% \text{subject to}
% & (A^{f})^{T} y + z^{f} = c^{f}, \\
% & (A^{l})^{T} y + z^{l} = c^{l}, \\
% & (A_{i}^{q})^{T} y + z_{i}^{q} = c_{i}^{q}, \\
% & (\mathcal{A}_{j}^{s})^{*} y + Z_{j}^{s} = C_{j}^{s},
% \end{array}$$
% where $z^{f} \in \{0\}^{n_{f}}$, $z^{l} \in \mathbb{R}^{n_{l}}_{+}$,
% $z_{i}^{q} \in \mathbb{L}^{q_i}$, $i = 1, \ldots, n_{q}$, and
% $Z_{j}^{s} \in \mathbb{S}^{s_{j}}_{+}$, $j = 1, \ldots, n_{s}$.
%
% The objective functions and equality constraints of the primal and dual
% problem are linear.  Thus conic programming can be seen an extension of linear
% programming with additional conic constraints.
%
% By definition the vector $x^{f}$ contains all unconstrained or free
% variables, whereas all other variables are bounded by conic constraints.
% In several applications some solvers (for example SDPA or CSDP) require that
% free variables are converted into the difference of nonnegative variables.
% Besides the major disadvantage that this transformation is numerical
% unstable, it also increases the number of variables of the particular
% problems.  In VSDP <https://vsdp.github.io/free_variables free variables>
% can be handled in a numerical stable manner.
%


%% Condensed form
%
% Occasionally, it is useful to represent the conic programming problem in a
% more compact form by using the symmetric vectorization operator.  This
% operator maps a symmetric matrix $X \in \mathbb{S}^{n \times n}$ to a
% $n(n + 1)/2$-dimensional vector
% $$svec(X, \alpha) :=
% \begin{pmatrix}
% X_{11} & \alpha X_{12} & X_{22} & \alpha X_{13} & \cdots &  X_{nn}
% \end{pmatrix}^{T},$$
% where $\alpha$ is a scaling factor for the off diagonal elements.  The
% inverse operation is denoted by $smat(x)$ such that $smat(svec(X)) = X$.
%
% By using $svec$ it is possible to map each symmetric matrix to a vector
% quantity.  Additionally, by using the scaling factor $\alpha = 1$ for all
% $A_{kj}^{s}$, $C_{j}^{s}$, and $Z_{j}^{s}$ and the scaling factor
% $\alpha = 2$ for all $X_{j}^{s}$, the inner product of symmetric matrices
% reduces to a simple scalar product of two vectors.  Thus all cones can be
% treated in exactly the same way.
%
% The condensed quantities $c$, $x$, and $z$ are $n \times 1$-vectors:
% $$c :=
% \begin{pmatrix}
% c^{f} \\ c^{l} \\ c_{1}^{q} \\ \vdots \\ c_{n_{q}}^{q} \\
% svec(C_{1}^{s},1) \\ \vdots \\ svec(C_{n_{s},1}^{s}) \\
% \end{pmatrix},
% x :=
% \begin{pmatrix}
% x^{f} \\ x^{l} \\ x_{1}^{q} \\ \vdots \\ x_{n_{q}}^{q} \\
% svec(X_{1}^{s},2) \\ \vdots \\ svec(X_{n_{s},2}^{s})
% \end{pmatrix},
% z :=
% \begin{pmatrix}
% z^{f} \\ z^{l} \\ z_{1}^{q} \\ \vdots \\ z_{n_{q}}^{q} \\
% svec(Z_{1}^{s},1) \\ \vdots \\ svec(Z_{n_{s}}^{s},1) \\
% \end{pmatrix},$$
% where
% $n = n_{f} + n_{l} + \sum_{i = 1}^{n_{q}} q_{i} + \sum_{j = 1}^{n_{s}} s_{i}(s_{i} + 1)/2$
% and $A^{T}$ becomes a $n \times m$ matrix
% $$A^{T} :=
% \begin{pmatrix}
% & A^{f} & \\
% & A^{l} & \\
% & A_{1}^{q} & \\
% & \vdots & \\
% & A_{n_{q}}^{q} & \\
% svec(A_{11}^{s},1) & \cdots & svec(A_{1m}^{s},1) \\
% \vdots & & \vdots \\
% svec(A_{n_{s}1}^{s},1) & \cdots & svec(A_{n_{s}m}^{s},1)
% \end{pmatrix}$$
%
% Let the constraint cone $K$ and its dual cone $K^{*}$ be
% $$\begin{align}
% \mathcal{K} &:=&
% \mathbb{R}^{n_{f}} &\times
% \mathbb{R}^{n_{l}}_{+} \times
% \mathbb{L}^{q_{1}} \times \ldots \times \mathbb{L}^{q_{n_{q}}} \times
% \mathbb{S}^{s_{1}}_{+} \times \ldots \times \mathbb{S}^{s_{n_{s}}}_{+}, \\
% \mathcal{K}^{*} &:=&
% \{0\}^{n_{f}} &\times
% \mathbb{R}^{n_{l}}_{+} \times
% \mathbb{L}^{q_{1}} \times \ldots \times \mathbb{L}^{q_{n_{q}}} \times
% \mathbb{S}^{s_{1}}_{+} \times \ldots \times \mathbb{S}^{s_{n_{s}}}_{+}.
% \end{align}$$
%
% With these abbreviations we obtain the following block form of the conic
% problem:
% $$\begin{array}{ll}
% \text{minimize}   & c^{T} x, \\
% \text{subject to} & Ax = b, \\
%                   & x \in \mathcal{K},
% \end{array}$$
% with optimal value $\hat{f_{p}}$ and the corresponding dual problem
% $$\begin{array}{ll}
% \text{maximize}   & b^{T} y, \\
% \text{subject to} & z = c - (A)^{T} y \in \mathcal{K}^{*},
% \end{array}$$
% with optimal value $\hat{f_{d}}$.  In VSDP each conic problem is fully
% described by the four variables |(A, b, c, K)|.  The first two quantities
% represent the affine constraints $Ax = b$.  The third is the primal objective
% vector |c|, and the last describes the underlying cone.  The cone |K| is a
% structure with four fields: |K.f|, |K.l|, |K.q|, and |K.s|.  The field |K.f|
% stores the number of free variables $n_{f}$, the field |K.l| stores the
% number of nonnegative variables $n_{l}$, the field |K.q| stores the
% dimensions $q_{1}, \ldots, q_{n_{q}}$ of the second order cones, and
% similarly |K.s| stores the dimensions $s_{1}, \ldots, s_{n_{s}}$ of the
% semidefinite cones.  If a component of |K| is empty, then it is assumed that
% the corresponding cone do not occur.
%
% It is well known that for linear programming problems strong duality
% $\hat{f_{p}} = \hat{f_{d}}$ holds without any constraint qualifications.
% General conic programs satisfy only the weak duality condition
% $\hat{f_{d}} \leq \hat{f_{p}}$.  Strong duality requires additional
% constraint qualifications, such as _Slater's constraint qualifications_ (see
% <https://vsdp.github.io/references.html#Vandenberghe1996 [Vandenberghe1996]>,
% <https://vsdp.github.io/references.html#BenTal2001 [BenTal2001]>).
%
% *Strong Duality Theorem*
%
% * If the primal problem is strictly feasible (i.e. there exists a primal
%   feasible point $x$ in the interior of $K$) and $\hat{f_{p}}$ is finite,
%   then $\hat{f_{p}} = \hat{f_{d}}$ and the dual supremum is attained.
% * If the dual problem is strictly feasible (i.e. there exists some $y$ such
%   that $z = c - (A)^{T} y$ is in the interior of $K^{*}$) and $\hat{f_{d}}$
%   is finite, then $\hat{f_{d}} = \hat{f_{p}}$, and the primal infimum is
%   attained.
%
% In general, the primal or dual problem formulation may have optimal solutions
% while its respective dual problem is infeasible, or the duality gap may be
% positive at optimality.
%
% Duality theory is central to the study of optimization.  Firstly, algorithms
% are frequently based on duality (like primal-dual interior point methods),
% secondly, they enable one to check whether or not a given feasible point is
% optimal, and thirdly, it allows one to compute verified results efficiently.
%

%% Interval arithmetic
%
% For the usage of VSDP a knowledge of interval arithmetic is not required.
% Intervals are only used to specify error bounds.  An interval vector or an
% interval matrix is defined as a set of vectors or matrices that vary between
% a lower and an upper vector or matrix, respectively.  In other words, these
% are quantities with interval components.  In
% <http://www.ti3.tu-harburg.de/rump/intlab/ INTLAB> these interval quantities
% can be initialized with the routine |infsup|.  Equivalently, these quantities
% can be defined by a midpoint-radius representation, using the routine
% |midrad|.
%