%% Conic Programming
%
%%

%% Primal and dual form
%
% Let $\mathbb{R}^{n}_{+}$ denote the nonnegative orthant, and let
% $\mathbb{L}^{n} := \{x \in \mathbb{R}^{n} \colon x_{1} \geq ||x_{2:n}||_{2}\}$
% be the Lorentz cone.  We denote by $\langle c, x \rangle := c^{T} x$ the
% usual Euclidean inner product of vectors in $\mathbb{R}^{n}$.  The set
% $$
% \mathbb{S}^{n}_{+} := \left\{ X \in \mathbb{R}^{n \times n} \colon
% X = X^{T}, v^{T} X v \geq 0, \forall v \in \mathbb{R}^{n} \right\},
% $$
% denotes the cone of symmetric positive semidefinite $n \times n$ matrices.
% For symmetric matrices $X$, $Y$ the inner product is given by
% $$
% \langle X,Y \rangle := \text{trace}(XY).
% $$
%
% Let $A^{f}$ and $A^{l}$ be a $m \times n_{f}$ and a $m \times n_{l}$ matrix,
% respectively, and let $A_{i}^{q}$ be $m \times q_{i}$ matrices for
% $i = 1,\ldots,n_{q}$.  Let $x^{f} \in \mathbb{R}^{n_{f}}$,
% $x^{l} \in \mathbb{R}^{n_{l}}$, $x_{i}^{q} \in \mathbb{R}^{q_i}$, and
% $b \in \mathbb{R}^{m}$.  Moreover, let $A_{1,j}^{s}, \ldots, A_{m,j}^{s}$,
% $C_{j}^{s}$, $X_{j}^{s}$ be symmetric $(s_{j} \times s_{j})$ matrices for
% $j = 1, \ldots, n_{s}$.
%
% Now we can define the conic semidefinite-quadratic-linear programming
% problem in primal standard form:
% $$
% \begin{array}{ll}
% \text{minimize} &
% \langle c^{f}, x^{f} \rangle + \langle c^{l}, x^{l} \rangle +
% \sum_{i=1}^{n_{q}} \langle c_{i}^{q}, x_{i}^{q} \rangle +
% \sum_{j=1}^{n_{s}} \langle C_{j}^{s}, X_{j}^{s} \rangle \\
% \text{subject to} &
% A^{f} x^{f} + A^{l} x^{l} + \sum_{i=1}^{n_{q}} A_{i}^{q} x_{i}^{q} +
% \sum_{j=1}^{n_{s}}\mathcal{A}_{j}^{s}(X_{j}^{s}) = b,
% \end{array}
% $$
% where $x^{f} \in \mathbb{R}^{n_{f}}$ are "free variables",
% $x^{l} \in \mathbb{R}^{n_{l}}_{+}$ are "non-negative variables",
% $x_{i}^{q} \in \mathbb{L}^{q_i}$, $i = 1, \ldots, n_{q}$, are "second-order
% cone (SOCP) variables", and finally $X_{j}^{s} \in \mathbb{S}^{s_{j}}_{+}$,
% $j = 1, \ldots, n_{s}$ "positive semidefinite (SDP) variables".  The linear
% operator
% $$
% \mathcal{A}_{j}^{s}(X_{j}^{s}) :=
% \begin{pmatrix}
% \langle A_{1j}^{s}, X_{j}^{s} \rangle \\
% \vdots \\
% \langle A_{mj}^{s}, X_{j}^{s} \rangle
% \end{pmatrix}
% $$
% maps the symmetric matrices $X_{j}^{s}$ to $\mathbb{R}^{m}$.  The adjoint
% linear operator is
% $$
% (\mathcal{A}_{j}^{s})^{*} y := \sum_{k=1}^{m} A_{kj}^{s} y_{k}.
% $$
%
% The dual problem associated with the primal standard form is
% $$
% \begin{array}{ll}
% \text{maximize} & b^{T} y \\
% \text{subject to}
% & (A^{f})^{T} y + z^{f} = c^{f}, \\
% & (A^{l})^{T} y + z^{l} = c^{l}, \\
% & (A_{i}^{q})^{T} y + z_{i}^{q} = c_{i}^{q}, \\
% & (\mathcal{A}_{j}^{s})^{*} y + Z_{j}^{s} = C_{j}^{s},
% \end{array}
% $$
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
% problems.  In VSDP </free_variables free variables> can be handled in a
% numerical stable manner.
%


%% Condensed form
%
% Occasionally, it is useful to represent the conic programming problem in a
% more compact form by using the symmetric vectorization operator.  This
% operator maps a symmetric matrix $X \in \mathbb{S}^{n \times n}$ to a
% $n(n + 1)/2$-dimensional vector
% $$
% svec(X, \alpha) :=
% \begin{pmatrix}
% X_{11} & \alpha X_{12} & X_{22} & \alpha X_{13} & \cdots &  X_{nn}
% \end{pmatrix}^{T},
% $$
% where $\alpha$ is a scaling factor for the off diagonal elements.  The
% inverse operation is denoted by $smat(x)$ such that $smat(svec(X)) = X$.
%
% By using $svec$ it is possible to map each symmetric matrix to a vector
% quantity and by using the scaling factor $\alpha = 1$ for all $A_{kj}^{s}$,
% $C_{j}^{s}$, and $Z_{j}^{s}$ and the scaling factor $\alpha = 2$ for all
% $X_{j}^{s}$, the inner product of symmetrix matrices reduces to a simple scalar product of
% two vectors.
We condense the above quantities as follows:
% $$
% \begin{equation}
% \begin{aligned}
% x^{s} &:= (\operatorname{vec}(X_{1}^{s}); \ldots;
%            \operatorname{vec}(X_{n_{s}}^{s})), &
% c^{s} &:= (\operatorname{vec}(C_{1}^{s}); \ldots;
%            \operatorname{vec}(C_{n_{s}}^{s})), \\
% x^{q} &:= (x_{1}^{q}; \ldots; x_{n_{q}}^{q}), &
% c^{q} &:= (c_{1}^{q}; \ldots; c_{n_{q}}^{q}), \\
% x     &:= (x^{f}; x^{l}; x^{q}; x^{s}), &
% c     &:= (c^{f}; c^{l}; c^{q}; c^{s}).
% \end{aligned}
% \label{condensedX}
% \end{equation}
% $$
%
% Here $c^{q}$ and $x^{q}$ consist of $\bar{q} = \sum_{i=1}^{n_{q}} q_i$
% components, and $c^{s}$, $x^{s}$ are vectors of length
% $\bar{s} = \sum_{j=1}^{n_{s}} s_{j}^{2}$.  The total length of $x$ and $c$
% is equal to $n = n_{f} + n_{l} + \bar{q} + \bar{s}$.  As in the syntax of
% Matlab the separating column denotes the vertical concatenation of the
% corresponding vectors.
%
% The matrices describing the linear equations are condensed as follows:
% $$
% \begin{equation}
% \begin{aligned}
% A^{s} &= (A_{1}^{s}, \ldots, A_{n_{s}}^{s}),
%          \text{ where } A_{j}^{s} =
%          (\text{vec}(A_{1,j}^{s}), \ldots, \text{vec}(A_{m,j}^{s}))^{T}, \\
% A^{q} &= (A_{1}^{q}, \ldots, A_{n_{q}}^{q}), \\
% A     &= (A^{f}, A^{l}, A^{q}, A^{s}).
% \end{aligned}
% \label{condensedA}
% \end{equation}
% $$
%
% Let the constraint cone $K$ and its dual cone $K^{*}$ be
% $$
% \begin{equation}
% \begin{aligned}
% K &:=&
% \mathbb{R}^{n_{f}} &\times
% \mathbb{R}^{n_{l}}_{+} \times
% \mathbb{L}^{q_{1}} \times \ldots \times \mathbb{L}^{q_{n_{q}}} \times
% \mathbb{S}^{s_{1}}_{+} \times \ldots \times \mathbb{S}^{s_{n_{s}}}_{+}, \\
% K^{*} &:=&
% \{0\}^{n_{f}} &\times
% \mathbb{R}^{n_{l}}_{+} \times
% \mathbb{L}^{q_{1}} \times \ldots \times \mathbb{L}^{q_{n_{q}}} \times
% \mathbb{S}^{s_{1}}_{+} \times \ldots \times \mathbb{S}^{s_{n_{s}}}_{+}.
% \end{aligned}
% \label{primalDualCone}
% \end{equation}
% $$
%
% With these abbreviations we obtain the following block form of the conic
% problem \eqref{stdPrim}:
% $$
% \label{cpPrim}
% \begin{array}{ll}
% \text{minimize}   & c^{T} x, \\
% \text{subject to} & Ax = b, \\
%                   & x \in K,
% \end{array}
% $$
% with optimal value $\hat{f}_{p}$ and the corresponding dual problem
% $$
% \label{cpDual}
% \begin{array}{ll}
% \text{maximize}   & b^{T} y, \\
% \text{subject to} & z = c - (A)^{T} y \in K^{*},
% \end{array}
% $$
% with optimal value $\hat{f}_{d}$.
%
% For a linear programming problem a vector $x \in \mathbb{R}^{n_{l}}$ is in
% the interior of the cone $K = \mathbb{R}^{n_{l}}_{+}$, if $x_{i} > 0$ for
% $i = 1,\ldots,n_{l}$.  For a vector $x \in \mathbb{R}^{n}$ let
% $\lambda_{\min}(x) := x_{1} - ||x_{:}||_{2}$ denote the smallest eigenvalue
% of $x$ (see </references#Alizadeh2003 [Alizadeh2003]>).  Then
% for second order cone programming problems a vector
% $x \in \mathbb{R}^{\bar{q}}$ is in the interior of the cone
% $K = \mathbb{L}^{q_{1}} \times, \ldots, \times \mathbb{L}^{q_{n_{q}}}$,
% if $\lambda_{\min}(x_{i}) > 0$ for $i = 1,\ldots,n_{q}$.  Furthermore, for
% a symmetric matrix $X \in \mathbb{S}^{n}$ let $\lambda_{\min}(X)$ denote the
% smallest eigenvalue of $X$.  Then for semidefinite programming problems a
% symmetric block matrix
% $$
% X = \begin{pmatrix}
% X_{1} & 0 & 0 \\
% 0 & \ddots & 0 \\
% 0 & 0 & X_{n_{s}}
% \end{pmatrix},
% $$
% is in the interior of the cone
% $K = \mathbb{S}^{s_{1}}_{+} \times,\ldots,\times \mathbb{S}^{s_{n_{s}}}_{+}$,
% if $\lambda_{\min}(X_{j}) > 0$ for $j = 1,\ldots,n_{s}$.
%
% It is well known that for linear programming problems strong duality
% $\hat{f}_{p} = \hat{f}_{d}$ holds without any constraint qualifications.
% General conic programs satisfy only the weak duality condition
% $\hat{f}_{d} \leq \hat{f}_{p}$.  Strong duality requires additional
% constraint qualifications, such as _Slater's constraint qualifications_ (see
% </references#Boyd1996 [Boyd1996]>,
% </references#NestNem [NestNem]>).
%
% *Strong Duality Theorem*
%
% * If the primal problem is strictly feasible (i.e. there exists a primal
%   feasible point $x$ in the interior of $K$) and $\hat{f}_{p}$ is finite,
%   then $\hat{f}_{p} = \hat{f}_{d}$ and the dual supremum is attained.
% * If the dual problem is strictly feasible (i.e. there exists some $y$ such
%   that $z = c - (A)^{T} y$ is in the interior of $K^{*}$) and $\hat{f}_{d}$
%   is finite, then $\hat{f}_{d} = \hat{f}_{p}$, and the primal infimum is
%   attained.
%
% In general, one of the problems \eqref{cpPrim} or \eqref{cpDual} may have
% optimal solutions and its dual problem is infeasible, or the duality gap may
% be positive at optimality.
%
% Duality theory is central to the study of optimization.  Firstly, algorithms
% are frequently based on duality (like primal-dual interior point methods),
% secondly, they enable one to check whether or not a given feasible point is
% optimal, and thirdly, it allows one to compute verified results efficiently.
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