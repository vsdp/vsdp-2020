%% Semidefinite Programming
%
% Let the SDP program be given in the standard form \eqref{cpPrim}
% $$
% \label{BlockDiagPSDP}
% \begin{array}{lll}
% \text{minimize}
% & \sum_{j=1}^{n_{s}} \langle C_{j}^{s}, X_{j}^{s} \rangle, & \\
% \text{subject to}
% & \sum_{j=1}^{n_{s}} \langle A_{i,j}^{s}, X_{j}^{s} \rangle = b_{i},
% & i = 1,\ldots,m, \\
% & X_{j}^{s} \in \mathbb{S}^{s_{j}}_{+},
% & j = 1,\ldots,n_{s}.
% \end{array}
% $$
% Its dual problem \eqref{cpDual} is
% $$
% \label{BlockDiagDSDP}
% \begin{array}{ll}
% \text{maximize} & b^{T} y, \\
% \text{subject to}
% & Z_{j}^{s} = C_{j}^{s} - \sum_{i=1}^{m} y_{i} A_{i,j}^{s}
%   \in \mathbb{S}^{s_{j}}_{+},\quad j = 1, \ldots, n_{s}.
% \end{array}
% $$
% The matrices $A_{i,j}^{s}, C_{j}^{s}, X_{j}^{s}$ are assumed to be
% symmetric $s_{j} \times s_{j}$ matrices. We store this problem in our
% condensed format \eqref{vec}, \eqref{condensedX}, and \eqref{condensedA}.
%
% Let us consider the example
% (see </references#Borchers2009 [Borchers2009]>):
% $$
% \begin{array}{lll}
% \text{minimize}
% & \sum_{j=1}^{3} \langle C_{j}^{s}, X_{j}^{s} \rangle, & \\
% \text{subject to}
% & \sum_{j=1}^{3} \langle A_{i,j}^{s}, X_{j}^{s} \rangle = b_{i},\quad
%      i = 1,2, \\
% & X_{1}^{s} \in \mathbb{S}^{2}_{+}, \\
% & X_{2}^{s} \in \mathbb{S}^{3}_{+}, \\
% & X_{3}^{s} \in \mathbb{S}^{2}_{+},
% \end{array}
% $$
% where
% $$
% \begin{array}{ccc}
%   A^{s}_{1,1} = \begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix}
% & A^{s}_{1,2} =
%   \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}
% & A^{s}_{1,3} = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} \\
%   A^{s}_{2,1} = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}
% & A^{s}_{2,2} =
%   \begin{pmatrix} 3 & 0 & 1 \\ 0 & 4 & 0 \\ 1 & 0 & 5 \end{pmatrix}
% & A^{s}_{2,3} = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix} \\
%   C^{s}_{1} = \begin{pmatrix} -2 & -1 \\ -1 & -2 \end{pmatrix}
% & C^{s}_{2} =
%   \begin{pmatrix} -3 & 0 & -1 \\ 0 & -2 & 0 \\ -1 & 0 & -3 \end{pmatrix}
% & C^{s}_{3} = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix} \\
%   b = \begin{pmatrix} 1 \\ 2 \end{pmatrix} & &
% \end{array}
% $$
%
% In the condensed format the corresponding coefficient matrix and the primal
% objective are
% $$
% \begin{aligned}
% A &= \begin{pmatrix}
%      3 & 1 & 1 & 3 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
%      0 & 0 & 0 & 0 & 3 & 0 & 1 & 0 & 4 & 0 & 1 & 0 & 5 & 0 & 0 & 0 & 1
%      \end{pmatrix},\\
% c &= \begin{pmatrix}
%   -2 & -1 & -1 & -2 & -3 & 0 & -1 & 0 & -2 & 0 & -1 & 0 & -3 & 0 & 0 & 0 & 0
%   \end{pmatrix}^{T}.
% \end{aligned}
% $$
%
% We enter the problem data
%

clear A b c K

A = [3, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
     0, 0, 0, 0, 3, 0, 1, 0, 4, 0, 1, 0, 5, 0, 0, 0, 1];
c = [-2; -1; -1; -2; -3; 0; -1; 0; -2; 0; -1; 0; -3; 0; 0; 0; 0];
b = [1; 2];

%%
% define the structure |K| for the PSD-cone
%

K.s = [2; 3; 2];

%%
% and call |mysdps|
%

vsdpinit('sdpt3');
[objt,xt,yt,zt,info] = mysdps(A,b,c,K);
objt

%%
% The other quantities are not displayed for brevity.
%
% By calling |vsdplow| and |vsdpup| we get verified error bounds
%

[fL,y,dl] = vsdplow(A,b,c,K,xt,yt,zt)

%%
%

[fU,x,lb] = vsdpup(A,b,c,K,xt,yt,zt)

%%
% The components |lb(j)| are lower bounds of the smallest eigenvalues
% $\lambda_{\min}([X_{j}^{s}])$ for |j = 1,2,3|.  Hence |lb > 0| proves that
% all real matrices $X_{j}^{s}$, that are contained in the corresponding
% interval matrices $[X_{j}^{s}]$, are in the interior of the cone
% $\mathbb{S}^{2}_{+} \times \mathbb{S}^{3}_{+} \times\mathbb{S}^{2}_{+}$,
% where the interval matrices $[X_{j}^{s}]$ are obtained by applying the |mat|
% operator to intval |x|.  Analogously, |dl(j)| are lower bounds for the
% smallest eigenvalues of the dual interval matrices $[Z_{j}^{s}]$ that
% correspond to the dual solution |y|.  Since, for this example, both |dl|
% and |lb| are positive, strict feasibility is proved for the primal and the
% dual problem, and strong duality holds valid.
%
% Now we consider the following example
% (see </references#Jansson2006 [Jansson2006]>):
% $$
% \label{SDPexample}
% \begin{array}{ll}
% \text{minimize} & \langle C_{1}, X \rangle, \\
% \text{subject to}
% & \langle A_{1,1}, X \rangle = 1, \\
% & \langle A_{2,1}, X \rangle = 2??, \\
% & \langle A_{3,1}, X \rangle = 0, \\
% & \langle A_{4,1}, X \rangle = 0, \\
% & X \in \mathbb{S}^{3}_{+},
% \end{array}
% $$
% where
% $$
% \begin{array}{cc}
% C_{1} = \begin{pmatrix}
%         0   & 0.5 & 0 \\
%         0.5 & ??   & 0 \\
%         0   & 0   & ??
%         \end{pmatrix},
% & b = \begin{pmatrix} 1 \\ 2?? \\ 0 \\ 0 \end{pmatrix}, \\
% A_{1,1} = \begin{pmatrix}
%            0   & -0.5 & 0 \\
%           -0.5 &  0   & 0 \\
%            0   &  0   & 0
%           \end{pmatrix},
% & A_{2,1} = \begin{pmatrix}
%             1 & 0 & 0 \\
%             0 & 0 & 0 \\
%             0 & 0 & 0
%             \end{pmatrix}, \\
% A_{3,1} = \begin{pmatrix}
%           0 & 0 & 1 \\
%           0 & 0 & 0 \\
%           1 & 0 & 0
%           \end{pmatrix},
% & A_{4,1} = \begin{pmatrix}
%             0 & 0 & 0 \\
%             0 & 0 & 1 \\
%             0 & 1 & 0
%             \end{pmatrix}.
% \end{array}
% $$
%
% It is easy to prove that for
%
% * $?? > 0$: the problem is feasible with
%   $\hat{f}_{p} = \hat{f}_{d} = -0.5$ (zero duality gap),
% * $?? = 0$: the problem is feasible but ill-posed with nonzero duality
%   gap,
% * $?? < 0$: the problem is infeasible.
%
% For $?? > 0$ the primal optimal solution of the problem is given by the
% matrix
% $$
% \label{OptSolSDPExp}
% X^{*} = \begin{pmatrix}
% 2?? & -1 & 0 \\ -1 & \dfrac{1}{2??} & 0 \\ 0 & 0 & 0
% \end{pmatrix}.
% $$
% The corresponding dual optimal vector is
% $y^{*} = \begin{pmatrix} 0 & -1/(4??) & 0 & 0 \end{pmatrix}^{T}$.
% We choose $?? = 10^{-4}$ and enter the problem.
%

% define delta parameter
d = 1e-4;
% define the constraint matrices
A1 = [ 0,   -0.5, 0;
      -0.5,  0,   0;
       0,    0,   0];
A2 = [1, 0, 0;
      0, 0, 0;
      0, 0, 0];
A3 = [0, 0, 1;
      0, 0, 0;
      1, 0, 0];
A4 = [0, 0, 0;
      0, 0, 1;
      0, 1, 0];

b = [1; 2*d; 0; 0];

C = [0,   0.5, 0;
     0.5, d,   0;
     0,   0,   d];

% define the cone structure K
K.s = 3;

% vectorize the matrices Ai and C
A = [A1(:), A2(:), A3(:), A4(:)];
c = C(:);

%%
% The SDPT3-solver provides the following results:
%

vsdpinit('sdpt3');
[objt,xt,yt,zt,info] = mysdps(A,b,c,K);
objt, xt, yt, info  % zt:  hidden for brevity

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
