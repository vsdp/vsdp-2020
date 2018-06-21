# VSDP: Verified SemiDefinite-quadratic-linear Programming (Version 2018)

VSDP is a software package that is designed for the computation of verified
results in conic programming.  It supports the constraint cone consisting of the
product of semidefinite cones, second-order cones and the nonnegative orthant.
It provides functions for computing rigorous error bounds of the true optimal
value, verified enclosures of epsilon-optimal  solutions, and verified
certificates of infeasibility.  All rounding errors due to floating-point
arithmetic are taken into account.

The software is completely written in [MATLAB](https://www.mathworks.com) /
[GNU Octave](https://www.gnu.org/software/octave) and requires the interval
toolbox [INTLAB](http://www.ti3.tuhh.de/rump/intlab).  Thus interval input is
supported as well.

Via its interface, VSDP provides easy access to the conic solvers:
- [CSDP](https://projects.coin-or.org/Csdp),
- [GLPK](https://www.gnu.org/software/glpk/),
- [LINPROG](https://www.mathworks.com/help/optim/ug/linprog.html),
- [lp_solve](https://lpsolve.sourceforge.io),
- [SeDuMi](https://github.com/sqlp/sedumi),
- [SDPA](https://sdpa.sourceforge.io), and
- [SDPT3](https://github.com/sqlp/sdpt3).

For more information, please visit <https://vsdp.github.io/>.
