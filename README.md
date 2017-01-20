# VSDP: A Matlab toolbox for verified semidefinite-quadratic-linear programming (Version 2012)

VSDP is a software package that is designed for the computation of verified
results in conic programming.  The current version of VSDP supports the
constraint cone consisting of the product of semidefinite cones, second-order
cones and the nonnegative orthant.  It provides functions for computing
rigorous error bounds of the true optimal value, verified enclosures of
ε-optimal  solutions, and verified certificates of infeasibility.  All rounding
errors due to floating point arithmetic are taken into account.

VSDP is completely written in [MATLAB](https://www.mathworks.com) /
[GNU Octave](https://www.gnu.org/software/octave).  It uses
[INTLAB](http://www.ti3.tuhh.de/rump/intlab), and thus interval input data are
supported as well.  Via its interface, VSDP provides an easy access to the
conic solvers [CSDP](https://projects.coin-or.org/Csdp),
[SeDuMi](https://github.com/sqlp/sedumi), [SDPA](https://sdpa.sourceforge.io),
[SDPT3](https://github.com/sqlp/sdpt3), as well as
[lp_solve](https://lpsolve.sourceforge.io) and
[LINPROG](https://www.mathworks.com/help/optim/ug/linprog.html).
Detailed numerical results for the
[NETLIB LP library](http://www.netlib.org/lp),
[SDPLIB](http://euler.nmt.edu/~brian/sdplib/sdplib.html),
[DIMACS](http://dimacs.rutgers.edu/Challenges/Seventh/Instances/) and
[Kocvara's test library](http://plato.asu.edu/ftp/kocvara) for problems from
structural optimization are presented in the second part of this manual.
Many of these test problems are challenging real-world problems of large scale.


## Prerequisites

VSDP requires one of the aforementioned approximate solvers and the
[interval toolbox "INTLAB"](http://www.ti3.tuhh.de/rump/intlab) in version 9
or higher.


## Contributors

Viktor Härter, Marko Lange, and Christian Jansson (jansson@tuhh.de).


## References

- [Jansson2009] C. Jansson, On verified numerical computations in convex
  programming, Japan J. Indust. Appl. Math. (2009) 26: 337,
  https://dx.doi.org/10.1007/BF03186539

- [Jansson2005] C. Jansson, Termination and Verification for Ill-posed
  Semidefinite Programming Problems,
  http://www.optimization-online.org/DB_HTML/2005/06/1150.html

- [Ogita2005] T. Ogita, S.M. Rump, and S. Oishi, Accurate Sum and Dot Product,
  SIAM J. Sci. Comput. 26(6), 1955-1988,
  https://dx.doi.org/10.1137/030601818

- [Jansson2004] C. Jansson, Rigorous Lower and Upper Bounds in Linear
  Programming, SIAM J. OPTIM. Vol.14, No.3, pp. 914-935,
  https://dx.doi.org/10.1137/S1052623402416839
