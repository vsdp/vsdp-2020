# VSDP: A Matlab toolbox for verified semidenite-quadratic-linear programming (Version 2012)

VSDP is a software package that is designed for the computation of verified
results in conic programming.  The current version of VSDP supports the
constraint cone consisting of the product of semidenite cones, second-order
cones and the nonnegative orthant.  It provides functions for computing
rigorous error bounds of the true optimal value, verified enclosures of
ε-optimal  solutions, and verified certicates of infeasibility.  All rounding
errors due to floating point arithmetic are taken into account.

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

VSDP requires one of the abovementioned approximate solvers and the
[interval toolbox "INTLAB"](http://www.ti3.tuhh.de/rump/intlab) in version 9
or higher.
