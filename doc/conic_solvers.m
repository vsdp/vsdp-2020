## Conic solvers
#
# In this section we provide some notes about our usage experience with VSDP
# and the following approximate conic solvers.
#
# For our tests we tried the following four combinations to run the conic
# solver:
#
#   +-------------------------------------------------------+
#   |                      <conic solver>                   |
#   +-------------------------------------------------------+
#   |                      VSDP & INTLAB                    |
#   +------------+--------------+------------+--------------+
#   |   MATLAB   |  GNU Octave  |   MATLAB   |  GNU Octave  |
#   +------------+--------------+------------+--------------+
#   |       MS Windows          |        GNU/Linux          |
#   +---------------------------+---------------------------+
#
##

## CSDP
#
# * Website: <https://projects.coin-or.org/Csdp>
# * GitHub: <https://github.com/coin-or/Csdp>
# * Documentation: <https://github.com/coin-or/Csdp/blob/master/doc/csdpuser.pdf>
# * *Cones*: (Free variables), LP, SDP
# * *Installation*: Binary distributions for Windows and Linux on the website.
#   Extract binary distribution to arbitrary location and use |addpath| within
#   Octave or MATLAB to add the |bin| (solver executables) and |matlab|
#   (interface routines) subdirectories.
# * *Invocation*: Call |csdp| from the Octave or MATLAB command prompt.
# * *Notes*: Free variables are only supported as difference of LP variables.
#   The resulting problem is ill-posed.
#

## GLPK
#
# * Website: <https://www.gnu.org/software/glpk>
# * Documentation: Part of the source code archive available from the website.
# * *Cones*: Free variables, LP
# * *Installation*: Built-in solver of GNU Octave.
# * *Invocation*: Call |glpk| from the Octave command prompt.
# * *Notes*: *No MATLAB.*
#

## LINPROG
#
# * Website: <https://www.mathworks.com/help/optim/ug/linprog.html>
# * Documentation: See website.
# * *Cones*: Free variables, LP
# * *Installation*: Built-in solver of MATLAB.
# * *Invocation*: Call |linprog| from the MATLAB command prompt.
# * *Notes*: *No Octave.*
#

## lp_solve
#
# * Website: <http://lpsolve.sourceforge.net/5.5/index.htm>
# * Documentation: See website.
# * *Cones*: Free variables, LP
# * *Installation*: Binary distributions for Windows and Linux are available
#   from <https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.5>.
#   For Linux it is easier to install the binary solver files from the
#   respective distribution package manager.  The Octave and MATLAB interface
#   is available from <https://github.com/vsdp/lp_solve> or
#   <https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.5>.
#   Inside the |extra/octave/lpsolve| or |extra/matlab/lpsolve| subdirectory
#   one has to follow the build instructions and has to use |addpath| within
#   Octave or MATLAB to make it work.
# * *Invocation*: Call |lp_solve| from the Octave or MATLAB command prompt.
#

## MOSEK
#
# * Website: <https://www.mosek.com/>
# * Documentation: <https://www.mosek.com/documentation/>
# * *Cones*: Free variables, LP, SOCP, SDP
# * *Installation*: Binary distributions for Windows and Linux on the website.
#   Very good description is given at
#   <https://docs.mosek.com/8.1/install/installation.html>
# * *Invocation*: Call |mosekopt| from the Octave or MATLAB command prompt.
# * *Notes*: One has to obtain a *license*, which is gratis for personal and
#   academic use.  *No Octave.*  There is some no longer maintained
#   <https://github.com/MOSEK/octmosek octmosek> Octave package, that does notes
#   work for recent Octave versions.
#

## SDPA
#
# * Website: <http://sdpa.sourceforge.net>
# * GitHub: <https://github.com/vsdp/sdpa>
# * Documentation:
#   <https://sourceforge.net/projects/sdpa/files/sdpa/sdpa.7.1.1.manual.20080618.pdf>
# * *Cones*: LP, SDP
# * *Installation*: For *Windows* binary distributions can be obtained from
#   <http://sdpa.sourceforge.net/download.html> (use the *SDPA-M* version).
#   Extract binary distribution to arbitrary location and use |addpath| within
#   MATLAB to add the solver and interface routines.
# * *Invocation*: Call |mexsdpa| from the Octave or MATLAB command prompt.
# * *Notes*: The "native" installation for *Windows and MATLAB* is the most
#   reliable.  The other three combinations worked partially or not at all.
#   Some of the Linux efforts are reflected in the GitHub repository.
#

## SDPT3
#
# * Website: <http://www.math.nus.edu.sg/~mattohkc/SDPT3.html>
# * GitHub: <https://github.com/sqlp/sdpt3> (preferred) or
#   <https://github.com/Kim-ChuanToh/SDPT3>
# * Documentation: <http://www.math.nus.edu.sg/~mattohkc/sdpt3/guide4-0-draft.pdf>
# * *Cones*: Free variables, LP, SOCP, SDP
# * *Installation*: Download the files from the preferred GitHub repository
#   extracted to an arbitrary location and run |install_sdpt3| from the Octave
#   or MATLAB command prompt.
# * *Invocation*: Call |sqlp| from the Octave or MATLAB command prompt.
# * *Notes:* In case of errors with the MEX-Interface, run
#   |install_sdpt3 -rebuild| from the Octave or MATLAB command prompt.
#   The SDPT3 function |randmat| is in conflict with the INTLAB function
#   |randmat|.  To resolve that conflict, we chose to rename the INTLAB
#   function.
#

## SeDuMi
#
# * Website: <http://sedumi.ie.lehigh.edu>
# * GitHub: <https://github.com/sqlp/sedumi> (preferred)
# * Documentation:
#   <http://sedumi.ie.lehigh.edu/sedumi/files/sedumi-downloads/SeDuMi_Guide_11.pdf>
#   and <http://sedumi.ie.lehigh.edu/sedumi/files/sedumi-downloads/usrguide.ps>
# * *Cones*: Free variables, LP, SOCP, SDP
# * *Installation*: Download the files from the preferred GitHub repository
#   extracted to an arbitrary location and run |install_sedumi| from the Octave
#   or MATLAB command prompt.
# * *Invocation*: Call |sedumi| from the Octave or MATLAB command prompt.
# * *Notes:* In case of errors with the MEX-Interface, run
#   |install_sedumi -rebuild| from the Octave or MATLAB command prompt.
#
