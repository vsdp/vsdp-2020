%% Installation
%
%%

%% Requirements
% To run VSDP, the following requirements have to be fulfilled:
%
% * A recent version of <http://www.octave.org GNU Octave> or
%   <http://www.mathworks.com/products/matlab MATLAB> has to be installed.
% * The interval toolbox <http://www.ti3.tu-harburg.de/rump/intlab INTLAB> is
%   required.
% * At least one of the following approximate solvers has to be installed:
%   <https://github.com/coin-or/Csdp CSDP>,
%   <https://www.gnu.org/software/glpk GLPK>,
%   <https://www.mathworks.com/help/optim/ug/linprog.html LINPROG>,
%   <https://lpsolve.sourceforge.io lp_solve>,
%   <https://www.mosek.com MOSEK>,
%   <https://sdpa.sourceforge.io SDPA>,
%   <https://github.com/sqlp/sdpt3 SDPT3>, or
%   <https://github.com/sqlp/sedumi SeDuMi>.
%

%% Obtaining VSDP
%
% *ZIP-File*
%
% The most recent version of VSDP and this manual are available at
% <https://vsdp.github.io>.  There you can download a ZIP-file
% |vsdp-2018-master.zip| and extract it to an arbitrary location.
%
% Legacy versions of VSDP are available from
% <http://www.ti3.tu-harburg.de/jansson/vsdp/>.
%
% *Using git*
%
% If you have <https://git-scm.com/ git> installed and about 700 MB of disk
% space available, you can easily obtain a full bundle of VSDP 2006, 2012, 2018,
% including some aforementioned approximate solvers, and some benchmark
% libraries by the command
%
%    git clone --recurse-submodules https://github.com/vsdp/vsdp.github.io
%
% In the cloned directory |vsdp.github.io/vsdp/2018| you find the latest version
% of VSDP.
%


%% Installing VSDP
%
% If all requirements are fulfilled, just call from the MATLAB or GNU Octave
% command prompt inside the VSDP directory
%
%   install_vsdp
%
% and all necessary paths are set and VSDP is fully functional.  To test the
% latter, you can run the small builtin test suite from MATLAB via
%
%   runtests ('testVSDP')
%
%  Totals:
%     5 Passed, 0 Failed, 0 Incomplete.
%     8.2712 seconds testing time.
%
% or from GNU Octave via
%
%   testVSDP
%
%  Test summary
%  ------------
%
%  testSINDEX      PASSED  0.036834
%  testSVEC_SMAT   PASSED  0.277806
%  testLP          PASSED  0.923324
%  testSOCP        PASSED  0.918687
%  testSDP         PASSED  1.242036
%
