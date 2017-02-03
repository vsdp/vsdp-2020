function opts = vsdpinit(opts,display)
% VSDPINIT  Initialization and defaults for VSDP.
%
%   opts = VSDPINIT() returns the persistently stored option structure for VSDP.
%      It contains following fields:
%
%               'SOLVER'    Select one of the supported solvers:
%                           'sedumi', 'sdpt3', 'sdpa', 'csdp', 'sdplr',
%                           'lp_solve', 'linprog', or 'glpk'.
%                           Default is 'sedumi'.
%
%   'USE_STARTING_POINT'    Decides whether initial starting point shall be
%                           used or not.  Default is false.
%
%             'ITER_MAX'    Maximum number of iterations that can be used to
%                           perturbate the problem and find a feasible solution.
%                           Default is 10.
%
%                'ALPHA'    Growing factor for problem perturbation.  Default
%                           is 0.5.
%
%  'FULL_EIGS_ENCLOSURE'    If true code for stronger complete eigenvalue
%                           enclosure will be applied.  Default is true.
%
%      'VERIFY_FULL_LSS'    If true the full non-symmetric matrix lss enclosure
%                           is applied within vuls-function.  Default is false.
%
%     'ALLOW_TRIANGULAR'    Some solvers like SeDuMi allow non-symmetric input.
%                           If true VSDP is increasing the speed of data
%                           conversion by processing only a triangular part of
%                           the SDP block matrices.  Default is false.
%
%      'MIN_SDPBLK_SIZE'    Minimum size of an SDP block, if block size is
%                           smaller than 'MIN_SDPBLK_SIZE' the blocks will be
%                           grouped when transforming into SDPT3 format.
%                           Default is 2500.
%
%       'SOLVER_OPTIONS'    Option structure for approximate solver that is
%                           empty by default.
%
%       'VERBOSE_OUTPUT'    If false, VSDP is trying to minimize verbose
%                           output of the supported solvers.  Otherwise the
%                           default solver messages are printed. Default is
%                           false.
%
%      The persistent settings can be cleared by calling 'clear vsdpinit'.
%
%   VSDPINIT('solver') optionally, pass the name of the solver to use.  See
%      'SOLVER' for supported solvers.
%
%   VSDPINIT(opts) optionally, update the default structure, by passing a
%      structure containing the same fields as the output structure.  Fields
%      not listed above are ignored.
%
%   VSDPINIT(opts,display) optionally, if 'display == true' the original VSDP
%      splash is printed to the command window.
%
%   Example:
%
%       vsdpinit('sdpt3')
%
%       opts.SOLVER = 'sdpt3'; % Same as above.
%       vsdpinit(opts)
%
%   See also mysdps.

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

% Quick return for calls like opts = vsdpinit()
persistent VSDP_OPTIONS;
% First time initialization?
if (~isstruct(VSDP_OPTIONS))
  VSDP_OPTIONS = getDefaultOptions();
  % Check whether path already exists, add path if necessary
  vsdp_path = fullfile(fileparts(which('vsdpinit')));
  if isempty(strfind(path,vsdp_path))
    addpath(vsdp_path);
    addpath(fullfile(vsdp_path, 'conversion'));
    if (exist('OCTAVE_VERSION','builtin'))
      addpath(fullfile(vsdp_path, 'octave'));
    end
    path(path); % refresh path
  end
  % Check for default solver
  defaultSolver = 'sedumi';
  defaultSolverAvailable = (exist('sedumi','file') == 2);
  if (~defaultSolverAvailable)
    warning('VSDP:VSDPINIT', ['Default solver ''%s'' is not available.\n', ...
      '\t Please choose another supported and available solver.\n'...
      '\t Type ''help vsdpinit'' for a list of supported solvers.'], ...
      defaultSolver);
  end
  % Check for INTLAB
  if (~(exist('INTLAB_Version_9.m', 'file') == 2))
    error('VSDP:VSDPINIT','%s.  %s\n\n\t%s\n\n%s.\n\n', ...
      'Interval toolbox "INTLAB" (>= 9) not found', ...
      'Get a recent version from', 'http://www.ti3.tuhh.de/rump/intlab', ...
      'and run "startintlab" from the root directory');
  end
end
if (nargin == 0)
  opts = VSDP_OPTIONS;
  return;
end

% update options using input opts
narginchk(1,2);
if (nargin < 2)
  display = false;
end

% solver
% for calls like vsdpinit('solver')
if (ischar(opts))
  opts = struct('SOLVER',opts);
end
if (isfield(opts,'SOLVER'))
  defaultSolver = 'sedumi';
  defaultSolverAvailable = (exist('sedumi','file') == 2);
  setSolver = false;
  switch (lower(opts.SOLVER))
    case 'sedumi'
      setSolver = defaultSolverAvailable;
    case 'sdpt3'
      setSolver = (exist('sqlp','file') == 2);
    case 'sdpa'
      setSolver = ((exist('sdpam','file') == 2) ...
        && ((exist('mexsdpa','file') == 3) ...
        || (exist('callSDPA','file') == 2)));
    case 'csdp'
      setSolver = (exist('csdp','file') == 2);
    case 'sdplr'
      setSolver = (exist('sdplr','file') == 2);
    case 'lp_solve'
      setSolver = (exist('lp_solve','file') == 2);
    case 'linprog'
      setSolver = (exist('linprog','file') == 2);
    case 'glpk'
      setSolver = (exist('glpk', 'file') == 2);
    otherwise
      warning('VSDP:VSDPINIT', 'Solver ''%s'' is not supported.', opts.SOLVER);
  end
  if (setSolver)
    VSDP_OPTIONS.SOLVER = lower(opts.SOLVER);
  elseif (~strcmp(VSDP_OPTIONS.SOLVER,defaultSolver))
    % Else if another solver was successully set
    warning('VSDP:VSDPINIT', ...
      'Solver ''%s'' is not available.  Keep solver: ''%s''.', opts.SOLVER, ...
      VSDP_OPTIONS.SOLVER);
  elseif (defaultSolverAvailable)
    warning('VSDP:VSDPINIT', ...
      'Solver ''%s'' is not available.  Use default solver: ''%s''.', ...
      opts.SOLVER, defaultSolver);
    VSDP_OPTIONS.SOLVER = 'sedumi';
  else
    error('VSDP:VSDPINIT', ...
      'Solver ''%s'' and the default solver ''%s'' are not available.', ...
      opts.SOLVER, defaultSolver);
  end
end

% use starting point
if (isfield(opts,'USE_STARTING_POINT'))
  try
    validateattributes(opts.USE_STARTING_POINT,{'logical'},{'scalar'});
    VSDP_OPTIONS.USE_STARTING_POINT = opts.USE_STARTING_POINT;
  catch
    error('VSDP:VSDPINIT', 'opts.USE_STARTING_POINT must be a logical scalar.');
  end
end

% maximum iteration in vsdplow and vsdpup
if (isfield(opts,'ITER_MAX'))
  try
    validateattributes(opts.ITER_MAX,{'numeric'},{'positive','scalar'});
    VSDP_OPTIONS.ITER_MAX = opts.ITER_MAX;
  catch
    error('VSDP:VSDPINIT', 'opts.ITER_MAX must be a numeric positive scalar.');
  end
end

% growing factor for perturbations
if (isfield(opts,'ALPHA'))
  try
    validateattributes(opts.ALPHA,{'numeric'},{'positive','scalar'});
    VSDP_OPTIONS.ALPHA = opts.ALPHA;
  catch
    error('VSDP:VSDPINIT', 'opts.ALPHA must be a numeric positive scalar.');
  end
end

% full eigenvalue enclosure
if (isfield(opts,'FULL_EIGS_ENCLOSURE'))
  try
    validateattributes(opts.FULL_EIGS_ENCLOSURE,{'logical'},{'scalar'});
    VSDP_OPTIONS.FULL_EIGS_ENCLOSURE = opts.FULL_EIGS_ENCLOSURE;
  catch
    error('VSDP:VSDPINIT', 'opts.FULL_EIGS_ENCLOSURE must be a logical scalar');
  end
end

% full lss verification
if (isfield(opts,'VERIFY_FULL_LSS'))
  try
    validateattributes(opts.VERIFY_FULL_LSS,{'logical'},{'scalar'});
    VSDP_OPTIONS.VERIFY_FULL_LSS = opts.VERIFY_FULL_LSS;
  catch
    error('VSDP:VSDPINIT', 'opts.VERIFY_FULL_LSS must be a logical scalar');
  end
end

% allow non-symmetric, upper triangular format conversion
if (isfield(opts,'ALLOW_TRIANGULAR'))
  try
    validateattributes(opts.ALLOW_TRIANGULAR,{'logical'},{'scalar'});
    VSDP_OPTIONS.ALLOW_TRIANGULAR = opts.ALLOW_TRIANGULAR;
  catch
    error('VSDP:VSDPINIT', 'opts.ALLOW_TRIANGULAR must be a logical scalar');
  end
end


% blocksize for sdp block groups
if (isfield(opts,'MIN_SDPBLK_SIZE'))
  try
    validateattributes(opts.MIN_SDPBLK_SIZE,{'numeric'},{'positive','scalar'});
    VSDP_OPTIONS.MIN_SDPBLK_SIZE = opts.MIN_SDPBLK_SIZE;
  catch
    error('VSDP:VSDPINIT', ...
      'opts.MIN_SDPBLK_SIZE must be a numeric positive scalar.');
  end
end

% solver options unchecked
if (isfield(opts,'SOLVER_OPTIONS'))
  VSDP_OPTIONS.SOLVER_OPTIONS = opts.SOLVER_OPTIONS;
end

% verbose output
if (isfield(opts,'VERBOSE_OUTPUT'))
  try
    validateattributes(opts.VERBOSE_OUTPUT,{'logical'},{'scalar'});
    VSDP_OPTIONS.VERBOSE_OUTPUT = opts.VERBOSE_OUTPUT;
  catch
    error('VSDP:VSDPINIT', 'opts.VERBOSE_OUTPUT must be a logical scalar');
  end
end

opts = VSDP_OPTIONS;

% message of successful start
if (logical(display))
  fprintf('\n___________________________________________________________________\n');
  fprintf('********************** Initialization of VSDP *********************\n');
  fprintf('*                                                                 *\n');
  fprintf('*     VSDP 2012:     Verified Conic Programming                   *\n');
  fprintf('*             Viktor Haerter, viktor.haerter@tu-harburg.de        *\n');
  fprintf('*             Christian Jansson, jansson@tu-harburg.de            *\n');
  fprintf('*             Marko Lange, marko.lange@tu-harburg.de              *\n');
  fprintf('*             Institute for Reliable Computing                    *\n');
  fprintf('*             Hamburg University of Technology,                   *\n');
  fprintf('*             Germany,  October 2012                              *\n');
  fprintf('*                                                                 *\n');
  fprintf('******************* VSDP is now ready for use *********************');
  fprintf('\n*******************************************************************\n\n');
end

end

function [opts] = getDefaultOptions()
opts.SOLVER = 'sedumi';
opts.USE_STARTING_POINT = false;
opts.ITER_MAX = 10;
opts.ALPHA = 0.5;
opts.FULL_EIGS_ENCLOSURE = true;
opts.VERIFY_FULL_LSS = false;
opts.MIN_SDPBLK_SIZE = 2500;
opts.ALLOW_TRIANGULAR = false;
opts.SOLVER_OPTIONS = [];
opts.VERBOSE_OUTPUT = false;
end
