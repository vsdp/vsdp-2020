function [stat opts] = vsdpinit(opts,display)
%% VSDPINIT - Initialization and defaults for VSDP3
%    [stat] = vsdpinit(par,val)
%
%% >> Description:
%     VSDPINIT can be used to initialize VSDP or setup option settings.
%
%     vsdpinit('clear') to clear all settings
%     vsdpinit(opts)  updates global setting structure with opts structur
%     vsdpinit('sdpt3')  same as  opts.SOLVER = 'sdpt3'; vsdpinit(opts)
%
%% >> Input:
%
%     - opts: option structure with the following fields
%
%               'SOLVER'    select one of the supported solvers:
%                           'sedumi','sdpt3','sdpa','csdp','sdplr', ...
%                           'lp_solve','linprog'
%
%   'USE_STARTING_POINT'    decides whether initial starting point shall be
%                           used or not
%
%             'ITER_MAX'    maximum number of iterations that can be used
%                           to perturbate the problem and find a feasible
%                           solution
%
%                'ALPHA'    growing factor for problem perturbation
%
%  'FULL_EIGS_ENCLOSURE'    if true code for stronger complete eigenvalue
%                           enclosure will be applied
%
%      'VERIFY_FULL_LSS'    if true the full non-symmetric matrix lss
%                           enclosure is applied within vuls-function
%
%     'ALLOW_TRIANGULAR'    some solvers like SeDuMi allow non-symmetric
%                           input, if 'ALLOW_TRIANGULAR' is set true VSDP
%                           is increasing the speed of data conversion by
%                           processing only a triangular part of the sdp
%                           block matrices
%
%        'SDPT3_VERSION'    version number of SDPT3 solver
%
%      'MIN_SDPBLK_SIZE'    minimum size of an sdp block, if block size is
%                           smaller the blocks will be grouped when
%                           transforming into SDPT3 format
%
%     '$SOLVER$_OPTIONS'    option structure for approximate solver
%                           the field for sdpt3 options is 'SDPT3_OPTIONS'
%
%	- display:  -1 - display only errors
%                0 - display warnings
%                1 - display warnings and additional information
%
%% >> Output:
%   - stat: number of option structure fields which have been detected in
%           opts
%
%   - opts: copy of global option structure
%

%% ********************************************************************* %%
%% This file is part of VSDP by V. Haerter, C. Jansson and M. Lange      %%
%% Copyright (c) 2012, C. Jansson                                        %%
%%                     Technical University of Hamburg (TUHH)            %%
%%                     Institute for Reliable Computing (IRC)            %%
%% VSDP can be freely used for private and academic purposes.            %%
%% Commercial use or use in conjunction with a commercial program which  %%
%% requires VSDP or part of it to function properly is prohibited.       %%
%% ********************************************************************* %%

%% Last modified:  
% 31/08/10    V. Haerter, comments added
% 09/09/10    V. Haerter, store Matlab version
% 09/09/12    M. Lange, rewrite for new options structure
%
%% ToDo:
%

%% check input

% global option settings
global VSDP_OPTIONS;

% check input parameter
if nargin<1
    opts = [];  display = 0;
elseif nargin<2 || isempty(display)
    display = 0;
end

% vsdpinit('solver') or vsdpinit('clear') call
if ischar(opts)
    if strcmpi(opts,'clear')
        clear GLOBAL VSDP_OPTIONS;
        stat = 0;
        opts = [];
        return;
    end
    opts = struct('SOLVER',opts);
end

wng = warning;  % old warning status
warning off;

% get current vsdp path
vsdppath = which('vsdpinit');
vsdppath(end-9:end) = [];

% check whether path already exists, add path if necessary
if isempty(strfind(path,vsdppath))
    % add vsdp to paths
    addpath(vsdppath);
    addpath([vsdppath 'conversion']);
    addpath([vsdppath 'examples']);
    % refresh path
    path(path);
end

% display warnings
if display<0
    warning off;
else
    warning on;
end


%% update options
stat = 0;  % number of read options fields
wfields = 0;  % number of fields that could not be imported

% solver
solverLIST = {'sedumi','sdpt3','sdpa','csdp','sdplr', ...
    'lp_solve','linprog'};
if isfield(opts,'SOLVER')
    stat = stat + 1;
    ofield = opts.SOLVER;
    if isnumeric(ofield) && ofield>0 && ofield<=length(solverLIST)
        VSDP_OPTIONS.SOLVER = solverLIST(ofield);
    elseif isempty(ofield) || any(strcmp(ofield,solverLIST))
        VSDP_OPTIONS.SOLVER = ofield;
    else
        warning('VSDP:VSDPINIT','selected solver is not supported');
        wfields = wfields + 1;
    end
elseif ~isfield(VSDP_OPTIONS,'SOLVER')
    warning('VSDP:VSDPINIT','no solver has been selected');
end

% use starting point
if isfield(opts,'USE_STARTING_POINT')
    stat = stat + 1;
    ofield = opts.USE_STARTING_POINT;
    if length(ofield)==1 && any(ofield==[0 1])
        VSDP_OPTIONS.USE_STARTING_POINT = ofield;
    else
        warning('VSDP:VSDPINIT','USE_STARTING_POINT has to be logical');
        wfields = wfields + 1;
    end
end

% maximum iteration in vsdplow and vsdpup
if isfield(opts,'ITER_MAX')
    stat = stat + 1;
    ofield = opts.MAX_ITER;
    if isnumeric(ofield) && length(ofield)==1 && ofield>=1
        VSDP_OPTIONS.MAX_ITER = ofield;
    else
        warning('VSDP:VSDPINIT','MAX_ITER not accepted');
        wfields = wfields + 1;
    end
end

% growing factor for perturbations
if isfield(opts,'ALPHA')
    stat = stat + 1;
    ofield = opts.ALPHA;
    if isnumeric(ofield) && length(ofield)==1 && ofield~=0
        VSDP_OPTIONS.ALPHA = ofield;
    else
        warning('VSDP:VSDPINIT','ALPHA has wrong format or is zero');
        wfields = wfields + 1;
    end
end

% full eigenvalue enclosure
if isfield(opts,'FULL_EIGS_ENCLOSURE')
    stat = stat + 1;
    ofield = opts.FULL_EIGS_ENCLOSURE;
    if length(ofield)==1 && any(ofield==[0 1])
        VSDP_OPTIONS.FULL_EIGS_ENCLOSURE = ofield;
    else
        warning('VSDP:VSDPINIT','FULL_EIGS_ENCLOSURE has to be logical');
        wfields = wfields + 1;
    end
end

% full lss verification
if isfield(opts,'VERIFY_FULL_LSS')
    stat = stat + 1;
    ofield = opts.VERIFY_FULL_LSS;
    if length(ofield)==1 && any(ofield==[0 1])
        VSDP_OPTIONS.VERIFY_FULL_LSS = ofield;
    else
        warning('VSDP:VSDPINIT','VERIFY_FULL_LSS has to be logical');
        wfields = wfields + 1;
    end
end

% allow non-symmetric, upper triangular format conversion
if isfield(opts,'ALLOW_TRIANGULAR')
    stat = stat + 1;
    ofield = opts.ALLOW_TRIANGULAR;
    if length(ofield)==1 && any(ofield==[0 1])
        VSDP_OPTIONS.ALLOW_TRIANGULAR = ofield;
    else
        warning('VSDP:VSDPINIT','ALLOW_TRIANGULAR has to be logical');
        wfields = wfields + 1;
    end
end

% SDPT3 version number
if isfield(opts,'SDPT3_VERSION')
    stat = stat + 1;
    ofield = opts.SDPT3_VERSION;
    if isnumeric(ofield) && length(ofield)==1
        VSDP_OPTIONS.SDPT3_VERSION = ofield;
    else
        warning('VSDP:VSDPINIT','wrong SDPT3_VERSION format');
        wfields = wfields + 1;
    end
end

% blocksize for sdp block groups
if isfield(opts,'MIN_SDPBLK_SIZE')
    stat = stat + 1;
    ofield = opts.MIN_SDPBLK_SIZE;
    if isnumeric(ofield) && length(ofield)==1
        VSDP_OPTIONS.MIN_SDPBLK_SIZE = ofield;
    else
        warning('VSDP:VSDPINIT','wrong MIN_SDPBLK_SIZE format');
        wfields = wfields + 1;
    end
end

% solver options  -  not checked
soLIST = cellfun(@(x) [upper(x) '_OPTIONS'],solverLIST,'UniformOutput',false);
for i = 1:length(soLIST)
    if isfield(opts,soLIST{i})
        stat = stat + 1;
        eval(['VSDP_OPTIONS.' soLIST{i} ' = opts.' soLIST{i} ';']);
    end
end

% output options + statistics
if display>=0 && ~isempty(opts)
    fnum = numel(structfun(@(x) 1,opts));
    if stat<fnum
        warning(['%d field(s) of "opts" has/have not been ' ...
            'detected as vsdp option'],fnum-stat);
    end
end
opts = VSDP_OPTIONS;


%% message of successful start
if display>0
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

% reset warning state
warning(wng);

%______________________________END OF VSDPINIT____________________________
