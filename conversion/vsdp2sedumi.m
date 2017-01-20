function [A,c,x,z] = vsdp2sedumi(A,c,x,z,K,opts)
%% VSDP2SEDUMI:  transforms problem data from VSDP (intern) format to SeDuMi format 
%    [A c x z] = vsdp2sedumi(A,c,x,z,K,opts)
%
%% >> Input:
% A: an M x nA3 matrix, or nA3 x M matrix in VSDP format 
%     whereas nA3 = dimf+diml+dimq+dims3 
%     dimf: number of free variables: dimf = sum(K.f>0)
%     diml: number of nonnegative variables: diml = sum(K.l>0)
%     dimq: sum of all socp variables: dimq = sum_i(K.q(i))
%     dims3: sum of sdp variables: dims3 = sum_i(K.s(i)*(K.s(i)+1)/2)
% c,x,z: nA3 x 1 vectors, whereas x has been created by svec(mu=2)
% K: a structure with following fields
%     - K.f stores the number of free variables
%     - K.l is the number of nonnegative components
%     - K.q lists the lengths of socp blocks
%     - K.s lists the dimensions of semidefinite blocks
% opts: structure for additional parameter settings:
%    regarded fields:
%     - 'ALLOW_TRIANGULAR'  Some solvers like SeDuMi allow non-symmetric
%                           input. If 'ALLOW_TRIANGULAR' is set true, VSDP
%                           is increasing the speed of data conversion by
%                           processing only a triangular part of the
%                           matrices.
%                               -> default: false
%
%% >> Output: 
% A: linear constraints matrix (nA x M or M x nA) in SeDuMi format
% c,x,z: nA x 1 vectors in SeDuMi format. The primal ojective vector may be
%        created such that lower triangulars of sdp parts are zero.
%
%%
% Note that the right hand side of the linear constraints (b) and the
% vector for primal objective function (c) as well as the cone structure
% (K) have the same format in VSDP and SeDuMi. Thus only At has to be
% converted.
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
% 31/07/10    V. Haerter, comments added
% 11/07/12    M. Lange, added opts parameter and adapted for new functions
% 27/07/12    M. Lange, adapted for new format and modified smat
%


%% input parameter
if nargin<5 || ~isstruct(K)
    error('VSDP:VSDP2SEDUMI','not enough input parameters\n');
elseif nargin<6
    opts = [];
end

global VSDP_OPTIONS;

sflag = 1;  % symmetry flag  -> default: symmetric output
if isfield(opts,'ALLOW_TRIANGULAR')
    sflag = (opts.ALLOW_TRIANGULAR==false);
elseif isfield(VSDP_OPTIONS,'ALLOW_TRIANGULAR')
    sflag = (VSDP_OPTIONS.ALLOW_TRIANGULAR==false);
end

% get problem data dimensions
fields = isfield(K,{'f','l','q','s'});
dim = 0;
if fields(1)
    dim = sum(K.f);
end
if fields(2)
    dim = dim + sum(K.l);
end
if fields(3)
    dim = dim + sum(K.q);
end
if fields(4)
    K.s = K.s(K.s>0);
    dim3 = dim + sum(K.s.*(K.s+1))/2;
    dim = dim + sum(K.s.*K.s);
else
    dim3 = dim;
end

% check dimension of input
adim = [0 dim3 dim];  % allowed dimensions
if all(length(A)~=adim) || all(length(c)~=adim) || all(length(x)~=adim) ...
        || all(length(z)~=adim)
    error('VSDP:VSDP2SEDUMI','cone structure does not match to matrix dimension');
end


%% transform data
% use mu=2 for triangular output else mu=1
mu = 2-sflag;
Imat = [];
adim = [0 dim];
if all(length(A)~=adim)
    [A Imat] = vsmat(A,K,mu,sflag);
end
if all(length(c)~=adim)
    [c Imat] = vsmat(c,K,mu,sflag,Imat);
end
if all(length(x)~=adim)
    [x Imat] = vsmat(x,K,0.5,1,Imat);  % x has been created by svec(mu=2)
end
if all(length(z)~=adim)
    z = vsmat(z,K,1,1,Imat);
end

%_____________________________End VSDP2SEDUMI___________________________