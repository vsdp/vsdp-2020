function vA = sscale(vA,K,mu)
%% SSCALE:  scale operator to scale the off-diagonal elements of sdp blocks
%    vA = vscale(vA,K,mu)
%
%% >> Description:
% For a symmetric matrix X it applies:
%   vx = vscale(X(:),K,mu)  =>  vx = vX(:), where vX = mu*X + (1-mu)*I
% The same holds valid for matrices in the SDPT3 format. (-> see svec)
%
%% >> Input:
% vA: an nA x M, M x nA, nA3 x M, or M x nA3 matrix, alternatively,
%     whereas nA = dimf+diml+dimq+dims,  nA3 = dimf+diml+dimq+dims3
%     dimf: number of free variables: dimf = sum_i(K.f(i)>0)
%     diml: number of nonnegative variables: diml = sum_i(K.l(i)>0)
%     dimq: sum of all socp variables: dimq = sum_i(K.q(i))
%     dims: sum of all sdp variables: dims = sum_i(K.s(i)^2)
%     dims3: sum of sdp variables: dims = sum_i(K.s(i)*(K.s(i)+1)/2)
% K: a structure with following fields
%     - K.f stores the number of free variables
%     - K.l is the number of nonnegative components
%     - K.q lists the lengths of socp blocks
%     - K.s lists the dimensions of semidefinite blocks
% mu: scaling factor for off-diagonal elements
%     for performance reasons and reversibilty mu=0 is not allowed
%
%% >> Output:
% vA: matrix of same dimension as input vA
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
% 10/07/12    M. Lange, rewrite for faster indexing and both formats
% 28/07/12    M. Lange, speed improvements + allow scaling along rows
%

%% check input parameter
if nargin~=3 || ~isstruct(K)
  error('VSDP:VSCALE','all input parameters have to be set\n');
end

% get problem data dimensions
nos = 0;  % number of variables that are not in SDP-cone
fields = isfield(K,{'f','l','q','s'});
if fields(1)
  nos = sum(K.f);
end
if fields(2)
  nos = nos + sum(K.l);
end
if fields(3)
  nos = nos + sum(K.q);
end
if fields(4)
  K.s = reshape(K.s(K.s>0),[],1);
  dim = nos + sum(K.s.*K.s);
  dim3 = nos + sum(K.s.*(K.s+1))/2;
else
  dim = nos;
  dim3 = nos;
end

% svec or vec formated, column or row-wise
[isize idim] = max(size(vA));
if all(isize~=[0 dim3 dim])
  error('VSDP:SSCALE','cone dimension does not fit to size of input');
end

% anything to do ?
if isize==0 || isempty(mu) || mu==1 || isempty(K.s)
  % warning('VSDP:SSCALE','nothing to do - no sdp block found or mu=1\n');
  return;
elseif mu==0
  error('VSDP:SSCALE','Input mu=0 is not allowed');
end


%% create index vector for diagonal entries
if isize==dim  % full matrix format
  I = zeros(sum(K.s),1);
  I(cumsum([2; K.s(1:end-1)])) = K.s;
  I(1) = 1;
else  % reduced vectorized format
  I = ones(sum(K.s),1);
end
I(cumsum(K.s(1:end-1))+1) = (isize==dim3) - K.s(1:end-1);
I = cumsum(I);
I(1) = nos + 1;
I = cumsum(I);


%% scale off-diagonal elements
transInt = false;  % transposed interval input
if isa(vA,'intval') && idim==1
  vA = vA';
  idim = 2;
  transInt = true;
end
switch 2*(nos~=0) + (idim==2)
  case 0  % nos==0 && idim==1
    vA = mu*vA + bsxfun(@times,vA,sparse(I,1,1-mu,isize,1));
  case 1  % nos==0 && idim==2
    tmp = vA(:,I);
    vA = mu*vA;
    vA(:,I) = tmp;
  case 2  % nos>0  && idim==1
    vA = [vA(1:nos,:); mu*vA(nos+1:end,:)] + ...
      bsxfun(@times,vA,sparse(I,1,1-mu,isize,1));
  otherwise  % nos>0 && idim==2
    tmp = vA(:,I);
    vA = [vA(:,1:nos) mu*vA(:,nos+1:end)];
    vA(:,I) = tmp;
end
if transInt
  vA = vA';
end

end
