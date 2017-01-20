function [A, I] = vsmat(vA,K,mu,sflag,I)
%% VSMAT:  smat operator for VSDP3
%    A = vsmat(vA,K,mu,sflag)
%    [A, I] = vsmat(vA,K,mu,sflag,I)
%
%% >> Description:
% For symmetric matrices X and Y it applies:
% <X,Y> := trace(X*Y) = svec(X,K,sqrt(2))'*svec(Y,K,sqrt(2)).
% The 'smat' operator is the converse operator of 'svec'.
% For vX = svec(X,mu=2) the vector mX = smat(vX,mu=0.5) equals the matrix X.
%
% For a symmetric matrix X let vA=X(:). Then for K.s=n  (n=length(X))
%   svec(vA,K,mu) = [x11 mu*x12 ... mu*x1n x22 mu*x23 ...mu*xn3 ... xnn]
%
% it is:    vA = smat(svec(vA,K,mu),K,1/mu)
%
%% >> Input:
% vA: an M x nA3 matrix, or nA3 x M matrix, alternatively
%     whereas nA3 = dimf+diml+dimq+dims3
%     dimf: number of free variables: dimf = sum_i(K.f(i)>0)
%     diml: number of nonnegative variables: diml = sum_i(K.l(i)>0)
%     dimq: sum of all socp variables: dimq = sum_i(K.q(i))
%     dims3: sum of all sdp variables: dims3 = sum_i(K.s(i)*(K.s(i)+1)/2)
% K: a structure with following fields
%     - K.f stores the number of free variables
%     - K.l is the number of nonnegative components
%     - K.q lists the lengths of socp blocks
%     - K.s lists the dimensions of semidefinite blocks
% mu: scaling factor for off-diagonal elements
%     for performance reasons and reversibilty mu=0 is not allowed
% sflag: symmetry flag: how to treat coefficient matrices
%     - 0 -> do not care about lower triangular part of matrix
%            the corresponding elements are filled with zeros
%     - 1 -> output will be symmetric
%     - default: 1
% I: index vector for matrix conversion. Setting the index vector "I" avoids
%    creation of the index vector and saves some computation time. If
%    "sflag" is changed the index vector has to be computed again.
%
%% >> Output:
% A: a M x nA matrix, or nA x M matrix, depending on the input
%    whereas nA = dimf+diml+dimq+dims
%    dims: sum of all sdp variables: dims = sum_i(K.s(i)^2)
% I: index vector for applied conversion. Can be used to speed up
%    following conversions of similar cone variables by smat.
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
% 09/07/12    M. Lange, rewrite for faster indexing and non-symmetric input
% 28/07/12    M. Lange, speed improvements + adaption for new format
%

%% check input parameter
if nargin<2
  error('VSDP:VSMAT','Not enough input parameter');
end
if nargin<3 || isempty(mu)
  mu = 1;
end
if nargin<4 || isempty(sflag)
  sflag = 1;  % assume symmetric case
end
if nargin<5
  I = [];
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
  K.s = K.s(K.s>0);
  dim = nos + sum(K.s.*K.s);
  dim3 = nos + sum(K.s.*(K.s+1))/2;
else
  dim = nos;
  dim3 = nos;
end


% column or row-wise
[isize,idim] = max(size(vA));
if isize>3 && all(isize~=[dim3 dim])
  error('VSDP:VSMAT','Cone dimension does not fit to size of input');
elseif isize<=3 || isize==dim
  % nothing to do
  A = vA;
  return;
end


%% mat-transformation index vector
ns = length(K.s);
if length(I)~=dim || (sflag==1 && ~isnumeric(I)) || (sflag~=1 && ~islogical(I))
  I = cell(ns+1,1);  % cell to hold index vectors for every sdp block
  if sflag==1  % symmetric output
    blks = nos + 1;
    I{1} = (1:nos)';
    for k = 2:ns+1
      nk = K.s(k-1);
      I{k} = tril(repmat((-1:nk-2)',1,nk),-1)+1;
      I{k}(1,:) = cumsum([blks 1:nk-1]);
      I{k} = reshape(cumsum(I{k},1),[],1);
      blks = blks + nk*(nk+1)/2;
    end
  else  % only upper triangular output
    I{1} = true(nos,1);
    for k = 1:ns
      I{k+1} = reshape(triu(true(K.s(k))),[],1);
    end
  end
  % create index vector from cell
  I = vertcat(I{:});
end


%% transform matrix
if mu~=1
  vA = sscale(vA,K,mu);
end

switch 2*(sflag==1)+idim
  case 1  % sflag==0 && idim==1
    A(I,:) = vA;
  case 2  % sflag==0 && idim==2
    A(:,I) = vA;
  case 3  % sflag==1  && idim==1
    A = vA';       % faster than A = vA(I,:)
    A = A(:,I)';
  otherwise  % sflag==1 && idim==2
    A = vA(:,I);
end

end
