function smat (obj, mu, isSymmetric)
% SMAT  Inverse operator of symmetric vectorization (svec).
%
%   obj.smat(mu, isSymmetric)
%
% Description:
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
% Input:
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
% Output:
% A: a M x nA matrix, or nA x M matrix, depending on the input
%    whereas nA = dimf+diml+dimq+dims
%    dims: sum of all sdp variables: dims = sum_i(K.s(i)^2)
% I: index vector for applied conversion. Can be used to speed up
%    following conversions of similar cone variables by smat.
%

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

narginchk(2, 3);
if ((nargin < 2) || isempty(mu))
  mu = 1;
end
if ((nargin < 3)|| isempty(isSymmetric))
  isSymmetric = true;
end

% Get problem data dimensions:
% a) The number of variables that are not in SDP-cone.
nos = sum(obj.K.f) + sum(obj.K.l) + sum(obj.K.q);
% b) The number of variables in SDP-cones.
dim = nos + sum(obj.K.s .* obj.K.s);
% c) The number of condensed variables in SDP-cones.
dim3 = nos + sum(obj.K.s .* (obj.K.s + 1)) / 2;

if (size(obj.At, 1) == dim)
  % nothing to do
  return;
end


% mat-transformation index vector
ns = length(K.s);
if length(I)~=dim || (isSymmetric==1 && ~isnumeric(I)) || (isSymmetric~=1 && ~islogical(I))
  I = cell(ns+1,1);  % cell to hold index vectors for every sdp block
  if isSymmetric==1  % symmetric output
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


% transform matrix
if mu~=1
  vA = sscale(vA,K,mu);
end

if (isSymmetric)
  A = vA'; % TODO: really faster than A = vA(I,:)?
  A = A(:,I)';
else
  A(I,:) = vA;
end

end
