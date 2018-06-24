function svec (obj, mu, isSymmetric)
% SVEC  Symmetric vectorization operator.
%
%   obj.svec(mu, isSymmetric)
%
%   For the trace inner product of symmetric matrices X and Y holds:
%
%     <X,Y> := trace(X*Y) = svec(X,K,sqrt(2))' * svec(Y,K,sqrt(2))
%                         = svec(X,K,1)'       * svec(Y,K,2).
%
%   Don't use mu = sqrt(2) for verified computations. Instead,
% to avoid rounding errors use: svec(X,K,1)'*svec(Y,K,2).
%
% For a symmetric matrix X let A=X(:). Then for K.s=n  (n=length(X))
%   svec(A,K,mu) = [x11 mu*x12 x22 mu*x13 mu*x23 x33 ... mu*x1n ... xnn]
%
% >> Input:
% A: nA x M matrix, alternatively
%     whereas nA = dimf+diml+dimq+dims
%     dimf: number of free variables: dimf = sum_i(K.f(i)>0)
%     diml: number of nonnegative variables: diml = sum_i(K.l(i)>0)
%     dimq: sum of all socp variables: dimq = sum_i(K.q(i))
%     dims: sum of all sdp variables: dims = sum_i(K.s(i)^2)
% K: a structure with following fields
%     - K.f stores the number of free variables
%     - K.l is the number of nonnegative components
%     - K.q lists the lengths of socp blocks
%     - K.s lists the dimensions of semidefinite blocks
% mu: scaling factor for off-diagonal elements
%     for performance reasons and reversibilty mu=0 is not allowed
% isSymmetric: how to treat coefficient matrices
%     - 0 -> matrices have no symmetry and both triangular parts has to be
%            processed  -  this is allowed for real matrices only
%     - 1 -> default symmetric case
% I: index cell for matrix conversion. Setting the index vector "I" avoids
%    creation of the index vector and saves some computation time.
%
% >> Output:
% vA: an M x nA3 matrix, or nA3 x M matrix (depends on input)
%     whereas nA3 = dimf+diml+dimq+dims3
%     dims3: sum of all sdp variables: dims3 = sum_i(K.s(i)*(K.s(i)+1)/2)
%

% check input parameter
if nargin<2 || ~isstruct(K)
  error('VSDP:VSVEC','not enough input parameters\n');
end
if nargin<3 || isempty(mu)
  mu = 1;
end
if nargin<4 || isempty(isSymmetric)
  isSymmetric = 1;  % assume symmetric case
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

isize = size(obj.A,1);
if isize>3 && all(isize~=[dim3 dim])
  error('VSDP:VSVEC','Cone dimension does not fit to size of input');
elseif isize<=3 || isize==dim3
  % nothing to do
  vA = obj.A;
  return;
end

% input index
if iscell(I) && length(I)==2 && length(I{1})==dim
  Ilow = I{2};  I = I{1};
else
  Ilow = [];  I = [];
end


% index creation  -  upper triangular part
ns = length(K.s);
if isempty(I)
  I = cell(ns+1,1);
  I{1} = true(nos,1);
  for k = 1:ns
    I{k+1} = reshape(triu(true(K.s(k))),[],1);
  end
  I = vertcat(I{:});
end


% index vector for lower triangular parts of sdp blocks
if isSymmetric==0 && isempty(Ilow)
  Ilow = cell(ns+1,1);
  Ilow{1} = (1:nos)';
  blks = nos + 1;
  for k = 2:ns+1
    nk = K.s(k-1);
    Ilow{k} = nk * ones(nk*(nk+1)/2,1);
    Ilow{k}(cumsum([1 1:nk-1])) = [blks 1:-nk:1-nk*(nk-2)];
    Ilow{k} = cumsum(Ilow{k});
    blks = blks + nk*nk;
  end
  Ilow = vertcat(Ilow{:});
end


% Remove all rows that are not indexed, regard `mu`.
if (isSymmetric)
  vA = 0.5 * (obj.A(I,:) + obj.A(Ilow,:));
else
  vA = obj.A(I,:);
end

% Scale if necessary.
if (mu ~= 1)
  vA = sscale (vA, K, mu);
end

% write index cell
obj.Ivec = {I, Ilow};
