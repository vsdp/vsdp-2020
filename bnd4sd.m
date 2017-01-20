function [lambdaMin,trN,dshift] = bnd4sd(A,sym,mode)
% BND4SD  Calculates a lower bound for the smallest eigenvalue and the
%         sum over all negative eigenvalues of A.
%
%           The matrix A is assumed to be nearly positive semidefinite.
%           BND4SD can be used to find rigorous bounds for non-positive
%           definiteness due to rounding errors. The function is designed
%           for the use in VSDP.
%
%           A may also be an interval matrix, full or sparse. This interval
%           matrix can be defined as:
%               1. a structure with the fields 'mid' and 'rad',
%               2. a cell array, where the first entry will be interpreted
%                  as 'mid' and the second as 'rad',
%               3. a class of type 'intval' [see INTLAB toolbox].
%
%   @params
%       A:      Real symmetric input matrix.
%
%       sym:    The function requires a symmetric input matrix.
%               In case that A is not symmetric (for memory resons etc.)
%               set 'sym=false'. Then S = 0.5*(A+A') will be used instead.
%                       -> default: true
%
%       mode:   inclusion mode:
%                   0 - fast verification by applying cholesky
%                       factorization, if cholesky fails, lower bounds is
%                       -Inf
%                   1 - if verification by cholesky decomposition failed
%                       function will call 'vsdpneig' for a stronger enclosure
%                       of all eigenvalues
%                   default: 1
%
%   @output
%       lambdaN:    A lower bound for the minimum eigenvalue of A.
%
%       trN:    A lower bound for the sum of all negative eigenvalues of A.
%
%       dshift:     Shift of diagonal elements. Possible permutation of the
%                   matrix: If positive semidefiniteness could not be
%                   verified, it may be possible for the slightly
%                   perturbated matrix A+diag(dshift).
%
%

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

% check input
if nargin<1 || isempty(A)
  error('BND4SD: No input matrix set.');
elseif ~(isnumeric(A) || isa(A,'intval') || all(isfield(A,{'mid','rad'})))
  error('BND4SD: Datatype of input matrix is not accepted.');
end   % size & dimension is not checked. Error will occur in the code.


% prepare input data

% get rounding mode and set round to nearest
rnd = getround();
setround(0);    % for preparation only

% prepare interval input
if ~isnumeric(A)  % A is an interval structure or of type 'intval'
  Arad = A.rad;
  A = A.mid;
else
  Arad = sparse(size(A,1),size(A,2));
end


% prepare A
if ~isreal(A)
  error('BND4SD: Complex input not yet supported')
elseif size(A,1)~=size(A,2) && size(A,2)~=1
  A = A(:);
  Arad = Arad(:);
end
n = size(A,1);

if size(A,2)==1  % compact vectorized input
  n = sqrt(2*n+0.25) - 0.5;
  if n~=round(n)
    error('BND4SD: Wrong dimension of vectorized input');
  end
  E(triu(true(n))) = full(A);  % E is place-holder, will be overwritten later
  A = reshape(E,n,n);
  A = A + triu(A,1)';
  if nnz(A)<max(n^1.5,n*n/20)
    A = sparse(A);
  end
  if any(Arad)
    L(triu(true(n))) = full(Arad);  % L is used as place-holder again
    Arad = reshape(L,n,n);
    Arad = Arad + triu(Arad,1)';
    if nnz(Arad)<n*n/5
      Arad = sparse(Arad);
    end
  else
    Arad = sparse(size(A,1),size(A,2));
  end
elseif nargin>=2 && ~isempty(sym) && ~sym  % non-symmetric input
  setround(1);
  E = 0.5 * A;  % E and L are used as place-holder
  L = E';
  A = E + L;  % sup(0.5*(A+A'))
  Arad = 0.5 * (Arad+Arad') + (A + ((-E) - L));  % Arad + (sup(A)-inf(A))
  setround(0);
end


% calculate expected error for precondition shift
E = abs(A);
d = sqrt(diag(E));
setround(-1);
L = E~=d*d';
setround(1);
L = L | E~=d*d';
E = E.*L;
dshift = 8*eps*sqrt(sum(max(E))*sum(E)) + sum(Arad) + realmin;
setround(0);


% remove zero rows/columns
index = any(A-diag(sparse(diag(A)))) | any(Arad - diag(sparse(diag(Arad))));
lambdaMin = inf;
trN = 0;

if ~all(index)
  setround(-1);
  lambdaMin = min(diag(A)-diag(Arad));
  trN = sum(min(diag(A)-diag(Arad),0));
  setround(0);
  A = A(index,index);
  Arad = Arad(index,index);
  n = length(A);
end

if isempty(A) % A was pure diagonal matrix
  dshift = min(lambdaMin/mean(dshift),0) * dshift;
  setround(rnd);
  return;
end


% approximative cholesky factorization with adaptive error shift
[L,p] = chol(A-diag(sparse(dshift(index))),'lower');
if p~=0  % not semidefinite
  [L,p] = chol(A-diag(sparse(0.4*dshift(index))),'lower');
end
if p~=0
  [L,p] = chol(A,'lower');
end
if p==0 && ~isempty(find(isnan(L),1))
  p = n;
end


% cholesky decomposition did not work
if nargin==2 && mode==0 && p~=0
  lambdaMin = -Inf;
  trN = -Inf;
  return;
elseif p~=0  % call vsdpneig
  A = struct('mid',A,'rad',Arad);
  clear Arad D E L;  % not needed anymore
  lambda = vsdpneig(A,1);  % full enclosure
  setround(-1);
  lambda = lambda.mid - lambda.rad;  % inf(lambda)
  lambdaMin = min(min(lambda),lambdaMin);
  trN = sum(min(lambda,0)) + trN;
  dshift = min((-lambdaMin)/mean(-dshift),0) * dshift;
  setround(rnd);
  return;
end


% calculate error matrix
setround(-1);
E = -(L*L'-A) - diag(sparse(inf(n,1)));  % -inf(L*L'-A) = sup(A-L*L'), without diagonal
setround(1);
E = max(E,L*L'-A)+Arad; % max(sup(A-L*L')-inf*I,sup(L*L'-A))


if find(isnan(E) | tril(E,-1)<0,1)
  disp('break');
end

% clear variables to save memory
clear L D;

% verify that E is symmetric
E = min(E,E');

% keep setround(1) as default rounding for verified norm bounds


% find lower bound for minimum eigenvalue

% apply adaptive shift of diagonal elements
d = -E(1:n+1:end);    % diagonal of E had opposite sign
maxd = max(d);
E(1:n+1:end) = maxd - d;

% we have Enew = abs(Eold - D) with D = maxd*I
% hence lambda_i(Eold)>= lambda_i(D)-norm(Enew,2)

% using collatz theorem for upper bound of 2-norm
y = sum(E,2);
norm2bnd = max(y);

while norm2bnd
  oldnorm = norm2bnd;
  x = y / norm2bnd;
  y = E*x;
  norm2bnd = max(y./x);
  if norm2bnd*1.005>=oldnorm
    break;
  end
end

% lambda_min(D+E) >= lambda_min(D) - norm(E,2)
%   => lambdaMin >= inf(maxd-norm(E,2)) = -sup(norm(E,2)-maxd)
lambdaMinO = lambdaMin;  % old upper bound for the eigenvalues
lambdaMin = min(-(norm2bnd-maxd), lambdaMinO);

if lambdaMin>=0
  trN = 0;
  return;
elseif nargin>2 && mode>0  % do full eigenvalue enclosure
  A = struct('mid',A,'rad',Arad);
  clear Arad E;
  lambda = vsdpneig(A,1);
  lambda = lambda.rad - lambda.mid;  % -inf(lambda)
  lambdaMin = max(min(-max(lambda),lambdaMinO),lambdaMin);
  trN = - (sum(max(lambda,0)) - trN);
  dshift = min(lambdaMin,0) * (dshift/mean(dshift));
  setround(rnd);
  return;
end


% calculate trace bound

% trN := trace of negative semidefinite part of A
% -> sum(|lambda(E)|) <= sqrt(n) * norm(E,'fro')
% -> sum(lambda_pos(E))+sum(lamda_neg(E)) = trace(E);
%  => 2*sum(|lamda_neg(E)|) <= sqrt(n) * norm(E,'fro') - trace(E)
lambda_sn = 0.5 * ( sqrtsup(n*(E(:)'*E(:))) + sum(-diag(E)) );

if maxd<0  % if maxd<0: trN = trace(D) - sum(|lamda_neg(E)|)
  trN = -((-n)*maxd + lambda_sn);
else
  % for a given lambda_sn = sum(|lamda_neg(E)|), with maxd>0 and
  % nneg := number of negative eigenvalues it holds:
  % trN = - lambda_sn - nneg*maxd
  %  => trN is maximized when nneg is minimized
  nneg = floor(lambda_sn/norm2bnd);  % number of negative eigenvalues
  trN = -max((-nneg)*lambdaMin,lambda_sn+(-nneg-1)*maxd);
end


% reset rounding mode
setround(rnd);
end



% the following is a small function to calculate a verified upper bound
% for the square root       -> setround(1) is assumed
% only used in case that the sqrt function does not regard the rounding mode
function res = sqrtsup(a)
res = sqrt(a);
res = res + realmin * (res.*(-res) > -a);
end
