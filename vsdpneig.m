function lambda = vsdpneig(A,hermit)
% VSDPNEIG  rigorous enclosure of all eigenvalues of a normal matrix.
%           Input matrix A has to be normal.
%           A may also be an interval matrix, full or sparse. This interval
%           matrix can be defined as:
%               1. a structure with the fields 'mid' and 'rad',
%               2. a cell array, where the first entry will be interpreted
%                  as 'mid' and the second as 'rad',
%               3. a class of type 'intval' [see INTLAB toolbox].
%
%   @params
%       A:  Normal input matrix. Normality is not checked. However, if the
%           condition of the approximative eigen-vector matrix is too bad
%           the function will throw a warning.
%
%       hermit: Can be used to explicitly force the function to
%               treat A as a hermitian matrix. If not set and A is not
%               equal to A' the matrix will be treated as an
%               arbitrary normal matrix.
%
%   @output
%       lambda: If A is of type 'intval' lambda will be of the same type,
%               else lambda is a structure with the fields 'mid' and 'rad'
%               which are describing the corresponding interval vector.
%
%   @dependencies: 'setround'-function [see Intlab]
%
%       Replaces "veigsym" by Christian Jansson
%

%
% written    19/03/12  Marko Lange
%  modified  20/04/12   - allow interval structure input
%                       - removed most dependencies, still need 'setround'
%                         function
%  modified  07/08/12   - allow compact vectorized input
%

%
% Example:
%
% A = [1 2 3; 2 1 4; 3 4 5];
% A = midrad(A, 0.01*ones(3));
% lambda = vsdpneig(A,true)
% intval lambda =
% [   -1.5125,   -1.4615]
% [   -0.6172,   -0.5679]
% [    9.0506,    9.1084]
%
% or:
%
% A.mid = [1 2 3; 2 1 4; 3 4 5];
% A.rad = 0.01*ones(3);
% lambda = vsdpneig(A,true)
% lambda =
%     mid: [3x1 double]
%     rad: [3x1 double]
% lambda.mid
% ans =
%    -1.4870
%    -0.5925
%     9.0795
% lambda.rad
% ans =
%     0.0254
%     0.0246
%     0.0289
%


%% check input data %%

% input matrix
if nargin<1 || isempty(A)
  error('VSDPNEIG: No input matrix set.');
elseif ~(isnumeric(A) || isa(A,'intval') || all(isfield(A,{'mid','rad'})))
  error('VSDPNEIG: Datatype of input matrix is not accepted.');
end   % size & dimension is not checked. Error will occur when calling eig.

if nargin<2 || ~islogical(hermit)
  hermit = false;
end


%% prepare input data %%

% get rounding mode and set round to nearest
rnd = getround();
setround(0);    % for preparation only

% is user using INTLAB?
isintval = isa(A,'intval');

% prepare interval input
if ~isnumeric(A)  % A is an interval structure or of type 'intval'
  Arad = A.rad;
  A = A.mid;
else
  Arad = sparse(size(A,1),size(A,2));
end

%% prepare A
if size(A,1)~=size(A,2) && size(A,2)~=1
  A = A(:);
  Arad = Arad(:);
end
n = size(A,1);

if size(A,2)==1  % vectorized input
  hermit = true;
  n = sqrt(2*n+0.25)-0.5;
  if n~=round(n)
    error('VSDPNEIG: Wrong dimension of vectorized input');
  end
  E(triu(true(n))) = full(A);  % E is place-holder, will be overwritten later
  A = reshape(E,n,n);
  A = A + triu(A,1)';
  if any(Arad)
    X(triu(true(n))) = full(Arad);  % L is used as place-holder again
    Arad = reshape(X,n,n);
    Arad = Arad + triu(Arad,1)';
    if nnz(Arad)<n*n/5
      Arad = sparse(Arad);
    end
  else
    Arad = sparse(size(A,1),size(A,2));
  end
end

% constants
e = ones(n,1);


%% approximate eigenvalue decomposition %%
[X,d] = eig(full(A));
d = diag(d);

% use Rayleigh quotient approximation instead
% d = ( sum(conj(X).*(A*X),1)./(sum(conj(X).*X,1).^0.5) ).';


%% calculate verified bounds by applying perturbation theorems %%

% default rounding mode for verification -> up
setround(1);

% E = mag(I-X'*X)
if isreal(X)
  setround(-1);
  E = -(X'*X-speye(n));  % -inf(X'*X-I)
  setround(1);
  % E = max(-inf(X'*X-I),sup(X'*X-I)) >= ||I-X'*X||
  E = max(E,X'*X-speye(n));
else
  realX = real(X);    % by using these variables the compiler will
  imagX = imag(X);    % recognize symmetry -> faster
  
  % sup(imag(X')*real(X)+real(X')*imag(X)) = sup(imag(X'*X))
  E = imagX'*(-realX) + realX'*imagX;
  % abs(imag(X'*X)) <= max(sup(imag(X'*X)),-inf(imag(X'*X)))
  % sup(imag(X'*X)) = -inf(imag(X'*X))'
  E = max(E,E');  % overwrite E for memory reasons, here E >= abs(imag(X'*X))
  
  % E >= abs(abs(real(X'*X)+j*abs(imag(X'*X)))
  % abs(real(X'*X)) <= max(sup(real(X'*X)),-inf(real(X'*X)))
  setround(-1);
  INF = realX'*realX + imagX'*imagX - speye(size(X,2));
  setround(1);
  SUP = realX'*realX + imagX'*imagX - speye(size(X,2));
  E = abs( max(SUP,-INF) + 1i*E);
  
  clear INF SUP realX imagX;
  
  % significantly easier code with Intlab, but slower + more dependencies
  % E = mag(eye(n)-intval(X)*X');
end

% kappa <= 1 - ||X'*X-I|| <= omin(X)^2
kappa = -(sqrtsup(norm22pos(E))-1);       % avoid switch of rounding mode

% warning if A may not be normal
if kappa<eps
  warning('VSDPNEIG:InclusionFailed','Eigenvector approximation is too bad.');
  if isintval
    lambda = midrad(d+1i,inf(size(d)));
  else
    lambda = struct('mid',d+1i,'rad',inf(size(d)));
  end
  return;
elseif kappa<1-eps*(n+1)^2
  warning('VSDPNEIG:WeakEnclosure','A may be not normal.');
end


% residual matrix R = magnitude(A*X-X*D)
if isreal(X) || isreal(A)
  SUP = A*X + X.*((-e)*real(d)') + (-1i*X).*(e*imag(d)');  % sup(A*X-X*D)
  INF = A*(-X) + X.*(e*real(d)') + (1i*X).*(e*imag(d)');  % -inf(A*X-X*D)
else   % A and X are complex
  SUP = real(A)*X + imag(A)*(1i*X) + ...
    X.*((-e)*real(d).') + (-1i*X).*(e*imag(d).');  % sup(A*X-X*D)
  INF = real(A)*(-X) + imag(A)*(-1i*X) + ...
    X.*(e*real(d).') + (1i*X).*(e*imag(d).');  % -inf(A*X-X*D)
end

if isreal(SUP) && isreal(INF)
  % abs(real(R)) <= max(sup(real(A*X-X*D)),-inf(real(A*X-X*D)))
  R = max(SUP,INF);
else
  % R >= abs(abs(real(A*X-X*D)) + j*abs(imag(A*X-X*D))) >= abs(A*X-X*D)
  R = abs( max(real(SUP),real(INF)) + 1i * max(imag(SUP),imag(INF)) );
end

if ~isempty(find(Arad,1))  % if interval radius~=0
  R = R + Arad*abs(X);
end

% free some memory
clear INF SUP Arad X


% r >= ||A*X-X*D|| / omin(X) %
if hermit
  r = sqrtsup( norm22pos(R) / kappa );
else  % else frobenius-norm
  r = sqrtsup( (R(:)'*R(:)) / kappa );
end

drad = r*e;  % radius of eigenvalue diagonal vector


%% using generalized perturbation theorem for tighter inclusion %%

% split into clusters of eigenvalues and update residual norms
if hermit  % Cao: sharp version of Kahan's perturbation theorem
  
  [d,index] = sort(real(d),1,'ascend');
  
  k = 1;  % index of first eigenvalue in current cluster
  
  for i = 1:n
    if i>=n || d(i+1)>d(i)+r  % if end of cluster
      % upper bound for matrix norm
      kappa = -( sqrtsup( norm22pos(E(index(k:i),index(k:i))) ) - 1 );
      drad(k:i) = sqrtsup( norm22pos(R(:,index(k:i))) / kappa );
      k = i+1;  % set new starting position for next cluster
    end
  end
  
else  % matrix is normal: use frobenius-norm for inclusion, see Kahan
  
  % build cluster graph  -> abs(e*d'-d*e')<r
  graph = abs( max(e*real(d)'+(-real(d))*e',real(d)*e'+e*(-real(d)')) + ...
    1i*max(e*imag(d)'+(-imag(d))*e',imag(d)*e'+e*(-imag(d)')) );
  graph = graph<r;
  
  % remove diagonal entries (self-connection) of graph
  graph(1:n+1:n*n) = 0;
  
  index = 1:n;  % indices of eigenvalues
  k = 1;  % first index of current eigenvalue cluster
  l = 1;  % last index of current eigenvalue cluster
  
  for i = 1:n
    
    % remove row entries for current index
    graph(index(i),:) = 0;
    
    % find neighbors of eigenvalue lambda_i
    cluster = find(graph(:,index(i)));
    
    if ~isempty(cluster)  % found clustered eigenvalues
      
      % extend cluster & switch corresponding indices
      index(cluster) = l+1:l+length(cluster);
      index(l+1:l+length(cluster)) = cluster;
      l = l + length(cluster);  % update l
      
      % remove row entries corresponding to added eigenvalues
      graph(cluster,:) = 0;
      
    elseif i>=l   % last position of current cluster ?
      
      kappa = -( sqrtsup(norm22pos(E(index(k:l),index(k:l)))) - 1 );
      Rp = reshape(R(:,index(k:l)),[],1);  % get residuum for current cluster
      drad(index(k:l)) = sqrtsup((Rp'*Rp) / kappa);
      k = i + 1;  % set position of next cluster
      l = k;
      
    end
    
  end
  
end  % is normal


% update interval vector for lambda inclusion
if isintval
  lambda = midrad(d,drad);
else
  lambda = struct('mid',d,'rad',drad);
end


% reset rounding mode
setround(rnd);
end



%%% modified normposA from INTLAB by Siegfried M. Rump %%%

function normbnd = norm22pos(A)
% upper bound for squared spectral norm of nonnegative matrix A
% applying collatz theorem to power iteration
%
% note: this is a helper function used by vsdpneig,
%       parameters are not checked

% setround(1)       % here not set, already done in calling function

x = full(sum(A,2)'*A)';
normbnd = max(x);
x = x / normbnd;

while normbnd
  oldnorm = normbnd;
  y = ((A*x)'*A)';
  normbnd = max(y./x);
  if normbnd*1.005>=oldnorm
    break;
  end
  x = y / normbnd;
end
end



% the following is a small function to calculate a verified upper bound
% for the square root       -> setround(1) is assumed
% only used in case the sqrt function does not regard the rounding mode
function res = sqrtsup(a)
res = sqrt(a);
res = res + realmin * (res*(-res) > -a);
end
