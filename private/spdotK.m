function [res, resrad, q] = spdotK(varargin)
%% SPDOTK   Dot products in K-fold precision
%
%   [res, resrad] = spdotK(A1,x1,...,An,xn,K)
%
% where res is approximation of A1*x1 + A2*x2 + ... ,
%   or A1.'*x1 + A2.'*x2 + ...  respectively (useful for sparse matrices)
% at least one of the factors has to be real
%
% the exact result is enclosured by [res-resrad, res+resrad]
% for K>2 the error vector of res(i) is saved in q{i}
%
% @parameters
%     Ai  - matrices, size(Ai,2) = length(res)
%     xi  - vectors, size(xi,1) = size(Ai,1)
%     K   - precision value: the greater K the more accurate the result
%           K has to be an integer >1
%               default: 2
%
% @dependencies
%     for the computation of 'resrad' the setround function is needed,
%          [for further information, see INTLAB]
%
% Note!
%   this function implements algorithms from
%   T. Ogita, S.M. Rump, S. Oishi: Accurate Sum and Dot Product,
%     SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005
%
%% last modified
%   27/08/12   M.Lange, written
%   05/09/12   M.Lange, added code for full matrices
%


% fast check for rounding mode
rnd = getround();
crnd = rnd;
if rnd~=0
  setround(0);
  crnd = 0;
end

% get precision value
if even(nargin)
  K = 2;
else
  K = varargin{end};
end

% check precision K
if K<2
  error('SPDOTK: invalid parameter K');
end

% initialization
factor = 134217729;  % splitting factor 2^27+1


%% prepare input

% check whether Ai*xi or Ai.'*xi shall be computed
szs = zeros(floor(nargin/2),4);  % [size(Ai) size(xi); ...]
for i = 2:2:nargin
  szs(i/2,:) = [size(varargin{i-1}) size(varargin{i})];
end
horz = all(szs(:,2)==szs(:,3)) && all(szs(:,1)==szs(1,1));
vert = ~horz && all(szs(:,1)==szs(:,3)) && range(szs(:,2))==0;

% check dimesion
if any(szs(:,4)~=1)
  error('SPDOTK: second factor must be vector');
elseif ~horz && ~vert
  error('SPDOTK: dimension of factors do not match');
end

% concat Ai and xi
if nargin>3 && horz
  A = cat(2,varargin{1:2:end-1});
  x = full(cat(1,varargin{2:2:end}));
elseif nargin>3
  A = cat(1,varargin{1:2:end-1});
  x = full(cat(1,varargin{2:2:end}));
else
  A = varargin{1};
  x = full(varargin{2});
end


if 10*nnz(A)<numel(A);
  %% two-product - sparse
  
  % vertical: A.'*x
  if horz
    A = A';
  end
  sA2 = size(A,2);
  
  % preallocation for output
  res = zeros(sA2,1);
  q = cell(sA2,1);
  
  % get non-zero elements
  [i,j,sa] = find(A);
  x = x(i);
  h = sa .* x;  % approximation
  
  % only zeros
  if ~any(h)
    resrad = res;  % zeros
    if crnd~=rnd
      setround(rnd)
    end
    return;
  end
  
  % index vector for columns
  i = cat(1,1,find(diff(j))+1);  % column-break indices
  j = j(i)';  % column indices
  i(end+1) = size(sa,1) + 1;  % last column-break
  
  % splitting
  au = factor * sa;
  au = au - (au - sa);  % upper part
  sa = sa - au;  % lower part
  
  xu = factor * x;
  xu = xu - (xu - x);
  x = x - xu;
  
  % sa.*x = h + r
  r = sa.*x - (((h - au.*xu) - sa.*xu) - au.*x);
  
  % free some memory
  clear A sa au x xu;
  
  % save one iteration if only res is needed
  K = K - (nargout<2);
  
  % walk through column indices
  nk = 0;
  for k = j
    % nonzero elements in column
    nk = nk + 1;
    ind = i(nk):i(nk+1)-1;
    hk = cat(1,h(ind));
    rk = r(ind);
    % first step of k-fold sum
    cK = K - 1;
    c = cumsum(hk);
    cl = cat(1,0,c);  cl(end) = [];
    z = c - cl;
    hk = (cl - (c - z)) + (hk - z);
    % finish if K==2 or no errors
    if cK<1 || ~(any(hk) || any(rk))
      res(k) = c(end);
      q{k} = nonzeros(cat(1,hk,rk));
      continue;
    end
    % k-fold sum loop
    hk = nonzeros(cat(1,hk,rk,c(end)));
    hlen = size(hk,1);
    while 1
      cK = cK - 1;
      c = cumsum(hk);
      cl = cat(1,0,c);  cl(end) = [];
      z = c - cl;
      hk = nonzeros((cl-(c-z))+(hk-z));
      nhlen = size(hk,1) + 1;
      if cK<1 || hlen==2 || nhlen<2 || ...
          (hlen==nhlen && isequal(hk,c(1:end-1)))
        break;
      end
      hk(end+1,1) = c(end);
      hlen = nhlen;
    end
    res(k) = c(end);
    q{k} = hk;
  end
  
else
  %% two-product - full
  
  % use full matrix
  if issparse(A)
    A = full(A);
  end
  % horizontal A*x
  if vert
    A = A';
  end
  x = x';
  sA1 = size(A,1);
  H = bsxfun(@times,A,x);  % approximation of A.*X
  
  % preallocation for output
  res = zeros(sA1,1);
  q = cell(sA1,1);
  
  % splitting
  Au = factor * A;
  Au = Au - (Au - A);  % upper part
  A = A - Au;  % lower part
  
  xu = factor * x;
  xu = xu - (xu - x);
  x = x - xu;
  
  % A.*X = H + R
  % here: A = [H; R]  with R = Al.*Xl-(((H-Au.*xu)-Al.*Xu)-Au.*Xl)
  R = bsxfun(@times,A,x) - (((H - bsxfun(@times,Au,xu)) - ...
    bsxfun(@times,A,xu)) - bsxfun(@times,Au,x));
  
  % free some memory
  clear A Au xu x;
  
  % first sumK step
  C = cumsum(H,2);
  Cl = cat(2,res,C);  Cl(:,end) = [];
  Z = C - Cl;
  H = (Cl - (C - Z)) + (H - Z);
  
  % finish if K==2
  if K==2
    % write result
    if nargout<2
      res = C(:,end) + sum(H+R,2);
    else
      setround(-1);
      inf = C(:,end) + sum(H+R,2);
      setround(1);
      sup = C(:,end) + sum(H+R,2);
      res = 0.5 * (inf + sup);  % mid
      resrad = res - inf;  % rad
      crnd = 1;
    end
    q = [];  % not computet for K==2
    % undo setround and return
    if crnd~=rnd
      setround(rnd)
    end
    return;
  end
  
  % second sumK step
  H = cat(2,H,R,C(:,end));  % add two-product error
  C = cumsum(H,2);
  Cl = cat(2,res,C);  Cl(:,end) = [];
  Z = C - Cl;
  H = (Cl - (C - Z)) + (H - Z);
  K = K - 2 - (nargout<2);  % save one iteration if only res is needed
  
  % finish if K==3 and only result is needed
  if K==0
    res = C(:,end) + sum(H,2); % false
    if crnd~=rnd
      setround(rnd);
    end
    return;
  end
  
  % faster processing along the columns
  H = cat(2,H,C(:,end))';
  
  % free some memory
  clear C Cl Z R;
  
  % main loop for K-fold sum
  for i = 1:sA1
    % read column
    h = nonzeros(H(:,i));
    % break if single element
    if size(h,1)<2
      res(i) = sum(h);
      continue;
    end
    % k-fold sum loop
    while 1
      K = K - 1;
      c = cumsum(h);
      cl = cat(1,0,c);  cl(end) = [];
      z = c - cl;
      h = nonzeros((cl-(c-z))+(h-z));
      if K<1 || size(c,1)==2 || size(h,1)==0 || isequal(h,c(1:end-1))
        break;
      end
      h(end+1,1) = c(end);
    end
    res(i) = c(end);
    q{i} = h;
  end
  
end


%% compute radius
if nargout<2
  % for-loop is faster than  res + cellfun(@(x) sum(x),q)
  for i = 1:length(q)
    res(i) = res(i) + sum(q{i});
  end
else
  setround(1);
  resrad = res;  % preallocation
  for i = length(q):-1:1
    resrad(i) = sum(abs(q{i}));
  end
  crnd = 1;
end

% reset rounding
if crnd~=rnd
  setround(rnd)
end

end
