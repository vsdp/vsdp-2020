function [K,A,c,x,z] = sdpt2vsdp(blk,A,c,x,z)
%% SDPT2VSDP - converts data in SDPT3 input format to VSDP format
%    [K,A,c,x,z] = sdpt2vsdp(blk,At,c,x,z)
%
% except for blk all input parameter are optional
%
%% >> Input:
% blk: a cell array describing the block diagonal structure of SQL data
% At: a cell array with At{1} = A(1:dimf,:), At{2} = A(dimf+1:dimf+diml,:)
%     and the same for socp cone, and At{i} = [svec(Ai1) ... svec(Aim)]
%     for positive semidefinite cone constraint matrices, where Aij is the
%     matrix of j-th block and i-th constraint.
% c: a cell array of matrices of dimensions given by K
%       -> correspond to primal objective function
% x: a cell array of matrices of dimensions given by K
%       -> approximate primal optimal solution
% z: a cell array of matrices of dimensions given by K
%
%% >> Output:
% K: a structure with following fields
%     - K.f stores the number of free variables
%     - K.l is the number of nonnegative components
%     - K.q lists the lengths of socp blocks
%     - K.s lists the dimensions of semidefinite blocks
% A: a nA3 x M Matrix,
%     whereas nA = dimf+diml+dimq+dims
%     dimf: number of free variables: dimf = sum(K.f>0)
%     diml: number of nonnegative variables: diml = sum(K.l>0)
%     dimq: sum of all socp variables: dimq = sum_i(K.q(i))
%     dims3: sum of sdp variables: dims3 = sum_i(K.s(i)*(K.s(i)+1)/2)
% c: an nA3 x 1 vector  - svec(mu=1)
% x: an nA3 x 1 vector   -> approximate optimal solution  - svec(mu=2)
% z: an nA3 x 1 vector   -> approximate LMI dual vector  - svec(mu=1)
%
%%
% Note that the right hand side of the linear constraints (b) and the
% dual optimal solution vector (y) have the same format in VSDP and SDPT3.
%

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

%% check input
if nargin<1 || isempty(blk)
  error('VSDP:SDPT2VSDP','not enough input parameter');
elseif nargin>1 && iscell(A) && all(min(size(A))~=[0 1])
  error('VSDP:SDPT2VSDP',['this is an old SDPT3 input format '...
    'not supported by this function']);
elseif size(blk,2)>2 && any([blk{:,3}])
  error('VSDP:SDPT2VSDP','unsupported SDPT3 low rank structure input');
end


%% permutation index for cells if necessary

% declare blk-vector
blkvec = cat(1,blk{:,1});

% indices to sort variables  f -> l -> q -> s
indf = find(blkvec=='u');
indl = find(blkvec=='l');
indq = find(blkvec=='q');
inds = find(blkvec=='s');

ind = [indf; indl; indq; inds];

if size(ind,1)~=size(blk,1)
  error('VSDP:SDPT2VSDP','unsupported cell format or empty cells');
end


%% create K
K = struct('f',sum([blk{indf,2}]),'l',sum([blk{indl,2}]),...
  'q',[blk{indq,2}]','s',[blk{inds,2}]');

% what has to be processed
doA = nargin>1 && iscell(A) && ~isempty(A{1}) && nargout>1;
doC = nargin>2 && iscell(c) && ~isempty(c{1}) && nargout>2;
doX = nargin>3 && iscell(x) && ~isempty(x{1}) && nargout>3;
doZ = nargin>4 && iscell(z) && ~isempty(z{1}) && nargout>4;


%% convert c,x,z
if doC || doX || doZ
  % vectorize sdp blocks
  for i = 1:length(inds)
    j = inds(i);  % index of current sdp block group
    js = blk{j,2}(:);  % index set for sdp blocks of j-th sdp group
    if length(js)==1
      blkI = triu(true(js));
    else  % extract block diagonal
      % block diagonal index computation to avoid additional for-loop
      blkI = ones((js'*(js+1))/2,1);
      ofI = ones(sum(js),1);
      ofV = -ofI;
      js(end) = [];
      ofI(cumsum([2; js])) = [0; 1-js];
      ofV(cumsum([1; js])) = [1; js-1];
      ofV(2) = length(ofV) - (js(1)>1);
      blkI(cumsum(cumsum(ofI(1:length(ofV))))) = cumsum(ofV);
      blkI = cumsum(blkI);
    end
    if doC
      c{j} = c{j}(blkI);
    end
    if doX
      x{j} = x{j}(blkI);
    end
    if doZ
      z{j} = z{j}(blkI);
    end
  end
end

% cell to mat
if doC
  c = cat(1,c{ind});
else
  c = [];
end
if doX
  x = sscale(cat(1,x{ind}),K,2);
else
  x = [];
end
if doZ
  z = cat(1,z{ind});
else
  z = [];
end


%% convert A
if doA
  A = sscale(cat(1,A{ind}),K,sqrt(0.5));
else
  A = [];
end

end
