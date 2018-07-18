function obj = solve_sdpa (obj)
% SOLVE_SDPA  Approximately solve conic problem instance with SDPA.
%
%   For more information on the SDPA-M format, see:
%
%     [1] https://sourceforge.net/projects/sdpa/files/sdpa-m/sdpamManual.pdf
%         Version 2005.
%
%   See also vsdp.solve.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% Check cones
if ((obj.K.f > 0) || (sum(obj.K.q) > 0))
  error('VSDP:solve_sdpa:unsupportedCones', ...
    ['solve_sdpa: free variables (K.f) second-order cone variables (K.q) ', ...
    'are not supported by SDPA.']);
end

% Note: for 'x0' and later 'x' see [1, p. 14] "mDIM -- All the letters after m
% through the end of the line are neglected".

% Should initial solution guess be taken into account?
if (obj.options.USE_STARTING_POINT)
  [x0, X0, Y0] = deal (obj.y, obj.z, obj.x);
  x0 = [x0; 0];  % expand to mDIM
  X0 = mat2cell (X0,  obj.K.dims, 1);
  Y0 = mat2cell (Y0,  obj.K.dims, 1);
  X0 = cellfun(@(x) vsdp.smat([], x, 1),   X0, 'UniformOutput', false);
  Y0 = cellfun(@(x) vsdp.smat([], x, 1/2), Y0, 'UniformOutput', false);
else
  [x0, X0, Y0] = deal ([], [], []);
end
% Clear previous solutions or initial points.
[obj.x, obj.y, obj.z] = deal ([], [], []);

if (~isempty (obj.options.SOLVER_OPTIONS))
  OPTIONS = obj.options.SOLVER_OPTIONS;
else
  OPTIONS = [];
end
if (~obj.options.VERBOSE_OUTPUT)
  OPTIONS.print = 'no';
end

% In case of interval data solve midpoint problem.
c = mid (-obj.b);
F = [ ...
  mat2cell(mid (-obj.c),  obj.K.dims, 1), ...
  mat2cell(mid (-obj.At), obj.K.dims, ones(1, obj.m))];
F = cellfun(@(x) vsdp.smat([], x, 1), F, 'UniformOutput', false);
[nBLOCK, mDIM] = size (F);
bLOCKsTRUCT = [-obj.K.l(obj.K.l > 0), obj.K.s'];

% Call solver.
if (exist ('mexsdpa', 'file') == 3)
  [~, x, X, Y, ~] = sdpam ...
    (mDIM, nBLOCK, bLOCKsTRUCT, c, F, x0, X0, Y0, OPTIONS);
elseif (exist('callSDPA','file') == 2)
  [x, X, Y] = callSDPA ...
    (mDIM, nBLOCK, bLOCKsTRUCT, c, F, x0, X0, Y0, OPTIONS);
else
  error('VSDP:MYSDPS', 'You need to compile the SDPA MEX-interface.');
end

% Store results.
obj.x = vsdp.svec (obj, vsdp.cell2mat (Y), 2);
obj.y = x(1:end-1);
obj.z = vsdp.svec (obj, vsdp.cell2mat (X), 1);

info = 0;

end

function [c, F, x0, X0, Y0] = convert (obj)
% CONVERT  Convert problem data to SDPA format.
%
%   [mDIM,nBLOCK,bLOCKsTRUCT,c,F,x0,X0,Y0] = vsdp2sdpam (obj)
%
%
% Output:
% mDIM: scalar - number of primal variables
% nBLOCK: scalar - number of blocks of F
% bLOCKsTRUCT: scalar - represents the block structure of F
% c: vector - coefficient vector
% F: cell array of coefficient matrices
% x0: vector for initial solution
% X0,Y0: cell arrays of initial points
%

% non-compact vectorized format
xdim = ~isempty(x);
zdim = ~isempty(z);
[c,Imat] = vsmat(c,K,1,1);
if zdim
  c = [vsmat(z,K,1,1,Imat) c];
end
if xdim
  c = [vsmat(x,K,0.5,1,Imat) c];
end

% A = [x z c A]
if size(A,1)<size(A,2)
  A = -[c vsmat(A,K,1,1,Imat)'];
else
  A = -[c vsmat(A,K,1,1,Imat)];
end

% initialize outputs
c = full(-b);  % full instead of sparse to catch a bug in sdpa
x0 = y;
F = cell(nBLOCK,mDIM+1);
X0 = cell(nBLOCK,1);
Y0 = cell(nBLOCK,1);

% transform linear part
if K.l
  At = A(1:K.l,:);
  Y0{1} = At(:,1:xdim);
  X0{1} = At(:,xdim+1:xdim+zdim);
  for i = 1:mDIM+1
    F{1,i} = At(:,xdim+zdim+i);
  end
end

% transform sdp part
blkj = K.l;
for j = (K.l>0)+1:nBLOCK
  nj = K.s(j-(K.l>0));
  At = reshape(A(blkj+1:blkj+nj*nj,:),nj,[]);
  blkj = blkj + nj*nj;
  Y0{j} = At(:,1:xdim*nj);
  X0{j} = At(:,xdim*nj+1:(xdim+zdim)*nj);
  blki = xdim + zdim;
  for i = 1:mDIM+1
    F{j,i} = At(:,blki+1:blki+nj);
    blki = blki + nj;
  end
end

end
