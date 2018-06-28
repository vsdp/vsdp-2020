function A = sscale(obj, A, mu)
% SSCALE  Scale off-diagonal elements of SDP blocks by 'mu' in 'A'.
%
%    A = obj.sscale (A, mu);
%
% In general this function performs the actual scaling part of the SVEC and SMAT
% operations in VSDP.  It is not recommended to use 'mu = sqrt(2)' for verified
% computations.
%
% For better comprehension of SVEC, SMAT, and SSCALE, assume a symmetric 3x3
% matrix 'A':
%
%                            [a]
%                            [b]
%                            [c]                                [a   ]
%       [a b c]              [b]                                [b*mu]
%   A = [b d e]  ==>  A(:) = [d]  -->  vsdp.SVEC(A(:),mu)  -->  [c*mu] = a
%       [c e f]              [e]  <--  vsdp.SMAT(a, 1/mu)  <--  [d   ]
%                            [c]                                [e*mu]
%                            [e]                                [f   ]
%                            [f]
%
% For a symmetric matrix X it applies:
%   vx = vscale(X(:),K,mu)  =>  vx = vX(:), where vX = mu*X + (1-mu)*I
% The same holds valid for matrices in the SDPT3 format. (-> see svec)
%
% mu: scaling factor for off-diagonal elements
%     for performance reasons and reversibilty mu=0 is not allowed
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% check input
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
[isize,idim] = max(size(vA));
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


% create index vector for diagonal entries
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


% scale off-diagonal elements
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
