function [blk, At, b, c] = to_sdpt3_fmt (obj)
% TO_SDPT3_FMT  Convert problem data to SDPT3 format.
%
%   [blk,At,ct,xt,zt] = VSDP2SDPT3(K,A,c,x0,z0,opts)
%
% blk: a cell array describing the block diagonal structure of problem data
% At: a cell array with At{1} = A(1:dimf,:), At{2} = A(dimf+1:dimf+diml,:)
%     and the same for socp cone, and At{i} = [svec(Ai1) ... svec(Aim)]
%     for positive semidefinite cone constraint matrices, where Aij is the
%     matrix of j-th block and i-th constraint
% ct: a cell array of matrices of dimensions given by K
% xt: a cell array of matrices of dimensions given by K
% zt: a cell array of matrices of dimensions given by K
% Is: index vector which describes the reorganization of the sdp blocks to
%     group small sdp blocks
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% convert to appropriate format
if length(A)~=dim3
  A = vsvec(A,K,sqrt(2),0);
else
  A = sscale(A,K,sqrt(2));
end
Imat = [];
if length(c)~=dim
  [c,Imat] = vsmat(c,K,2,0);
end
if all(length(x0)~=[0 dim])  % vx0 = svec(x0(:),K,2), !! mu=2 !!
  [x0,Imat] = vsmat(x0,K,1,0,Imat);
end
if all(length(z0)~=[0 dim])
  z0 = vsmat(z0,K,2,0,Imat);
end

% row- or column-wise
idim = 1 + (size(A,2)==dim3);

% create cell arrays for output
blk = cell(n,2);
At = cell(n,1);
ct = cell(n,1);
xt = cell(n,1);
zt = cell(n,1);

% reducing the size of the input data after building the corresponding
% cells saves memory and increases the speed for the remaining data
% access
% the reduction, however, can have a significant effect on cpu-time
% the following is a tradeoff for cpu and memory efficiency
mem3tol = dim3 / min(n+1,3+sqrt(n));

% cell index
bli = 1;

% position indices
blk2e = 0;
blk3e = 0;


% transform unconstrained variables
if K.f>0
  blk{bli,1} = 'u';
  blk{bli,2} = K.f;
  if ~isempty(A) && idim==1
    At{bli} = A(blk3e+1:blk3e+K.f,:);
    blk3e = blk3e+K.f;
  elseif ~isempty(A)
    At{bli} = A(:,blk3e+1:blk3e+K.f)';
    blk3e = blk3e+K.f;
  end
  if ~isempty(c)
    ct{bli} = c(blk2e+1:blk2e+K.f);
  end
  if ~isempty(x0)
    xt{bli} = x0(blk2e+1:blk2e+K.f);
  end
  if ~isempty(z0)
    zt{bli} = z0(blk2e+1:blk2e+K.f);
  end
  blk2e = blk2e + K.f;
  bli = bli + 1;
end


% transform linear variables
if K.l>0
  blk{bli,1} = 'l';
  blk{bli,2} = K.l;
  if ~isempty(A) && idim==1
    At{bli} = A(blk3e+1:blk3e+K.l,:);
    blk3e = blk3e+K.l;
  elseif ~isempty(A)
    At{bli} = A(:,blk3e+1:blk3e+K.l)';
    blk3e = blk3e+K.l;
  end
  if ~isempty(c)
    ct{bli} = c(blk2e+1:blk2e+K.l);
  end
  if ~isempty(x0)
    xt{bli} = x0(blk2e+1:blk2e+K.l);
  end
  if ~isempty(z0)
    zt{bli} = z0(blk2e+1:blk2e+K.l);
  end
  blk2e = blk2e + K.l;
  bli = bli + 1;
end


% transform socp part
if ~isempty(K.q)
  dimq = sum(K.q);
  blk{bli,1} = 'q';
  blk{bli,2} = K.q(:)';
  if ~isempty(A) && idim
    At{bli} = A(blk3e+1:blk3e+dimq,:);
    blk3e = blk3e+dimq;
  elseif ~isempty(A)
    At{bli} = A(:,blk3e+1:blk3e+dimq)';
    blk3e = blk3e+dimq;
  end
  if ~isempty(c)
    ct{bli} = c(blk2e+1:blk2e+dimq);
  end
  if ~isempty(x0)
    xt{bli} = x0(blk2e+1:blk2e+dimq);
  end
  if ~isempty(z0)
    zt{bli} = z0(blk2e+1:blk2e+dimq);
  end
  blk2e = blk2e + dimq;
  bli = bli + 1;
end


% transform sdp part, consider small blocks
if ~isempty(K.s)
  
  % save index positions before processing sdp cells
  bls = bli;
  blk2s = blk2e;
  
  % starting index for current group
  kstart = 1;
  
  % create blk
  VSDP_OPTIONS = vsdpinit(opts);
  blke = 0;
  for k = 1:length(K.s)
    blke = blke + K.s(k)*(K.s(k)+1)/2;
    
    % group great enough or end of sdp data?
    if ((blke > VSDP_OPTIONS.MIN_SDPBLK_SIZE) || (k == length(K.s)))
      blk{bli,1} = 's';
      blk{bli,2} = reshape(K.s(kstart:k),1,[]);
      blke = 0;
      kstart = k + 1;
      bli = bli + 1;
    end
  end
  
  
  % create At cells for sdp part
  if ~isempty(A)
    % keep sdp part only
    if blk3e>0 && idim
      A(1:blk3e,:) = [];
      blk3e = 0;
    elseif blk3e>0
      A(:,1:blk3e) = [];
      blk3e = 0;
    end
    % create cell blocks
    if idim==1
      for k = bls:bli-1
        nk3 = sum(blk{k,2}.*(blk{k,2}+1))/2;
        At{k} = A(blk3e+1:blk3e+nk3,:);
        blk3e = blk3e + nk3;
        if blk3e>mem3tol
          A = A(blk3e+1:end,:);
          blk3e = 0;
        end
      end
    else
      for k = bls:bli-1
        nk3 = sum(blk{k,2}.*(blk{k,2}+1))/2;
        At{k} = A(:,blk3e+1:blk3e+nk3)';
        blk3e = blk3e + nk3;
        if blk3e>mem3tol
          A = A(:,blk3e+1:end);
          blk3e = 0;
        end
      end
    end
  end  % transformation of A
  
  
  % block - cell transformation for sdp part of c
  if ~isempty(c)
    blk2e = blk2s;  % reset index position
    for k = bls:bli-1
      nk = blk{k,2};
      lk = length(nk);
      if lk==1  % single sdp block
        ct{k} = reshape(c(blk2e+1:blk2e+nk*nk),nk,nk);
        ct{k} = 0.5 * (ct{k} + ct{k}');
        blk2e = blk2e + nk*nk;
      else  % block diagonal form of sdp group
        ct{k} = cell(1,lk);
        for i = 1:lk
          nki = nk(i);
          ct{k}{i} = reshape(c(blk2e+1:blk2e+nki*nki),nki,nki);
          ct{k}{i} = 0.5 * sparse(ct{k}{i} + ct{k}{i}');
          blk2e = blk2e + nki*nki;
        end
        ct{k} = blkdiag(ct{k}{:});
      end
    end
  end  % transformation of c
  
  
  % block - cell transformation for sdp part of x0
  if ~isempty(x0)
    blk2e = blk2s;  % reset index position
    for k = bls:bli-1
      nk = blk{k,2};
      lk = length(nk);
      if lk==1  % single sdp block
        xt{k} = reshape(x0(blk2e+1:blk2e+nk*nk),nk,nk);
        xt{k} = 0.5 * (xt{k} + xt{k}');
        blk2e = blk2e + nk*nk;
      else  % block diagonal form of sdp group
        xt{k} = cell(1,lk);
        for i = 1:lk
          nki = nk(i);
          xt{k}{i} = reshape(x0(blk2e+1:blk2e+nki*nki),nki,nki);
          xt{k}{i} = 0.5 * sparse(xt{k}{i} + xt{k}{i}');
          blk2e = blk2e + nki*nki;
        end
        nki = nk(end);
        xt{k} = blkdiag(xt{k}{:});
        blk2e = blk2e + nki*nki;
      end
    end
  end  % transformation of x0
  
  
  % block - cell transformation for sdp part of z0
  if ~isempty(z0)
    blk2e = blk2s;  % reset index position
    for k = bls:bli-1
      nk = blk{k,2};
      lk = length(nk);
      if lk==1  % single sdp block
        zt{k} = reshape(z0(blk2e+1:blk2e+nk*nk),nk,nk);
        zt{k} = 0.5 * (zt{k} + zt{k}');
        blk2e = blk2e + nk*nk;
      else  % block diagonal form of sdp group
        zt{k} = cell(1,lk);
        for i = 1:lk
          nki = nk(i);
          zt{k}{i} = reshape(z0(blk2e+1:blk2e+nki*nki),nki,nki);
          zt{k}{i} = 0.5 * sparse(zt{k}{i} + zt{k}{i}');
          blk2e = blk2e + nki*nki;
        end
        nki = nk(end);
        zt{k} = blkdiag(zt{k}{:});
        blk2e = blk2e + nki*nki;
      end
    end
  end  % transformation of z0
  
end  % transform sdp part


% remove unused cells
blk(bli:end,:) = [];
At(bli:end) = [];
ct(bli:end) = [];
xt(bli:end) = [];
zt(bli:end) = [];

if isempty(At{1})
  At = [];
end
if isempty(c)
  ct = [];
end
if isempty(x0)
  xt = [];
end
if isempty(z0)
  zt = [];
end

end
