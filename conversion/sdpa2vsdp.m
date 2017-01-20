function [A,b,c,K] = sdpa2vsdp(filename,blksize)
%% SDPA2VSDP: reads data from a SDPA file and converts them to VSDP
%    [A,b,c,K] = sdpa2vsdp(sparse-problem-data.dat-s)
%    [A,b,c,K] = sdpa2vsdp(dense-problem-data.dat)
%    [x,y,INFO] = sdpa2vsdp(optimization-result.out)
%    [x0,y0,z0,K] = sdpa2vsdp(sparse-initial-data.ini-s,blksize)
%    [x0,y0,z0,K] = sdpa2vsdp(dense-initial-data.ini,blksize)
%
%% >> Input:
% filename: string with the path + filname to the problem
%
%% >> Output:
% A: a nA3 x M Matrix,
%     whereas nA3 = dimf+diml+dimq+dims3
%     dimf: number of free variables: dimf = sum(K.f>0)
%     diml: number of nonnegative variables: diml = sum(K.l>0)
%     dimq: sum of all socp variables: dimq = sum_i(K.q(i))
%     dims3: sum of all sdp variables: dims3 = sum_i(K.s(i)*(K.s(i)+1)/2)
% b: a M x 1 vector
% c: a nC x 1 vector,
%     whereas nC = dimf+diml+dimq+dims
%     dimf, diml, dimq are just the same as for nA3
%     dims = sum_i(K.s(i)^2)
% K: a structure with following fields
%     - K.f stores the number of free variables
%     - K.l is the number of nonnegative components
%     - K.q lists the lengths of socp blocks
%     - K.s lists the dimensions of semidefinite blocks
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
% 25/07/12    M. Lange, rewrite for improved performance
% 04/08/12    M. Lange, add support for initial and full data
%
%%
% TODO: (error when reading large scale problems)
%

%% open file and import data
if nargin~=1 || ~ischar(filename) || length(filename)<4
  error('VSDP:SDPA2VSDP','input argument must be a filename with extension');
end

% type of sdpa-data
extLIST = {'at-s','.dat','.out','ni-s','.ini'};
dtype = find(strncmpi(filename(end-3:end),extLIST,4),1);
if isempty(dtype)
  error('VSDP:SDPA2VSDP','extension "%s" is not supported',filename(end-3:end));
end

% open file
fid = fopen(filename,'r');
if fid==-1
  error('VSDP:SDPA2VSDP','Cannot open %s',filename);
end

switch dtype
  case 1  % '.dat-s'  sparse problem data
    
    % skip comments
    str = fgetl(fid);
    while str(1)=='*' || str(1)=='"'
      str = fgetl(fid);
    end
    
    m = sscanf(str,'%d',1);  % number of decision variables
    nblocks = fscanf(fid,'%d',1);  % number of blocks
    blksize = fscanf(fid,'%*[^0-9+-]%d',nblocks);  % size of each block
    
    % right hand side vector (b)
    b = -fscanf(fid,'%*[^0-9+-]%lf',m);
    
    %% read data of A and c
    [data,cnt] = fscanf(fid,'%*[^0-9+-]%d %d %d %d %lf',[5 inf]);
    fclose(fid);
    if cnt==0 || mod(cnt,5)~=0
      error('VSDP:SDPA2VSDP','Could not read SDPA data');
    end
    [col,blk,i,j,data] = deal(data(1,:), data(2,:), data(3,:), data(4,:), -data(5,:));
    
    % any i>j ?
    idx = find(i>j);
    if ~isempty(idx)
      warning('VSDP:SDPA2VSDP','Lower triangular elements will be ignored.');
      col(idx) = [];  blk(idx) = [];  i(idx) = [];  j(idx) = [];  data(idx) = [];
    end
    
    %% block information
    idx = blksize>1;  % index vector for sdp blocks
    dims = abs(blksize(~idx));  % dimensions of lp blocks
    
    % cone structure
    K.l = sum(dims);  % number of all linear variables
    K.s = reshape(blksize(idx),[],1);  % dimensions of sdp blocks
    
    %% calculate offset positions
    ofs = zeros(1,length(blksize));  % block offsets
    ofs(~idx) = cumsum([0 dims(1:end-1)]);  % linear block offsets
    ofs(idx) = cumsum([K.l; K.s(1:end-1).*(K.s(1:end-1)+1)/2]);
    % offsets = i + block-offsets + column-offsets
    idx = 0.5*idx';  % factor 0.5 for sdp part column offset:
    i = i + ofs(blk) + idx(blk).*j.*(j-1);  % overwrite i by offsets
    
    %% create output
    dim3 = K.l + (K.s'*(K.s+1))/2;
    j = col>0;  % overwrite j by index vector for A
    A = sparse(i(j),col(j),data(j),dim3,m);
    c = sparse(i(~j),1,data(~j),dim3,1);
    
  case 2  % '.dat'  dense problem data
    
    % skip comments
    str = fgetl(fid);
    while str(1)=='*' || str(1)=='"'
      str = fgetl(fid);
    end
    
    m = sscanf(str,'%d',1);  % number of decision variables
    nblocks = fscanf(fid,'%d',1);  % number of blocks
    blksize = fscanf(fid,'%*[^0-9+-]%d',nblocks);  % size of each block
    
    % right hand side vector (b)
    b = -fscanf(fid,'%*[^0-9+-]%lf',m);
    
    %% block information
    idx = blksize>1;  % index vector for sdp blocks
    dims = abs(blksize(~idx));  % dimensions of lp blocks
    
    % cone structure
    K.l = sum(dims);  % number of all linear variables
    K.s = reshape(blksize(idx),[],1);  % dimensions of sdp blocks
    
    dim = K.l + K.s'*K.s;
    
    %% read data of A and c
    [c,cnt] = fscanf(fid,'%*[^0-9+-]%f',[dim 1]);
    if cnt~=dim
      fclose(fid);
      error('VSDP:SDPA2VSDP','Could not read SDPA data');
    end
    [A,cnt] = fscanf(fid,'%*[^0-9+-]%f',[dim m]);
    if cnt~=dim*m
      fclose(fid);
      error('VSDP:SDPA2VSDP','Could not read SDPA data');
    end
    fclose(fid);
    
    %% sort data: lp first
    if find(idx,1)<=find(~idx,1,'last')
      idxs = cell(length(idx),1);
      for j = find(idx)'  % sdp part
        idxs{j} = true(blksize(j)^2,1);
      end
      for j = find(~idx)'  % lp part
        idxs{j} = false(abs(blksize(j)),1);
      end
      idxs = vertcat(idxs{:});  % sdp index
      idx = ~idxs;  % lp index
      c = -[c(idx); c(idxs)];
      A = -[A(idx,:); A(idxs,:)];
    else
      c = -c;
      A = -A;
    end
    
  case 3  % '.out'  solver output
    
    error('output data not yet supported');
    
  case 4  % '.ini-s'  sparse initial data
    
    % cone information from blksize
    if nargin<2 || isempty(blksize)
      error('VSDP:SDPA2VSDP','block sizes "blksize" has to be set to read initial data');
    end
    
    % dual initial vector y0
    str = fgetl(fid);
    b = sscanf(str,'%*[^0-9+-]%lf',[inf 1]);
    
    %% read data of x0 and z0
    [data,cnt] = fscanf(fid,'%*[^0-9+-]%d %d %d %d %lf',[5 inf]);
    fclose(fid);
    if cnt==0 || mod(cnt,5)~=0
      error('VSDP:SDPA2VSDP','Could not read SDPA data');
    end
    [col,blk,i,j,data] = deal(data(1,:), data(2,:), data(3,:), data(4,:), data(5,:));
    
    % any i>j ?
    idx = i>j;
    if any(idx)
      warning('VSDP:SDPA2VSDP','Lower triangular elements will be ignored.');
      col(idx) = [];  blk(idx) = [];  i(idx) = [];  j(idx) = [];  data(idx) = [];
    end
    
    %% block information
    idx = blksize>1;  % index vector for sdp blocks
    dims = abs(blksize(~idx));  % dimensions of lp blocks
    
    % cone structure
    K.l = sum(dims);  % number of all linear variables
    K.s = reshape(blksize(idx),[],1);  % dimensions of sdp blocks
    
    %% calculate offset positions
    ofs = zeros(1,length(blksize));  % block offsets
    ofs(~idx) = cumsum([0 dims(1:end-1)]);  % linear block offsets
    ofs(idx) = cumsum([K.l; K.s(1:end-1).*(K.s(1:end-1)+1)/2]);
    % offsets = i + block-offsets + column-offsets
    idx = 0.5*idx';  % factor 0.5 for sdp part column offset:
    blk = i + ofs(blk) + idx(blk).*j.*(j-1);  % overwrite blk by offsets
    
    %% create output
    dim3 = K.l + (K.s'*(K.s+1))/2;
    col = col<2;  % overwrite col by index vector for z0
    A = sparse(blk(col),1,data(col),dim3,1);  % z0
    data = data .* (1 + (i<j));  % mu=2, for x0
    c = sparse(blk(~col),1,data(~col),dim3,1);  % x0
    
  case 5  % '.ini'  dense initial data
    
    % cone information from blksize
    if (nargin<2 || isempty(blksize))
      error('VSDP:SDPA2VSDP','block sizes "blksize" has to be set to read initial data');
    end
    
    % dual initial vector y0
    str = fgetl(fid);
    b = -sscanf(str,'%*[^0-9+-]%lf',[inf 1]);
    
    %% block information
    idx = blksize>1;  % index vector for sdp blocks
    dims = abs(blksize(~idx));  % dimensions of lp blocks
    
    % cone structure
    K.l = sum(dims);  % number of all linear variables
    K.s = reshape(blksize(idx),[],1);  % dimensions of sdp blocks
    
    dim = K.l + K.s'*K.s;
    
    %% read data of z0 and x0
    [c,cnt] = fscanf(fid,'%*[^0-9+-]%f',[dim 1]);
    if (cnt~=dim)
      fclose(fid);
      error('VSDP:SDPA2VSDP','Could not read SDPA data');
    end
    [A,cnt] = fscanf(fid,'%*[^0-9+-]%f',[dim 1]);
    if (cnt~=dim)
      fclose(fid);
      error('VSDP:SDPA2VSDP','Could not read SDPA data');
    end
    fclose(fid);
    
    %% sort data: lp first
    if (find(idx,1)<=find(~idx,1,'last'))
      idxs = cell(length(idx),1);
      for j = find(idx)'  % sdp part
        idxs{j} = true(blksize(j)^2,1);
      end
      for j = find(~idx)'  % lp part
        idxs{j} = false(abs(blksize(j)),1);
      end
      idxs = vertcat(idxs{:});  % sdp index
      idx = ~idxs;  % lp index
      c = [c(idx); c(idxs)];
      A = [A(idx,:); A(idxs,:)];
    end
end

end
