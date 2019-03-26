function obj = from_sdpa_file (fname, blksize)
% FROM_SDPA_FILE  Import conic problem from a SDPA file 'fname'.
%
%   The parameter 'blksize' is used for initial data import.
%
%   Example:
%
%   obj = vsdp.FROM_SDPA_FILE ('C:\path\to\sparse-problem-data.dat-s');
%   obj = vsdp.FROM_SDPA_FILE ('C:\path\to\dense-problem-data.dat');
%   obj = vsdp.FROM_SDPA_FILE ('C:\path\to\sparse-initial-data.ini-s', blksize);
%   obj = vsdp.FROM_SDPA_FILE ('C:\path\to\dense-initial-data.ini',    blksize);
%
%   See also vsdp.
%

% Copyright 2004-2019 Christian Jansson (jansson@tuhh.de)

narginchk(1, 2);
if (exist (fname, 'file') ~= 2)
  error ('VSDP:from_sdpa_file:file_not_exists', ...
    'from_sdpa_file: Input file does not exist.');
end

fid = fopen (fname, 'r');
if (fid == -1)
  error ('VSDP:from_sdpa_file:cannot_open_file', ...
    'from_sdpa_file: Cannot open %s', fname);
end

[~,~,fext] = fileparts (fname);
switch (fext)
  case '.dat-s'  % Sparse problem data
    
    % Skip comment lines
    str = fgetl(fid);
    while ((str(1) == '*') || (str(1) == '"'))
      str = fgetl(fid);
    end
    
    m = sscanf (str, '%d', 1);  % number of decision variables
    nblocks = fscanf (fid, '%d', 1);  % number of blocks
    blksize = fscanf (fid, '%*[^0-9+-]%d', nblocks);  % size of each block
    
    % Right-hand side vector.
    b = -fscanf(fid, '%*[^0-9+-]%lf', m);
    
    % read data of A and c
    [data, cnt] = fscanf(fid, '%*[^0-9+-]%d %d %d %d %lf', [5, inf]);
    fclose (fid);
    if ((cnt == 0) || (mod(cnt,5) ~= 0))
      error ('VSDP:FROM_SDPA_FILE:data_import_failed', ...
        'from_sdpa_file: Could not read SDPA data');
    end
    [col, blk, i, j, data] = deal (data(1,:), data(2,:), data(3,:), ...
      data(4,:), -data(5,:));
    
    % any i>j ?
    idx = find(i>j);
    if ~isempty(idx)
      warning ('VSDP:from_sdpa_file:readError', ...
        'Lower triangular elements will be ignored.');
      col(idx) = [];
      blk(idx) = [];
      i(idx) = [];
      j(idx) = [];
      data(idx) = [];
    end
    
    % block information
    idx = blksize>1;  % index vector for sdp blocks
    dims = abs(blksize(~idx));  % dimensions of lp blocks
    
    % cone structure
    K.l = sum(dims);  % number of all linear variables
    K.s = reshape(blksize(idx),[],1);  % dimensions of sdp blocks
    
    % calculate offset positions
    ofs = zeros(1,length(blksize));  % block offsets
    ofs(~idx) = cumsum([0 dims(1:end-1)]);  % linear block offsets
    ofs(idx) = cumsum([K.l; K.s(1:end-1).*(K.s(1:end-1)+1)/2]);
    % offsets = i + block-offsets + column-offsets
    idx = 0.5*idx';  % factor 0.5 for sdp part column offset:
    i = i + ofs(blk) + idx(blk).*j.*(j-1);  % overwrite i by offsets
    
    % create output
    dim3 = K.l + (K.s'*(K.s+1))/2;
    j = col>0;  % overwrite j by index vector for A
    A = sparse(i(j),col(j),data(j),dim3,m);
    c = sparse(i(~j),1,data(~j),dim3,1);
    
  case '.dat'  % Dense problem data
    
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
    
    % block information
    idx = blksize>1;  % index vector for sdp blocks
    dims = abs(blksize(~idx));  % dimensions of lp blocks
    
    % cone structure
    K.l = sum(dims);  % number of all linear variables
    K.s = reshape(blksize(idx),[],1);  % dimensions of sdp blocks
    
    dim = K.l + K.s'*K.s;
    
    % read data of A and c
    [c,cnt] = fscanf(fid,'%*[^0-9+-]%f',[dim 1]);
    if cnt~=dim
      fclose(fid);
      error('VSDP:from_sdpa_file:readError', ...
        'from_sdpa_file: Could not read SDPA data');
    end
    [A,cnt] = fscanf(fid,'%*[^0-9+-]%f',[dim m]);
    if cnt~=dim*m
      fclose(fid);
      error('VSDP:from_sdpa_file:readError', ...
        'from_sdpa_file: Could not read SDPA data');
    end
    fclose(fid);
    
    % sort data: lp first
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
    
  case '.ini-s'  % Sparse initial data.
    
    % Cone information from blksize
    if nargin < 2 || isempty(blksize)
      error ('VSDP:from_sdpa_file:readError', ...
        ['from_sdpa_file: Block sizes "blksize" has to be set to read ', ...
        'initial data']);
    end
    
    % dual initial vector y0
    str = fgetl(fid);
    b = sscanf(str,'%*[^0-9+-]%lf',[inf 1]);
    
    % read data of x0 and z0
    [data,cnt] = fscanf(fid,'%*[^0-9+-]%d %d %d %d %lf',[5 inf]);
    fclose(fid);
    if cnt==0 || mod(cnt,5)~=0
      error ('VSDP:from_sdpa_file:readError', ...
        'from_sdpa_file: Could not read SDPA data');
    end
    [col, blk, i, j, data] = deal(data(1,:), data(2,:), data(3,:), ...
      data(4,:), data(5,:));
    
    % Any i > j?
    idx = i>j;
    if any(idx)
      warning ('VSDP:from_sdpa_file:readError', ...
        'from_sdpa_file: Lower triangular elements will be ignored.');
      col(idx) = [];
      blk(idx) = [];
      i(idx) = [];
      j(idx) = [];
      data(idx) = [];
    end
    
    % block information
    idx = blksize>1;  % index vector for sdp blocks
    dims = abs(blksize(~idx));  % dimensions of lp blocks
    
    % cone structure
    K.l = sum(dims);  % number of all linear variables
    K.s = reshape(blksize(idx),[],1);  % dimensions of sdp blocks
    
    % calculate offset positions
    ofs = zeros(1,length(blksize));  % block offsets
    ofs(~idx) = cumsum([0 dims(1:end-1)]);  % linear block offsets
    ofs(idx) = cumsum([K.l; K.s(1:end-1).*(K.s(1:end-1)+1)/2]);
    % offsets = i + block-offsets + column-offsets
    idx = 0.5*idx';  % factor 0.5 for sdp part column offset:
    blk = i + ofs(blk) + idx(blk).*j.*(j-1);  % overwrite blk by offsets
    
    % create output
    dim3 = K.l + (K.s'*(K.s+1))/2;
    col = col<2;  % overwrite col by index vector for z0
    A = sparse(blk(col),1,data(col),dim3,1);  % z0
    data = data .* (1 + (i<j));  % mu=2, for x0
    c = sparse(blk(~col),1,data(~col),dim3,1);  % x0
    
  case '.ini'  % Dense initial data
    
    % cone information from blksize
    if (nargin<2 || isempty(blksize))
      error ('VSDP:from_sdpa_file:readError', ...
        'block sizes "blksize" has to be set to read initial data');
    end
    
    % dual initial vector y0
    str = fgetl(fid);
    b = -sscanf(str,'%*[^0-9+-]%lf',[inf 1]);
    
    % block information
    idx = blksize>1;  % index vector for sdp blocks
    dims = abs(blksize(~idx));  % dimensions of lp blocks
    
    % cone structure
    K.l = sum(dims);  % number of all linear variables
    K.s = reshape(blksize(idx),[],1);  % dimensions of sdp blocks
    
    dim = K.l + K.s'*K.s;
    
    % read data of z0 and x0
    [c,cnt] = fscanf(fid,'%*[^0-9+-]%f',[dim 1]);
    if (cnt~=dim)
      fclose(fid);
      error ('VSDP:from_sdpa_file:readError', ...
        'Could not read SDPA data');
    end
    [A,cnt] = fscanf(fid,'%*[^0-9+-]%f',[dim 1]);
    if (cnt~=dim)
      fclose(fid);
      error ('VSDP:from_sdpa_file:readError', ...
        'Could not read SDPA data');
    end
    fclose(fid);
    
    % sort data: lp first
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
    
  otherwise
    error ('VSDP:from_sdpa_file:unsupported_file_extension', ...
      'from_sdpa_file: Unsupported file extension ''%s''.', fext);
end

obj = vsdp (A, b, c, K);

end
