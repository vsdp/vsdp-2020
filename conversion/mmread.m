function  [A,rows,cols,entries,rep,field,symm] = mmread (filename, cache)
% MMREAD  Read Matrix from MatrixMarket collection.
%
%    MMREAD(matrix_name) Reads the MatrixMarket matrix named 'matrix_name'
%       into the current working directory.  A comprehensive list of available
%       matrices is given at: http://math.nist.gov/MatrixMarket/matrices.html.
%       If 'matrix_name' is already in the current working directory, no
%       internet access is performed.  Additionally 'matrix_name' can also be a
%       local file in MatrixMarket format, that is not compressed.
%
%    [A,rows,cols,entries,rep,field,symm] = MMREAD(matrix_name) Same as before,
%       but reads the MatrixMarket matrix into the workspace, where 'A' is the
%       resulting matrix.  'A' will be either sparse or full, depending on the
%       MatrixMarket format indicated by 'coordinate' (coordinate sparse
%       storage), or 'array' (dense array storage).  The data will be
%       duplicated as appropriate if symmetry is indicated in the header.
%
%       Optionally, size information about the matrix can be obtained by using
%       the return values 'rows', 'cols', and 'entries', where 'entries' is the
%       number of nonzero entries in the final matrix.  Type information can
%       also be retrieved using the optional return values 'rep'
%       (representation), 'field', and 'symm' (symmetry).
%
%    [...] = MMREAD(..., cache) Instead of the current working directory, an
%       arbitrary path 'cache' can be set to dowload the matrices to.
%
%   This file is a modified version of the mmread function taken from
%   http://math.nist.gov/MatrixMarket/mmio/matlab/mmread.m.
%
%   Example:
%
%       A = mmread('494_bus');
%

% Copyright 2016-2017 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

if (nargin == 2)
  cache = char(cache);  % Ensure input type
  if (~exist(cache,'dir'))
    warning('Create cache directory');
    mkdir(cache);
  end
else
  cache = pwd ();
end

filename = char(filename);  % Ensure input type
if (~(exist(filename,'file') == 2))
  % Identify that matrix is from MatrixMarket
  if (~exist(fullfile(cache,'MMmatrices.mat'),'file'))
    matrices = getMMmatrices ();
    if (nargin == 2)
      save(fullfile(cache,'MMmatrices.mat'),'matrices');
    end
  else
    load(fullfile(cache,'MMmatrices.mat'));
  end
  matIdx = find(strcmpi(matrices(:,3),filename),1);
  if (isempty(matIdx))
    error('''%s'' is not a local file or part of the %s.\n', ...
      filename, 'MatrixMarket (http://math.nist.gov/MatrixMarket/)');
  end
  % Download, if not in cache
  if (~exist(fullfile(cache,[filename, '.mtx']),'file'))
    url = sprintf('ftp://math.nist.gov/pub/MatrixMarket2/%s/%s/%s.mtx.gz', ...
      matrices{matIdx,1}, matrices{matIdx,2}, filename);
    gunzip(url, cache);
  end
  filename = [filename, '.mtx'];
end
filename = fullfile(cache, filename);

% Return if no matrix is wanted.
if (nargout == 0)
  return;
end

fid = fopen(filename,'r');
if (fid == -1)
  error('Cannot open ''%s''', filename);
end

% Read header line
header = fgets(fid);
if (header == -1)
  error('Empty file.')
end
header = strsplit(header);
if ((length(header) < 5) || (~strcmp(header{1},'%%MatrixMarket')))
  error('MatrixMarket:badHeader', ...
    'Bad header line in file ''%s'' expected format is:\n\n\t%s\n\n', ...
    filename, '%%MatrixMarket matrix representation field symmetry');
end
if (~strcmpi(header{2},'matrix'))
  error('MatrixMarket:badHeader', ...
    '''%s'' seems to be a MatrixMarket %s file.\n%s\n', filename, header{2}, ...
    'This function only knows how to read MatrixMarket matrix files.');
end
rep   = lower(header{3});
field = lower(header{4});
symm  = lower(header{5});
validatestring(symm, {'symmetric','hermitian','skew-symmetric','general'}, ...
  'mmread','MatrixMarket header symmetry');

% Skip comments
commentline = fgets(fid);
while ((~isempty(commentline)) && (commentline(1) == '%'))
  commentline = fgets(fid);
end

% Read size information, then branch according to sparse or dense format
if (strcmp(rep,'coordinate'))
  % Read matrix given in sparse coordinate matrix format
  [sizeinfo,count] = sscanf(commentline,'%d%d%d');
  while (count == 0)
    commentline =  fgets(fid);
    if (commentline == -1 )
      error('End-of-file reached before size information was found.')
    end
    [sizeinfo,count] = sscanf(commentline,'%d%d%d');
    if ((count > 0) && (count ~= 3))
      error('Invalid size specification line.')
    end
  end
  rows = sizeinfo(1);
  cols = sizeinfo(2);
  entries = sizeinfo(3);
  
  if (strcmp(field,'real'))  % real valued entries
    T = fscanf(fid,'%f',3);
    T = [T; fscanf(fid,'%f')];
    if (size(T) ~= 3*entries)
      error('MatrixMarket:invalidData', '%s\n%s\n%s\n', ...
        'Data file does not contain expected amount of data.', ...
        'Check that number of data lines matches nonzero count.', ...
        'Invalid data.');
    end
    T = reshape(T,3,entries)';
    A = sparse(T(:,1), T(:,2), T(:,3), rows , cols);
  elseif (strcmp(field,'complex'))  % complex valued entries
    T = fscanf(fid,'%f',4);
    T = [T; fscanf(fid,'%f')];
    if (size(T) ~= 4*entries)
      error('MatrixMarket:invalidData', '%s\n%s\n%s\n', ...
        'Data file does not contain expected amount of data.',...
        'Check that number of data lines matches nonzero count.', ...
        'Invalid data.');
    end
    T = reshape(T,4,entries)';
    A = sparse(T(:,1), T(:,2), T(:,3) + T(:,4)*sqrt(-1), rows , cols);
  elseif (strcmp(field,'pattern'))  % pattern matrix (no values given)
    T = fscanf(fid,'%f',2);
    T = [T; fscanf(fid,'%f')];
    if (size(T) ~= 2*entries)
      error('MatrixMarket:invalidData', '%s\n%s\n%s\n', ...
        'Data file does not contain expected amount of data.',...
        'Check that number of data lines matches nonzero count.', ...
        'Invalid data.');
    end
    T = reshape(T,2,entries)';
    A = sparse(T(:,1), T(:,2), ones(entries,1) , rows , cols);
  end
  
elseif (strcmp(rep,'array'))
  % read matrix given in dense array (column major) format
  [sizeinfo,count] = sscanf(commentline,'%d%d');
  while (count == 0)
    commentline =  fgets(fid);
    if (commentline == -1)
      error('End-of-file reached before size information was found.')
    end
    [sizeinfo,count] = sscanf(commentline,'%d%d');
    if ((count > 0) && (count ~= 2))
      error('Invalid size specification line.')
    end
  end
  rows = sizeinfo(1);
  cols = sizeinfo(2);
  entries = rows*cols;
  if (strcmp(field,'real'))  % real valued entries
    A = fscanf(fid,'%f',1);
    A = [A; fscanf(fid,'%f')];
    if (strcmp(symm,'symmetric') || strcmp(symm,'hermitian') ...
        || strcmp(symm,'skew-symmetric') )
      for j = 1:cols-1
        currenti = j*rows;
        A = [A(1:currenti); zeros(j,1);A(currenti+1:length(A))];
      end
    end
    A = reshape(A,rows,cols);
  elseif (strcmp(field,'complex'))  % complex valued entries
    tmpr = fscanf(fid,'%f',1);
    tmpi = fscanf(fid,'%f',1);
    A  = tmpr+tmpi*1i;
    for j=1:entries-1
      tmpr = fscanf(fid,'%f',1);
      tmpi = fscanf(fid,'%f',1);
      A = [A; tmpr + tmpi*1i];
    end
    if (strcmp(symm,'symmetric') || strcmp(symm,'hermitian') ...
        || strcmp(symm,'skew-symmetric'))
      for j = 1:cols-1
        currenti = j*rows;
        A = [A(1:currenti); zeros(j,1);A(currenti+1:length(A))];
      end
    end
    A = reshape(A,rows,cols);
  elseif (strcmp(field,'pattern'))  % pattern (makes no sense for dense)
    error('mmread:invalidMatrixTypeSpecification', ...
      'Pattern matrix type ''%s'' invalid for array storage format.', field);
  else  % Unknown matrix type
    error('mmread:invalidMatrixTypeSpecification', ...
      'Invalid matrix type specification ''%s''. %s.\n', ...
      field, 'Check header against MatrixMarket documentation');
  end
end

% If symmetric, skew-symmetric or Hermitian, duplicate lower
% triangular part and modify entries as appropriate:
if (strcmp(symm,'symmetric') || strcmp(symm,'hermitian'))
  A = A + A' - diag(diag(A));
  entries = nnz(A);
elseif ( strcmp(symm,'skew-symmetric') )
  A = A - A';
  entries = nnz(A);
end

fclose(fid);
end

function matrices = getMMmatrices ()
% GETMMMATRICES  Create local database from MatrixMarket URL.
%
%   http://math.nist.gov/MatrixMarket/matrices.html
%

s = urlread('http://math.nist.gov/MatrixMarket/matrices.html');
s = strsplit(s,'\n');
% The URL format is "/MatrixMarket/data/(collection)/(set)/(matrix).html"
s = regexp(s,'<A HREF="\/MatrixMarket\/data\/(.*)\/(.*)\/(.*)\.html"', ...
  'tokens');
s = s(~cellfun(@isempty,s));
s = [s{:}];
s = [s{:}];
%           collection   set          matrix
matrices = [s(1:3:end)', s(2:3:end)', s(3:3:end)'];
end
