function [obj, pd] = from_mps_file (fname)
% FROM_MPS_FILE  Import LP problem data from MPS file.
%
%   [obj, pd] = vsdp.FROM_MPS_FILE ('C:\path\to\problem.mps');
%
%   The output is identical to from_lp_solve_fmt.
%
%   For further references on the MPS-format, see:
%
%      [1] http://plato.asu.edu/ftp/mps_format.txt
%      [2] http://lpsolve.sourceforge.net/5.5/mps-format.htm
%
%   See also from_lp_solve_fmt.
%

% Copyright 2004-2020 Christian Jansson (jansson@tuhh.de)

narginchk(1, 1);
if (exist (fname, 'file') ~= 2)
  error ('VSDP:FROM_MPS_FILE:file_not_exists', ...
    'from_mps_file: Input MPS-file does not exist.');
end
problem = read_mps (fname);
if (isempty (problem) || ~all (isfield (problem, {'A', 'rowtypes'})))
  error ('VSDP:FROM_MPS_FILE:import_error', ...
    'insufficient input data');
end

% The MPS file format distinguishes the following row types:
%
%   'N'        objective / no restriction
%   'E'  '=='  equality
%   'L'  '<='  less than or equal
%   'G'  '>='  greater than or equal

% read A and c.
idx = (problem.rowtypes ~= 'N');
c = -sum (problem.A(~idx,:), 1)';  % (-) because of maximization
A = problem.A(idx,:);


% read e
e = problem.rowtypes(idx);
e = (e == 'G') - (e == 'L');  % Create +1 for `>` and -1 for `<`.


% read b
if isfield(problem,'rhs') && ~isempty(problem.rhs)
  b = problem.rhs(idx,1);
else
  b = zeros(size(A,1),1);
end
if isfield(problem,'ranges') && ~isempty(problem.ranges)
  ranges = problem.ranges(idx);
  % set equalities: ranges==0
  e(ranges==0) = 0;
  % find two-sided inequalities
  idx = find(ranges~=0 & ~isinf(ranges));
  ranges = ranges(idx);
  rhsH = b(idx);
  eR = e(idx);
  rhsL = rhsH - (eR<0).*abs(ranges) + (eR==0).*min(ranges,0);
  rhsH = rhsH + (eR>0).*abs(ranges) + (eR==0).*max(ranges,0);
  % extend problem by ranges
  b(idx) = rhsL;
  b = [b; rhsH];
  e(idx) = 1;
  e = [e; -ones(size(rhsH))];
  A = [A; A(idx,:)];
end


% read lb and ub
if isfield(problem,'lowerbounds') && ~isempty(problem.lowerbounds)
  lb = problem.lowerbounds(:,1);
else
  lb = zeros(size(A,2),1);
end
if isfield(problem,'upperbounds') && ~isempty(problem.upperbounds)
  ub = problem.upperbounds(:,1);
else
  ub = inf(size(A,2),1);
end

[obj, pd] = vsdp.from_lp_solve_fmt (A, b, c, e, lb, ub);
end

function problem = read_mps(filename)
% READ_MPS  Read a file in MPS format and creates a problem structure
%
%   @output  problem with the fields:
%       - name:  name of the problem
%       - objsense:  'MIN','MINIMIZE','MAX','MAXIMIZE'
%       - objname:  name of the objective function rw
%       - refrow:  name of reference row for SOS
%       - rownames:  rowcnt x 1 cell array of row names
%       - rowtypes:  rowcnt x 1 vector of row types ('L','G','N','E')
%       - columnnames:  colcnt x 1 cell array of column names
%       - rhsnames:  cell array of names of right hand sides
%       - rangenames:  cell array of names of ranges
%       - boundnames:  cell array of names of bounds
%       - rhs:  rowcnt x length(rhsnames) matrix of right hand sides
%       - ranges:  rowcnt x length(rangenames) matrix of ranges
%       - lowerbounds:  colcnt x length(boundnames) matrix of lower bounds
%       - upperbounds:  colcnt x length(boundnames) matrix of upper bounds
%       - bintflags:  colcnt x length(boundnames) matrix of bintflags
%                  bintflag(i,j)=1 for every integer column i in boundset j
%       - intflags:  numbers of interger columns
%       - sosflags:  cell of SOS sets,  sosflags{1}=[2 i j k] means that
%                  columns i, j and k are in an SOS2 set
%

% Copyright 2004-2020 Christian Jansson (jansson@tuhh.de)

% initialize default output
problem.name = '';
problem.objsense = 'MINIMIZE';
problem.objname = '';
problem.refrow = '';

% initialize variables
[acnt,rowcnt,colcnt,bndcnt,rngcnt,intflag,setsos] = ...
  deal(0,0,0,0,0,0,0);
premem = 100;

% prepare internal variables for output data, preallocate memory
A = zeros(3,2*premem);
colnames = cell(premem,1);
rownames = cell(premem,1);
rowtypes = zeros(premem,1);
rhsnames = cell(0);
rngnames = cell(0);
bndnames = cell(0);
rhs = [];
ranges = [];
intflags = [];
sosflags = cell(0);
lbnds = sparse([]);
ubnds = [];
bintflags = sparse([]);

% open file
fid = fopen(filename,'r');
if fid==-1
  error('READ_MPS: cannot open file');
end

% first line.
[line,f,lf] = getfields(fid);


% main loop
while 1
  
  % check for missing section card.
  if line(1)==' '
    warning('VSDP:READ_MPS','Missing Section Card');
    [line,f,lf] = getfields(fid);
    continue;
  end
  
  % handle simple section cards
  section = find(strcmp(f{1},{'ENDATA','NAME','OBJSENSE','REFROW','OBJNAME'}),1);
  
  if ~isempty(section)
    switch section
      case 1  % ENDATA
        problem.rhsnames = rhsnames;
        problem.rhs = rhs;
        problem.rangenames = rngnames;
        problem.ranges = ranges;
        problem.intflags = intflags;
        problem.sosflags = sosflags;
        problem.boundnames = bndnames;
        problem.lowerbounds = lbnds;
        problem.upperbounds = ubnds;
        problem.bintflags = bintflags;
        return;
      case 2  % NAME
        if lf>1
          problem.name = f{2};
        end
      case 3  % OBJSENSE
        [~,f] = getfields(fid);
        problem.objsense = f{1};
      case 4  % REFROW
        [~,f] = getfields(fid);
        problem.refrow = f{1};
      case 5  % OBJNAME
        [~,f] = getfields(fid);
        problem.objname = f{1};
    end
    
    [line,f,lf] = getfields(fid);
    continue;
  end
  
  
  % Handle the ROWS Section.
  if strcmp(f{1},'ROWS')
    [line,f,lf] = getfields(fid);
    precnt = premem;  % preallocation count
    while line(1)==' '
      rowcnt = rowcnt + 1;
      rownames{rowcnt} = f{2};
      rowtypes(rowcnt) = upper(f{1}(1));
      if rowcnt==precnt  % preallocate memory
        precnt = 2 * precnt;
        rownames{precnt} = '';
        rowtypes(precnt) = 0;
      end
      [line,f,lf] = getfields(fid);
    end
    problem.rownames = rownames(1:rowcnt);
    problem.rowtypes = char(rowtypes(1:rowcnt));
    
    % create row hash-table
    rown = 2 * rowcnt;
    rowtbl = cell(rown,1);
    for i = 1:rowcnt
      rowtbl{hash(rownames{i},rown)}(end+1) = i;
    end
    
    continue;  % move on to the next section
  end
  
  
  % handle COLUMNS section.
  if strcmp(f{1},'COLUMNS')
    colname = '';
    [line,f,lf] = getfields(fid);
    precnt = 1;  % preallocation count
    preacnt = 1;  % preallocation count for A
    while line(1)==' '
      olf = lf ~= 2*floor(0.5*lf);
      
      % check number of fields
      if lf<2 || lf>5
        error('READ_MPS: wrong number of fields in record');
      end
      
      % check for 'MARKER' records
      idx = 0;
      if lf>2 && strcmp('''MARKER''',f{2})
        idx = 3;
      elseif lf>3 && strcmp('''MARKER''',f{3})
        idx = 4;
      end
      if idx>0
        intflag = (intflag + ...
          strcmp(f{idx},{'''INTORG''','''INTEND'''})*[1;-1]) > 0;
        if strcmp(f{idx},'''SOSORG''')  % SOSORG -> find sosflag
          setsos = true;
          sosflags{end+1}(1) = 1;
          if f{1}(1)=='S'
            sosflags{end}(1) = str2double(f{1}(2:end));
          end
        elseif strcmp(f{idx},'''SOSEND''')
          setsos = false;
        end
        
        [line,f,lf] = getfields(fid);
        continue;
      end
      
      % new column
      if olf && ~strcmp(f{1},colname)
        colname = f{1};
        colcnt = colcnt + 1;
        colnames{colcnt} = f{1};
        if intflag
          intflags{end}(end+1) = colcnt;
        end
        if setsos
          sosflags{end}(end+1) = colcnt;
        end
        if colcnt==precnt  % preallocate memory
          precnt = 2 * precnt;
          colnames{precnt} = '';
        end
      end
      
      % read entries in record
      for i = 1+olf:2:lf-1
        rowno = rowtbl{hash(f{i},rown)};
        if length(rowno)>1
          rowno = rowno(strcmp(f{i},rownames(rowno)));
        end
        if length(rowno)~=1
          error('READ_MPS: COLUMNS entry specifies no unique row');
        end
        acnt = acnt + 1;
        A(:,acnt) = [rowno; colcnt; str2double(f{i+1})];
        if acnt==preacnt  % preallocation
          preacnt = 2 * preacnt;
          A(3,preacnt) = 0;
        end
      end
      
      [line,f,lf] = getfields(fid);  % move on to the next record
    end
    
    % all columns are read, add to problem
    problem.columnnames = colnames(1:colcnt);
    A = A(:,1:acnt)';
    problem.A = sparse(A(:,1),A(:,2),A(:,3),rowcnt,colcnt);
    clear A;
    
    % create column hash-table
    coln = 2 * colcnt;
    coltbl = cell(coln,1);
    for i = 1:colcnt
      coltbl{hash(colnames{i},coln)}(end+1) = i;
    end
    
    continue;  % to the next section
  end
  
  
  % handle RHS section
  if strcmp(f{1},'RHS')
    [line,f,lf] = getfields(fid);
    
    % initial rhs
    rhsname = 'DRHS';  % a default RHS name
    if lf ~= 2*floor(0.5*lf)
      rhsname = f{1};
    end
    rhscnt = 1;
    rhsno = 1;
    rhsnames{rhscnt} = rhsname;
    rhs = zeros(rowcnt,rhscnt);
    
    % read rhs
    while line(1)==' '
      olf = lf ~= 2*floor(0.5*lf);
      
      if olf && ~strcmp(f{1},rhsname)  % new rhs
        rhsname = f{1};
        rhsno = find(strcmp(rhsname,rhsnames),1);
        if isempty(rhsno)
          rhscnt = rhscnt + 1;
          rhsnames{rhscnt} = rhsname;
          rhsno = rhscnt;
          rhs(rowcnt,rhscnt) = 0;
        end
      end
      
      for i = 1+olf:2:lf-1  % read entries
        rowno = rowtbl{hash(f{i},rown)};
        if length(rowno)>1
          rowno = rowno(strcmp(f{i},rownames(rowno)));
        end
        if length(rowno)~=1
          error('READ_MPS: RHS entry specifies row that does not exist');
        end
        rhs(rowno,rhsno) = str2double(f{i+1});
      end
      
      [line,f,lf] = getfields(fid);  % move on to the next record
    end
    
    continue;  % to the next section
  end
  
  
  % handle the RANGES section.
  if strcmp(f{1},'RANGES')
    [line,f,lf] = getfields(fid);
    rngname = 'DRNG';
    % default range values for N, E, L, G constraints
    rngVAL = zeros(rowcnt,1);
    rngVAL(rowtypes=='L' | rowtypes=='G') = Inf;
    while line(1)==' '
      olf = lf ~= 2*floor(0.5*lf);
      
      if isempty(rngnames) || (olf && ~strcmp(f{1},rngname))  % new range
        rngname = f{1};
        rngno = find(strcmp(rngname,rngnames),1);
        if isempty(rngno)
          rngcnt = rngcnt + 1;
          rngnames{rngcnt} = rngname;
          rngno = rngcnt;
          ranges(:,rngno) = rngVAL;  % initial range values
        end
      end
      
      for i = 1+olf:2:lf-1
        rowno = rowtbl{hash(f{i},rown)};
        if length(rowno)>1
          rowno = rowno(strcmp(f{i},rownames(rowno)));
        end
        if length(rowno)~=1
          error('READ_MPS: RANGES entry specifies row that does not exist');
        end
        ranges(rowno,rngno) = str2double(f{i+1});
      end
      
      [line,f,lf] = getfields(fid);  % move on to the next record
    end
    
    continue;  % to the next section
  end
  
  
  % handle the BOUNDS section.
  if strcmp(f{1},'BOUNDS')
    [line,f,lf] = getfields(fid);
    bndname = 'DBOUND';
    lastbndname = '';
    while line(1)==' '
      
      switch lf
        case 4
          bndname = f{2};
          colname = f{3};
          entry = str2double(f{4});
        case 2
          % only legal if no boundname is given and type is FR or BV
          if any(strcmpi(f{1},{'BV','FR'}))
            colname = f{2};
            entry = 0;
          else
            disp(line);
            error('READ_MPS: invalid entry in BOUNDS section.');
          end
        case 3
          % a BV, FR, MI or PL bound with no entry
          if any(strcmpi(f{1},{'BV','FR','MI','PL'}))
            bndname = f{2};
            colname = f{3};
            entry = 0;
          else  % no bound name given
            colname = f{2};
            entry = str2double(f{3});
          end
        otherwise
          disp(line);
          error('READ_MPS: wrong number of fields in BOUNDS record');
      end
      
      if ~strcmp(bndname,lastbndname)  % new bounds
        lastbndname = bndname;
        bndno = find(strcmp(bndname,bndnames),1);
        if isempty(bndno)
          bndcnt = bndcnt + 1;
          bndnames{bndcnt} = bndname;
          bndno = bndcnt;
          % default bounds
          lbnds(colcnt,bndno) = 0;
          ubnds(1:colcnt,bndno) = Inf;
          ubnds(intflags,bndno) = 1;
          bintflags(colcnt,bndno) = 0;
        end
      end
      
      colno = coltbl{hash(colname,coln)};
      if length(colno)>1
        colno = colno(strcmp(colname,colnames(colno)));
      end
      try
        if length(colno)~=1
          disp(line);
          error('READ_MPS: invalid column in BOUNDS section');
        end
      catch
        disp('hallo');
      end
      
      % handle bound types
      btype = strfind(['LO','UP','FX','FR','MI','PL','LI',...
        'UI','BV'],upper(f{1})) / 2 - 0.5;
      
      if isempty(btype) || btype~=round(btype) || length(f{1})~=2
        disp(line);
        error('unhandled bound type');
      end
      
      switch btype
        case 0  % LO
          lbnds(colno,bndno) = entry;
        case 1  % UP
          ubnds(colno,bndno) = entry;
        case 2  % FX
          lbnds(colno,bndno) = entry;
          ubnds(colno,bndno) = entry;
        case 3  % FR
          lbnds(colno,bndno) = -Inf;
          ubnds(colno,bndno) = +Inf;
        case 4  % MI
          lbnds(colno,bndno) = -Inf;
        case 5  % PL
          ubnds(colno,bndno) = +Inf;
        case 6  % LI
          lbnds(colno,bndno) = entry;
          bintflags(colno,bndno) = 1;
        case 7  % UI
          ubnds(colno,bndno) = entry;
          bintflags(colno,bndno) = 1;
        case 8  % BV
          ubnds(colno,bndno) = 1;
          bintflags(colno,bndno) = 1;
      end
      
      [line,f,lf] = getfields(fid);  % move to next record
    end
    
    continue  % to the next section
  end
  
  % section that cannot be handled
  disp('unhandled section');
  disp(line);
  [line,f,lf] = getfields(fid);
  while line(1)==' '
    [line,f,lf] = getfields(fid);
  end
end

end



function [line,fields,lf] = getfields(fid)
% GETFIELDS
%       [line,fields,lf] = getfields(fid)
% reads the next nonempty, non-comment line of the MPS input file and
% returns the fields of this line in a cell array

while 1
  line = fgets(fid);
  if line==-1
    error('READ_MPS:GETFIELDS','unexpected end of file');
  elseif ~isempty(line) && line(1)~='*'
    fields = textscan(line,'%s');
    fields = fields{1};
    lf = length(fields);
    if lf>0
      break;
    end
  end
end
end



function h = hash(key,n)
% HASH
%       h = hash(key,n)
% computes a hash value for the 'key' string between 1 and n

h = sum(key./(27.3:27+length(key)));
h = round((h-floor(h))*(n-1)) + 1;
end
