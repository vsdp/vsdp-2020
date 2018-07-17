function info (obj)
% INFO  Display detailed information about VSDP object.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

fprintf ('\n\n');
fprintf ('  VSDP Conic programming problem\n');
fprintf ('  ------------------------------\n\n');

dig1 = num2str(digits(obj.N));
fprintf (['    m = %', dig1,'d  [Number of constraints]\n'], obj.m);
fprintf (['    N = %', dig1,'d  [Uncondensed cone dimension]\n'], obj.N);
fprintf (['    n = %', dig1,'d  [  Condensed cone dimension]\n'], obj.n);

[At, b, c, x, y, z] = deal (obj.At, obj.b, obj.c, obj.x, obj.y, obj.z);
S = whos ('At', 'b', 'c', 'x', 'y', 'z');
clear ('At', 'b', 'c', 'x', 'y', 'z');

fprintf ('\n\n  Parameter details:\n\n');
% Display details about variables.
names = sprintf('%2s:\n', S.name);
names = strsplit (names, '\n');
names = char (names(1:end-1)');
sizes = [S.size];
sizes = [sizes(1:2:end)', sizes(2:2:end)'];
sizes = as_dim (sizes);
classes = char ({S.class}');
sparsity = {' ', 'sparse'};
sparsity = char (sparsity(1 + [S.sparse]')');
byte_size = repmat ({'Bytes'}, 6, 1);
bytes = [S.bytes]';
byte_size(bytes > 1024^2) = {'MB'};
bytes(bytes > 1024^2) = bytes(bytes > 1024^2) ./ 1024^2;
byte_size(bytes > 1024^1) = {'KB'};
bytes(bytes > 1024^1) = bytes(bytes > 1024^1) ./ 1024^1;
bytes = num2str (bytes, '%.1f');
byte_size = char (byte_size);
space = repmat (' ', length (S), 1);
disp ([space, space, space, names, ...
  space, sizes, space, space, sparsity, space, classes, ...
  space, space, bytes, space, byte_size]);


fprintf (['\n\n  Cone structure of ''K''', ...
  '  [dimensions, condensed dims., index ranges]:\n\n']);
% Cone types
ctype = [ ...
  repmat('f', obj.K.f > 0, 1);
  repmat('l', obj.K.l > 0, 1);
  repmat('q', sum(obj.K.q > 0), 1);
  repmat('s', sum(obj.K.s > 0), 1)];
idxs = [obj.K.idx.f; obj.K.idx.l; obj.K.idx.q; obj.K.idx.s];
% Uncondensed cone dimensions.
dims = [obj.K.f(obj.K.f > 0), ones(sum (obj.K.f > 0), 1);
  obj.K.l(obj.K.l > 0), ones(sum (obj.K.l > 0), 1);
  obj.K.q(obj.K.q > 0), ones(sum (obj.K.q > 0), 1);
  obj.K.s(obj.K.s > 0), obj.K.s(obj.K.s > 0)];
% Condensed cone dimensions.
cdims = dims;
cdims(ctype == 's', :) = [...
  cdims(ctype == 's', 1) .* (cdims(ctype == 's', 1) + 1) / 2, ...
  ones(sum (obj.K.s > 0), 1)];
% Numbers of the individual cones, 1,2,3,...
cnums = [ ...
  ones(obj.K.f > 0);
  ones(obj.K.l > 0);
  (1:length(obj.K.q))';
  (1:length(obj.K.s))'];
% Cone strings: "K.s(12)"
cstrs = [num2cell(ctype)'; num2cell(cnums)'];
cstrs = sprintf(['K.%c(%', num2str(digits (max (length (obj.K.q), ...
  length (obj.K.s)))), 'd) =\n'], cstrs{:});
cstrs = strsplit (cstrs, '\n');
cstrs = char (cstrs(1:end-1)');

space = repmat (' ', size (idxs, 1), 1);

outstr = [space, space, space, space, cstrs, ...
  space, space, as_dim(dims), ...
  space, space, space, as_dim(cdims), ...
  space, space, space, as_range(idxs)];

% Add spaces for quadratic and semidefinite cone
if (~isempty (obj.K.s))
  idx = find (ctype == 's', 1);
  outstr = [outstr(1:(idx-1),:);
    repmat(' ', 1, size(outstr, 2));
    outstr(idx:end,:)];
end

if (~isempty (obj.K.q))
  idx = find (ctype == 'q', 1);
  outstr = [outstr(1:(idx-1),:);
    repmat(' ', 1, size(outstr, 2));
    outstr(idx:end,:)];
end

disp (outstr);
disp (' ')
end

function d = digits (vec)
d = numel (num2str (max (vec)));
end

function str = as_dim (d)
str = sprintf(['[ %', num2str(digits (d(:,1))), 'd x %', ...
  num2str(digits (d(:,2))), 'd ]\n'], d');
str = strsplit (str, '\n');
str = char (str(1:end-1)');
end

function str = as_range (d)
str = sprintf(['[ %', num2str(digits (d(:,1))), 'd:%-', ...
  num2str(digits (d(:,2))), 'd ]\n'], d');
str = strsplit (str, '\n');
str = char (str(1:end-1)');
end
