function info (obj)
% INFO  Display detailed information about VSDP object.
%
%   See also vsdp.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% Short conic programming theory.
fprintf ('\n  VSDP conic programming problem in primal (P) dual (D) form:\n\n');
fprintf ('       (P)  min   c''*x          (D)  max  b''*y\n');
fprintf ('            s.t. At''*x = b           s.t. z := c - At*y\n');
fprintf ('                     x in K               z in K^*\n\n');

[At, b, c] = deal (obj.At, obj.b, obj.c);
S = whos ('At', 'b', 'c');
clear ('At', 'b', 'c');

fprintf ('  Parameter details:\n\n');
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
byte_size = repmat ({'Bytes'}, length (S), 1);
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

dig1 = num2str(digits(obj.N));
fprintf ('\n');
fprintf (['    m = %', dig1,'d  [     Number of constraints]\n'], obj.m);
fprintf (['    N = %', dig1,'d  [Uncondensed cone dimension]\n'], obj.N);
fprintf (['    n = %', dig1,'d  [  Condensed cone dimension]\n'], obj.n);

fprintf ('\n\n  Solution details:\n\n');
if (~isempty (obj.solutions.initial))
  fprintf ('\n  - Initial:\n\n');
  obj.solutions.initial.mem_info ();
end
if (~isempty (obj.solutions.approximate))
  fprintf ('\n  - obj.solutions.approximate for (P) and (D):\n\n');
  disp (obj.solutions.approximate);
  obj.solutions.approximate.mem_info ();
end
if (~isempty (obj.solutions.rigorous_lower_bound))
  fprintf ('\n  - obj.solutions.rigorous_lower_bound  fL <= c''*x   for (P):\n\n');
  disp (obj.solutions.rigorous_lower_bound);
  obj.solutions.rigorous_lower_bound.mem_info ();
end
if (~isempty (obj.solutions.rigorous_upper_bound))
  fprintf ('\n  - obj.solutions.rigorous_upper_bound  b''*y <= fU   for (D):\n\n');
  disp (obj.solutions.rigorous_upper_bound);
  obj.solutions.rigorous_upper_bound.mem_info ();
end
if (~isempty (obj.solutions.certificate_primal_infeasibility))
  fprintf ('\n  - obj.solutions.certificate_primal_infeasibility:\n\n');
  disp (obj.solutions.certificate_primal_infeasibility);
  obj.solutions.certificate_primal_infeasibility.mem_info ();
end
if (~isempty (obj.solutions.certificate_dual_infeasibility))
  fprintf ('\n  - obj.solutions.certificate_dual_infeasibility:\n\n');
  disp (obj.solutions.certificate_dual_infeasibility);
  obj.solutions.certificate_dual_infeasibility.mem_info ();
end

fprintf (['\n\n  Cone structure of ''K''', ...
  '  [dimensions, condensed dims., index ranges]:\n\n']);
% Cone types
ctype = [ ...
  repmat('f', obj.K.f > 0, 1);
  repmat('l', obj.K.l > 0, 1);
  repmat('q', sum(obj.K.q > 0), 1);
  repmat('s', sum(obj.K.s > 0), 1)];
idxs = [];
if (obj.K.f > 0)
  idxs = obj.K.idx.f;
end
if (obj.K.l > 0)
  idxs = [idxs; obj.K.idx.l];
end
if (~isempty (obj.K.q))
  idxs = [idxs; obj.K.idx.q];
end
if (~isempty (obj.K.s))
  idxs = [idxs; obj.K.idx.s];
end
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

space = repmat (' ', size (cnums, 1), 1);

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
