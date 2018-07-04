function disp (obj)
% DISP  Display a VSDP object.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

fprintf('VSDP Conic programming problem\n\n');
[At, b, c, x, y, z] = deal (obj.At, obj.b, obj.c, obj.x, obj.y, obj.z);
S = whos ('At', 'b', 'c', 'x', 'y', 'z');
clear ('At', 'b', 'c', 'x', 'y', 'z');
space = repmat (' ', 6, 1);
names = sprintf('  %2s: \n', S.name);
names = strsplit (names, '\n');
names = char (names(1:end-1)');
sizes = [S.size];
sizes = [sizes(1:2:end)', sizes(2:2:end)'];
digits1 = num2str (ceil (log10 (max(sizes(:,1)))) + 1);
digits2 = num2str (ceil (log10 (max(sizes(:,2)))) + 1);
sizes = sprintf(['[%', digits1, 'd x %', digits2, 'd]\n'], sizes');
sizes = strsplit (sizes, '\n');
sizes = char (sizes(1:end-1)');
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
disp ([names, sizes, space, space, sparsity, space, classes, ...
  space, space, bytes, space, byte_size]);
fprintf('\n\n  Cone structure:\n\n');
vsep = ('+||||+')';
cone_str = [space, space, space, space, space, vsep];
if (obj.K.f > 0)
  l{2,1} = '     free';
  l{3,1} = ' ';
  l{4,1} = sprintf('   K.f = %d ', obj.K.f);
  l{5,1} = sprintf(' Range = %d:%d ', obj.K.idx.f);
  l{6,1} = repmat ('-', 1, length (l{5,1}));
  l{1,1} = l{6,1};
  cone_str = [cone_str, char(l), vsep];
end
if (obj.K.l > 0)
  l{2,1} = '     linear';
  l{3,1} = ' ';
  l{4,1} = sprintf('   K.l = %d ', obj.K.l);
  l{5,1} = sprintf(' Range = %d:%d ', obj.K.idx.l);
  l{6,1} = repmat ('-', 1, length (l{5,1}));
  l{1,1} = l{6,1};
  cone_str = [cone_str, char(l), vsep];
end
disp(cone_str);

end
