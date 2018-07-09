function disp (obj)
% DISP  Display short information about VSDP object.
%
%   See also vsdp.info.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% Short conic programming theory.
fprintf ('  VSDP conic programming problem in primal (P) dual (D) form:\n\n');
fprintf ('       (P)  min   c''*x          (D)  max  b''*y\n');
fprintf ('            s.t. At''*x = b           s.t. z := c - At*y\n');
fprintf ('                     x in K               z in K^*\n\n');

% Display problem dimensions.
d = num2str (numel (num2str (obj.n)));
fprintf (['    Constraints    m = length(%s.b) = %', d, 'd\n'], ...
  inputname(1), obj.m);
fprintf (['    Cone dimension n = length(%s.c) = %', d, 'd'], ...
  inputname(1), obj.n);
if (obj.n < obj.N)
  fprintf (' (condensed)');
end
fprintf ('\n\n');

% Display short cone structure.
if (obj.K.f > 0)
  fprintf('          K.f = %d\n', obj.K.f);
end
if (obj.K.l > 0)
  fprintf('          K.l = %d\n', obj.K.l);
end
if (~isempty (obj.K.q))
  fprintf('          K.q = [ ');
  fstr = ['%', num2str(numel (num2str (max (obj.K.q)))), 'd'];
  for i = 1:length (obj.K.q)
    fprintf (fstr, obj.K.q(i));
    if (i < length (obj.K.q))
      fprintf (', ');
    end
    if (mod (i, 5) == 0)
      fprintf ('\n                  ');
    end
  end
  fprintf(' ]\n');
end
if (~isempty (obj.K.s))
  fprintf('          K.s = [ ');
  fstr = ['%', num2str(numel (num2str (max (obj.K.s)))), 'd'];
  for i = 1:length (obj.K.s)
    fprintf (fstr, obj.K.s(i));
    if (i < length (obj.K.s))
      fprintf (', ');
    end
    if (mod (i, 5) == 0)
      fprintf ('\n                  ');
    end
  end
  fprintf(' ]\n');
end

% Display short VSDP reference, what can be done now.
fprintf ('\n');
fprintf ('    Common operations for VSDP problems:\n\n');
fprintf ('      - Detailed object information:           %s.info()\n', ...
  inputname(1));
fprintf (['      - Approximate solution:          ', ...
  '[x,y] = %s.solve()\n'], inputname(1));
fprintf (['      - Rigorous lower bound for (P):   ', ...
  'lbnd = %s.lower_bound()\n'], inputname(1));
fprintf (['      - Rigorous upper bound for (D):   ', ...
  'ubnd = %s.upper_bound()\n'], inputname(1));
fprintf (['      - Check if (P) is infeasible:     ', ...
  'cert = %s.check_prim_infeasible()\n'], inputname(1));
fprintf (['      - Check if (D) is infeasible:     ', ...
  'cert = %s.check_dual_infeasible()\n'], inputname(1));
fprintf ('\n');
fprintf ('    See also  help vsdp.solve\n');
fprintf ('              help vsdp.lower_bound\n');
fprintf ('              help vsdp.upper_bound\n');
fprintf ('              help vsdp.check_prim_infeasible\n');
fprintf ('              help vsdp.check_dual_infeasible\n\n');

end
