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
fprintf ( '  with dimensions  [n,m] = size(At)\n');
d = num2str (numel (num2str (obj.n)));
fprintf (['                    n    = %', d, 'd variables\n'], obj.n);
fprintf (['                      m  = %', d, 'd constraints\n\n'], obj.m);

% Display short cone structure.
if (obj.K.f > 0)
  fprintf('        K.f = %d\n', obj.K.f);
end
if (obj.K.l > 0)
  fprintf('        K.l = %d\n', obj.K.l);
end
if (~isempty (obj.K.q))
  fprintf('        K.q = [ ');
  fstr = ['%', num2str(numel (num2str (max (obj.K.q)))), 'd'];
  for i = 1:length (obj.K.q)
    fprintf (fstr, obj.K.q(i));
    if (i < length (obj.K.q))
      fprintf (', ');
    end
    if (mod (i, 5) == 0)
      fprintf ('\n                ');
    end
  end
  fprintf(' ]\n');
end
if (~isempty (obj.K.s))
  fprintf('        K.s = [ ');
  fstr = ['%', num2str(numel (num2str (max (obj.K.s)))), 'd'];
  for i = 1:length (obj.K.s)
    fprintf (fstr, obj.K.s(i));
    if (i < length (obj.K.s))
      fprintf (', ');
    end
    if (mod (i, 5) == 0)
      fprintf ('\n                ');
    end
  end
  fprintf(' ]\n');
end

% Display short VSDP reference, what can be done now.
fprintf ('\n');
fprintf ('  %s.solutions(''Approximate solution'')  for (P) and (D).\n', ...
  inputname(1));
if (isempty (obj.solutions('Approximate solution')))
  fprintf ('\n      None.  Compute with ''%s = %s.solve()''\n\n', ...
    inputname(1), inputname(1));
else
  disp (obj.solutions('Approximate solution'))
end

fprintf (['  %s.solutions(''Rigorous lower bound'')  fL <= c''*x   ', ...
  'for (P):\n'], inputname(1));
if (isempty (obj.solutions('Rigorous lower bound')))
  fprintf ( ...
    '\n      None.  Compute with ''%s = %s.rigorous_lower_bound()''\n\n', ...
    inputname(1), inputname(1));
else
  disp (obj.solutions('Rigorous lower bound'))
end

fprintf (['  %s.solutions(''Rigorous upper bound'')  b''*y <= fU   ', ...
  'for (D):\n'], inputname(1));
if (isempty (obj.solutions('Rigorous upper bound')))
  fprintf ( ...
    '\n      None.  Compute with ''%s = %s.rigorous_upper_bound()''\n\n', ...
    inputname(1), inputname(1));
else
  disp (obj.solutions('Rigorous upper bound'))
end

fprintf ('  %s.solutions(''Certificate primal infeasibility''):\n', ...
  inputname(1));
if (isempty (obj.solutions('Certificate primal infeasibility')))
  fprintf ( ...
    '\n      None.  Check with ''%s = %s.check_primal_infeasible()''\n\n', ...
    inputname(1), inputname(1));
else
  disp (obj.solutions('Certificate primal infeasibility'))
end

fprintf ('  %s.solutions(''Certificate dual infeasibility''):\n', ...
  inputname(1));
if (isempty (obj.solutions('Certificate dual infeasibility')))
  fprintf ( ...
    '\n      None.  Check with ''%s = %s.check_dual_infeasible()''\n\n', ...
    inputname(1), inputname(1));
else
  disp (obj.solutions('Certificate dual infeasibility'))
end

fprintf (' For more information type:  %s.info()\n\n', inputname(1));

end
