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
obj_name = inputname(1);
if (isempty (obj_name))
  obj_name = 'obj';
end
fprintf ('\n');
fprintf ('  %s.solutions.approximate  for (P) and (D):\n', obj_name);
if (isempty (obj.solutions.approximate))
  fprintf ('\n      None.  Compute with ''%s = %s.solve()''\n\n', ...
    obj_name, obj_name);
else
  disp (obj.solutions.approximate)
end

fprintf (['  %s.solutions.rigorous_lower_bound  fL <= c''*x   ', ...
  'for (P):\n'], obj_name);
if (isempty (obj.solutions.rigorous_lower_bound))
  fprintf ( ...
    '\n      None.  Compute with ''%s = %s.rigorous_lower_bound()''\n\n', ...
    obj_name, obj_name);
else
  disp (obj.solutions.rigorous_lower_bound)
end

fprintf (['  %s.solutions.rigorous_upper_bound  b''*y <= fU   ', ...
  'for (D):\n'], obj_name);
if (isempty (obj.solutions.rigorous_upper_bound))
  fprintf ( ...
    '\n      None.  Compute with ''%s = %s.rigorous_upper_bound()''\n\n', ...
    obj_name, obj_name);
else
  disp (obj.solutions.rigorous_upper_bound)
end

fprintf ('  %s.solutions.certificate_primal_infeasibility:\n', ...
  obj_name);
if (isempty (obj.solutions.certificate_primal_infeasibility))
  fprintf ( ...
    '\n      None.  Check with ''%s = %s.check_primal_infeasible()''\n\n', ...
    obj_name, obj_name);
else
  disp (obj.solutions.certificate_primal_infeasibility)
end

fprintf ('  %s.solutions.certificate_dual_infeasibility:\n', ...
  obj_name);
if (isempty (obj.solutions.certificate_dual_infeasibility))
  fprintf ( ...
    '\n      None.  Check with ''%s = %s.check_dual_infeasible()''\n\n', ...
    obj_name, obj_name);
else
  disp (obj.solutions.certificate_dual_infeasibility)
end

fprintf (' For more information type:  %s.info()\n\n', obj_name);

end
