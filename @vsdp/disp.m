function disp (obj)
% DISP  Display short information about VSDP object.
%
%   See also vsdp.info.

% Copyright 2004-2019 Christian Jansson (jansson@tuhh.de)

% Determine object name.
obj_name = inputname(1);
if (isempty (obj_name))
  obj_name = 'obj';
end

% Display problem dimensions.
fprintf ('  VSDP conic programming problem with dimensions:\n\n');
fprintf ('    [n,m] = size(%s.At)\n', obj_name);
d = num2str (numel (num2str (obj.n)));
fprintf (['     n    = %', d, 'd variables\n'], obj.n);
fprintf (['       m  = %', d, 'd constraints\n\n'], obj.m);

% Display short cone structure.
fprintf ('  and cones:\n\n');
if (obj.K.f > 0)
  fprintf('     K.f = %d\n', obj.K.f);
end
if (obj.K.l > 0)
  fprintf('     K.l = %d\n', obj.K.l);
end
if (~isempty (obj.K.q))
  if (length (obj.K.q) <= 20)
    fprintf('     K.q = [ ');
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
  else
    fprintf('     K.q = [ %d cones (length = %d) ]\n', length (obj.K.q), ...
      sum (obj.K.q));
  end
end
if (~isempty (obj.K.s))
  if (length (obj.K.s) <= 20)
    fprintf('     K.s = [ ');
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
  else
    fprintf('     K.s = [ %d cones (diag. sum = %d) ]\n', length (obj.K.s), ...
      sum (obj.K.s));
  end
end

% Display short VSDP reference, what can be done now.
if (isempty (obj.solutions.approximate))
  fprintf ('\n  Compute an approximate solution:\n\n');
  fprintf ('    ''%s = %s.solve()''\n', obj_name, obj_name);
else
  fprintf ('\n  %s.solutions.approximate:\n\n', obj_name);
  disp (obj.solutions.approximate)
end

if (isempty (obj.solutions.rigorous_lower_bound))
  if (~isempty (obj.solutions.approximate))
    fprintf ('\n  Compute a rigorous lower bound:\n\n');
    fprintf ('    ''%s = %s.rigorous_lower_bound()''\n', obj_name, obj_name);
  end
else
  fprintf ('\n  %s.solutions.rigorous_lower_bound:\n\n', obj_name);
  disp (obj.solutions.rigorous_lower_bound)
end

if (isempty (obj.solutions.rigorous_upper_bound))
  if (~isempty (obj.solutions.approximate))
    fprintf ('\n  Compute a rigorous upper bound:\n\n');
    fprintf ('    ''%s = %s.rigorous_upper_bound()''\n', obj_name, obj_name);
  end
else
  fprintf ('  %s.solutions.rigorous_upper_bound:\n\n', obj_name);
  disp (obj.solutions.rigorous_upper_bound)
end

if (isempty (obj.solutions.certificate_primal_infeasibility))
  if (~isempty (obj.solutions.rigorous_upper_bound) ...
      && isinf (obj.solutions.rigorous_upper_bound.f_objective(2)))
    fprintf ('\n  The rigorous upper bound is infinite, ');
    fprintf ('check primal infeasibility:\n\n');
    fprintf ('    ''%s = %s.check_primal_infeasible()''\n', obj_name, obj_name);
  end
else
  fprintf ('\n  %s.solutions.certificate_primal_infeasibility:\n\n', obj_name);
  disp (obj.solutions.certificate_primal_infeasibility)
end

if (isempty (obj.solutions.certificate_dual_infeasibility))
  if (~isempty (obj.solutions.rigorous_lower_bound) ...
      && isinf (obj.solutions.rigorous_lower_bound.f_objective(1)))
    fprintf ('\n  The rigorous lower bound is infinite, ');
    fprintf ('check dual infeasibility:\n\n');
    fprintf ('    ''%s = %s.check_dual_infeasible()''\n', obj_name, obj_name);
  end
else
  fprintf ('\n  %s.solutions.certificate_dual_infeasibility:\n\n', obj_name);
  disp (obj.solutions.certificate_dual_infeasibility)
end

fprintf ('\n\n  Detailed information:  ''%s.info()''\n\n', obj_name);

end
