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
fprintf ('\n  State:\n\n');
state = {' ', 'x'};

appr_sol = (~isempty (obj.x) && ~isempty (obj.y) && ~isempty (obj.z));
fprintf ('     [%c]  Approximate solution:\n\n', state{appr_sol + 1});
if (~appr_sol)
  fprintf ('            [x,y,z] = %s.solve();\n\n', inputname(1));
else
  fprintf ('               c''*x = %.15e\n',   obj.c' * obj.x);
  fprintf ('               b''*y = %.15e\n\n', obj.b' * obj.y);
end

rig_lb = false;
fprintf ('     [%c]  Rigorous lower bound   lbnd <= c''*x   for (P):\n\n', ...
  state{rig_lb + 1});
if (~rig_lb)
  fprintf ('               lbnd = %s.rigorous_lower_bound();\n\n', ...
    inputname(1));
else
  fprintf ('               lbnd = %f\n\n', obj.lbnd);
end

rig_ub = false;
fprintf ('     [%c]  Rigorous upper bound   b''*y <= ubnd   for (D):\n\n', ...
  state{rig_ub + 1});
if (~rig_ub)
  fprintf ('               ubnd = %s.rigorous_upper_bound();\n\n', ...
    inputname(1));
else
  fprintf ('               ubnd = %f\n\n', obj.ubnd);
end

p_inf = false;
fprintf ('     [%c]  Check if (P) is infeasible:\n\n', state{p_inf + 1});
if (~p_inf)
  fprintf ('               cert = %s.check_prim_infeasible();\n\n', ...
    inputname(1));
else
end

d_inf = false;
fprintf ('     [%c]  Check if (D) is infeasible:\n\n', state{d_inf + 1});
if (~d_inf)
  fprintf ('               cert = %s.check_dual_infeasible()\n\n', inputname(1));
else
end
fprintf (' For more information type:  %s.info()\n\n', inputname(1));

end
