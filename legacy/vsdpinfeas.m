function [isinfeas,X,Y] = vsdpinfeas(varargin)
% VSDPINFEAS  Infeasibility-check for the block-diagonal problem:
%
%         min  sum(j=1:n| <  C{j}, X{j}>)
%         s.t. sum(j=1:n| <A{i,j}, X{j}>) = b(i)  for i = 1:m
%              X{j} must be positive semidefinite for j = 1:n
%
%   as well as its dual:
%
%         max  b'*y
%         s.t. C{j} - sum(i=1:m| y{i}*A{i,j}) must be positive semidefinite
%                                             for j = 1:n.
%
%   [isinfeas,X,Y] = VSDPINFEAS(blk,A,C,b,choose)
%      The block-diagonal format (blk,A,C,b) is explained in 'mysdps.m'.
%
%         'choose'    If the character is 'p', primal infeasibility should be
%                     verified.  If 'd', dual infeasibility should be verified.
%
%      The output is:
%
%         'isinfeas'  Returns 1 if the primal or dual problem is proved to
%                     be infeasible and 0 if infeasibility cannot be verified.
%
%         'X'         Contains a rigorous certificate (improving ray) of dual
%                     infeasibility, if it is not equal to NaN.
%
%         'Y'         Is a rigorous certificate (improving ray) of primal
%                     infeasibility, if it is not equal to NaN.
%
%   VSDPINFEAS(...,Xt,yt,Zt) optionally provide already by 'mysdps' computed
%      approximate solutions (Xt,yt,Zt).  This avoids calling 'mysdps' from
%      within this function, if approximate solutions are already present.
%
%   Example:
%
%       EPS = -0.01;
%       DELTA = 0.1;
%       blk(1,:) = {'s'; 2};
%       A{1,1} = [1 0; 0 0];
%       A{2,1} = [0 1; 1 DELTA];
%         C{1} = [0 0; 0 0];
%            b = [EPS; 1];
%       [objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);
%       choose = 'p';
%       [isinfeas,X,Y] = vsdpinfeas(blk,A,C,b,choose,Xt,yt,Zt);
%
%   See also mysdps.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

obj = vsdp (varargin{1:4});
obj.options.SOLVER = 'sdpt3';
obj.options.VERBOSE_OUTPUT = false;

% If no approximate certificate is given.
if (nargin <= 5)
  obj.solve ();
end

switch (varargin{5})
  case 'p'
    if (nargin >= 7)
      obj.check_primal_infeasible (varargin{7});
    else
      obj.check_primal_infeasible ();
    end
    isinfeas = obj.solutions.certificate_primal_infeasibility.f_objective(1);
    X = nan;
    if (isinfeas)
      Y = obj.solutions.certificate_primal_infeasibility.y;
    else
      Y = nan;
    end
  case 'd'
    if (nargin >= 6)
      x = vsdp.svec (obj, vsdp.cell2mat (varargin{6}), 2);
      obj.check_dual_infeasible (x);
    else
      obj.check_dual_infeasible ();
    end
    isinfeas = obj.solutions.certificate_dual_infeasibility.f_objective(2);
    if (isinfeas)
      x = obj.solutions.certificate_dual_infeasibility.x;
      X = cell (length (obj.K.s), 1);
      x = vsdp_indexable (x, obj);
      for j = 1:length(obj.K.s)
        X{j} = vsdp.smat ([], x.s(j), 1/2);
      end
    else
      X = nan;
    end
    Y = nan;
  otherwise
    error ('Choose "p" or "d".');
end

end
