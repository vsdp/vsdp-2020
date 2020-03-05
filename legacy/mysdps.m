function [obj,X,y,Z,info] = mysdps(varargin)
% MYSDPS  My Semidefinite Programming solver for the block-diagonal problem:
%
%         min  sum(j=1:n| <  C{j}, X{j}>)
%         s.t. sum(j=1:n| <A{i,j}, X{j}>) = b(i)  for i = 1:m
%              X{j} must be positive semidefinite for j = 1:n
%
%   Moreover, a certificate of feasibility for LMI's is provided.  The routine
%   uses the semidefinite solver SDPT3 by default.
%
%   [obj,X,y,Z,info] = MYSDPS(blk,A,C,b)
%      The problem data of the block-diagonal structure:
%
%         'blk'  cell(n,2)
%         'A'    cell(m,n)
%         'C'    cell(n,1)
%         'b'  double(m,1)
%
%      The j-th block C{j} and the blocks A{i,j}, for i = 1:m, are real
%      symmetric matrices of common size s_j, and blk(j,:) = {'s'; s_j}.
%
%      The blocks C{j} and A{i,j} must be stored as individual matrices in
%      dense or sparse format.
%
%      The function returns:
%
%         'obj'    Primal and dual objective value [<C,X> <b,y>].
%
%         'X,y,Z'  An approximate optimal solution or a primal or dual
%                  infeasibility certificate.
%
%         'info'   Termination-code with
%                   0: indication of optimality (normal termination),
%                   1: indication of primal infeasibility,
%                   2: indication of dual infeasibility,
%                   3: SDPT3 and SDPA-M: indication of both primal and dual
%                      infeasibilities,
%                  -1: otherwise.
%
%   [...] = MYSDPS(...,X0,y0,Z0) optionally provide an initial guess (X0,y0,Z0)
%      to the approximate solver.
%
%   Example:
%
%       blk(1,:) = {'s'; 2};
%       A{1,1} = [0 1; 1 0];
%       A{2,1} = [1 1; 1 1];
%         C{1} = [1 0; 0 1];
%            b = [1; 2.0001];
%       [objt,X,y,Z,info] = mysdps(blk,A,C,b);
%
%   See also vsdpcheck, vsdplow, vsdpup, vsdpinfeas.
%

% Copyright 2004-2020 Christian Jansson (jansson@tuhh.de)

vobj = vsdp (varargin{:});
vobj.options.VERBOSE_OUTPUT = false;
vobj.solve ('sdpt3');

obj = vobj.solutions.approximate.f_objective.';
[X, y, Z] = vobj.solutions.approximate.to_2006_fmt (vobj);
switch (vobj.solutions.approximate.solver_info.termination)
  case 'Normal termination'
    info = 0;
  case 'Primal infeasible'
    info = 1;
  case 'Dual infeasible'
    info = 2;
  otherwise
    info = -1;
end

end
