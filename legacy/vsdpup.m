function [fU,X,lb] = vsdpup(varargin)
% VSDPUP  Rigorous upper bound for the min. value of the block-diagonal problem:
%
%         min  sum(j=1:n| <  C{j}, X{j}>)
%         s.t. sum(j=1:n| <A{i,j}, X{j}>) = b(i)  for i = 1:m
%              X{j} must be positive semidefinite for j = 1:n
%
%   Moreover, a rigorous certificate of primal feasibility is provided.
%
%   [fU,X,lb] = VSDPUP(blk,A,C,b,Xt,yt,Zt)
%      The problem data (blk,A,C,b) of the block-diagonal problem is described
%      in 'mysdps.m'.  A, C, and b may be floating-point or interval quantities.
%
%      The inputs Xt,yt, and Zt are an approximate floating-point solution or
%      initial starting point for the primal (this is Xt) and the dual problem
%      (this is yt and Zt) that are computed with 'mysdps'.
%
%      The function returns:
%
%         fU   A rigorous upper bound of the minimum value for all real input
%              data (A,C,b) within the interval input data.  fU = inf, if no
%              finite rigorous upper bound can be computed.
%
%          X   =NAN and lb=NaN(n,1), if primal feasibility is not verified.
%              Otherwise, X is an interval quantity which contains a primal
%              feasible solution (certificate) of the block-diagonal problem
%              for all real input data within the interval data.
%
%         lb   An n-vector, where lb(j) is a rigorous lower bound of the
%              smallest eigenvalue of block X{j}.  lb > 0 means that all
%              symmetric matrices within the interval matrix X are rigorously
%              certified as positiv definite.  In this case, the existence of
%              strictly feasible solutions and strong duality is proved.
%
%   [...] = VSDPUP(...,yu)
%      Optionally, finite upper bounds yu in form of a nonnegative m-vector can
%      be provided.  The following dual boundedness assumption is assumed:
%      an optimal dual solution y satisfies
%
%         -yu(i) <= y(i) <= yu(i),  i = 1:m.
%
%      We recommend to use infinite bounds yu(i) instead of unreasonable large
%      bounds yu(i).  This improves the quality of the lower bound in many cases,
%      but may increase the computational time.
%
%   Example:
%
%       blk(1,:) = {'s'; 2};
%       A{1,1} = [0 1; 1 0];
%       A{2,1} = [1 1; 1 1];
%         C{1} = [1 0; 0 1];
%            b = [1; 2.0001];
%       [objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b); % Computes approximations
%       [fU,X,lb] = vsdpup(blk,A,C,b,Xt,yt,Zt);
%
%   See also mysdps.

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

obj = vsdp (varargin{1:7});
obj.options.SOLVER = 'sdpt3';
obj.options.VERBOSE_OUTPUT = false;
x = obj.solutions.initial.x;
y = obj.solutions.initial.y;
z = obj.solutions.initial.z;
obj.add_solution ('Approximate', x, y, z);
if (nargin < 8)
  obj.rigorous_upper_bound ();
else
  obj.rigorous_upper_bound (varargin{8});
end

fU = obj.solutions.rigorous_upper_bound.f_objective(2);
if ((isfinite (fU)) && (nargin < 8))
  x = obj.solutions.rigorous_upper_bound.x;
  if (~isempty (x))
    X = cell (length (obj.K.s), 1);
    x = vsdp_indexable (x, obj);
    for j = 1:length(obj.K.s)
      X{j} = vsdp.smat ([], x.s(j), 1/2);
    end
  end
  lb = obj.solutions.rigorous_upper_bound.z;
else
  X = nan;
  lb = nan;
end
end
