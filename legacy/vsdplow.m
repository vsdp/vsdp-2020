function [fL,Y,dl] = vsdplow(varargin)
% VSDPLOW  Rigorous upper bound for the min. value of the block-diagonal problem:
%
%         min  sum(j=1:n| <  C{j}, X{j}>)
%         s.t. sum(j=1:n| <A{i,j}, X{j}>) = b(i)  for i = 1:m
%              X{j} must be positive semidefinite for j = 1:n
%
%   Moreover, a rigorous certificate of feasibility for LMI's is provided.
%
%   [fU,X,lb] = VSDPLOW(blk,A,C,b,Xt,yt,Zt)
%      The problem data (blk,A,C,b) of the block-diagonal problem is described
%      in 'mysdps.m'.  A, C, and b may be floating-point or interval quantities.
%
%      The inputs Xt,yt, and Zt are an approximate floating-point solution or
%      initial starting point for the primal (this is Xt) and the dual problem
%      (this is yt and Zt) that are computed with 'mysdps'.
%
%      The function returns:
%
%         fL   A rigorous lower bound of the minimum value for all real input
%              data (A,C,b) within the interval input data.  fL = -inf, if no
%              finite rigorous lower bound can be computed.
%
%          Y   =NAN and dl=NaN(n,1), if dual feasibility is not verified.
%              Otherwise, Y is a floating-point vector which is for all
%              real input data (C,A,b) within the interval input data a dual
%              feasible solution (i.e. LMI certificate of feasibility)
%
%         dl   is a n-vector, where dl(j)is a rigorous lower bound of
%              the smallest eigenvalue of C(j)-sum(i=1:m| Y(i)*A{j,i}).
%              dl > 0 implies the existence of strictly dual feasible solutions
%              and strong duality.
%
%   [...] = VSDPLOW(...,xu)
%      Optionally, finite upper bounds xu in form of a nonnegative n-vector can
%      be provided for the case that upper bounds xu(j), j = 1:n, for the
%      maximal eigenvalues of some optimal X(j) are known.
%
%      We recommend to use infinite bounds xu(j) instead of unreasonable large
%      bounds xu(j).  This improves the quality of the lower bound in many
%      cases, but may increase the computational time.
%
%   Example:
%
%       blk(1,:) = {'s'; 2};
%       A{1,1} = [0 1; 1 0];
%       A{2,1} = [1 1; 1 1];
%         C{1} = [1 0; 0 1];
%            b = [1; 2.0001];
%       [objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b); % Computes approximations
%       [fL,Y,dl] = vsdplow(blk,A,C,b,Xt,yt,Zt);
%
%   See also mysdps.
%

% Copyright 2004-2020 Christian Jansson (jansson@tuhh.de)

obj = vsdp (varargin{1:7});
obj.options.SOLVER = 'sdpt3';
obj.options.VERBOSE_OUTPUT = false;
x = obj.solutions.initial.x;
y = obj.solutions.initial.y;
z = obj.solutions.initial.z;
obj.add_solution ('Approximate', x, y, z);
if (nargin < 8)
  obj.rigorous_lower_bound ();
else
  obj.rigorous_lower_bound (varargin{8});
end

fL = obj.solutions.rigorous_lower_bound.f_objective(1);
if (isfinite (fL))
  Y = obj.solutions.rigorous_lower_bound.y;
  dl = obj.solutions.rigorous_lower_bound.z;
else
  Y = nan;
  dl = nan;
end

end
