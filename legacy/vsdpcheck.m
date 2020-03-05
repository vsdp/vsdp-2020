function [m,n] = vsdpcheck(varargin)
% VSDPCHECK  Routine for checking the VSDP-2006 format.
%
%   [m,n] = VSDPCHECK(blk,A,C,b)
%      The block-diagonal format is explained in 'mysdps.m'.
%
%   [...] = VSDPCHECK(...,X0,y0,Z0) optionally the format of an initial guess
%      (X0,y0,Z0) is checked as well.
%
%   Example:
%
%       blk(1,:) = {'s'; 2};
%       A{1,1} = [0 1; 1 0];
%       A{2,1} = [1 1; 1 1];
%         C{1} = [1 0; 0 1];
%            b = [1; 2.0001];
%       [m,n] = vsdpcheck(blk,A,C,b);
%
%   See also mysdps.
%

% Copyright 2004-2020 Christian Jansson (jansson@tuhh.de)

obj = vsdp (varargin{:});
m = obj.m;
n = length (obj.K.s);

end
