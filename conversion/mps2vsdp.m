function [A,b,c,K,pd] = mps2vsdp(problem)
%% MPS2VSDP:  reads/transform problem from MPS format to VSDP3 format
%       [A,b,c,K,pd] = mps2vsdp(problem)
%
%% >> Input:
% problem: can be the filename of a text file in MPS format or a problem
%          structure as created by the read_mps function [see VSDP/read_mps]
%
%% >> Output: 
% A: matrix of linear equations
% b, c: - coefficients of primal/dual objective functions
% K: a structure with following fields
%     - K.f stores the number of free variables
%     - K.l is the number of nonnegative components
%     - K.q lists the lengths of socp blocks
%     - K.s lists the dimensions of semidefinite blocks
% pd: a scalar with one of 'p' or 'd'. If 'p' problem was not dualized.
%     If 'd' problem was dualized and optimal value must be negated. 
%

%% ********************************************************************* %%
%% This file is part of VSDP by V. Haerter, C. Jansson and M. Lange      %%
%% Copyright (c) 2012, C. Jansson                                        %%
%%                     Technical University of Hamburg (TUHH)            %%
%%                     Institute for Reliable Computing (IRC)            %%
%% VSDP can be freely used for private and academic purposes.            %%
%% Commercial use or use in conjunction with a commercial program which  %%
%% requires VSDP or part of it to function properly is prohibited.       %%
%% ********************************************************************* %%

%% Last modified:
% 29/08/12    M. Lange, written
%
%% ToDo:
%       - support ranges
%       - check of unsupported flags etc
%

%% check input
if ischar(problem)
    problem = read_mps(problem);
elseif isempty(problem) || ~all(isfield(problem,{'A','rowtypes'}))
    error('VSDP:MPS2VSDP','insufficient input data');
end


%% read A and c
idx = problem.rowtypes~='N';
c = -sum(problem.A(~idx,:),1)';  % (-) because of maximization
A = problem.A(idx,:);


%% read e
e = problem.rowtypes(idx);
e = (e=='G') - (e=='L');


%% read b
if isfield(problem,'rhs') && ~isempty(problem.rhs)
    b = problem.rhs(idx,1);
else
    b = zeros(size(A,1),1);
end
if isfield(problem,'ranges') && ~isempty(problem.ranges)
    ranges = problem.ranges(idx);
    % set equalities: ranges==0
    e(ranges==0) = 0;
    % find two-sided inequalities
    idx = find(ranges~=0 & ~isinf(ranges));
    ranges = ranges(idx);
    rhsH = b(idx);
    eR = e(idx);
    rhsL = rhsH - (eR<0).*abs(ranges) + (eR==0).*min(ranges,0);
    rhsH = rhsH + (eR>0).*abs(ranges) + (eR==0).*max(ranges,0);
    % extend problem by ranges
    b(idx) = rhsL;
    b = [b; rhsH];
    e(idx) = 1;
    e = [e; -ones(size(rhsH))];
    A = [A; A(idx,:)];
end


%% read lb and ub
if isfield(problem,'lowerbounds') && ~isempty(problem.lowerbounds)
    lb = problem.lowerbounds(:,1);
else
    lb = zeros(size(A,2),1);
end
if isfield(problem,'upperbounds') && ~isempty(problem.upperbounds)
    ub = problem.upperbounds(:,1);
else
    ub = inf(size(A,2),1);
end


%% transform LP format to VSDP format 
[A,b,c,K,pd] = lp2vsdp(A,b,c,e,lb,ub);

%___________________________End of MPS2VSDP____________________________