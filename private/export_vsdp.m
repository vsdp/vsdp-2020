function [x z] = export_vsdp(IF,K,x,z)
%% EXPORT_VSDP - transform internal format back to imported problem format
%    [x z] = read_vsdp(IF,K,x,z)
%
%% >> Input:
% IF: 'SEDUMI', 'VSDP01' or [] for SeDuMi, old VSDP and new internal VSDP
%     format, respectively
% K: a structure with following fields
%     - K.f stores the number of free variables
%     - K.l is the number of nonnegative components
%     - K.q lists the lengths of socp blocks
%     - K.s lists the dimensions of semidefinite blocks
% x,z: primal/dual solution vector in SeDuMi or VSDP internal format
%
%% >> Output:
% x,y: in chosen format
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
% 17/08/12    M. Lange, written for common data export
%
%%
% TODO: rotated quadratic cones
%


%% check input
if nargin<3 || ~isstruct(K)
    error('VSDP:EXPORT_VSDP','not enough input parameter or wrong data format');
elseif nargin<4
    z = [];
end

%% conversion
if isempty(IF)
    %% export internal VSDP format
    % export x
    if ~isempty(x)
        x = vsvec(x,K,2,0);  % mu=2, does nothing if already compact vec format
    end
    % export z
    if ~isempty(z)
        z = vsvec(z,K,1,0);
    end     
elseif strcmpi(IF,'VSDP01')
    %% export old vsdp format
    % check supported cones
    if isfield(K,'f') && any(K.f>0) || isfield(K,'l') && any(K.l>0) || ...
            isfield(K,'q') && any(K.q>0)
        error('VSDP:EXPORT_VSDP','old vsdp format supports only semidefinite cone');
    end
    % convert x
    if ~isempty(x)
        x = vsvec(x,K,2,0);  % to verify that internal VSDP format is given
        blke = length(x);
        for j = length(K.s):-1:1
            nj = K.s(j);
            blks = blke - nj*(nj+1)/2 + 1;
            Xt{j}(triu(true(nj))) = 0.5 * x(blks:blke);  % mu==2
            Xt{j} = reshape(Xt{j},nj,nj);
            Xt{j} = Xt{j} + Xt{j}';
            blke = blks - 1;
        end
        x = Xt;
    end
    % convert z
    if ~isempty(z)
        z = vsvec(z,K,1,0);  % to verify that internal VSDP format is given
        blke = length(z);
        for j = length(K.s):-1:1
            nj = K.s(j);
            blks = blke - nj*(nj+1)/2 + 1;
            Zt{j}(triu(true(nj))) = z(blks:blke);
            Xt{j} = reshape(Xt{j},nj,nj);
            Zt{j} = Xt{j} + Xt{j}' - diag(sparse(diag(Xt{j})));
            blke = blks - 1;
        end
        z = Zt;
    end 
elseif strcmpi(IF,'SEDUMI')
    %% export sedumi format
    Imat = [];  % index vector to improve speed of smat
    % convert x
    if ~isempty(x)
        [x Imat] = vsmat(x,K,0.5,1);
    end
    % convert z
    if ~isempty(z)
        z = vsmat(z,K,1,1,Imat);
    end
else
    %% format not supported or detected
    error('VSDP:EXPORT_VSDP','tries to export unsupported format');
end

%_____________________________End EXPORT_VSDP____________________________