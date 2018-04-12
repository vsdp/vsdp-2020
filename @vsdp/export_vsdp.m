function [x,z] = export_vsdp(IF,K,x,z)
% EXPORT_VSDP  Export internal format back to imported problem format.
%
%   [x,z] = export_vsdp(IF,K,x,z)
%
% Input:
% IF: 'SEDUMI', 'VSDP01' or [] for SeDuMi, old VSDP and new internal VSDP
%     format, respectively
% K: a structure with following fields
%     - K.f stores the number of free variables
%     - K.l is the number of nonnegative components
%     - K.q lists the lengths of socp blocks
%     - K.s lists the dimensions of semidefinite blocks
% x,z: primal/dual solution vector in SeDuMi or VSDP internal format
%
% Output:
% x,y: in chosen format
%

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

% check input
if nargin<3 || ~isstruct(K)
  error('VSDP:EXPORT_VSDP','not enough input parameter or wrong data format');
elseif nargin<4
  z = [];
end

% conversion
if isempty(IF)
  % export internal VSDP format
  % export x
  if ~isempty(x)
    x = vsvec(x,K,2,0);  % mu=2, does nothing if already compact vec format
  end
  % export z
  if ~isempty(z)
    z = vsvec(z,K,1,0);
  end
elseif strcmpi(IF,'VSDP01')

elseif strcmpi(IF,'SEDUMI')
  % export sedumi format
  Imat = [];  % index vector to improve speed of smat
  % convert x
  if ~isempty(x)
    [x,Imat] = vsmat(x,K,0.5,1);
  end
  % convert z
  if ~isempty(z)
    z = vsmat(z,K,1,1,Imat);
  end
else
  error('VSDP:EXPORT_VSDP','tries to export unsupported format');
end

end
