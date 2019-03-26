function install_vsdp ()
% INSTALL_VSDP  Ensures VSDP to work properly.
%
%   See also vsdp.
%

% Copyright 2004-2019 Christian Jansson (jansson@tuhh.de)

% Check whether path already exists and add path if necessary.
vsdp_path = fullfile (fileparts (which ('install_vsdp')));
addpath (vsdp_path);
addpath (fullfile (vsdp_path, 'test'));
if (exist ('OCTAVE_VERSION', 'builtin'))
  addpath (fullfile (vsdp_path, 'octave'));
end
path (path);  % Refresh path

% Check for available solvers.
solver.registry.status ();

end
