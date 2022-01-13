function install_vsdp (settings_file)
% INSTALL_VSDP  Ensures VSDP to work properly.
%
%   For loading a custom settings file use:
%
%      install_vsdp ('my_vsdp_settings.yaml')
%
%   See also vsdp.
%

% Copyright 2004-2022 Christian Jansson (jansson@tuhh.de)

% Check whether path already exists and add path if necessary.
vsdp_path = fullfile (fileparts (which ('install_vsdp')));
addpath (vsdp_path);
addpath (fullfile (vsdp_path, 'test'));
if (exist ('OCTAVE_VERSION', 'builtin'))
  warning ('off', 'Octave:shadowed-function', 'local');
  addpath (fullfile (vsdp_path, 'octave'));
end
path (path);  % Refresh path

% Check for available solvers.
if (nargin)
  vsdp.settings ('settings_file', 'path', settings_file);
end

solver.registry.status ();

end
