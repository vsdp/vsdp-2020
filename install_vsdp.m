function install_vsdp ()
% INSTALL_VSDP  Ensures VSDP to work properly.
%
%   See also vsdp.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% Check whether path already exists and add path if necessary.
vsdp_path = fullfile (fileparts (which ('install_vsdp')));
addpath (vsdp_path);
addpath (fullfile (vsdp_path, 'test'));
if (exist ('OCTAVE_VERSION', 'builtin'))
  addpath (fullfile (vsdp_path, 'octave'));
end
path (path);  % Refresh path

% Check for Interval-toolbox.
if (exist ('startintlab.m', 'file') ~= 2)
  error ('VSDP:vsdp_options:INTLAB', ...
    '%s.  %s\n\n\t%s\n\n%s.\n\n', ...
    'Interval toolbox "INTLAB" not found', ...
    'Get a recent version from', ...
    'http://www.ti3.tuhh.de/rump/intlab', ...
    'and run "startintlab" from the root directory');
end
end