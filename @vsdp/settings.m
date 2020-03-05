function val = settings (id, key, val)
% SETTINGS  Function to store an load global VSDP settings.
%
%   Store solver setting:
%
%      vsdp.settings ('id', 'key', 'val');
%      vsdp.settings ('vsdp', 'path', '/path/to/vsdp');
%
%   Load solver setting:
%
%      val = vsdp.settings ('id', 'key');
%

% Copyright 2004-2020 Christian Jansson (jansson@tuhh.de)

persistent vsdp_settings;  % In memory copy of VSDP settings file.
persistent settings_file;

narginchk (2, 3);

% Assume first initialization.
if (isempty (settings_file))
  vsdp_path = fullfile (fileparts (which ('install_vsdp')));
  settings_file = fullfile (vsdp_path, ...
    sprintf ('vsdp_settings_%d.mat', sys_hash ()));
  if (exist (settings_file, 'file') == 2)
    load (settings_file, 'vsdp_settings');
  else
    vsdp_settings.vsdp.path = vsdp_path;
    save (settings_file, 'vsdp_settings', '-v7');
  end
end

switch (nargin)
  case 2  % Load setting
    val = [];
    if (isfield (vsdp_settings, id) ...
        && isfield (getfield (vsdp_settings, id), key))
      val = getfield (getfield (vsdp_settings, id), key);
    end
  case 3  % Store setting
    vsdp_settings = setfield (vsdp_settings, id, struct (key, val));
    save (settings_file, 'vsdp_settings', '-v7');
  otherwise
    error ('VSDP:SOLVER:UTIL:settings', 'solver.util.settings: Bad usage.');
end
end


function val = sys_hash ()
% SYS_HASH  Compute a system unique hash value.
%

val = sum (double ([computer(), '_vsdp_2020_', version()]));
end
