function val = settings (id, key, val)
% SETTINGS  Function to store an load global VSDP settings.
%
%   Store solver setting:
%
%      vsdp.settings ('id', 'key', 'val');
%      vsdp.settings ('vsdp', 'path', '/path/to/vsdp');
%      vsdp.settings ('settings_file', 'path', '/path/to/vsdp_settings_file');
%
%   Load solver setting:
%
%      val = vsdp.settings ('id', 'key');
%

% Copyright 2004-2022 Christian Jansson (jansson@tuhh.de)

persistent vsdp_settings;  % In memory copy of VSDP settings file.
persistent settings_file;

narginchk (2, 3);

vsdp_path = fullfile (fileparts (which ('install_vsdp')));

% Set the settings_file string itself.
if (strcmp (id, 'settings_file'))
  if (nargin == 2)
    val = settings_file;
    return;
  else
    settings_file = val;
    vsdp_settings = [];
  end
elseif (isempty (settings_file))  % Ensure default value.
  settings_file = fullfile (vsdp_path, 'vsdp_settings.yaml');
end

% vsdp_settings initialization.
if (isempty (vsdp_settings))
  if (exist (settings_file, 'file') == 2)
    vsdp_settings = load_vsdp_settings (settings_file);
  else
    vsdp_settings.vsdp.path = vsdp_path;
    save_vsdp_settings (settings_file, vsdp_settings);
  end
end

% Finished when changing the settings_file.
if (strcmp (id, 'settings_file'))
  return;
end

% Load setting
if (nargin == 2)
  val = [];
  if (isfield (vsdp_settings, id) ...
      && isfield (getfield (vsdp_settings, id), key))
    val = getfield (getfield (vsdp_settings, id), key);
  end
else  % Store setting
  vsdp_settings = setfield (vsdp_settings, id, struct (key, val));
  save_vsdp_settings (settings_file, vsdp_settings);
end
end


function save_vsdp_settings (settings_file, vsdp_settings)
  items = fieldnames (vsdp_settings);
  str = sprintf ('%s\n', '---');
  for i = 1:length (items)
    item = vsdp_settings.(items{i});
    keys = fieldnames (item);
    str = sprintf ('%s%s\n', str, ['- "', items{i}, '":']);
    for j = 1:length (keys)
      key = keys{j};
      val = item.(key);
      str = sprintf ('%s%s\n', str, ['  - ', key, ': "', val, '"']);
    end
  end
  fid = fopen (settings_file, 'w');
  fprintf (fid, '%s', str);
  fclose (fid);
end


function vsdp_settings = load_vsdp_settings (settings_file)
  vsdp_settings = struct ();
  str = fileread (settings_file);
  str = strsplit (str, {'\n', '\r'});
  idx = cellfun (@isempty, str) | cellfun (@(s) strcmp (s, '---'), str);
  str = str(~idx);
  items = cellfun (@(s) (s(1) == '-'), str);
  items = cumsum (items);  % 1 1 2 2 ...
  for i = 1:max(items)
    % [1,1] = - "vsdp":
    % [1,2] =   - path: "/path/to/vsdp/"
    % [1,3] =   - key: "val"
    item = str(items == i);
    id = get_quot_str_content (item{1});
    item = item(2:end);
    attrs = cell (1, length (item));
    for j = 1:length (item)
      attr = strsplit (item{j}, {':'});
      attrs{2*j - 1} = get_quot_str_content (attr{1});
      if (length (attr) > 1)
        attrs{2*j} = get_quot_str_content (attr{2});
      end
    end
    vsdp_settings = setfield (vsdp_settings, id, struct (attrs{:}));
  end
end


function str = get_quot_str_content (str)
  rstr = regexp (str, '.*"(.*)"', 'tokens');
  if (~ isempty (rstr))
    str = rstr{1};
    str = str{1};
  else
    % Heuristic to extract relevant string from YAML.
    str = strtrim (str);
    if (str(1) == '-')
      str = str(2:end);
    end
    if (str(end) == ':')
      str = str(1:end-1);
    end
    str = strtrim (str);
  end
end
