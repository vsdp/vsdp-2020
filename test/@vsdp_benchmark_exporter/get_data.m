function obj = get_data (obj)
% GET_DATA  Get the data for the listed test cases.

% Copyright 2018-2019 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

len = size (obj.cdata, 1);
for i = 1:len
  fprintf ('  (%3d/%3d) %-10s %-10s %-10s\n', i, len, obj.cdata{i,1:3});
  vsdp_obj = get_vsdp_obj (obj, obj.cdata{i,1:3});
  get_vsdp_solutions (obj, obj.cdata(i,1:3), vsdp_obj);
  obj.cdata(i,4:5) = {vsdp_obj.m, vsdp_obj.n};
  obj.cdata(i,6:9) = {vsdp_obj.K.f > 0, vsdp_obj.K.l > 0, ...
    ~isempty(vsdp_obj.K.q), ~isempty(vsdp_obj.K.s)};
  if (~isempty (vsdp_obj.solutions.approximate))
    obj.cdata(i,10:12) = { ...
      vsdp_obj.solutions.approximate.f_objective(1), ...
      vsdp_obj.solutions.approximate.f_objective(2), ...
      vsdp_obj.solutions.approximate.solver_info.elapsed_time};
  end
  if (~isempty (vsdp_obj.solutions.rigorous_lower_bound))
    obj.cdata(i,13:14) = { ...
      vsdp_obj.solutions.rigorous_lower_bound.f_objective(1), ...
      vsdp_obj.solutions.rigorous_lower_bound.solver_info.elapsed_time};
  end
  if (~isempty (vsdp_obj.solutions.rigorous_upper_bound))
    obj.cdata(i,15:16) = { ...
      vsdp_obj.solutions.rigorous_upper_bound.f_objective(2), ...
      vsdp_obj.solutions.rigorous_upper_bound.solver_info.elapsed_time};
  end
end
% Add header line.
obj.cdata = [{'lib', 'name' 'sname', 'm', 'n', ...
  'K_f', 'K_l', 'K_q', 'K_s', 'fp', 'fd', 'ts', 'fL', 'tL', ...
  'fU', 'tU'}; obj.cdata];
end

function obj = get_vsdp_obj (obj, lib, name, solver)
persistent last_obj;
persistent last_lib;
persistent last_name;
persistent last_solver;
if (strcmp (lib, last_lib) && strcmp (name, last_name) ...
    && strcmp (solver, last_solver))
  obj = last_obj;  % Just return the already constructed object.
  return;
elseif (strcmp (lib, last_lib) && strcmp (name, last_name))
  obj = vsdp (last_obj);  % Make a clean copy.
else
  switch (lib)
    case 'DIMACS'
      src_file = dir (fullfile (lib, 'data', '**',[name, '.mat.gz']));
      if (length (src_file) > 1)
        error ('VSDP:vsdp_benchmark_exporter', ...
          'vsdp_benchmark_exporter: ''%s'' and ''%s'' is not unique.', ...
          lib, name);
      end
      tmp_file = extract_gz_file (obj, ...
        fullfile (src_file.folder, src_file.name));
      load (tmp_file, 'A*', 'b', 'c', 'K');
      if (exist ('A', 'var') == 1)
        obj = vsdp (A, b, c, K);
        clear ('A', 'b', 'c', 'K');
      else
        obj = vsdp (At, b, c, K);
        clear ('At', 'b', 'c', 'K');
      end
    case {'ESC', 'RDM', 'SPARSE_SDP'}
      if (strcmp (lib, 'ESC'))
        src_file = dir (fullfile (lib, 'data', [name, '_*.dat-s.gz']));
      elseif (strcmp (lib, 'RDM'))
        src_file = dir (fullfile (lib, 'data', [name, '.*.dat-s.gz']));
      else
        src_file = dir (fullfile (lib, 'data', [name, '*.dat-s.gz']));
      end
      if (length (src_file) > 1)
        error ('VSDP:vsdp_benchmark_exporter', ...
          'vsdp_benchmark_exporter: ''%s'' and ''%s'' is not unique.', ...
          lib, name);
      end
      tmp_file = extract_gz_file (obj, ...
        fullfile (src_file.folder, src_file.name));
      obj = vsdp.from_sdpa_file (tmp_file);
    case 'SDPLIB'
      obj = vsdp.from_sdpa_file (fullfile (lib, 'data', [name, '.dat-s']));
  end
  
  % Finally, delete temporary files.
  if (exist ('tmp_file', 'var') == 1)
    delete (sprintf('%s*', tmp_file));
  end
  
  % Optimize problem structure automatically, no output.
  obj = obj.analyze (true, false);
end
last_obj = obj;
last_lib = lib;
last_name = name;
last_solver = solver;
end

function get_vsdp_solutions (obj, lib_name_solver, vsdp_obj)
% Try to add approximate solution.
try
  load (fullfile (obj.data_dir, [strjoin([lib_name_solver, ...
    {'approximate_solution'}], '__'), '.mat']), 'app_sol');
  vsdp_obj.add_solution (app_sol.sol_type, app_sol.x, app_sol.y, ...
    app_sol.z, app_sol.f_objective, app_sol.solver_info);
catch
  % Ignore.
end
% Try to add rigorous lower bound.
try
  load (fullfile (obj.data_dir, [strjoin([lib_name_solver, ...
    {'rigorous_lower_bound'}], '__'), '.mat']), 'rig_lbd');
  vsdp_obj.add_solution (rig_lbd.sol_type, rig_lbd.x, rig_lbd.y, ...
    rig_lbd.z, rig_lbd.f_objective, rig_lbd.solver_info);
catch
  % Ignore.
end
% Try to add rigorous upper bound.
try
  load (fullfile (obj.data_dir, [strjoin([lib_name_solver, ...
    {'rigorous_upper_bound'}], '__'), '.mat']), 'rig_ubd');
  vsdp_obj.add_solution (rig_ubd.sol_type, rig_ubd.x, rig_ubd.y, ...
    rig_ubd.z, rig_ubd.f_objective, rig_ubd.solver_info);
catch
  % Ignore.
end
end

function tmp_file = extract_gz_file (obj, fpath)
[~, fname, fext] = fileparts (fpath);
% Copy file to temporary directory.
tmp_file = fullfile (obj.tmp_dir, [fname, fext]);
copyfile (fpath, tmp_file);
% Extract.
gunzip (tmp_file);
% Update data for working copy.
tmp_file((end - length ('.gz') + 1):end) = [];
end
