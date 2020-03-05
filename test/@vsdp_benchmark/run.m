function obj = run (obj, options)
% RUN  Run the VSDP benchmark.
%
%   obj.run ()  Runs all benchmarks, specified in obj.BENCHMARK with all
%               solvers from obj.SOLVER.  Computed are an approximate solution
%               and rigorous lower and upper bounds.
%
%   obj.run (options)  Same as before, but specify a structure 'options' with
%                      fields:
%
%       - 'log' :    true (Default)  Write to log file     to 'log'  directory.
%       - 'save':    true (Default)  Save computed results to 'data' directory.
%       - 'compute': true (Default)  Compute any missing values.
%       - 'skip_log_exists': true (Default)  Skip tests, if log file for
%                                            approx. solution is already given.
%
%   NOTE [*]: This benchmark relaxes the problems of the ESC and RDM library.
%
%     The original ESC/RDM problems have equality constraints of the form
%
%       trace(gamma) - N = 0.
%
%     In the ESC/RDM data those are reformulated to two inequality conditions
%     with a small tolerance of epsilon (ESC) or epsilon = 0 (RDM):
%
%       trace(gamma) <=  N + epsilon
%      -trace(gamma) <= -N + epsilon
%
%     Those are in the last diagonal semidefinite block of 'vsdp_obj.At'.
%     By applying 'vsdp_obj.analyze(true);' this last diagonal semidefinite
%     block is turned into a linear block.  Furthermore, we remove any existing
%     tolerances 'epsilon' of that linear block by rounding and add
%     'epsilon = 1e-7' to 'vsdp_obj.c(1:vsdp_obj.K.l)'.
%
%   See also vsdp_benchmark, vsdp_benchmark.filter.
%

% Copyright 2004-2020 Christian Jansson (jansson@tuhh.de)

if (isempty (obj.RESULT_DIR))
  error ('VSDP_BENCHMARK:run:resultDirUnknown', ...
    'add_benchmark: Please set ''obj.RESULT_DIR''.');
end

if (nargin > 1)
  options = check_options (options);
else
  options = check_options ([]);
end

% Define helper functions for diary display.
print_header = @(str) fprintf ('\n%s\n%s\n%s\n\n', repmat ('-', 1, 50), ...
  str, repmat ('-', 1, 50));

% Solve selected test cases.
for j = 1:size (obj.BENCHMARK, 1)
  fprintf ('(%3d/%3d) %s/%s\n', j, size (obj.BENCHMARK, 1), ...
    obj.BENCHMARK{j,1}, obj.BENCHMARK{j,2});
  
  % Call all selected solvers.
  for i = 1:length (obj.SOLVER)
    
    % Specify file names.
    prefix = sprintf('%s__%s__%s', obj.BENCHMARK{j,1}, obj.BENCHMARK{j,2}, ...
      obj.SOLVER{i});
    data_file_prefix = fullfile (obj.RESULT_DIR, 'data', prefix);
    app_sol_file = sprintf ('%s__approximate_solution.mat', data_file_prefix);
    rig_lbd_file = sprintf ('%s__rigorous_lower_bound.mat', data_file_prefix);
    rig_ubd_file = sprintf ('%s__rigorous_upper_bound.mat', data_file_prefix);
    log_file = fullfile (obj.RESULT_DIR, 'log', sprintf('%s.log', prefix));
    
    % Check if other worker is already on this test case.
    if (options.skip_log_exists && (exist (log_file, 'file') == 2))
      fprintf ('\n  >>> Skipping, log file exists.\n\n');
      continue;
    end
    
    % Start logging.
    if (options.log)
      diary (log_file);
    end
    
    % Display comprehensive header.
    fprintf ('\n\n');
    print_header (sprintf ('Test : %s/%s', obj.BENCHMARK{j,1}, ...
      obj.BENCHMARK{j,2}));
    obj.disp ('system');
    
    try
      [fpath, fname, fext] = fileparts (obj.BENCHMARK{j,3});
      
      % Extract *.gz-archive if necessary.
      if (strcmp (fext, '.gz'))
        % Copy file to temporary directory.
        tmp_file = fullfile (obj.TMP_DIR, [fname, fext]);
        copyfile (obj.BENCHMARK{j,3}, tmp_file);
        % Extract.
        gunzip (tmp_file);
        % Update data for working copy.
        tmp_file((end - length ('.gz') + 1):end) = [];
        [fpath, fname, fext] = fileparts (tmp_file);
      end
      
      % Import data to VSDP object 'obj' depending on the file type.
      dfile = fullfile (fpath, [fname, fext]);
      switch (fext)
        case '.mat'   % MAT-file.
          load (dfile, 'A*', 'b', 'c', 'K');
          if (exist ('A', 'var') == 1)
            vsdp_obj = vsdp (A, b, c, K);
            clear ('A', 'b', 'c', 'K');
          else
            vsdp_obj = vsdp (At, b, c, K);
            clear ('At', 'b', 'c', 'K');
          end
        case '.dat-s'  % Sparse SDPA data.
          vsdp_obj = vsdp.from_sdpa_file (dfile);
        case '.SIF'    % MPS data.
          vsdp_obj = vsdp.from_mps_file (dfile);
        otherwise
          warning ('VSDP_BENCHMARK:run:unsupportedData', ...
            'run: Unsupported file ''%s''.', obj.BENCHMARK{j,3});
          continue;
      end
      
      % Finally, delete temporary files.
      if (exist ('tmp_file', 'var') == 1)
        delete (sprintf('%s*', tmp_file));
      end
    catch err
      fprintf (2, '\n\n%s\n\n', err.message);
      continue;
    end
    
    % Optimize problem structure.
    vsdp_obj = vsdp_obj.analyze (true);
    
    % Relax ESC and RDM problems (see note above [*]).
    if (strcmp (obj.BENCHMARK{j,1}, 'ESC') ...
        || strcmp (obj.BENCHMARK{j,1}, 'RDM'))
      [K.f, K.l, K.q, K.s] = deal (vsdp_obj.K.f, vsdp_obj.K.l, vsdp_obj.K.q, ...
        vsdp_obj.K.s);
      c = vsdp_obj.c;
      c(1:K.l) = round(c(1:K.l)) + 1e-7;
      data = {vsdp_obj.At, vsdp_obj.b, c, K};
      vsdp_obj = vsdp (data{:});
      print_header ('Relax ESC and RDM problem parameter ''c''.');
    end
    
    try
      print_header (sprintf ('Start: %s   Solver: ''%s''', ...
        datestr (now ()), obj.SOLVER{i}));
      
      % Make a clean copy and set the solver to be used.
      vsdp_obj = vsdp (vsdp_obj);
      vsdp_obj.options.SOLVER = obj.SOLVER{i};
      
      % Compute approximate solution, if not already computed.
      print_header ('Approximate solution');
      if (exist (app_sol_file, 'file') ~= 2)
        if (options.compute)
          vsdp_obj.solve (obj.SOLVER{i});
        else
          fprintf ('>>> skipped <<<');
        end
        if (options.save)
          app_sol = get_solution_as_struct (vsdp_obj.solutions.approximate);
          save (app_sol_file, 'app_sol', '-v7');
        end
      else  % ... or load from file.
        load (app_sol_file, 'app_sol')
        vsdp_obj.add_solution (app_sol.sol_type, app_sol.x, app_sol.y, ...
          app_sol.z, app_sol.f_objective, app_sol.solver_info);
        fprintf ('>>> cached <<<');
      end
      fprintf ('\n>>> done   <<<\n');
      
      % Compute rigorous lower bound, if not already computed.
      print_header ('Rigorous lower bound');
      if (exist (rig_lbd_file, 'file') ~= 2)
        if (options.compute)
          vsdp_obj.rigorous_lower_bound ();
        else
          fprintf ('>>> skipped <<<');
        end
        if (options.save)
          rig_lbd = get_solution_as_struct ( ...
            vsdp_obj.solutions.rigorous_lower_bound);
          save (rig_lbd_file, 'rig_lbd', '-v7');
        end
      else  % ... or load from file.
        load (rig_lbd_file, 'rig_lbd')
        vsdp_obj.add_solution (rig_lbd.sol_type, rig_lbd.x, rig_lbd.y, ...
          rig_lbd.z, rig_lbd.f_objective, rig_lbd.solver_info);
        fprintf ('>>> cached <<<');
      end
      fprintf ('\n>>> done   <<<\n');
      
      % Compute rigorous upper bound, if not already computed.
      print_header ('Rigorous upper bound');
      if (exist (rig_ubd_file, 'file') ~= 2)
        if (options.compute)
          vsdp_obj.rigorous_upper_bound ();
        else
          fprintf ('>>> skipped <<<');
        end
        if (options.save)
          rig_ubd = get_solution_as_struct ( ...
            vsdp_obj.solutions.rigorous_upper_bound);
          save (rig_ubd_file, 'rig_ubd', '-v7');
        end
      else  % ... or load from file.
        load (rig_ubd_file, 'rig_ubd')
        vsdp_obj.add_solution (rig_ubd.sol_type, rig_ubd.x, rig_ubd.y, ...
          rig_ubd.z, rig_ubd.f_objective, rig_ubd.solver_info);
        fprintf ('>>> cached <<<');
      end
      fprintf ('\n>>> done   <<<\n');
      
      % Compute rigorous upper bound with a priori bounds, if not already
      % computed.
      if ((strcmp (obj.BENCHMARK{j,1}, 'ESC') ...
          || strcmp (obj.BENCHMARK{j,1}, 'RDM')) ...
          && (exist ('esc_rdm_a_priori_bounds', 'file') == 2))
        esc_rdm_a_priori_bounds (vsdp_obj, obj.BENCHMARK{j,1}, ...
          obj.BENCHMARK{j,2}, obj.SOLVER{i}, options);
      end
    catch err
      fprintf (2, '\n\n%s\n\n', err.message);
      continue;
    end
    print_header (sprintf ('End: %s', datestr (now ())));
    if (options.log)
      diary ('off');
    end
  end
end
end


function sol = get_solution_as_struct (vsdp_sol)
% GET_SOLUTION_AS_STRUCT  Convert VSDP solution class to struct without warning.

if (exist ('OCTAVE_VERSION', 'builtin'))
  S = warning ('off', 'Octave:classdef-to-struct');
else
  S = warning ('off', 'MATLAB:structOnObject');
end
sol = struct (vsdp_sol);
warning (S);
end


function opts = check_options (options)
% CHECK_OPTIONS  Check given structure options or return default values.
opts.log = true;
opts.save = true;
opts.compute = true;
opts.skip_log_exists = true;

if (~isempty (options))
  if (isfield (options, 'log'))
    opts.log = options.log;
  end
  if (isfield (options, 'save'))
    opts.save = options.save;
  end
  if (isfield (options, 'compute'))
    opts.compute = options.compute;
  end
  if (isfield (options, 'skip_log_exists'))
    opts.skip_log_exists = options.skip_log_exists;
  end
end
end
