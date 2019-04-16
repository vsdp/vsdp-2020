function export_default_tables (obj, outdir)
% EXPORT_DEFAULT_TABLES  Generates some standard tables for publication.
%
%   obj.export_default_tables():  Export a set of default tables to the
%                                 folder 'result_tables' in the current
%                                 working directory.
%
%   export_default_tables (outdir):  Export a set of default tables to the
%                                    folder 'outdir'.
%
%   See also vsdp_benchmark.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

if (nargin < 2)
  outdir = 'result_tables';
end
mkdir (outdir);

% For all benchmark libraries.
libs = unique ({obj.BENCHMARK.lib});
for i = 1:length (libs)
  disp (libs{i});
  % Export huge tables for the website publication.
  fmts = {'cell', 'csv', 'html'};
  fnames = fullfile (outdir, strcat (libs{i}, {'.mat', '.csv', '.html'}));
  use_columns = { ...
    'name',     'Name',    '%s';
    'm',        '$m$',     '%d';
    'n',        '$n$',     '%d';
    'K_f',      '$K_f$',   '%d';
    'K_l',      '$K_l$',   '%d';
    'K_q',      '$K_q$',   '%d';
    'K_s',      '$K_s$',   '%d';
    'sname',    'Solver',  '%s';
    'fp',       '$f_p$',   '%e';
    'fd',       '$f_d$',   '%e';
    'fL',       '$\underline{f_p}$',   '%e';
    'fU',       '$\overline{f_d}$',    '%e';
    'fU-fL',    '$\overline{f_d} - \underline{f_p}$',     '%e';
    'mu_fU_fL', '$\mu(\overline{f_d}, \underline{f_p})$', '%e';
    'ts',       '$t_s$',   '%.2f';
    'tL',       '$\underline{t}$',     '%.2f';
    'tL/ts',    '$\underline{t}/t_s$', '%.2e';
    'tU',       '$\overline{t}$',      '%.2f';
    'tU/ts',    '$\overline{t}/t_s$',  '%.2e';};
  obj.export (fmts, fnames, obj.filter (libs{i}), use_columns);
  
  % Export solver statistics.
  fmts = {'markdown'};
  use_columns = { ...
    'mu_fU_fL', '$\mu(\overline{f_d}, \underline{f_p})$', '%e';
    'ts',       '$t_s$',   '%.2f';
    'tL/ts',    '$\underline{t}/t_s$', '%.2e';
    'tU/ts',    '$\overline{t}/t_s$',  '%.2e';};
  stat_funs = {'min', 'mean', 'median', 'max'};
  solvers = unique ({obj.SOLVER(3:end).name});
  for j = 1:length (solvers)
    fnames = {fullfile(outdir, sprintf ('%s_%s.md', libs{i}, solvers{j}))};
    obj.export (fmts, fnames, obj.filter (libs{i}, [], solvers{j}), ...
      use_columns, stat_funs);
  end
  
  % Export tables for thesis.
  fmts = {'latex'};
  fnames = {fullfile(outdir, sprintf ('%s.tex', libs{i}))};
  use_columns = { ...
    'name',     'Name',    '%s';
    'm',        '$m$',     '%d';
    'n',        '$n$',     '%d';
    'sname',    'Solver',  '%s';
    'fp',       '$f_p$',   '%e';
    'fd',       '$f_d$',   '%e';
    'fL',       '$\underline{f_p}$',   '%e';
    'fU',       '$\overline{f_d}$',    '%e';
    'fU-fL',    '$\overline{f_d} - \underline{f_p}$',     '%e';
    'mu_fU_fL', '$\mu(\overline{f_d}, \underline{f_p})$', '%e';
    'ts',       '$t_s$',   '%.2f';
    'tL',       '$\underline{t}$',     '%.2f';
    'tL/ts',    '$\underline{t}/t_s$', '%.2e';
    'tU',       '$\overline{t}$',      '%.2f';
    'tU/ts',    '$\overline{t}/t_s$',  '%.2e';};
  obj.export (fmts, fnames, obj.filter (libs{i}), use_columns);
end
end
