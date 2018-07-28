function obj = solve_sdpa (obj)
% SOLVE_SDPA  Approximately solve conic problem instance with SDPA.
%
%   For more information on the SDPA-M format, see:
%
%     [1] https://sourceforge.net/projects/sdpa/files/sdpa-m/sdpamManual.pdf
%         Version 2005.
%
%   See also vsdp.solve.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% Check solver availability.
if ((exist ('mexsdpa', 'file') ~= 3) && (exist ('callSDPA', 'file') ~= 2))
  error ('VSDP:solve_sdpa:notAvailable', ...
    ['solve_sdpa: SDPA does not seem to be ready.\n\n', ...
    'To select another solver, run:  %s.solve()'], inputname(1));
end

% Check cones
if ((obj.K.f > 0) || (sum(obj.K.q) > 0))
  error('VSDP:solve_sdpa:unsupportedCones', ...
    ['solve_sdpa: free variables (K.f) second-order cone variables (K.q) ', ...
    'are not supported by SDPA.']);
end

% Note: for 'x0' and later 'x' see [1, p. 14] "mDIM -- All the letters after m
% through the end of the line are neglected".

% Should initial solution guess be taken into account?
if ((obj.options.USE_STARTING_POINT) && (~isempty (obj.solution)))
  [x0, X0, Y0] = deal (obj.solution.y, obj.solution.z, obj.solution.x);
  x0 = [x0; 0];  % expand to mDIM
  X0 = mat2cell (X0,  obj.K.dims, 1);
  Y0 = mat2cell (Y0,  obj.K.dims, 1);
  X0 = cellfun(@(x) vsdp.smat([], x, 1),   X0, 'UniformOutput', false);
  Y0 = cellfun(@(x) vsdp.smat([], x, 1/2), Y0, 'UniformOutput', false);
else
  [x0, X0, Y0] = deal ([], [], []);
end

% Should special solver options be taken into account?
if (~isempty (obj.options.SOLVER_OPTIONS))
  OPTIONS = obj.options.SOLVER_OPTIONS;
else
  OPTIONS = [];
end

% Adapt output verbosity.
if (~obj.options.VERBOSE_OUTPUT)
  OPTIONS.print = 'no';
end

% Prepare data for solver.
[At, c, b] = obj.get_perturbed_midpoint_problem ();  % Note b <--> c!
c = -c;
F = [ ...
  mat2cell(-b,  obj.K.dims, 1), ...
  mat2cell(-At, obj.K.dims, ones(1, obj.m))];
F = cellfun(@(x) vsdp.smat([], x, 1), F, 'UniformOutput', false);
[nBLOCK, mDIM] = size (F);
bLOCKsTRUCT = [-obj.K.l(obj.K.l > 0), obj.K.s'];

% Call solver.
tic;
if (exist ('mexsdpa', 'file') == 3)
  [~, x, X, Y, INFO] = sdpam ...
    (mDIM, nBLOCK, bLOCKsTRUCT, c, F, x0, X0, Y0, OPTIONS);
elseif (exist('callSDPA','file') == 2)
  [x, X, Y] = callSDPA ...
    (mDIM, nBLOCK, bLOCKsTRUCT, c, F, x0, X0, Y0, OPTIONS);
  INFO = 0;
else
  error('VSDP:MYSDPS', 'You need to compile the SDPA MEX-interface.');
end
elapsed_time = toc;

% Store solution.
x = vsdp.svec (obj, vsdp.cell2mat (Y), 2);
y = x(1:end-1);
z = vsdp.svec (obj, vsdp.cell2mat (X), 1);
f_objective = [obj.c'*x; obj.b'*y];
if (isstruct (INFO))
  info = INFO.phase.value;
else
  info = INFO;
end
obj.add_solution(x, y, z, f_objective, obj.options.SOLVER, info, elapsed_time);

end
