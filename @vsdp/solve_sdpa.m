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

% Check cones
if ((obj.K.f > 0) || (sum(obj.K.q) > 0))
  error('VSDP:solve_sdpa:unsupportedCones', ...
    ['solve_sdpa: free variables (K.f) second-order cone variables (K.q) ', ...
    'are not supported by SDPA.']);
end

% Note: for 'x0' and later 'x' see [1, p. 14] "mDIM -- All the letters after m
% through the end of the line are neglected".

% Should initial solution guess be taken into account?
if (obj.options.USE_STARTING_POINT)
  [x0, X0, Y0] = deal (obj.y, obj.z, obj.x);
  x0 = [x0; 0];  % expand to mDIM
  X0 = mat2cell (X0,  obj.K.dims, 1);
  Y0 = mat2cell (Y0,  obj.K.dims, 1);
  X0 = cellfun(@(x) vsdp.smat([], x, 1),   X0, 'UniformOutput', false);
  Y0 = cellfun(@(x) vsdp.smat([], x, 1/2), Y0, 'UniformOutput', false);
else
  [x0, X0, Y0] = deal ([], [], []);
end
% Clear previous solutions or initial points.
[obj.x, obj.y, obj.z] = deal ([], [], []);

if (~isempty (obj.options.SOLVER_OPTIONS))
  OPTIONS = obj.options.SOLVER_OPTIONS;
else
  OPTIONS = [];
end
if (~obj.options.VERBOSE_OUTPUT)
  OPTIONS.print = 'no';
end

% In case of interval data solve midpoint problem.
c = mid (-obj.b);
F = [ ...
  mat2cell(mid (-obj.c),  obj.K.dims, 1), ...
  mat2cell(mid (-obj.At), obj.K.dims, ones(1, obj.m))];
F = cellfun(@(x) vsdp.smat([], x, 1), F, 'UniformOutput', false);
[nBLOCK, mDIM] = size (F);
bLOCKsTRUCT = [-obj.K.l(obj.K.l > 0), obj.K.s'];

% Call solver.
if (exist ('mexsdpa', 'file') == 3)
  [~, x, X, Y, ~] = sdpam ...
    (mDIM, nBLOCK, bLOCKsTRUCT, c, F, x0, X0, Y0, OPTIONS);
elseif (exist('callSDPA','file') == 2)
  [x, X, Y] = callSDPA ...
    (mDIM, nBLOCK, bLOCKsTRUCT, c, F, x0, X0, Y0, OPTIONS);
else
  error('VSDP:MYSDPS', 'You need to compile the SDPA MEX-interface.');
end

% Store results.
obj.x = vsdp.svec (obj, vsdp.cell2mat (Y), 2);
obj.y = x(1:end-1);
obj.z = vsdp.svec (obj, vsdp.cell2mat (X), 1);

info = 0;

end
