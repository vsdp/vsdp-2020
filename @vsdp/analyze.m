function obj = analyze (obj, yes_to_all)
% ANALYZE  Analyze the conic program for typical problems.
%
%   obj = obj.ANALYZE ()       Interactive mode.  If a typical problem is
%                              recognized for the VSDP conic program 'obj',
%                              ask the user if it should be fixed.
%
%   obj = obj.ANALYZE (true)   Automatically fix all typical problems.
%   obj = obj.ANALYZE (false)  Only report typical problems.
%
%   Typical problems:
%
%     #1 Check for diagonal SDP blocks:  If the constraint vector 'c.s(j)'
%        and the constraint matrix 'A.s(j)' only have entries on the main
%        diagonal, it is recommended to convert them into a LP block.
%        This saves computation time and memory: 'n*(n+1)/2' vs. 'n' entries.
%        Additionally, VSDP can easier compute rigorous cone bounds, as no
%        eigenvalue computations are involved.
%
%   See also vsdp.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (1, 2);
if (nargout ~= 1)
  % Determine object name.
  obj_name = inputname(1);
  if (isempty (obj_name))
    obj_name = 'obj';
  end
  error ('VDSP:analyze:updateHandle', ...
    ['analyze: Call this function by ''%s = %s.analyse()'' to avoid ', ...
    'stale object handles.'], obj_name, obj_name);
end
if (nargin == 1)
  interactive = true;
  yes_to_all = false;
else
  interactive = false;
end
obj = pattern1 (obj, yes_to_all, interactive);  % Check for diagonal SDP blocks.
end


function obj = pattern1 (obj, yes_to_all, interactive)
% PATTERN1  Check for diagonal only SDP blocks ==> LP blocks.

At = vsdp_indexable (obj.At, obj);
c  = vsdp_indexable (obj.c , obj);
A_linear = [];
c_linear = [];
drop_idx = [];
for i = 1:length(obj.K.s)
  cs = c.s(i);   % Extract the i-th semidefinite cone.
  As = At.s(i);
  idx = cumsum(1:obj.K.s(i));  % Compute diagonal indices
  Al = As(idx,:);  % Extract the diagonal elements of SDP cone.
  cl = cs(idx,:);
  % If the number of non-zero elements of the matrix and the diagonal match,
  % it must be a diagonal matrix, thus a LP matrix.
  if ((nnz (Al) == nnz (As)) && (nnz (cl) == nnz (cs)))
    warning ('VDSP:analyze:possibleLpCone', ...
      'analyze: K.s(%d) seems to only have diagonal elements.  ', i);
    if (yes_to_all)
      fprintf(' --> Convert it to LP block.\n');
    elseif (interactive)
      answer = input ('Convert it to LP block? [y/n] ', 's');
      if (isempty (answer) || (answer ~= 'y'))
        continue;
      end
    else
      continue;
    end
    A_linear = [A_linear; Al];
    c_linear = [c_linear; cl];
    drop_idx = [drop_idx; i];
  end
end
% Perform the SDP to LP update.
if (~isempty (drop_idx))
  % First delete the SDP entries.
  for i = flipud(drop_idx)
    idx = obj.K.idx.s(i,:);
    obj.At(idx(1):idx(end),:) = [];
    obj.c (idx(1):idx(end),:) = [];
    obj.K.s(i) = [];
  end
  % Then add the LP entries.
  obj.At = [obj.At(1:(obj.K.f + obj.K.l),:); ...
    A_linear; obj.At((obj.K.f + obj.K.l + 1):end,:)];
  obj.c = [obj.c(1:(obj.K.f + obj.K.l),:); ...
    c_linear; obj.c((obj.K.f + obj.K.l + 1):end,:)];
  obj.K.l = obj.K.l + size (A_linear, 1);
  % Finally, create a new fresh VSDP object.
  obj = vsdp(obj);
end
end
