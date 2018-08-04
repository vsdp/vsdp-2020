function obj = analyze (obj, yes_to_all)
% ANALYZE  Analyze the conic program for pathological patterns.
%
%   1) Check for diagonal SDP blocks:  If the constraint vector 'c.s(j)' and the 
%      constraint matrix 'A.s(j)' only contain entries on the main diagonal,
%      it is recommended to convert them into a LP block.  This saves memory
%      'n*(n+1)/2' vs. 'n' entries and computation time.
%
%   See also vsdp.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

% Check for diagonal only SDP blocks ==> LP blocks.
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
  % If the number of non-zero elements of the matrix and the diagonal match, it
  % must be a diagonal matrix, thus a LP matrix.
  if ((nnz (Al) == nnz (As)) && (nnz (cl) == nnz (cs)))
    warning ('VDSP:analyze:possibleLpCone', ...
      'analyze: K.s(%d) seems to only have diagonal elements.  ', i);
    if (nargout == 1)  % Properly assign the handle.
      if ((nargin == 2) && (yes_to_all))
        fprintf(' --> Convert it to LP block.\n');
      else
        answer = input ('Convert it to LP block? [y/n] ', 's');
        if (isempty (answer) || (answer ~= 'y'))
          continue;
        end
      end
      A_linear = [A_linear; Al];
      c_linear = [c_linear; cl];
      drop_idx = [drop_idx; i];
    end
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
  % Finally, create a new fresh VSDP object and analyse it recursively.
  obj = vsdp(obj).analyze();
end
end
