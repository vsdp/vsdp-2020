function [slist, plist] = registry (cmd)
% REGISTRY  Helpter to handle by VSDP supported solvers.
%
%   For a solver to be supported by VSDP, it must implement a proxy class in the
%   +solver directory.
%
%   Possible calls:
%
%       slist = solver.registry ('list')  List all supported solver and return
%                                         them as cell array of strings.
%                                                 
%      [slist, plist] = solver.registry ('install')  Like above.  Additionally
%                                                    check their availability
%                       and return the installed solver location as cell array
%                       of strings.
%
%       str = solver.registry ('status')  Displays an human readable overview
%                                         about supported and available solver.
%

narginchk (0, 1);
if (nargin == 0)
  cmd = 'status';
end

switch (cmd)
  case 'list'
    nargoutchk (0, 1);
    % Change to VSDP root directory.
    old_dir = cd (vsdp.settings ('vsdp', 'path'));
    
    % Generate a list of supported solvers.
    slist = dir ('+solver');
    cd (old_dir);
    
    slist = {slist.name};
    idx = (cellfun (@length, slist) > 2);      % Strip items like '.' and '..'.
    idx = idx & ~strcmp ('registry.m', slist); % Strip this function.
    slist = slist(idx);
    % Strip *.m suffix.
    [~, slist] = cellfun (@fileparts, slist, 'UniformOutput', false);
    slist = slist(:);
  case 'install'
    slist = solver.registry ('list');
    plist = cell (size (slist));
    for i = 1:length (slist)
      plist{i} = eval (sprintf ('solver.%s.install ();', slist{i}));
    end
  case 'status'
    nargoutchk (0, 1);
    [slist, plist] = solver.registry ('install');
    idx = cellfun (@isempty, plist);
    plist(idx) = {'-- not available --'};
    output = [slist, plist]';
    slist = sprintf ('    %-10s (%s)\n', output{:});
    slist = sprintf ('\n\n  Solver detected by VSDP:\n\n%s', slist);
    disp (slist);
  otherwise
    error ('VSDP:SOLVER:registry:unknownCommand', ...
      'solver.registry: unknown command ''%s''', cmd);
end
end
