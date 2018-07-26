function val = cache (obj, val)
% CACHE  Getter and setter for cached entries.
%
%   Example:
%
%      A = rand(20);
%      obj.cache(A);           % Stores the variable A in cache.
%      B = obj.cache('A');     % Retrieves previous variable A from cache.
%      obj.cache('clear');     % Clears the cache if the current scope.
%      obj.cache('clearall');  % Clears the entire cache.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

narginchk (2, 2);

% If cache is disabled.
if (~obj.options.CACHE)
  val = [];
  return;
elseif (strcmpi (val, 'clearall'))
  obj.cache_memory = [];
  val = [];
  return;
end

% Get the calling function as scope.
ST = dbstack();
if (length (ST) == 2)
  caller = ST(2).name;
else
  caller = 'global';
end

% Extract or clear the relevant sub structure.
% Small note on 'clear': If field 'caller' does not exist, there is nothing
% to clear and the program flow leads to proper termination.
if (isfield (obj.cache_memory, caller))
  if (strcmpi (val, 'clear'))
    obj.cache_memory = rmfield (obj.cache_memory, caller);
    val = [];
    return;
  else
    s = getfield (obj.cache_memory, caller);
  end
else
  s = [];
end

if (ischar (val)) % Request variable from sub structure.
  if (isfield (s, val))
    val = getfield (s, val);
  else
    val = [];
  end
else % Store sub structure in cache.
  s = setfield (s, inputname(2), val);
  obj.cache_memory = setfield (obj.cache_memory, caller, s);
end

end
