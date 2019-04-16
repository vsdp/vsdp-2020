function save_state (obj)
% SAVE_STATE  Restores the state of a VSDP benchmark.
%
%   Save the VSDP benchmark state to the file 'benchmark_state.mat' in
%   the object's RESULT_DIR.
%
%   See also vsdp_benchmark.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

bm_state_file = fullfile (obj.RESULT_DIR, 'benchmark_state.mat');
bm_data  = obj.BENCHMARK;
sol_data = obj.SOLVER;
save (bm_state_file, 'bm_data', 'sol_data', '-v7');

end
