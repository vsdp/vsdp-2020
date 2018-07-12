function tests = testVSDP()
% TESTVSDP  Runs a testsuite for VSDP (version 2018).
%
%   Example:
%
%       clc; table (runtests ('testVSDP')) % Matlab
%       clc; testVSDP;                     % Octave
%
%   See also demovsdp.

% Copyright 2016-2018 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

if (exist ('OCTAVE_VERSION', 'builtin'))
  addpath (fullfile (pwd (), '..', 'octave'));
end

tests = functiontests(localfunctions);
end

function testSINDEX(testCase)
testVSDP_sindex (testCase);
end

function testSVEC(testCase)
%testVSDP_svec (testCase);
end
