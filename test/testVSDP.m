function tests = testVSDP()
% TESTVSDP  Runs a testsuite for VSDP (version 2018).
%
%   Example:
%
%       clc; table (runtests ('testVSDP')) % Matlab
%       clc; testVSDP;                     % Octave
%
%   See also vsdp.
%

% Copyright 2016-2018 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

tests = functiontests (localfunctions);
end

function testSINDEX (testCase)
testVSDP_sindex (testCase);
end

function testSVEC_SMAT (testCase)
testVSDP_svec_smat (testCase);
end

function testLP (testCase)
testVSDP_LP (testCase);
end

function testSOCP (testCase)
testVSDP_SOCP (testCase);
end

function testSDP (testCase)
testVSDP_SDP (testCase);
end
