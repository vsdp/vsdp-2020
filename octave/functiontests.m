function tests = functiontests(localfcns)
% FUNCTIONTESTS  Run all functions prefixed by "test" and print a test summary.
%
%   This MATLAB function is not present in GNU Octave (<= 5.1.0).
%
%   tests = functiontests(localfcns) runs all function handles from the cell
%                                    array localfcns.
%
%   Example:
%
%       tests = functiontests(localfunctions())
%
%   See also testVSDP, testVSDP2006.
%

% Copyright 2016-2019 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

j = 1;
for i = 1:length(localfcns)
  test_case = localfcns{i};
  test_case_name = func2str (test_case);
  if (! strncmp (test_case_name, "test", 4))
    continue;
  endif
  tests(j).Name = test_case_name;
  passed = true;
  t = tic;
  try
    test_case([]);
  catch (err)
    disp (sprintf ("\n--> %s\n", test_case_name))
    disp (err.message)
    passed = false;
  end_try_catch
  tests(j).Passed = passed;
  tests(j).Failed = ~passed;
  tests(j).Duration = toc(t);
  j = j + 1;
  fprintf (stdout, ".");
endfor
disp(" ")
disp("Test summary")
disp("------------")
disp(" ")
bool = {"FAIL", "PASSED"};
name_len = 0;
for i = 1:length(tests)
  name_len = max (name_len, length (tests(i).Name));
endfor
for i = 1:length(tests)
  disp(sprintf("%s%s\t%s\t%f", tests(i).Name, ...
    blanks (name_len - length (tests(i).Name)), ...
    bool{tests(i).Passed + 1}, tests(i).Duration));
endfor
endfunction
