function testVSDP_sindex (testCase)

vidx{1} = logical ([1 0]);
midx{1} = logical ([1 0]);
lidx{1} = [];
lidx{1} = lidx{1}(:);

vidx{2} = logical ([ ...
  1 0 1;
  0 1 0]');
midx{2} = logical ([ ...
  1 0 0 1;
  0 0 1 0]');
lidx{2} = 2;

vidx{3} = logical ([ ...
  1 0 1 0 0 1;
  0 1 0 1 1 0]');
midx{3} = logical ([ ...
  1 0 0 0 1 0 0 0 1;
  0 0 0 1 0 0 1 1 0]');
lidx{3} = [2 3 6]';

vidx{4} = logical ([ ...
  1 0 1 0 0 1 0 0 0 1;
  0 1 0 1 1 0 1 1 1 0]');
midx{4} = logical ([ ...
  1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1;
  0 0 0 0 1 0 0 0 1 1 0 0 1 1 1 0]');
lidx{4} = [2 3 7 4 8 12]';

for i = 1:4
  [v, m, l] = vsdp.sindex (i);
  verifyEqual (testCase, v, vidx{i})
  verifyEqual (testCase, m, midx{i})
  verifyEqual (testCase, l, lidx{i})
end

% Combined dimension test.
d = [2 2 4 1 3 3 4 2];
vv = vidx(d);
vv = cell2mat(vv(:));
mm = midx(d);
mm = cell2mat(mm(:));
ll = lidx(d);
offset = cumsum ([0, d(1:end-1).^2]);
for i = 1:length(offset)
  ll{i} = ll{i} + offset(i);
end
ll = cell2mat(ll(:));

[v, m, l] = vsdp.sindex (d);
verifyEqual (testCase, v, vv)
verifyEqual (testCase, m, mm)
verifyEqual (testCase, l, ll)

% Test for bad input.
verifyError(testCase, @() vsdp.sindex (-1),       'VSDP:sindex:badD');
verifyError(testCase, @() vsdp.sindex ([1 2 -1]), 'VSDP:sindex:badD');
verifyError(testCase, @() vsdp.sindex ([1 0 12]), 'VSDP:sindex:badD');

end
