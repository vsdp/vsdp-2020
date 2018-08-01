function testVSDP_sindex (testCase)

vidx{1} = true;
midx{1} = logical ([1 0]);
lidx{1} = [];
lidx{1} = lidx{1}(:);

vidx{2} = logical ([1 0 1]');
midx{2} = logical ([ ...
  1 0 0 1;
  0 0 1 0]');
lidx{2} = 2;

vidx{3} = logical ([1 0 1 0 0 1]');
midx{3} = logical ([ ...
  1 0 0 0 1 0 0 0 1;
  0 0 0 1 0 0 1 1 0]');
lidx{3} = [2 3 6]';

vidx{4} = logical ([1 0 1 0 0 1 0 0 0 1]');
midx{4} = logical ([ ...
  1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1;
  0 0 0 0 1 0 0 0 1 1 0 0 1 1 1 0]');
lidx{4} = [2 3 7 4 8 12]';

for i = 1:4
  K.s = i;
  [v, m, l] = vsdp.sindex (K);
  verifyEqual (testCase, v, vidx{i})
  verifyEqual (testCase, m, midx{i})
  verifyEqual (testCase, l, lidx{i})
end

% Combined dimension test.
K.s = [2 2 4 1 3 3 4 2];
vv = vidx(K.s);
vv = cell2mat(vv(:));
mm = midx(K.s);
mm = cell2mat(mm(:));
ll = lidx(K.s);
offset = cumsum ([0, K.s(1:end-1).^2]);
for i = 1:length(offset)
  ll{i} = ll{i} + offset(i);
end
ll = cell2mat(ll(:));

[v, m, l] = vsdp.sindex (K);
verifyEqual (testCase, v, vv)
verifyEqual (testCase, m, mm)
verifyEqual (testCase, l, ll)

% Mixed cone test.
K.f = 2;
K.l = 3;
K.q = [2 3 4];
offset = K.f + K.l + sum (K.q);

[v, m, l] = vsdp.sindex (K);
verifyEqual (testCase, v, [false(offset, 1); vv])
verifyEqual (testCase, m, [false(offset, 2); mm])
verifyEqual (testCase, l, ll + offset)

end
