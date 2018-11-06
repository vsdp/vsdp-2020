function testVSDP_sindex (testCase)

vidx{1} = logical ([1 0]);
midx{1} = logical ([1 0]);
mlidx{1} = [];
mlidx{1} = mlidx{1}(:);
vlidx{1} = 1;

vidx{2} = logical ([ ...
  1 0 1;
  0 1 0]');
midx{2} = logical ([ ...
  1 0 0 1;
  0 0 1 0]');
mlidx{2} = 2;
vlidx{2} = (1:3)';

vidx{3} = logical ([ ...
  1 0 1 0 0 1;
  0 1 0 1 1 0]');
midx{3} = logical ([ ...
  1 0 0 0 1 0 0 0 1;
  0 0 0 1 0 0 1 1 0]');
mlidx{3} = [2 3 6]';
vlidx{3} = [1 2 4 3 5 6]';

vidx{4} = logical ([ ...
  1 0 1 0 0 1 0 0 0 1;
  0 1 0 1 1 0 1 1 1 0]');
midx{4} = logical ([ ...
  1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1;
  0 0 0 0 1 0 0 0 1 1 0 0 1 1 1 0]');
mlidx{4} = [2 3 7 4 8 12]';
vlidx{4} = [1 2 5 3 6 8 4 7 9 10]';

for i = 1:4
  K.s = i;
  [v, m, l] = vsdp.sindex (K);
  verifyEqual (testCase, v, vidx{i})
  verifyEqual (testCase, m, midx{i})
  verifyEqual (testCase, l, mlidx{i})
  [v, m, l, vl] = vsdp.sindex (K);
  verifyEqual (testCase, v,  vidx{i})
  verifyEqual (testCase, m,  [])
  verifyEqual (testCase, l,  [])
  verifyEqual (testCase, vl, vlidx{i})
end

% Combined dimension test.
K.s = [2 2 4 1 3 3 4 2];
vv = vidx(K.s);
vv = cell2mat(vv(:));
mm = midx(K.s);
mm = cell2mat(mm(:));
ll = mlidx(K.s);
offset = cumsum ([0, K.s(1:end-1).^2]);
for i = 1:length(offset)
  ll{i} = ll{i} + offset(i);
end
ll = cell2mat(ll(:));
vvll = vlidx(K.s);
offset = cumsum ([0, K.s(1:end-1) .* (K.s(1:end-1) + 1) / 2]);
for i = 1:length(offset)
  vvll{i} = vvll{i} + offset(i);
end
vvll = cell2mat(vvll(:));

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
verifyEqual (testCase, v, [false(offset, 2); vv])
verifyEqual (testCase, m, [false(offset, 2); mm])
verifyEqual (testCase, l, ll + offset)
[v, m, l, vl] = vsdp.sindex (K);
verifyEqual (testCase, v, [false(offset, 2); vv])
verifyEqual (testCase, m,  [])
verifyEqual (testCase, l,  [])
verifyEqual (testCase, vl, vvll + offset)

end
