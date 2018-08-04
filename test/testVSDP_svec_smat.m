function testVSDP_svec_smat (testCase)

% Unsymmetric matrix.
A = [...
  11, 12, 13, 14; ...
  21, 22, 23, 24; ...
  31, 32, 33, 34; ...
  41, 42, 43, 44];
A_sym = triu (A) + triu (A, 1)';
A_svec_1_sym   = A([1, 5:6, 9:11, 13:16])';
A_svec_1_unsym = diag (diag (A)) + 0.5 * (triu(A, 1) + tril(A, -1)');
A_svec_1_unsym = A_svec_1_unsym([1, 5:6, 9:11, 13:16])';

% Symmetric matrix.
B = [...
  1, 2, 4; ...
  2, 3, 5; ...
  4, 5, 6];
B_svec = @(mu) [1 2*mu 3 4*mu 5*mu 6]';
B_svec_1     = B_svec(1);
B_svec_2     = B_svec(2);
B_svec_sqrt2 = B_svec(sqrt(2));

% Try multiple vectorized matrices.
K.s   = [4, 3];
C     = [A(:)    ; B(:)];
C_sym = [A_sym(:); B(:)];
C_svec_1_sym   = [A_svec_1_sym;   B_svec_1];
C_svec_1_unsym = [A_svec_1_unsym; B_svec_1];

D     = repmat (C,     1, 3);
D_sym = repmat (C_sym, 1, 3);
D_svec_1_sym   = repmat (C_svec_1_sym,   1, 3);
D_svec_1_unsym = repmat (C_svec_1_unsym, 1, 3);

fun = {@(x) x, @(x) sparse(x)};
if (exist ('intval', 'file') == 2)
  fun{end + 1} = @(x) intval (x);
else
  warning ('VSDP:testVSDP_svec_smat:noIntval', ...
    'testVSDP_svec_smat: Skip intval tests.');
end

% SVEC
for i = 1:length(fun)
  f = fun{i};
  verifyEqual (testCase, vsdp.svec ([], f(A)   , 1, 'sym'), f(A_svec_1_sym))
  verifyEqual (testCase, vsdp.svec ([], f(A(:)), 1, 'sym'), f(A_svec_1_sym))
  % Skip these test for intervals 'fun{3}', due to sqrt(2).
  if (i ~= 3)
    verifyEqual (testCase, vsdp.svec ([], f(B)   ), f(B_svec_sqrt2))
    verifyEqual (testCase, vsdp.svec ([], f(B(:))), f(B_svec_sqrt2))
  end
  verifyEqual (testCase, vsdp.svec ([], f(B)   , 1, 'sym'), f(B_svec_1))
  verifyEqual (testCase, vsdp.svec ([], f(B(:)), 1, 'sym'), f(B_svec_1))
  verifyEqual (testCase, vsdp.svec ([], f(B)   , 2, 'sym'), f(B_svec_2))
  verifyEqual (testCase, vsdp.svec ([], f(B(:)), 2, 'sym'), f(B_svec_2))
  verifyEqual (testCase, vsdp.svec (K,  f(C)   , 1, 'sym'), f(C_svec_1_sym))
  verifyEqual (testCase, vsdp.svec (K,  f(D)   , 1, 'sym'), f(D_svec_1_sym))
  verifyEqual (testCase, vsdp.svec ([], f(A)   , 1, 'unsym'), f(A_svec_1_unsym))
  verifyEqual (testCase, vsdp.svec ([], f(A(:)), 1, 'unsym'), f(A_svec_1_unsym))
  verifyEqual (testCase, vsdp.svec ([], f(B)   , 1, 'unsym'), f(B_svec_1))
  verifyEqual (testCase, vsdp.svec ([], f(B(:)), 1, 'unsym'), f(B_svec_1))
  verifyEqual (testCase, vsdp.svec ([], f(B)   , 2, 'unsym'), f(B_svec_2))
  verifyEqual (testCase, vsdp.svec ([], f(B(:)), 2, 'unsym'), f(B_svec_2))
  verifyEqual (testCase, vsdp.svec (K,  f(C)   , 1, 'unsym'), f(C_svec_1_unsym))
  verifyEqual (testCase, vsdp.svec (K,  f(D)   , 1, 'unsym'), f(D_svec_1_unsym))
  
  % Test rescaling.
  warning ('off', 'VSDP:svec:justScale');
  verifyEqual (testCase, vsdp.svec ([], f(B_svec_1),   1), f(B_svec_1))
  verifyEqual (testCase, vsdp.svec ([], f(B_svec_2),   1), f(B_svec_2))
  verifyEqual (testCase, vsdp.svec ([], f(B_svec_1),   2), f(B_svec_2))
  verifyEqual (testCase, vsdp.svec ([], f(B_svec_2), 1/2), f(B_svec_1))
  warning ('on', 'VSDP:svec:justScale');
end

% SMAT
for i = 1:length(fun)
  f = fun{i};
  verifyEqual (testCase, vsdp.smat ([], f(A_svec_1_sym), 1  ), f(A_sym))
  % Skip these test for intervals 'fun{3}', due to sqrt(2).
  if (i ~= 3)
    verifyEqual (testCase, vsdp.smat ([], f(B_svec_sqrt2)   ), f(B), ...
      'AbsTol', 8*eps ())
  end
  verifyEqual (testCase, vsdp.smat ([], f(B_svec_1),     1  ), f(B))
  verifyEqual (testCase, vsdp.smat ([], f(B_svec_2),     1/2), f(B))
  verifyEqual (testCase, vsdp.smat (K,  f(C_svec_1_sym), 1  ), f(C_sym))
  verifyEqual (testCase, vsdp.smat (K,  f(D_svec_1_sym), 1  ), f(D_sym))
end

% Test matrices without semidefinite cones.
clear K
K.l = size (D, 1);
for i = 1:length(fun)
  verifyEqual (testCase, vsdp.svec (K, f(D)   ), f(D))
  verifyEqual (testCase, vsdp.svec (K, f(D), 1), f(D))
  verifyEqual (testCase, vsdp.svec (K, f(D), 1, 'sym'  ), f(D))
  verifyEqual (testCase, vsdp.svec (K, f(D), 1, 'unsym'), f(D))
  verifyEqual (testCase, vsdp.smat (K, f(D)   ), f(D))
  verifyEqual (testCase, vsdp.smat (K, f(D), 1), f(D))
end

end
