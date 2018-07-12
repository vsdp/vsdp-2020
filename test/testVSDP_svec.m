function testVSDP_svec (testCase)

% Unsymmetric matrix.
A = [...
  11, 12, 13, 14; ...
  21, 22, 23, 24; ...
  31, 32, 33, 34; ...
  41, 42, 43, 44];
A_svec_1_sym   = A([1, 5:6, 9:11, 13:16])';
A_svec_1_unsym = diag (diag (A)) + 0.5 * (triu(A, 1) + tril(A, -1)');
A_svec_1_unsym = A_svec_1_unsym([1, 5:6, 9:11, 13:16])';

% Symmetric matrix.
B = [...
  1, 2, 4; ...
  2, 3, 5; ...
  4, 5, 6];
B_svec = @(mu) [1 2*mu 3 4*mu 5*mu 6]';
B_svec_1 = B_svec(1);
B_svec_2 = B_svec(2);
B_svec_sqrt2 = B_svec(sqrt(2));

% Try multiple vectorized matrices.
K.s = [4, 3];
C = [A(:); B(:)];
C_svec_1_sym   = [A_svec_1_sym;   B_svec_1];
C_svec_1_unsym = [A_svec_1_unsym; B_svec_1];

D = repmat (C, 1, 3);
D_svec_1_sym   = repmat (C_svec_1_sym,   1, 3);
D_svec_1_unsym = repmat (C_svec_1_unsym, 1, 3);

verifyEqual (testCase, vsdp.svec ([], A   , 1, 'sym'), A_svec_1_sym)
verifyEqual (testCase, vsdp.svec ([], A(:), 1, 'sym'), A_svec_1_sym)
verifyEqual (testCase, vsdp.svec ([], B   ), B_svec_sqrt2)
verifyEqual (testCase, vsdp.svec ([], B(:)), B_svec_sqrt2)
verifyEqual (testCase, vsdp.svec ([], B   , 1, 'sym'), B_svec_1)
verifyEqual (testCase, vsdp.svec ([], B(:), 1, 'sym'), B_svec_1)
verifyEqual (testCase, vsdp.svec ([], B   , 2, 'sym'), B_svec_2)
verifyEqual (testCase, vsdp.svec ([], B(:), 2, 'sym'), B_svec_2)
verifyEqual (testCase, vsdp.svec (K,  C   , 1, 'sym'), C_svec_1_sym)
verifyEqual (testCase, vsdp.svec (K,  D   , 1, 'sym'), D_svec_1_sym)
verifyEqual (testCase, vsdp.svec ([], A   , 1, 'unsym'), A_svec_1_unsym)
verifyEqual (testCase, vsdp.svec ([], A(:), 1, 'unsym'), A_svec_1_unsym)
verifyEqual (testCase, vsdp.svec ([], B   , 1, 'unsym'), B_svec_1)
verifyEqual (testCase, vsdp.svec ([], B(:), 1, 'unsym'), B_svec_1)
verifyEqual (testCase, vsdp.svec ([], B   , 2, 'unsym'), B_svec_2)
verifyEqual (testCase, vsdp.svec ([], B(:), 2, 'unsym'), B_svec_2)
verifyEqual (testCase, vsdp.svec (K,  C   , 1, 'unsym'), C_svec_1_unsym)
verifyEqual (testCase, vsdp.svec (K,  D   , 1, 'unsym'), D_svec_1_unsym)

% Try preserving sparsity.
As = sparse (A);
Bs = sparse (B);
Cs = sparse (C);
Ds = sparse (D);

verifyEqual (testCase, vsdp.svec ([], As   , 1, 'sym'), sparse (A_svec_1_sym))
verifyEqual (testCase, vsdp.svec ([], As(:), 1, 'sym'), sparse (A_svec_1_sym))
verifyEqual (testCase, vsdp.svec ([], Bs   ), sparse (B_svec_sqrt2))
verifyEqual (testCase, vsdp.svec ([], Bs(:)), sparse (B_svec_sqrt2))
verifyEqual (testCase, vsdp.svec ([], Bs   , 1, 'sym'), sparse (B_svec_1))
verifyEqual (testCase, vsdp.svec ([], Bs(:), 1, 'sym'), sparse (B_svec_1))
verifyEqual (testCase, vsdp.svec ([], Bs   , 2, 'sym'), sparse (B_svec_2))
verifyEqual (testCase, vsdp.svec ([], Bs(:), 2, 'sym'), sparse (B_svec_2))
verifyEqual (testCase, vsdp.svec (K,  Cs   , 1, 'sym'), sparse (C_svec_1_sym))
verifyEqual (testCase, vsdp.svec (K,  Ds   , 1, 'sym'), sparse (D_svec_1_sym))
verifyEqual (testCase, vsdp.svec ([], As   , 1, 'unsym'), sparse (A_svec_1_unsym))
verifyEqual (testCase, vsdp.svec ([], As(:), 1, 'unsym'), sparse (A_svec_1_unsym))
verifyEqual (testCase, vsdp.svec ([], Bs   , 1, 'unsym'), sparse (B_svec_1))
verifyEqual (testCase, vsdp.svec ([], Bs(:), 1, 'unsym'), sparse (B_svec_1))
verifyEqual (testCase, vsdp.svec ([], Bs   , 2, 'unsym'), sparse (B_svec_2))
verifyEqual (testCase, vsdp.svec ([], Bs(:), 2, 'unsym'), sparse (B_svec_2))
verifyEqual (testCase, vsdp.svec (K,  Cs   , 1, 'unsym'), sparse (C_svec_1_unsym))
verifyEqual (testCase, vsdp.svec (K,  Ds   , 1, 'unsym'), sparse (D_svec_1_unsym))

% Try preserving interval datatype.
if (exist ('intval', 'file') == 2)
  Ai = intval (A);
  Bi = intval (B);
  Ci = intval (C);
  Di = intval (D);
  
  verifyEqual (testCase, vsdp.svec ([], Ai   , 1, 'sym'), intval (A_svec_1_sym))
  verifyEqual (testCase, vsdp.svec ([], Ai(:), 1, 'sym'), intval (A_svec_1_sym))
  % Due to interval arithmetic, no test for B_svec_sqrt2.
  verifyEqual (testCase, vsdp.svec ([], Bi   , 1, 'sym'), intval (B_svec_1))
  verifyEqual (testCase, vsdp.svec ([], Bi(:), 1, 'sym'), intval (B_svec_1))
  verifyEqual (testCase, vsdp.svec ([], Bi   , 2, 'sym'), intval (B_svec_2))
  verifyEqual (testCase, vsdp.svec ([], Bi(:), 2, 'sym'), intval (B_svec_2))
  verifyEqual (testCase, vsdp.svec (K,  Ci   , 1, 'sym'), intval (C_svec_1_sym))
  verifyEqual (testCase, vsdp.svec (K,  Di   , 1, 'sym'), intval (D_svec_1_sym))
  verifyEqual (testCase, vsdp.svec ([], Ai   , 1, 'unsym'), intval (A_svec_1_unsym))
  verifyEqual (testCase, vsdp.svec ([], Ai(:), 1, 'unsym'), intval (A_svec_1_unsym))
  verifyEqual (testCase, vsdp.svec ([], Bi   , 1, 'unsym'), intval (B_svec_1))
  verifyEqual (testCase, vsdp.svec ([], Bi(:), 1, 'unsym'), intval (B_svec_1))
  verifyEqual (testCase, vsdp.svec ([], Bi   , 2, 'unsym'), intval (B_svec_2))
  verifyEqual (testCase, vsdp.svec ([], Bi(:), 2, 'unsym'), intval (B_svec_2))
  verifyEqual (testCase, vsdp.svec (K,  Ci   , 1, 'unsym'), intval (C_svec_1_unsym))
  verifyEqual (testCase, vsdp.svec (K,  Di   , 1, 'unsym'), intval (D_svec_1_unsym))
else
  warning ('VSDP:testVSDP_svec:noIntval', ...
    'testVSDP_svec: Skip intval tests.');
end

end
