-----------------
-- polygen.m2 ---
-----------------
--
-- A script to compute the defining polynomials for varieties of
-- real (t,d,n)-designs, and to compute sum-of-squares decompositions
-- using the SumsOfSquares package.
--

needsPackage "SumsOfSquares"
changeSolver("MOSEK","/Users/aelz176/Desktop/mosek/9.2/tools/platform/osx64x86/bin/mosek")

kk = QQ
kkDim = 1

-- compute coefficient
computeCoefficient =
  (n,d,t,reDim) -> product toList apply (
              (0..(t-1)),
              i -> (reDim + 2*i)/(reDim*d + 2*i) )

-- parameters to list (t,d,n); note d = 1 is trivial
params = sort toList((set (1..10)**set(2..10)**set({1}))/splice)

for i in params do (
  t = i_0; d = i_1; n = i_2;
  << "(" << t << ", " << d << ", " << n << "): " << flush;
  c = computeCoefficient (n, d, t, kkDim);

  symbols = (x_1 .. x_(d*n));
  R = kk[symbols];
  pointMatrix = matrix pack (generators R, n);

  dotProduct =
    (i,j) -> sum apply (entries pointMatrix_{i,j}, product);

  q = sum flatten table ( toList (0..n-1), toList (0..n-1),
                            (i,j) -> (dotProduct(i,j))^(2*t) );
  p = sum apply         ( toList (0..n-1),
                             i    -> (dotProduct(i,i))^t );

  definingPolynomial = q - c*p^2;

  result = solveSOS(definingPolynomial);
  cleanedSOS = clean(1e-5,sosPoly result);
  << result#Status << ", #terms = " << length cleanedSOS;
  << ", degrees = " << unique ((gens cleanedSOS)/degree/first);
  << "\n" << flush;
)

exit
