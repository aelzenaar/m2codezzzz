-----------------
-- polygen.m2 ---
-----------------
--
-- A script to compute the defining polynomials for varieties of
-- real (t,d,n)-designs, and to compute sum-of-squares decompositions
-- using the SumsOfSquares package.
--

needsPackage "SumsOfSquares"

kk = QQ
kkDim = 1

-- compute coefficient
computeCoefficient =
  (n,d,t,reDim) -> product toList apply (
              (0..(t-1)),
              i -> (reDim + 2*i)/(reDim*d + 2*i) )

-- parameters to list (t,d,n)
params = [(1,2,2),(1,3,3),(2,2,3),(2,3,6),(2,4,12),(2,5,20),(2,6,24),(2,7,28),
          (3,2,4),(3,3,16)]

for i in params do (
  t = i_0; d = i_1; n = i_2;
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
  << "(" << t << ", " << d << ", " << n << "): " << result#Status;
  cleanedSOS = clean(1e-5,sosPoly result);
  << ", #terms = " << length cleanedSOS;
  << ", #terms/term = " << (for term in gens cleanedSOS list length terms term);
  << ", degree = " << first degree (gens cleanedSOS)_1;
  << "\n" << flush;
)

exit
