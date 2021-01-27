-------------------
-- factoring.m2 ---
-------------------
--
-- A script to compute factorisations of the defining polynomials for
-- varieties of real (t,d,n)-designs.
--

kk = QQ
kkDim = 1

-- compute coefficient
computeCoefficient =
  (n,d,t,reDim) -> product toList apply (
              (0..(t-1)),
              i -> (reDim + 2*i)/(reDim*d + 2*i) )

-- parameters to list (t,d,n)
T = 10
N = 10
D = 10

for t from 1 to T do
  for d from 1 to D do
    for n from 1 to N do (
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

      factors = factor definingPolynomial;
      factorCount =
        if definingPolynomial == 0 then 0
          else if first degree value last factors == 0 then #factors - 1
            else #factors;
      << t << ", " << d << ", " << n << ", " << factorCount << "\n" << flush;
)

exit
