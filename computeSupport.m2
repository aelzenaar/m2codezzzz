------------------------
-- computeSupport.m2 ---
------------------------
--
-- A script to compute "supp f_t,d,n" for various parameters, placing the output
-- in the file "computeSupport.out".
--
-- Usage: M2 --script computeSupport.m2 <t> <d> <n>
--

needsPackage "SpechtModule" -- for permutationMatrix

t = value scriptCommandLine_1
d = value scriptCommandLine_2
n = value scriptCommandLine_3

assert(d>1)

-- all the cyclic column permutations
cycle = mutableMatrix matrix id_(QQ^(n))
columnPermute(cycle,0,toList append(1..(n-1),0))
cyclePowers = accumulate(times,(n+1):(matrix cycle))

-- Generate matrices with the single non-zero columns consisting of compositions of 2t
generatingMatrices = compositions(d,2*t) /
                      (k -> join(toList ((n*d - d):0) ,k)) /
                        (m -> pack(m,d)) /
                          matrix /
                            transpose
largeCompositions = flatten (table(generatingMatrices,cyclePowers,times)/toList)

-- Generate matrices with the single non-zero columns consisting of compositions of t
generatingMatrices = compositions(d,t) /
                      (k -> join(toList ((n*d - d):0) ,k)) /
                        (m -> pack(m,d)) /
                          matrix /
                            transpose
smallCompositions = flatten (table(generatingMatrices,cyclePowers,times)/toList)


-- First `type' of pairs are sums A + f(A) for A a composition of 2t, f a column permutation
permutationMatrices = permutations(n)/permutationMatrix
sumsType1 = flatten (table(largeCompositions, permutationMatrices, (A,f) -> A + A*f)/toList)

-- Second `type' of pairs are sums 2A + 2B for A, B compositions of t
sumsType2 = unique flatten (table(smallCompositions, smallCompositions, (A,B)->(2*A + 2*B)) /toList)

-- Take the union to compute supp f_t,d,n for d > 1.
exponentVectors = unique join (sumsType1, sumsType2)

-- Make into a format we can do text processing on
f = "computeSupport.out" << ""
scan(unique (exponentVectors/entries/flatten), e -> (f << e << endl) )
f << close

quit
