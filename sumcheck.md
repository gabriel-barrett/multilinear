CLAIM: 𝚺 p(x1, x2, x3) = C

If p is multilinear, then this is pretty easy.
Prover sends:
p1(X) = a1 * x + b1
Verifier sends:
r1

Prover sends:
p2(X) = a2 * x + b2
Verifier sends:
r2

Prover sends:
p3(X) = a3 * x + b3
Verifier sends:
r3

Verifier checks:
1) p3(r3) = p(r1, r2, r3)
2) p2(r2) = p3(0) + p3(1)
3) p1(r1) = p2(0) + p2(1)
4) C      = p1(0) + p1(1)

In terms of the variables:
1) a3 * r3 + b3 = p(r1, r2, r3)
2) a2 * r2 + b2 = a3 + 2 * b3
3) a1 * r1 + b1 = a2 + 2 * b2
4) C            = a1 + 2 * b1

Now, how does the prover efficiently compute these polynomials,
specially if the number of variables is high but it is sparse.

Suppose you have a sparse representation of a multilinear.
Say there are 10 variables but only 3 coefficients that are non-zero.
(MASK1, A1), (MASK2, A2), (MASK3, A3)
where the mask values take the form MASK1 = 0b1001100001, etc

In the case of our SNARK, we have to do sumcheck on this expression
𝚺{y,z} A(r1,..., y1,..., z1,...) * W(y1,...) * W(z1,...) = 0
where A(x1,..., y1,..., z1,...) is sparse. The polynomial we want to
do sumcheck A(r1,..., y1,..., z1,...) * W(y1,...) * W(z1,...) is
sparse still, but it is quadratic.
