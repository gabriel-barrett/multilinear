I want to show that
T(x1,...) = 𝚺{y,z} A(x1,..., y1,..., z1,...) * W(y1,...) * W(z1,...) = 0
in the hypercube which means T is identically zero (being multilinear)

I can try to prove that T(r1,...) = 0 for random r. I can then use the sumcheck
protocol to reduce the claim to an evaluation of
A(r1,..., s1,..., t1,...) * W(s1,...) * W(t1,...)

A polynomial commitment scheme will give me W(s1,...) and W(t1,...)
without revealing W, although it must also be shown that W(0,...) = 1

To evaluate A we can do it directly since the verifier knows it.
This is not succinct, except if we know A has a structure. For instance,
if A is a repetition of the same constraints for the same section of W,
then you can actually compute it succinctly.

Let me show how this works. Suppose A is a single "step". To extend A to get
two steps of computation, we can do

Aext(r1,..., s1,..., t1,... | a, b, c) = A(r1,..., s1,..., t1,...) * (abc + (1-a)(1-b)(1-c))

More generally, this is the shape the extended A
Aext(I1, I2, I3 | J1, J2, J3) = A(I1, I2, I3) * 𝚷{j1, j2, j3} (j1 * j2 * j3 + (1-j1) * (1-j2) * (1-j3))

Questions:
1) How do you introduce GKR?
2) How do you introduce switches? (this one is not so important)
3) Can we overlap constraints (for input/output)?

Actually we can extend this constraint system to a more general
𝚺{j,k} (Aijk * Wj * Wk + Bijk * Wj) = Ci
and this will have basically the same properties, but there's no need to fix a column
to be equal to 1.
