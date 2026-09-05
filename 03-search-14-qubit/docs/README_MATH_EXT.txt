=========================================================================================
stab14_ext -- mathematics behind the seven-extra-generator optimiser
=========================================================================================

PROBLEM.  S0 = <g1,g2,g3,g4> is fixed (weight 8, commuting, rank 4).  U is the set of
Pauli errors of qubit weight 1..4 with zero S0-syndrome; the program rebuilds U from the
four generators and refuses to run unless |U| = 5430 with weights 12/69/528/4821.
Find h1..h7, each of Pauli weight exactly 10 (Problem A) or in {10,12} (Problem B), all
commuting with S0 and with each other, with rank(H0;H1) = 11, maximising
        M = #{ E in U : <h_a,E> = 1 for at least one a }.

-----------------------------------------------------------------------------------------
1.  EVERYTHING LIVES IN A 20-DIMENSIONAL SYMPLECTIC SPACE
-----------------------------------------------------------------------------------------
Both the candidates h and the targets E lie in C(S0) = S0^perp (dim 24).  For h,E in C(S0)
and s in S0,  <h+s,E> = <h,E> + <s,E> = <h,E>, and symmetrically in E.  So coverage,
mutual commutation and independence are all functions of the CLASSES in
        V = C(S0)/S0 ,   dim V = 24 - 4 = 20 ,
on which the symplectic form is non-degenerate (the radical of the form on C(S0) is
exactly S0 -- the program constructs the basis and verifies this).

The program builds a symplectic basis u_1..u_10, w_1..w_10 of V by symplectic
Gram-Schmidt and encodes a class as a 20-bit integer
        class(v) = ( <w_i,v> )_{i<10}  |  ( <u_i,v> )_{i<10} << 10 ,
so that the form becomes  <v,v'> = parity( (a & b') ^ (b & a') )  -- two ANDs and a
popcount.  Candidate generation, commutation and independence are then 20-bit operations.

2.  THE OBJECTIVE DEPENDS ONLY ON THE SUBSPACE, AND IS A 127-TERM SUM
-----------------------------------------------------------------------------------------
Mutual commutation + independence = the seven classes span a 7-dimensional TOTALLY
ISOTROPIC subspace D <= V (rank 11 for the whole stabilizer is exactly dim D = 7).  A
target E is undetected iff its class lies in D^perp.  Writing
        cov1(v) = #{ E in U : <v,E> = 1 } ,
Poisson summation over V gives, for EVERY 7-dimensional D,

        M(D)  =  (1/64) * SUM_{v in D, v != 0} cov1(v) .                     (*)

  Proof.  #undetected = sum_{p in D^perp} f(p) where f is the multiplicity function of the
  target classes.  With sigma the half-swap that turns the symplectic form into the dot
  product, D^perp = sigma(D)^{perp,dot}, so sum_{p in D^perp} f(p) = 2^-7 sum_{v in
  sigma(D)} fhat(v).  And fhat(sigma(w)) = sum_E (-1)^{<w,E>} = 5430 - 2 cov1(w).  Summing
  the 128 terms gives #undetected = 5430 - (1/64) sum_{v in D\0} cov1(v).

Consequences used everywhere in the program:
  * a complete candidate set is scored with 127 lookups in a 2^20 table -- no 5430-bit
    union and no symplectic products;
  * the greedy gain of adding v to a d-dimensional D is EXACTLY G(v) = sum_{u in D}
    cov1(v^u), so greedy/beam ranking optimises the true objective, not a surrogate;
  * cov1 is obtained for all 2^20 classes at once by one Walsh-Hadamard transform
    (20 * 2^20 additions, milliseconds).

The 5430-bit bitset engine required by the specification is implemented as well (coverage
masks are linear in the class, so mask(v) is the XOR of 20 basis masks) and is used for
INDEPENDENT verification: on every reported solution the union popcount, the brute-force
recount over all 91770 errors, and identity (*) are printed side by side and must agree.

3.  CODE-THEORETIC READING -- THIS IS THE [[14,3,5]] QUESTION
-----------------------------------------------------------------------------------------
L = S0 + span(h_1..h_7) is an 11-dimensional isotropic subspace, i.e. a [[14,3]] stabilizer
code containing S0, and the undetected targets are exactly the weight-<=4 Paulis of N(L).
The same Krawtchouk identity as in the first-stage program gives
        #undetected = (1/2048) [ 91770 + SUM_{h in L, h != 0} P(wt h) ] ,
        P(m) = sum_{w=1..4} sum_k C(m,k)(-1)^k C(14-m,w-k) 3^{w-k}.
Sanity check: for L = S0 (dim 4) the same formula reads (1/16)[91770 + 15*(-326)] = 5430,
the target count itself.  Therefore

        M = 5430  <=>  N(L) contains no non-zero Pauli of weight <= 4
                  <=>  L is a PURE [[14,3,5]] stabilizer code containing S0.

Two structural corollaries the program uses and prints:
  * B_1 (undetected weight-1 targets) = 3*n0 + n1, where n_d counts the physical qubits
    whose column span in the 11-generator matrix has dimension d.  Full coverage forces
    every one of the 14 columns to be 2-dimensional (--full-columns searches that region).
  * B_2 = 0 additionally forces every pair of columns to intersect trivially, B_3 = 0 every
    triple to be independent, and so on: the classical "any 4 columns direct" condition.

4.  SEARCH
-----------------------------------------------------------------------------------------
Candidates.  A Gray-code sweep of all 2^24 elements of C(S0) records, for every class, the
lexicographically smallest weight-10 and weight-12 representative.  Result: 3 693 519
weight-10 Paulis reaching 514 153 classes and 3 031 911 weight-12 Paulis reaching 516 615
classes (their union is 524 796, so the {10,12} relaxation adds only ~10 600 new classes).
Generation is centralizer-aware by construction -- no weight-10 Pauli outside C(S0) is ever
built.

Heuristics.  (i) randomized greedy on the exact incremental gain G(v); (ii) steepest-descent
1-opt where a cached per-candidate "commutes with which of the seven" byte mask is updated
one column at a time, so a sweep is seven linear passes; (iii) ruin-and-recreate LNS that
drops 1-3 generators and rebuilds them; (iv) a beam search over subspaces de-duplicated by a
hash of the sorted span (this removes the many bases of the same subspace).  Every candidate
set passes a hard gate (pairwise orthogonality + rank 7 by Gaussian elimination) before it
can become the incumbent.

Exact branch and bound.  Increasing candidate indices (the h_a are unordered, so this is
sound and removes the 7! orderings), incremental span sums, allowed-set filtering by
orthogonality and span membership, and the admissible bound
        M  <=  ( S_d + max_{v allowed} G(v) + (128 - 2^{d+1}) * COVMAX ) / 64 ,
which relaxes isotropy, independence and weight-realisability, hence is valid.

HONEST LIMITATION.  That bound can prune only while COVMAX <= 64*incumbent/127.  With the
incumbent 5423 the threshold is 2733.0 while the actual COVMAX = 2808, so the relaxation
cannot prune anywhere near the root: the exhaustive search over 514 153 candidate classes is
out of reach and the program says GLOBAL OPTIMUM NOT PROVEN.  The only bound that is
genuinely proved is the trivial M <= 5430 (the two computed relaxations give 5572).

5.  RESULT
-----------------------------------------------------------------------------------------
        M_10    >= 5423 / 5430   (99.871087 %),  7 targets left undetected
        M_10/12 >= 5423 / 5430   -- the relaxation buys nothing, and the tie-break then
                                   prefers the all-weight-10 solution
5423 is reached from every seed tried (1, 777, 12345, 20250830), by plain greedy restarts,
by LNS and by beam search, and also bounds the --full-columns region (5422 there), so it is
a very robust local optimum -- but it is NOT proven optimal.

The seven leftovers form a group-like family, as they must: N(L) is a group, so a product of
two undetected targets is undetected whenever it still has weight <= 4.  In the reported
solution they are X0 ; X1X4 ; X0X1X4 ; X3X7X8 ; X0X3X7X8 ; and two more weight-4 Paulis.
