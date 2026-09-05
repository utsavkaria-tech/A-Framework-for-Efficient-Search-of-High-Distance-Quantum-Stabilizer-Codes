=========================================================================================
stab14_2stage -- mathematics of the two-stage optimiser
    STAGE 1: three weight-8 generators        (exhaustive, PROVEN optimal)
    STAGE 2: eight weight-10 generators on the frozen Stage-1 stabilizer
=========================================================================================
The two problems are kept strictly separate.  Stage 1 maximises detection over all 91770
weight-1..4 Pauli errors.  Stage 2 freezes that solution and maximises coverage of ONLY the
set U_3 that Stage 1 missed.  The final 11-generator numbers are reported, never optimised.

-----------------------------------------------------------------------------------------
0.  THE IDENTITY BOTH STAGES ARE BUILT ON
-----------------------------------------------------------------------------------------
For an r-dimensional isotropic subspace L of F_2^28 (a stabilizer group of rank r), the
Walsh transform over the syndrome space of the per-qubit error-counting function is
F_j(h) = 3 if h acts trivially on qubit j and -1 otherwise.  The number of weight-w errors
with zero syndrome is then the elementary symmetric polynomial e_w(F_1,...,F_14), which
depends on h only through wt(h):

    #{E : wt(E)=w, s(E)=0} = (1/2^r) SUM_{h in L} Ew[wt h][w],
    Ew[m][w] = SUM_k C(m,k)(-1)^k C(14-m,w-k) 3^{w-k},   P(m) = SUM_{w=1..4} Ew[m][w].

    #undetected(weight 1..4) = (1/2^r) [ 91770 + SUM_{h in L, h != 0} P(wt h) ].

P(m), m = 0..14:  91770 57914 34154 18314 8474 2970 394 -406 -326 -6 170 74 -166 -166 714.
Both stages therefore reduce to shaping the WEIGHT ENUMERATOR of the stabilizer group.

=========================================================================================
STAGE 1 -- exhaustive, and the optimum is proven
=========================================================================================
1.  COLUMN-TYPE REDUCTION OVER F_2^3.  Qubit j contributes (x_j,z_j) in F_2^3 x F_2^3.  A
    single-qubit Clifford acts as GL(2,2) = S_3 on that pair.  It
      * preserves each row's Pauli weight ((x[i],z[i]) != 0 iff the image is != 0),
      * preserves the commutation form exactly, since
            (ax+bz)(cx+dz)^T + (cx+dz)(ax+bz)^T = (ad+bc)(x z^T + z x^T) = x z^T + z x^T,
      * permutes the three single-qubit syndromes {z, x, x+z} setwise, and the objective
        sums over all three Paulis of every support qubit, so it is invariant.
    Hence the complete invariant of a column is the SUBSPACE W_j = span{x_j,z_j} <= F_2^3
    of dimension 0, 1 or 2: only 1 + 7 + 7 = 15 column types instead of 4^14 Pauli strings.

2.  CONSTRAINTS AND OBJECTIVE IN COLUMN LANGUAGE.  With U_j = W_j^perp and
    c(chi) = #{j : chi in U_j}, the group element g_chi has weight m(chi) = 14 - c(chi), so
      weight 8 for the three generators  <=>  c(e_i) = 6 for i = 1,2,3   (3 linear equations)
      mutual commutation                 <=>  XOR_j plucker(W_j) = 0 in F_2^3, where
                                              plucker(W) = x z^T + z x^T is basis-independent
      rank 3                             <=>  c(chi) <= 13 for every chi != 0
      objective                          <=>  minimise SUM_{chi != 0} P(14 - c(chi)).

3.  SYMMETRIES REMOVED.  S_14 qubit permutations (the state is a MULTISET of column types,
    built in non-decreasing type order so each multiset is visited once); (S_3)^14 local
    Cliffords (columns are subspaces); S_3 generator permutations (lex-minimality test on
    the sorted prefix -- valid because for any permutation pi and prefix p of the sorted
    multiset S, the k smallest entries of sorted(pi(S)) are <= sorted(pi(p)) entrywise).

4.  BOUND.  At a node with degree vector c, r columns left, all of type index >= t0:
    hard feasibility on the three weight equations, the rank cap c <= 13, plus the maximum
    of (a) SUM over the four non-generator points of min_{0<=d<=r} P(14-c-d) and (b) a
    knapsack bound using the fact that every remaining column adds exactly |U|-1 in {7,3,1}
    to SUM_chi c(chi).  Both relax an exact minimisation, so pruning on bound >= incumbent
    never discards a better solution.

5.  RESULT.  Multisets of 14 out of 15 types number C(28,14) = 40 116 600; the constraints
    and the bound cut the tree to 19 553 nodes and the search finishes in 8 ms.  The
    optimum is therefore PROVEN:

        M_3 = 80584 / 91770 = 87.81083143 %,   |U_3| = 11186
        weights 1/2/3/4 undetected: 14 / 119 / 1148 / 9905

    and it is attained by a group in which ALL SEVEN non-identity elements have weight 8
    (SUM P = 7*(-326) = -2282, so #undetected = (91778 - 2282)/8 - 1 = 11186).  Ties are
    broken in favour of a representative that leaves no qubit idle, then S_3-lex-minimum.

=========================================================================================
STAGE 2 -- the Stage-1 optimum is frozen
=========================================================================================
6.  QUOTIENT SPACE.  Candidates and targets both lie in C(S3) = S3^perp, dim 28-3 = 25, and
    <h,E> is unchanged by shifts of h or E by an element of S3.  Everything -- coverage,
    mutual commutation, independence -- is therefore a function of the class in

        V = C(S3)/S3 ,   dim V = 25 - 3 = 22   (4 194 304 classes),

    on which the form is non-degenerate; the program builds a symplectic basis
    u_1..u_11, w_1..w_11 by symplectic Gram-Schmidt and VERIFIES that the radical of the
    form on C(S3) is exactly S3.  A class is a 22-bit integer
        class(v) = ( <w_i,v> )_i | ( <u_i,v> )_i << 11
    and the symplectic form becomes parity( (a & b') ^ (b & a') ).

7.  ISOTROPIC-SUBSPACE VIEW.  Mutual commutation plus independence means the eight classes
    span an 8-dimensional TOTALLY ISOTROPIC subspace D <= V, and rank 11 for the complete
    stabilizer is exactly dim D = 8.  A target is undetected iff its class lies in D^perp,
    so the objective depends on the SUBSPACE, not on which basis of it we report.

8.  THE COVERAGE IDENTITY (normalisation re-derived for d = 8, not copied).  With
    cov1(v) = #{E in U_3 : <v,E> = 1},

        #covered  =  (1/2^{d-1}) SUM_{v in D\{0}} cov1(v)      for dim D = d
                  =  (1/128) SUM over the 255 non-zero elements  for d = 8.

    Derivation: #uncovered = sum_{p in D^perp} f(p), where f is the multiplicity function of
    the target classes.  D^perp is the dot-product dual of sigma(D) (sigma = swap of the two
    11-bit halves), so Poisson summation gives #uncovered = 2^{-d} sum_{v in sigma(D)}
    fhat(v), and fhat(sigma(w)) = |U_3| - 2 cov1(w).  Summing the 2^d terms leaves
    |U_3| - 2^{1-d} sum_{w != 0} cov1(w).  (For the earlier d = 7 problem this constant was
    1/64; here it is 1/128.)  The identity is checked against the explicit |U_3|-bit union
    on real 8-generator sets in the self-test, and again on every reported solution.

    Consequences: scoring a complete 8-set is 255 lookups in a 2^22 table; the greedy gain
    of adding v to a d-dimensional D is exactly SUM_{u in D} cov1(v^u), so greedy and beam
    ranking optimise the true objective rather than a surrogate; and cov1 for all 2^22
    classes comes from ONE Walsh-Hadamard transform (22 * 2^22 additions).

9.  CANDIDATES.  A Gray-code sweep of all 2^25 elements of C(S3) records, for each class,
    the lexicographically smallest weight-10 representative: 7 387 821 weight-10 Paulis
    reaching 2 016 879 of the 4 194 303 non-zero classes.  Generation is centralizer-aware
    by construction -- no weight-10 Pauli outside C(S3) is ever built.

10. SEARCH.  Randomized greedy on the exact gain; steepest-descent 1-opt with a cached
    per-candidate byte mask recording which of the eight current generators it commutes with
    (only the swapped column is recomputed); ruin-and-recreate LNS dropping 1-3 generators;
    optional beam search over subspaces de-duplicated by a hash of the sorted span.  Every
    candidate set passes a hard gate (pairwise orthogonality plus rank 8 by Gaussian
    elimination) before it can become the incumbent.  Branch and bound enumerates increasing
    candidate indices (h_1 < ... < h_8 -- the eight are unordered, so this removes 8!
    orderings) with the admissible bound

        coverage <= ( S_d + max_{v allowed} G(v) + (256 - 2^{d+1}) * COVMAX ) / 128 ,

    which relaxes isotropy, independence and weight-10 realisability.

11. HONEST LIMITATION.  That bound prunes only while COVMAX <= 128 * incumbent / 255.  With
    an incumbent near 11 179 the threshold is about 5610 while the actual COVMAX = 5796, so
    the relaxation cannot prune near the root: exhausting 2 016 879 candidate classes at
    depth 8 is out of reach and Stage 2 reports GLOBAL OPTIMUM NOT PROVEN.  The only proved
    bound is the trivial M_8 <= |U_3| = 11186 (the Fourier relaxation gives 11 549, worse).

12. WHY FULL COVERAGE IS HARD.  M_8 = |U_3| would mean the 11-dimensional stabilizer L has
    no non-zero Pauli of weight <= 4 in N(L), i.e. L is a PURE [[14,3,5]] code.  The
    still-undetected errors necessarily form a group-like family, because N(L) is a group:
    a product of two undetected errors is undetected whenever it still has weight <= 4.

13. RESULTS.
        M_3     = 80584 / 91770   PROVEN GLOBAL OPTIMUM
        |U_3|   = 11186
        M_8     >= 11179 / 11186  BEST KNOWN, optimum not proven (7 leftovers, see 14)
        M_final =  91763 / 91770  = 99.99237223 %
        M_final = M_3 + M_8       (verified independently by re-enumerating all 91770 errors)

14. WHAT THE SEVEN LEFTOVERS ACTUALLY ARE.  In the reported solution the still-undetected
    set is exactly the group generated by the three weight-1 errors X_4, X_9, X_13:
        X4 ; X9 ; X13 ; X4X9 ; X4X13 ; X9X13 ; X4X9X13     (2^3 - 1 = 7 elements)
    i.e. three physical qubits whose column in the 11-generator matrix spans only one
    dimension.  This is forced by the algebra: B_1 = 3*n0 + n1 counts those qubits, and
    N(L) being a group closes the set under products of weight <= 4.  Steering the search to
    kill all deficient columns (--full-columns) does reach 0 of them, but only 11178 coverage
    (8 leftovers of a different shape), so the unconstrained objective remains better.
