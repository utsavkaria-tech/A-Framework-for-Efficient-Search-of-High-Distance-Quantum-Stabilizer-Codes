=========================================================================================
stab14 -- mathematical derivation behind the search-space reduction
=========================================================================================

PROBLEM.  H = [X|Z] in F_2^{4x28}, rows = stabilizer generators g_1..g_4 on 14 qubits.
Constraints: wt(g_i) = w (default 8, Pauli weight, wt(Y)=1), [g_i,g_j] = 0, rank(H) = 4.
Objective:  maximise  N_det = #{ e : 1 <= wt(e) <= 4 , s(e) != 0 },  |universe| = 91770.

-----------------------------------------------------------------------------------------
1.  COLUMN VIEW.  Qubit j contributes (x_j, z_j) in F_2^4 x F_2^4, where x_j[i], z_j[i] are
    the X/Z bits of generator i on qubit j.  The syndrome of a single-qubit error at j is
        X -> z_j ,      Z -> x_j ,      Y -> x_j + z_j
    (from  <g,e> = x_g.z_e + z_g.x_e).  The syndrome of a general error is the XOR of the
    per-qubit contributions.

2.  LOCAL CLIFFORD.  A single-qubit Clifford acts on (x_j,z_j) as GL(2,2) = S_3:
        (x,z) -> (a x + b z, c x + d z),   ad + bc = 1.
    * per-row Pauli weight is preserved: (x[i],z[i]) != (0,0) <=> image != (0,0);
    * the commutation form is preserved exactly:
        (ax+bz)(cx+dz)^T + (cx+dz)(ax+bz)^T = (ad+bc)(x z^T + z x^T) = x z^T + z x^T ;
    * the set {z, x, x+z} of the three single-qubit syndromes is preserved setwise, and the
      objective sums over all three Paulis on every support qubit, so it is invariant.
    Therefore the only invariant that matters is the SUBSPACE
        W_j = span{x_j, z_j} <= F_2^4,      dim W_j in {0,1,2}.
    Number of column types drops from 4^14 patterns / 256 column bit-patterns to 51
    subspaces (1 of dim 0, 15 of dim 1, 35 of dim 2).

3.  GROUP WEIGHTS.  For chi in F_2^4 let g_chi = prod_i g_i^{chi_i}.  g_chi acts trivially
    on qubit j iff chi . x_j = chi . z_j = 0, i.e. iff chi in U_j := W_j^perp.  Hence
        m(chi) := wt(g_chi) = 14 - c(chi),        c(chi) := #{ j : chi in U_j }.
    c is the point-degree function of a multiset of 14 subspaces of PG(3,2).

4.  EXACT OBJECTIVE FROM THE WEIGHT ENUMERATOR.  Give qubit j the counting function
    f_j(v) = #{ P in {X,Y,Z} : syndrome of P at j equals v }.  The number of weight-w
    zero-syndrome errors is  sum_{|S|=w} (XOR-convolution of f_j, j in S)(0), and the
    Walsh-Hadamard transform of f_j is
        F_j(chi) = (-1)^{chi.z_j} + (-1)^{chi.(x_j+z_j)} + (-1)^{chi.x_j}
                 = 3  if g_chi is trivial on qubit j,   -1 otherwise.
    Consequently, with e_w the elementary symmetric polynomial,
        #{ e : wt(e)=w, s(e)=0 } = (1/16) sum_chi e_w(F_1(chi),...,F_14(chi))
                                 = (1/16) [ C(14,w)3^w + sum_{chi!=0} Ew[m(chi)][w] ],
        Ew[m][w] = sum_k C(m,k)(-1)^k C(14-m, w-k) 3^{w-k}.
    Writing P(m) = sum_{w=1..4} Ew[m][w], the total number of undetected errors (identity
    included) is
        U = (1/16) ( 91786 + sum_{chi != 0} P(m(chi)) )
    and N_det = 91770 - (U - 1).  So the objective depends ONLY on the weight enumerator of
    the 16-element stabilizer group.  Scoring costs 15 table lookups (~70 ns) instead of
    91770 symplectic products (~950 us): a measured 13000x speed-up.

    P(m), m = 0..14:
      91770 57914 34154 18314 8474 2970 394 -406 -326 -6 170 74 -166 -166 714
    (minimised by m = 7, then m = 8; note it is NOT monotone -- weight 10 and 14 elements
    are much worse than weight 12 or 13 ones.)

5.  CONSTRAINTS IN COLUMN LANGUAGE.
      weight:       c(e_i) = 14 - w  for the four points e_1..e_4  (4 linear equations);
      commutation:  XOR_j plucker(W_j) = 0 in F_2^6, where plucker(W) = x z^T + z x^T for
                    any basis (x,z) of W (basis-independent by step 2) and 0 if dim W < 2;
      rank 4:       c(chi) <= 13 for every chi != 0  (c(chi) = 14 means g_chi = I).

6.  SYMMETRIES REMOVED.
      qubit permutations S_14  -> the state is a MULTISET of column types (columns are added
                                  in non-decreasing type index, so each multiset appears once);
      local Cliffords (S_3)^14 -> columns are subspaces, not (x,z) pairs;
      generator permutations S_4 -> lex-minimality test on the sorted prefix.  Valid because
                                  for any permutation pi and prefix p of the sorted multiset S,
                                  the k smallest entries of sorted(pi(S)) are <= sorted(pi(p))
                                  entrywise; so sorted(pi(p)) < p implies pi(S) < S and S is
                                  not the canonical representative of its orbit.  Measured
                                  effect: 662 M nodes -> 45 M nodes (14.6x).
    Not used: full GL(4,2) basis changes, because fixing the four generators to be e_1..e_4
    is exactly the weight constraint; GL(4,2) would only be legitimate if the constraint were
    "some basis has weight 8", which is a different (weaker) problem.

7.  BRANCH AND BOUND.  At a node with degree vector c, r columns still to place and all of
    them of type index >= t0 (types are sorted by dim W ascending, so degU = |U|-1 is
    non-increasing in the type index):
      (a) feasibility: need_i = TARGET_C - c(e_i) must satisfy 0 <= need_i <= r, and some
          remaining type must still touch e_i;  sum_i need_i <= r * max popcount(sig);
      (b) rank: c(chi) <= 13;
      (c) bound 1: sum over the 11 non-generator points of  min_{0<=d<=r} P(14-c-d), plus the
          exact contribution 4*P(w) of the four generator points;
      (d) bound 2 (budget/knapsack): every remaining column increases sum_chi c(chi) by
          exactly degU in {3,7,15}, so the budget still available to the 11 non-generator
          points is at most degMax(t0)*r - sum_i need_i.  With
          ratio(c,r) = max_{1<=d<=min(r,13-c)} (P(14-c) - P(14-c-d))/d  (rounded up),
          the total achievable improvement is at most budget * max_chi ratio(c(chi),r).
      Both bounds are admissible (they relax an exact minimisation), so pruning on
      bound >= incumbent never removes a strictly better solution: the completed search is
      an optimality proof.

8.  COMPLEXITY.  Node cost is O(15) integer ops for the bound plus O(degU) for the
    incremental degree update; the S_4 test costs 24 * O(k log k) and is applied only for
    prefix length <= --sym-depth.  The raw multiset space is C(64,14) ~ 4.4e14; the search
    actually visits ~4.5e7 nodes for w=8 (0.5 s on 12 threads).

9.  RESULT (w = 8).  N_det = 86340 / 91770 = 94.08303367 %, undetected 5430, and the optimum
    is attained by a stabilizer group in which ALL 15 non-identity elements have weight
    exactly 8; the syndrome distribution is perfectly flat (5756 errors on each of the 15
    non-zero syndromes).  For w = 4 the optimum is N_det = 83968.
