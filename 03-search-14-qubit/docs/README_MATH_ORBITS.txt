=========================================================================================
stab14_2stage (orbit edition) -- what this experiment does, and what it does NOT do
=========================================================================================

THIS PROGRAM DOES NOT OPTIMISE THE 11-GENERATOR CODE DIRECTLY.

It is a strictly two-level experiment:

  level 1   find the GLOBALLY OPTIMAL Stage-1 family
                { S3 : three weight-8 commuting independent generators, M_3(S3) = 80584 }
            and enumerate that family completely, up to the symmetries that provably
            preserve the whole problem;

  level 2   for EACH member of that family, freeze it and solve the Stage-2 problem --
            eight weight-10 generators commuting with S3 and with each other, rank 11
            overall, maximising coverage of ONLY the errors U_3 that S3 missed.

The central quantity is therefore

        max over { S3 : M_3(S3) = 80584 }   of   max over H8   M_8(S3, H8),

and the scientific question is whether the residual errors left by the previous single-
solution run were an artefact of that one Stage-1 choice or a structural obstruction shared
by the entire optimal family.  A code that maximised detection over all 11 generators
jointly could well beat every number here; that is a different problem and is not asked.

-----------------------------------------------------------------------------------------
1.  WHY THE STAGE-1 ENUMERATION IS COMPLETE AND THE QUOTIENT IS LEGITIMATE
-----------------------------------------------------------------------------------------
Stage 1 is searched in the column-type representation: qubit j contributes the SUBSPACE
W_j = span{x_j, z_j} <= F_2^3, of which there are only 1 + 7 + 7 = 15.  Three group actions
are quotiented out, and each is proved to preserve generator weight, mutual commutation,
rank, AND the complete weight-<=4 detection objective:

  * S_14 qubit permutations -- the search state is a MULTISET of column types, and columns
    are added in non-decreasing type index, so each multiset is generated exactly once;
  * (S_3)^14 local Cliffords -- a single-qubit Clifford acts as GL(2,2) on (x_j,z_j); it
    preserves the row weights, preserves the commutation form exactly, since
        (ax+bz)(cx+dz)^T + (cx+dz)(ax+bz)^T = (ad+bc)(x z^T + z x^T) = x z^T + z x^T,
    and permutes the three single-qubit syndromes {z, x, x+z} setwise -- and the objective
    sums over all three Paulis on every support qubit, so it is invariant.  Hence only the
    subspace W_j matters, not the pair (x_j,z_j);
  * S_3 generator permutations -- a lex-minimality test on the sorted prefix.  Valid
    because for any permutation pi and any prefix p of the sorted multiset S, the k smallest
    entries of sorted(pi(S)) are <= sorted(pi(p)) entrywise, so sorted(pi(p)) < p implies
    pi(S) < S and S is not the canonical representative of its orbit.

Completeness of the optimal family: the DFS prunes a subtree only when its ADMISSIBLE
lower bound on the final objective is strictly worse than the current incumbent.  A pruned
subtree therefore cannot contain any optimal leaf, whatever the incumbent was at the time.
Collecting every leaf whose score equals the final optimum (and restarting the collection
whenever the incumbent improves) yields the entire optimal family in one pass.

The program reports two counts:
  "raw optimal Stage-1 solutions"          = distinct optimal column multisets, i.e. already
                                             modulo qubit permutations and local Cliffords;
  "inequivalent optimal Stage-1 solutions" = the same modulo S_3 as well.
(The number of literal 3 x 28 binary matrices is astronomically larger and is not useful.)

Two orbit representatives may in principle still be related by a symmetry OUTSIDE this
group -- e.g. a global (non-local) Clifford -- but such a map need not preserve Pauli
weight, so it is not usable here.  Running Stage 2 separately for each orbit is safe; the
only cost of a missed identification is duplicated work, never a missed solution.

-----------------------------------------------------------------------------------------
2.  WHY EACH ORBIT NEEDS ITS OWN STAGE-2 MACHINERY
-----------------------------------------------------------------------------------------
The Stage-2 quotient space V = C(S3)/S3 depends on S3, so class indices from one orbit are
meaningless in another.  Every orbit therefore rebuilds, from scratch:
  * the symplectic basis of C(S3)/S3 (dim 25 - 3 = 22), with a run-time check that the
    radical of the form on C(S3) really is S3;
  * the class map and the cov1 table for all 2^22 classes (one Walsh-Hadamard transform);
  * the weight-10 class representatives (Gray-code sweep of all 2^25 elements of C(S3));
  * the candidate list and the |U_3|-bit verification bitsets.
Only genuinely objective-independent tables are shared: binomials, the Krawtchouk/P table,
the error enumerator and the popcount helpers.

-----------------------------------------------------------------------------------------
3.  THE STAGE-2 COVERAGE IDENTITY (verified, not assumed)
-----------------------------------------------------------------------------------------
For an 8-dimensional totally isotropic D <= V, with cov1(v) = #{E in U_3 : <v,E> = 1},

        M_8(D) = (1/2^{d-1}) SUM_{v in D\{0}} cov1(v)  with d = 8
               = (1/128) SUM over the 255 non-zero elements of D.

Derivation: #uncovered = sum_{p in D^perp} f(p) where f is the multiplicity function of the
target classes; D^perp is the dot-product dual of sigma(D) (sigma swaps the two 11-bit
halves), so Poisson summation gives 2^{-d} sum_{v in sigma(D)} fhat(v), and
fhat(sigma(w)) = |U_3| - 2 cov1(w).  The constant is 1/2^{d-1}: it is 1/128 for d = 8, and
was 1/64 in the earlier d = 7 problem -- it is re-derived, not copied.  The program checks
it against the explicit |U_3|-bit bitset union in the self-test AND again on every reported
solution, where the bitset popcount, a from-scratch brute-force recount over all 91770
errors, and the algebraic value must all agree.

-----------------------------------------------------------------------------------------
4.  THE THREE SUCCESS CONDITIONS, TESTED SEPARATELY
-----------------------------------------------------------------------------------------
  A  full target coverage        M_8 = |U_3| = 11186  (=> M_final = 91770)
  B  all weight-4 errors detected  -- tested independently, by re-scanning every weight-4
     error against all 11 generators; the weight-4 total is COMPUTED (C(14,4)*3^4 = 81081),
     never hard-coded
  C  full weight-<=4 protection   M_final = 91770
A and C coincide here because Stage 2 is measured against exactly the Stage-1 misses, so
M_final = M_3 + M_8 -- an identity the program verifies independently.  B is strictly
weaker and can hold while A and C fail.

-----------------------------------------------------------------------------------------
5.  RESIDUAL-STRUCTURE ANALYSIS
-----------------------------------------------------------------------------------------
For every orbit the still-undetected set R is printed as Pauli strings AND in compact
support notation (e.g. XIIIIIIIIIIIII = X_0), and analysed for:
  * size and weight distribution;
  * whether R together with the identity is closed under multiplication, i.e. a subgroup
    (it always is when the closure test passes, since N(L) is a group and the closure only
    fails when a product leaves weight <= 4);
  * the dimension of span(R) and an explicit basis;
  * whether every residual error is X-only;
  * which physical qubits are involved.
The structure < X_i, X_j, X_k > is exactly the case "three physical qubits whose column in
the 11-generator matrix spans only one dimension": the Krawtchouk identity at w = 1 gives
B_1 = 3*n0 + n1 = that count, and N(L) being a group then forces the 2^3 - 1 = 7 products.

-----------------------------------------------------------------------------------------
6.  HONESTY ABOUT OPTIMALITY
-----------------------------------------------------------------------------------------
Stage 1: PROVEN globally optimal, and the optimal family is enumerated exhaustively.
Stage 2: heuristic per orbit (randomized greedy on the exact gain, steepest-descent 1-opt,
ruin-and-recreate LNS, optional beam search), with an exact branch-and-bound available.
Its admissible bound

    M_8 <= ( S_d + max_{v allowed} G(v) + (256 - 2^{d+1}) * COVMAX ) / 128

prunes only while COVMAX <= 128 * incumbent / 255; with COVMAX ~ 5796 and an incumbent near
11179 the threshold is ~5610, so the relaxation cannot prune near the root and exhausting
~2 x 10^6 candidate classes at depth 8 is out of reach.  Stage-2 results are therefore
BEST KNOWN unless a row of the comparison table says PROVEN.  A "no full coverage" verdict
is a statement about what was found, NOT a proof of impossibility.

-----------------------------------------------------------------------------------------
7.  RESULTS OF THE EXPERIMENT
-----------------------------------------------------------------------------------------
Stage-1 optimal family, enumerated exhaustively:
    148  optimal column multisets (already modulo qubit permutations and local Cliffords)
     45  inequivalent orbits (additionally modulo the S_3 generator permutations)
    every one of them has M_3 = 80584 / 91770 and |U_3| = 11186, as it must.

Stage 2 was then run independently on all 45 orbits (25 s of heuristic each, 19.4 min in
total, one fresh 22-dimensional quotient / cov1 table / candidate set per orbit):

    M_8 ranges from 11175 to 11179 across the family.

    best  M_8 = 11179  (7 residual errors) on orbits 0002, 0006, 0012, 0023, 0031, 0036, 0040
    worst M_8 = 11175  (11 residual errors) on orbits 0004, 0007, 0008, 0039

So the Stage-1 choice DOES matter, but only within a band of 5 errors, and no member of the
optimal family reached full coverage.

WEIGHT-4 ERRORS.  The choice matters much more for the weight-4 sub-question.  The number of
undetected weight-4 errors ranges from 9 down to 0, and five orbits -- 0002, 0012, 0023,
0031, 0036 -- detect ALL 81081 weight-4 errors (the total is computed, never hard-coded).
For those orbits the residual is a 3-dimensional subgroup lying entirely in weights 1..3:
    orbits 0002, 0023            residual weights 3/3/1/0, X-only, e.g. < X_1, X_2, X_3 >
    orbit  0012                  residual weights 3/3/1/0, not X-only
    orbits 0031, 0036            residual weights 1/3/3/0, not X-only
This is Condition B holding while Conditions A and C fail.

VERDICT.
    Condition A (M_8 = 11186)   NOT achieved on any of the 45 orbits
    Condition B (all weight-4)  ACHIEVED, on 5 of the 45 orbits
    Condition C (M_final=91770) NOT achieved; best M_final = 91763 / 91770 = 99.99237223 %

INTERPRETATION.  The previously observed seven-error residual is BOTH: its size (7) is
shared by the best members of the whole optimal family and never improved on, which looks
structural, while its position and weight profile are an artefact of the Stage-1 choice --
moving to orbit 0002 keeps 7 residual errors but moves all of them out of weight 4.
Because a residual of 0 would require the 11-dimensional stabilizer to have no non-zero
Pauli of weight <= 4 in its normaliser -- i.e. a PURE [[14,3,5]] code containing the frozen
S3 -- and because every Stage-2 search here is heuristic, this is evidence for a structural
obstruction but NOT a proof of one.

DEEP SECOND PASS.  The seven orbits attaining M_8 = 11179 were re-run with six times the
heuristic budget (150 s) and a different seed: all seven reproduced exactly 11179, none
improved, and orbit 0040 additionally moved its residual out of weight 4 (so six of the
forty-five orbits now achieve Condition B).  See deep_pass_refinement.txt.
