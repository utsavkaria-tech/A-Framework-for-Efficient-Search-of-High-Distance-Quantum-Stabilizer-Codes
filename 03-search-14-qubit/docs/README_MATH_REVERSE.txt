=========================================================================================
REVERSE TWO-STAGE OPTIMISATION OF A 14-QUBIT STABILIZER CODE
        STAGE 1 :  8 generators of Pauli weight 10
        STAGE 2 :  3 generators of Pauli weight 8
=========================================================================================
SELF-CONTAINED.  Everything produced by this experiment lives in reverse_8w10_then_3w8/.
Nothing in stab14*.cpp, stage1_orbit_*/, targeted_3opt4opt/, targeted_6opt/,
overnight_orbit0006/ or overnight_refine/ is read for input or written to.  The earlier
experiments did the two stages the other way round (3 x weight-8 first); this one reverses
the order to see whether the residual structure is an artefact of that ordering.

-----------------------------------------------------------------------------------------
0.  THE HEADLINE, BEFORE ANY RESULT
-----------------------------------------------------------------------------------------
Stage 1 CANNOT be exhausted, and Stage 2 CAN.  Concretely, all figures computed by the
program itself (--spacesize):

    8-dimensional isotropic subspaces of F_2^28          1.878e+40
    ... divided by the full 14! * 6^14 = 6.83e+21 group  2.749e+18
    reduced column types (subspaces of F_2^8, dim <= 2)      11051
    multisets of 14 such types                           4.686e+45

    for comparison, the ranks that WERE proven exhaustively in this project:
      rank 3:     15 types,   4.012e+07 multisets   (branch and bound: 0.4 s)
      rank 4:     51 types,   4.786e+13 multisets   (branch and bound: 224 s)
      rank 5:    187 types,   1.180e+21 multisets   (does not finish)

So every Stage-1 number this program prints is BEST KNOWN, never OPTIMUM.  A rigorous
lower bound on the residual |U_8| is computed separately (section 7) so the distance
between "best known" and "provable" is always explicit.  Stage 2, by contrast, is swept
completely for each frozen Stage-1 code and is a PROVEN OPTIMUM conditional on it.

-----------------------------------------------------------------------------------------
1.  REPRESENTATION AND THE EXACT SEARCH SPACE
-----------------------------------------------------------------------------------------
A Pauli modulo phase is (x|z) in F_2^14 x F_2^14, packed as one uint32  x | (z << 16).
Pauli weight is wt = popcount(x | z), so wt(X) = wt(Y) = wt(Z) = 1; this is NOT the binary
Hamming weight of the 28-bit vector.  Commutation is the symplectic form
<a,b> = x_a.z_b + z_a.x_b, and the syndrome bit of error e under generator g is exactly
<g,e>.  An error is DETECTED iff some generator anticommutes with it.

The objective set is every non-identity Pauli of weight 1..4.  Its size is computed, never
hard-coded: SUM_{w=1..4} C(14,w) 3^w = 91770, of which 81081 have weight 4.

Stage 1 searches the 8-dimensional totally isotropic subspaces L <= F_2^28 that admit a
basis of eight weight-10 elements.  Stage 2, with L frozen, searches the 3-dimensional
isotropic subspaces of the quotient L^perp / L that admit weight-8 lifts.

-----------------------------------------------------------------------------------------
2.  THE OBJECTIVE IN CLOSED FORM  (the engine of the whole program)
-----------------------------------------------------------------------------------------
By Fourier inversion over the group L (dim r),

    #{e : wt(e) = w, s(e) = 0}  =  (1/2^r) SUM_{h in L} SUM_{wt(e)=w} (-1)^{<h,e>} .

The inner sum factorises over qubits.  Where h is trivial, all three non-identity Paulis
commute and contribute +1 each, a factor 3; where h is non-trivial exactly one of the three
commutes, giving 1 - 2 = -1.  With m = wt(h) the generating function is
(1 - y)^m (1 + 3y)^{14-m}, so

    EW[m][w]  =  SUM_k C(m,k) (-1)^k C(14-m, w-k) 3^{w-k} ,
    P[m]      =  SUM_{w=1..4} EW[m][w] ,
    U(L)      =  (1/2^r) SUM_{h in L} P(wt(h)) ,      M = 91770 - U(L) .

P[m] for m = 0..14, recomputed at run time and printed by --selftest:

    91770  57914  34154  18314  8474  2970  394  -406  -326  -6  170  74  -166  -166  714

Scoring an 8-generator code is therefore 256 table lookups instead of 8 x 91770 symplectic
products.  Two properties matter for everything below:

  (a) it is the EXACT objective, not a proxy -- verified against direct enumeration on
      hundreds of random subspaces at every rank (--selftest);
  (b) it depends only on the SUBSPACE L, not on the chosen basis.  That is what makes the
      k-opt reformulation of section 5 legitimate.

-----------------------------------------------------------------------------------------
3.  COLUMN TYPES: THE QUBIT-STRUCTURED FORMULATION AND ITS CANONICALISATION
-----------------------------------------------------------------------------------------
Write H = [X|Z] with X,Z in F_2^{8x14}.  Qubit j contributes (x_j, z_j) in F_2^8 x F_2^8;
set W_j = span{x_j, z_j} <= F_2^8, so dim W_j <= 2.  For h = SUM_i a_i g_i the Pauli of h on
qubit j is (a.x_j, a.z_j), which is the identity exactly when a is orthogonal to W_j.  Hence

        wt(h_a)  =  14 - #{ j : a _|_ W_j } ,

so the whole weight function -- and therefore the whole objective -- is determined by the
MULTISET { W_1, ..., W_14 }.  This is the canonical form, and it absorbs two symmetries at
once:

  LOCAL CLIFFORD.  Per qubit, GL(2,2) = S_3 acts as (x,z) -> (ax+bz, cx+dz) with ad+bc = 1.
  It fixes W_j, hence preserves every element weight and the objective; and since the
  qubit-j contribution to <g_i,g_k> is multiplied by det = 1 it preserves all commutators.
  A column is thus completely described by its subspace: 65536 raw types collapse to
  1 + 255 + 10795 = 11051.

  QUBIT PERMUTATION.  S_14 permutes the multiset, so only the multiset matters.  The branch
  and bound enforces this by adding columns in non-decreasing type index.

  GENERATOR RELABELLING.  S_8 acts on F_2^8 and hence on the types.  The k-opt engine does
  not need it at all, because it searches SUBSPACES and never enumerates bases; the branch
  and bound uses it explicitly (section 6).

The three constraints all become statements about the multiset, and the program checks each
of them against the Pauli-level definition on every reported solution (columns.txt):

  WEIGHT       generator i has weight 10  <=>  exactly 4 columns satisfy W_j <= H_i,
               where H_i = {v : v_i = 0}.  These are 8 linear conditions on the multiset,
               propagated incrementally rather than tested at the end.
  COMMUTATION  SUM_j plucker(W_j) = 0, where plucker(W) = x z^T + z x^T is symmetric with
               zero diagonal, vanishes when dim W <= 1, and is basis-independent because a
               change of basis multiplies it by det = 1.  For r = 8 this is 28 bits.
  RANK         rank H = 8  <=>  W_1 + ... + W_14 = F_2^8  <=>  wt(h_a) > 0 for all a != 0.

-----------------------------------------------------------------------------------------
4.  WHY THE MULTISET FORMULATION DOES NOT SOLVE RANK 8
-----------------------------------------------------------------------------------------
At rank 3 and rank 4 this representation is decisive: 15 and 51 types, and the trees close.
At rank 8 there are 11051 types and 4.686e+45 multisets.  Pruning in the earlier
experiments bought roughly seven orders of magnitude; forty-five orders are needed.  The
same conclusion follows from the subspace side: 1.878e+40 isotropic 8-subspaces, and even
after dividing by the entire 6.83e+21 local-Clifford-and-permutation group, 2.75e+18
remain.  (That quotient is an accounting of scale, not an exact orbit count: it ignores
that only some isotropic subspaces admit a weight-10 basis, and orbits need not be free.)

Exhaustive Stage 1 is therefore out of reach, and the program says so rather than running
something it calls exhaustive.  --mode exact refuses to start without --force-exact.

-----------------------------------------------------------------------------------------
5.  WHAT THE STAGE-1 SEARCH ACTUALLY DOES:  k-OPT AS A SUBSPACE SEARCH
-----------------------------------------------------------------------------------------
Retain t = 8-k generators, K = their span (dim t, isotropic).  Every valid completion is an
8-dimensional isotropic D' with K <= D' <= K^perp, i.e. exactly a k-dimensional totally
isotropic subspace of

        Q = K^perp / K ,      dim Q = 28 - 2t = 12 + 2k .

With phi(q) = SUM_{u in K} P(wt(lift(q) + u)) the objective telescopes exactly:

        S(D')  =  SUM_{v in span(q_1..q_k)} phi(v) ,        U = S / 256 .

    k = 1 : dim Q = 14,  16384 classes   -- swept over all 8 drops
    k = 2 : dim Q = 16,  65536 classes   -- swept over all 28 retained sextuples
    k = 3 : dim Q = 18, 262144 classes   -- swept over all 56 retained quintuples

ADMISSIBILITY IS APPLIED TO THE SUBSPACE, NOT TO THE REPLACEMENTS.  D' is admissible iff
its weight-10 elements span it (rank 8), which is precisely "D' has a basis of eight
weight-10 generators".  This is a strictly larger and more correct neighbourhood than the
usual shortcut of demanding that each replaced class individually carry a weight-10 lift.

PRUNING (exact).  Sort the non-zero classes by phi ascending and let pre[k] be the sum of
the k smallest.  Any k DISTINCT non-zero classes have phi-sum at least pre[k] -- true
whatever has already been chosen, because excluding classes can only raise the sum.  At
each level the still-unknown terms are a known number of distinct classes, giving for 2-opt

    c1 : phi(c1) + pre[2] >= need        c2 : phi(c2) + pre[1] >= need - phi(c1)

Because the list is sorted ascending these are hard breaks, and nothing that could beat the
incumbent is discarded.  This replaces the cruder k*phimin bound and is what makes an
exhaustive 2-opt cost a fraction of a second instead of sweeping 2.1e9 pairs (measured:
about 1e5 pairs survive).

3-OPT NEEDS ONE MORE IDEA.  At k=3 the quotient has 262144 classes and a plain triple loop
is hopeless.  Telescoping one level fixes it: with psi(w) = phi(w) + phi(w ^ c1) the seven
non-zero classes of span(c1,c2,c3) regroup as

    S = sK + phi(c1) + psi(c2) + psi(c3) + psi(c2 ^ c3) ,

so once c1 is fixed the rest is exactly the cheap three-term pair problem that 2-opt
already solves.  Let pre2[k] be the sum of the k smallest psi over the classes other than 0
and c1.  If a triple beats the incumbent then psi(c2)+psi(c3)+psi(c2^c3) < need2, and since
any TWO of those three distinct classes sum to at least pre2[2], each one individually
satisfies psi(x) < need2 - pre2[2].  So c2, c3 AND c2^c3 all lie in the single candidate
list L = { x : psi(x) + pre2[2] < need2 }, built once per c1 and normally tiny; sorting L by
psi then gives exact breaks at both remaining levels.  Nothing that could beat the incumbent
is discarded.

3-opt is parallelised over the 56 drop-triples.  Measured on this machine: a complete 3-opt
of one solution is about 1h 05m to 1h 20m on 12 threads (two measurements, extrapolated from
the fraction of c1 prefixes retired), against 0.5 s for a complete 2-opt.  It is
therefore OFF by default (--threeopt turns it on) and is meant as a deep certification of a
champion rather than something to run on every family member.  It earns its keep: in the
calibration run a 120-second slice of it alone improved |U_8| from 249 to 228.

CONSTRUCTION.  Each restart builds a code greedily: at every step it samples candidate
weight-10 Paulis from the centralizer of what is already chosen (a random element of that
subspace has weight 10 about 22 % of the time, so rejection sampling is cheap) and keeps
the one minimising the partial span sum.  Then 1-opt runs to a local optimum.

-----------------------------------------------------------------------------------------
6.  THE BRANCH AND BOUND, AND ITS BOUND
-----------------------------------------------------------------------------------------
A state is a non-decreasing sequence of column-type indices, i.e. a multiset, so S_14 is
already quotiented out.  Incremental state is c(chi) = number of chosen columns on which
chi acts, for each chi in F_2^r.

BOUND (admissible).  chi = 0 contributes P(0), and each of the r generator points is pinned
at exactly weight w.  Every remaining column raises SUM_{chi != 0} c(chi) by
delta(W) = 2^r - 2^{r-dim W}, at most deltaMax = 3*2^{r-2}; the generator points already
claim needGen = SUM_i (w - c(e_i)) of that, and the rest is a BUDGET the other characters
share.  A character at count c cannot fall below PMIN[c][rem] and cannot fall faster than
RATIO[c][rem] per unit of budget, so the total decrease is bounded by both sums and the
smaller is used.  Additionally, characters still at zero each need at least one increment
(the rank constraint), which prunes when their number exceeds the budget.  Every ingredient
relaxes an exact minimisation, so pruning on "bound >= incumbent" cannot discard a strictly
better solution and a completed tree is a proof.

S_r SYMMETRY.  Relabelling generators permutes F_2^r and hence the types.  A node is pruned
when some permutation makes the sorted image of the partial multiset lexicographically
smaller.  Validity: for a canonical complete solution S, the k smallest entries of
sorted(pi(S)) are entrywise <= sorted(pi(p)) for any k-prefix p, so sorted(pi(p)) < p would
force sorted(pi(S)) < S and contradict canonicity.  The test is applied for the first
--sym-depth columns only, which is still valid because it is a subset of the legal prunes.
Measured at rank 4, weight 8: depth 3 does not close in 600 s; depth 6 closes in 224 s;
depth 14 cuts the tree 6.5x further but costs more per node than it saves.  Default 6.

INCUMBENT.  The tree is seeded with a greedy solution before it starts.  Any ACHIEVABLE
value is a legitimate upper bound, so seeding cannot change which optimum a completed tree
reports -- only how much of the tree is walked.  Because the default prune is
"bound >= incumbent", the seed is stored one unit loose so the leaf achieving the seed
value is not itself pruned.  --bb-all switches to "bound > incumbent", which is slower but
enumerates every optimal leaf.

VALIDATION.  Run by --selftest, and the reason the engine can be trusted:
    rank 3, weight 8  ->  M = 80584   matches the optimum proven earlier in this project
    rank 4, weight 8  ->  M = 86340   matches the optimum proven earlier in this project
Both trees are completed, not cut short.

-----------------------------------------------------------------------------------------
7.  A CERTIFIED LOWER BOUND ON |U_8|   (the only rigorous Stage-1 statement available)
-----------------------------------------------------------------------------------------
L is an additive self-orthogonal code over the four-element alphabet {I,X,Y,Z}, length 14,
256 codewords.  Let A_m be its weight distribution and B_w that of its dual (2^20 words).
With K_w(m) = SUM_k C(m,k)(-1)^k C(14-m,w-k) 3^{w-k},

    B_w = (1/256) SUM_m A_m K_w(m),   B_w >= A_w  (because L <= L^perp),
    A_m >= 0,   A_0 = 1,   SUM_m A_m = 256,   A_10 >= 8  (eight weight-10 generators),
    |U_8| = SUM_{w=1..4} B_w = (1/256) SUM_m A_m P(m).

For ANY multipliers y_w, z_w, s >= 0 define
    d_m = P(m)/256 - (1/256) SUM_w (y_w + z_w) K_w(m) + z_m - s*[m = 10] .
Then, since SUM_{m>=1} A_m = 255 with A_m >= 0,

        |U_8|  >=  8 s  +  d_0  +  255 * min_{m>=1} d_m .

This is weak duality: EVERY choice of multipliers yields a valid bound, so the hill-climb
that searches for good multipliers cannot produce a wrong answer, only a weak one.  No LP
solver is trusted anywhere in the argument, and the certificate is re-evaluated
independently of the search loop before being printed.

The bound this currently finds is |U_8| >= 36.  Against a best-known |U_8| in the low 200s
from short runs, that is a wide gap -- it says the LP relaxation is weak here, not that
|U_8| = 36 is reachable.  It is nevertheless the only statement in this experiment that
holds over ALL 8-generator weight-10 codes rather than over the ones that were visited.

-----------------------------------------------------------------------------------------
8.  STAGE 2 -- EXHAUSTIVE, HENCE A PROVEN OPTIMUM ONCE STAGE 1 IS FROZEN
-----------------------------------------------------------------------------------------
A Stage-2 generator must lie in L8^perp, have weight 8, and be independent of L8; the three
must commute.  The form descends to Q2 = L8^perp / L8, of dimension 28 - 16 = 12, so only
4096 classes exist and the admissible additions are exactly the 3-dimensional totally
isotropic subspaces of Q2 spanned by classes carrying a weight-8 lift.  That is small
enough to enumerate completely.  (Measured: about half the classes are weight-8-usable --
2079 of 4095 in the calibration run.)

COVERAGE IDENTITY.  For e in U_8 the bit <h,e> depends only on the class of h, so with
cov1(q) = #{e in U_8 : <q,e> = 1},

    #{e : <q,e> = 0 for all q in D} = (1/2^d) SUM_{q in D} (|U_8| - 2 cov1(q))
    =>  #covered by a d-dimensional D  =  (1 / 2^{d-1}) SUM_{q in D} cov1(q) ,

a quarter of the sum over the seven non-zero classes when d = 3.  cov1 is obtained by one
Walsh-Hadamard transform over the 4096 classes.  The identity is checked against brute
force on random 3-spaces in --selftest.

OBJECTIVE.  M_3(2) = #{e in U_8 detected by h1, h2 or h3}, and M_final = M_8 + M_3(2).
PRUNE (admissible): with base = cov1(a) + cov1(b) + cov1(a^b) the four remaining terms are
each at most COVMAX, so base + 4*COVMAX below the incumbent kills the pair outright.

STAGE-2 SYMMETRY.  Only the subspace matters, so the three generators are chosen as a
subspace of Q2 and never as an ordered triple, and one weight-8 representative per class is
fixed in advance.  No symmetry that moves the frozen Stage-1 code is used, so two
inequivalent Stage-1 residual structures are never accidentally identified.

-----------------------------------------------------------------------------------------
9.  THE OPTIMAL STAGE-1 FAMILY, AND THE DIRECTION OF ITS ERROR
-----------------------------------------------------------------------------------------
Every solution attaining the best known |U_8| is kept.  Distinct subspaces are recognised
by the reduced row echelon form of the 8x28 matrix, which is a true canonical form for the
subspace.  Solutions are then grouped by a fingerprint invariant under qubit permutation
and local Clifford:

    (a) the weight enumerator A[0..14] of the 256-element span;
    (b) for each qubit j the vector n_j[m] = #{h : wt(h) = m, h acts on j}, 14 vectors
        sorted (this kills the qubit permutation);
    (c) for each qubit pair the analogous vector, 91 vectors sorted.

Both operations preserve every element's weight and exact support, so equal-class implies
equal fingerprint.  THE CONVERSE DOES NOT HOLD.  Different fingerprints PROVE inequivalence;
equal fingerprints only suggest equivalence.  The reported class count is therefore a LOWER
bound on the number of inequivalent solutions found.  To guard against over-merging, several
distinct members of each fingerprint class are carried independently through Stage 2: if two
members produce different Stage-2 optima or different residual structure they are provably
inequivalent, and the comparison table shows it.

-----------------------------------------------------------------------------------------
10.  COMPUTATIONAL COMPLEXITY, MEASURED
-----------------------------------------------------------------------------------------
See CALIBRATION.txt for the machine-generated numbers; the shape is:

    scoring one 8-generator code        256 table lookups
    one restart (build + 1-opt descent) ~0.15 s single-threaded, reaching |U_8| ~ 265 on
                                        average and ~245 at best from a single descent
    one exhaustive 2-opt (28 sextuples) ~0.4 s single-threaded (about 1e5 pairs survive
                                        the filter out of 2.1e9)
    one exhaustive 3-opt (56 triples)   ~1h 05m-1h 20m on 12 threads -- opt-in, section 5
    one exhaustive Stage-2              ~0.4-1.8 s on 12 threads
    branch and bound, rank 3 weight 8   0.4 s, tree closed
    branch and bound, rank 4 weight 8   ~540 s, tree closed (sym-depth 6)
    branch and bound, rank 5+           does not close

    multi-start budget on 12 threads:   1e3 restarts 10 s, 1e4 ~1m45s, 1e5 ~17m, 1e6 ~2h54m

Phase 2 (multi-start) dominates and is the only phase with a tunable budget.  Phase 3
(2-opt) and phase 4 (Stage 2) are negligible next to it; phase 3b (3-opt) is not, which is
why it is opt-in.

-----------------------------------------------------------------------------------------
11.  INDEPENDENT VERIFICATION
-----------------------------------------------------------------------------------------
No reported solution is accepted on the strength of the algebraic score.  For every one,
all 91770 errors are regenerated from scratch and their syndromes recomputed by direct
symplectic products, then:

    1. eight Stage-1 generators of weight exactly 10
    2. three Stage-2 generators of weight exactly 8
    3. all 55 generator pairs commute
    4. rank of the 11-generator stabilizer is 11
    5. |U_8| rebuilt by brute force equals the algebraic score
    6. the closed-form identity SUM P(wt)/256 == |U_8| holds
    7. Stage-2 coverage recomputed independently equals the searched value
    8. the residual is listed explicitly, with its weight distribution, span dimension,
       subgroup test, X-only test and the qubits involved

--selftest additionally checks the error count against the formula, the P table against
brute force on hundreds of random subspaces at every rank, the column-multiset weight
function and objective against the Pauli-level ones, the quotient construction
(class(lift(c)) == c, the descended form, the phi table against a direct sum), the Stage-2
coverage identity against brute force, and the branch and bound against the two optima
proven earlier in this project.

-----------------------------------------------------------------------------------------
12.  WHAT MAY AND MAY NOT BE CLAIMED
-----------------------------------------------------------------------------------------
    PROVEN GLOBAL OPTIMUM   the branch and bound at rank 3 and rank 4 (trees closed)
    PROVEN LOWER BOUND      the |U_8| certificate of section 7
    PROVEN OPTIMUM,         each Stage-2 result, conditional on its frozen Stage-1 code
      CONDITIONAL
    BEST KNOWN              every Stage-1 result at rank 8, and hence every M_final
    HEURISTIC               nothing is reported under this label; the construction is
                            heuristic but only ever feeds the exact machinery above

If no complete [[14,3,5]] code is found, that is NOT an impossibility result.  Stage 2 is
exhaustive only for the Stage-1 codes actually visited, and the Stage-1 search covers a
vanishing fraction of a 1.878e+40-element space.  The correct statement is:

    "Within the best-known Stage-1 family found here, no Stage-2 completion detects every
     weight-<=4 error; the Stage-2 search was exhaustive for each frozen Stage-1 code."

The scientific question this experiment exists to answer is whether reversing the order
changes the residual structure, i.e. whether optimising the eight weight-10 generators
first leads anywhere different from optimising the three weight-8 generators first.  That
comparison lives in stage2_comparison.txt (Condition C) and is meaningful even though
Stage 1 is not exhaustive, because the whole best-known family is compared, not one
representative.

-----------------------------------------------------------------------------------------
13.  HOW TO RUN IT
-----------------------------------------------------------------------------------------
    cd "C:\Users\utsav\OneDrive\ドキュメント\Sidon Set Search\stabilizer_14q\reverse_8w10_then_3w8"
    .\reverse_8w10_3w8.exe --restarts 200000 --threads 12

    --time-limit S instead of --restarts to run for a fixed wall time
    --resume            continue from checkpoint.txt after any interruption
    --selftest          the full correctness suite (~10 minutes; it closes the rank-4 tree)
    --threeopt          add the exhaustive 3-opt phase (~1.1 h per member; see section 5)
    --calibrate         measure the engine on this machine and project runtimes
    --spacesize         the search-space accounting of section 0
    --bb R W            run only the branch and bound at rank R, generator weight W
    --threeopt          add the exhaustive 3-opt certification (expensive, see CALIBRATION)
    --members N         members retained per fingerprint class (default 4)
    --classes N         maximum fingerprint classes retained (default 512)

Output: reverse_search.log, checkpoint.txt (rewritten continuously), stage1_summary.txt,
stage1_optimal_family.txt, stage2_comparison.txt, final_summary.txt, and one directory
stage1_orbit_NNNN/ per class member holding its representative, its undetected set, its
Stage-2 optimum, the full 11-generator matrix, the remaining errors, the verification block
and the column-type report.
