=========================================================================================
STAGED SEARCH   3 x weight-8  ->  3 x weight-10  ->  3 x weight-10  ->  2 x weight-10
target: an 11-generator stabilizer on 14 qubits with |U_{<=4}| = 0, i.e. [[14,3,5]]
=========================================================================================
SELF-CONTAINED.  Everything lives in staged_3w8_3w10_3w10_2w10/.  No file, binary, log,
checkpoint or result of any earlier experiment in this project is read or written.

-----------------------------------------------------------------------------------------
0.  THE HEADLINE, BEFORE ANY DESIGN DETAIL
-----------------------------------------------------------------------------------------
Two findings, both measured, decide the shape of the whole thing:

  (A) NAIVE STAGED PREFIX-OPTIMISATION IS INVALID, and here it is not a theoretical
      worry: Stage 1 has 1057 inequivalent classes, of which only 10 attain the best
      |U_3|.  Keeping "the optimal prefix" would discard 1047 classes that no exact
      argument permits discarding.

  (B) EXHAUSTIVE STAGE 2 IS OUT OF REACH BY ~11 ORDERS OF MAGNITUDE.  Measured on this
      machine: the 3-dimensional isotropic subspaces of Q_3 number 5.49e16 per Stage-1
      class, the search runs at ~1.0e7 nodes/s, giving ~1.5e6 hours (about 170 years)
      PER CLASS, and there are 1057 classes.  The exact feasibility bound (F2) does not
      rescue this: it passes with only a 1.6-2.6 % margin, i.e. it is nearly vacuous.

So the honest deliverable is: Stage 1 solved exactly and completely, the exact conditions
under which the staged decomposition is valid, and measured proof of where it stops.

-----------------------------------------------------------------------------------------
1.  THE OBJECTIVE IN CLOSED FORM
-----------------------------------------------------------------------------------------
A Pauli mod phase is (x|z) in F_2^14 x F_2^14; weight is popcount(x|z), so
wt(X)=wt(Y)=wt(Z)=1.  Commutation is the symplectic form <a,b> = x_a.z_b + z_a.x_b and the
syndrome bit of e under g is <g,e>; an error is UNDETECTED iff it lies in L^perp, L being
the span of the generators.  Fourier inversion over L (dim r) and factorisation over qubits
(a trivial qubit contributes 3, a non-trivial one contributes 1-2 = -1) give

    EW[m][w] = SUM_k C(m,k)(-1)^k C(14-m,w-k) 3^{w-k},   P[m] = SUM_{w=1..4} EW[m][w],
    |U(L)|   = (1/2^r) SUM_{h in L} P(wt h).

P[0] = 91770 is the total number of weight-<=4 errors.  Scoring an r-generator prefix is
2^r table lookups.  P[] is computed at run time and checked against brute force on 400
random isotropic prefixes of every rank (0 mismatches).  P[m], m = 0..14:

    91770 57914 34154 18314 8474 2970 394 -406 -326 -6 170 74 -166 -166 714

-----------------------------------------------------------------------------------------
2.  IS STAGED PREFIX-OPTIMISATION VALID?   NO -- AND HERE IS THE EXACT REPLACEMENT
-----------------------------------------------------------------------------------------
Freeze a prefix L_r.  Every later generator commutes with L_r, so it lies in L_r^perp, and
shifting it by an element of L_r changes no syndrome (for g in L_r^perp and u in L_r,
<g,e+u> = <g,e>).  The entire remaining problem therefore lives in

    Q_r = L_r^perp / L_r ,     dim Q_r = 28 - 2r ,

a non-degenerate symplectic space.  U_r maps into Q_r and, for a completion spanning an
isotropic D <= Q_r of dimension k,

    e in U_r survives  <=>  [e] _|_ D ,
    detected(D) = (1 / 2^{k-1}) SUM_{d in D \ 0} cov(d),   cov(d) = #{e in U_r : <d,[e]>=1}.
                                                                                      (COV)

(COV) is proved by (1/2^k) SUM_{d in D} SUM_e (-1)^{<d,[e]>} and is checked against brute
force on random isotropic triples at run time.

WHY OPTIMALITY DOES NOT COMPOSE.  What a prefix can still become is governed by how the
classes [e] are DISTRIBUTED inside Q_r -- by the shape of cov -- and not by |U_r|.  Two
prefixes with |U_r| = 11186 and 11446 are not ordered by their futures; the worse one may
have a better-spread class distribution and admit strictly better completions.  There is no
exchange argument and hence no greedy or matroid structure, so "keep the optimal prefix" is
not a valid reduction.  This program never does it.

WHAT IS VALID.  The target is |U_11| = 0, not "minimise", so prefix OPTIMALITY is not what
we need -- prefix FEASIBILITY is, and that admits exact tests:

  (F1) PURITY OF EVERY PREFIX.  If h in L_r is non-identity with wt(h) <= 4 then for every
       completion L_r <= L_11 <= L_11^perp, so h stays undetected forever and |U_11| >= 1.
       Hence every prefix of a solution satisfies  min non-zero weight of L_r >= 5.
       Cheap and exact.  (Valid only for the pure target; see section 6.)

  (F2) COVERAGE BUDGET -- an admissible upper bound on ALL future improvement.  By (COV)
       with k = 11 - r generators still to come, reaching zero requires
             SUM_{d in D \ 0} cov(d) = 2^{k-1} |U_r| ,
       a sum of exactly 2^k - 1 terms.  So with S_top(n) = the sum of the n largest values
       of cov over the whole of Q_r,
             S_top(2^k - 1)  >=  2^{k-1} |U_r|                                       (F2)
       is NECESSARY.  It relaxes an exact maximisation -- it ignores that D must be a
       subspace, must be isotropic, and must be spanned by weight-10-liftable classes -- so
       it can never reject a prefix that could still reach zero.  Weaker corollary, useful
       as a quick reject: cov_max >= 2^{k-1}|U_r|/(2^k - 1) > |U_r|/2, i.e. SOME single
       class must detect more than half of everything still undetected.

  (F3) MONOTONICITY.  U_{r+1} is contained in U_r for every extension.  Asserted at run
       time (100 random chains, 0 violations) rather than assumed.

The staged structure is kept exactly as specified, 3 -> 6 -> 9 -> 11, but every stage
RETAINS EVERY PREFIX PASSING (F1) AND (F2), never only the optimal ones.  That is the
strongest form of the decomposition that provably preserves every globally optimal code.

-----------------------------------------------------------------------------------------
3.  STAGE 1 IS EXACTLY SOLVABLE -- THE COLUMN-MULTISET FORMULATION
-----------------------------------------------------------------------------------------
Read the check matrix column-wise: qubit j contributes x_j, z_j in F_2^r and we set
W_j = span{x_j,z_j}, of dimension at most 2.  The element h_a = SUM a_i g_i acts on qubit j
as (a.x_j, a.z_j), the identity exactly when a _|_ W_j.  So with c(a) = #{j : a _|_ W_j},

    wt(h_a) = 14 - c(a),

and every weight, the objective and purity are functions of the MULTISET {W_1..W_14} alone.
Two symmetries are absorbed exactly and for free:

  * LOCAL CLIFFORD: re-choosing the basis inside W_j is GL(2,2) = S_3 permuting X,Y,Z on
    that qubit; it fixes W_j, so it changes nothing.  Factor 6^14 = 7.8e10.
  * QUBIT PERMUTATION: only the multiset matters.  Factor up to 14! = 8.7e10.

For r = 3 the subspaces of F_2^3 of dim <= 2 number 1 + 7 + 7 = 15, so a Stage-1
configuration is a vector of 15 counts summing to 14: C(28,14) = 4.01e7 states, exhaustible
outright.  The conditions become:

    weight       c(a) = 6 for three independent a          (wt = 8)
    rank 3       c(a) <= 13 for every a != 0
    commuting    SUM_j pl(W_j) = 0 in Lambda^2(F_2^3) = F_2^3, pl(W) = x z^T + z x^T
                 (basis-independent: a change of basis multiplies it by det = 1; zero when
                 dim W <= 1)
    purity (F1)  c(a) <= 9 for every a != 0
    objective    |U_3| = (1/8) SUM_{a in F_2^3} P(14 - c(a))

The residual symmetry is the change of generator basis GL(3,2) (order 168) acting on F_2^3
and hence permuting the 15 types; the canonical form is the lexicographically smallest count
vector over that action.  GL(3,2) is legitimate here because the objective depends only on
the SUBSPACE L_3; the weight condition is imposed as "the set {a : c(a) = 6} contains a
basis", which is GL(3,2)-invariant, rather than as a condition on three named rows.

A CONSEQUENCE THAT BIT.  Because the canonical representative is GL(3,2)-canonical, its
weight-8 elements are generally NOT e_1,e_2,e_3.  Emitting rows 1..3 blindly produces
generators of the wrong weight.  The realisation routine therefore finds three INDEPENDENT
a with c(a) = 6 and emits h_a for those.  The run-time weight check caught this during
calibration, which is exactly what it is there for.

RESULT (exhaustive, 0.26 s):
    search nodes                       26453950
    complete multisets                  2769506
    surviving weight-8                   102759
    surviving commutation                 13493
    surviving rank 3                      13421
    surviving purity (F1)                 10341
    RAW optimal-feasible solutions        10341
    INEQUIVALENT classes                   1057
    distinct |U_3| values                    58
    best |U_3|                            11186   (raw 148, inequivalent classes 10)

-----------------------------------------------------------------------------------------
4.  STAGES 2-4, AND THE SYMMETRY THAT MAY BE USED THERE
-----------------------------------------------------------------------------------------
With L_r frozen, a completion is a k-dimensional isotropic subspace of Q_r spanned by
classes carrying a weight-10 lift, and the objective follows from (COV) with no
approximation.  The equivalence available for identifying two completions is NOT the full
group: it is the subgroup FIXING THE FROZEN PREFIX -- the stabiliser of the Stage-1 column
multiset inside (qubit permutations) x (local Cliffords), together with the changes of
generator basis preserving L_r setwise.  Canonicalising each new block on its own would
merge genuinely inequivalent branches, and is not done.

MEASURED STAGE-2 COST (4 classes sampled evenly across the 1057):
    dim Q_3 = 22,  |Q_3| = 4194304,  weight-10-liftable classes ~2.06e6  (49 %)
    cov_max                5796 - 5872        threshold |U_3|/2 = 5593 - 5723
    S_top(255)             1.47e6 - 1.49e6    need 2^7|U_3| = 1.43e6 - 1.47e6
    (F2)                   PASSES with a margin of only 1.6 - 2.6 %
    (F2) applied to ALL 1057 Stage-1 classes  ->  1057 pass, 0 rejected
    3-dim isotropic subspaces of Q_3          5.490e+16   per Stage-1 class
    measured DFS rate                         ~1.0e7 nodes/s
    projected exhaustive Stage 2, ONE class   ~1.5e6 hours  (~170 years)
    ... times 1057 classes                    ~1.8e5 years

WHY (F2) CANNOT SAVE IT.  A random class already covers about half of U_3, so the required
sum over the 2^k - 1 elements of D is only marginally above the average of an arbitrary
choice.  The bound is exact and correctly implemented, but the quantity it bounds simply has
no headroom, so it rejects nothing at all: run over the complete Stage-1 family it retained
1057 classes out of 1057.  That is a property of the problem, not of the
implementation, and no amount of tightening of THIS bound will change it.

-----------------------------------------------------------------------------------------
5.  WHAT THIS MEANS FOR THE STAGED STRATEGY
-----------------------------------------------------------------------------------------
Per-parent costs, computed by --analyze:

    stage   frozen r   dim Q_r   k   k-dim isotropic subspaces of Q_r (per parent)
      1         0         28     3   column multisets C(28,14) = 4.012e+07   EXHAUSTED
      2         3         22     3   5.490e+16                               WALL
      3         6         16     3   2.094e+11                               reachable
      4         9         10     2   8.696e+04                               trivial

So the wall is Stage 2 alone.  Stage 3 is about 5.8 hours per parent on one thread and
Stage 4 is instant; if a 6-generator prefix could be handed to the program the rest of the
ladder would run.  The obstruction is structural: freezing only 3 generators leaves
dim Q_3 = 22, and the 3-dimensional isotropic subspaces of a 22-dimensional symplectic space
number 5.5e16.

An important caveat about "just use smaller blocks": splitting Stage 2 into 1+1+1 does NOT
reduce the total, because 5.49e16 counts the 3-dimensional subspaces themselves, however one
walks to them; it only redistributes the work.  Any real gain has to come from a pruning
rule strictly stronger than (F1)/(F2), and the measurement above shows the coverage bound
cannot be it.  Options that preserve exactness, in decreasing order of how much they help:

  (a) FREEZE MORE FIRST.  Every extra frozen generator removes 2 from dim Q, dividing the
      Stage-2 count by roughly 2^(3*2) = 64 per generator.  A first block of 5 or 6
      generators would bring the second stage into range, at the cost of a much larger first
      stage (the column-multiset alphabet grows from 15 types at r = 3 to 155 at r = 5 and
      651 at r = 6, and the multiset count from 4.0e7 to 1.2e21) -- so this trades one hard
      stage for another.
  (b) SMALLER BLOCKS.  1 + 1 + ... instead of 3 + 3 + 3 + 2 makes each stage a sweep over
      |Q_r| classes rather than over k-dimensional subspaces, i.e. 4.2e6 instead of 5.5e16
      at the first extension.  The cost reappears as branching across stages, but it is a
      genuine depth-first tree with (F1)/(F2) applied at every level, and it is the natural
      exact continuation of this framework.
  (c) ACCEPT NON-EXHAUSTIVENESS at Stage 2 and label it BEST KNOWN.  Not done here without
      instruction, because the request was explicitly for correctness over a small space.

-----------------------------------------------------------------------------------------
6.  DEGENERACY
-----------------------------------------------------------------------------------------
The stated target |U_{<=4}| = 0 forbids weight-<=4 elements anywhere in L_11^perp, in
particular inside L_11 itself, so it FORCES a pure code.  A degenerate code is instead
allowed weight-<=4 stabilizer elements and its figure of merit is |U \ S|.  Both are
monotone under adding generators (if e is in U_new \ S_new then e is in U_old, and it cannot
be in S_old since S_old <= S_new), so the staged framework applies unchanged.  Filter (F1)
is valid only for the pure target and is switched off by --degenerate.

-----------------------------------------------------------------------------------------
7.  CORRECTNESS
-----------------------------------------------------------------------------------------
The program aborts loudly on any failed consistency check (check() -> exit 2), and did so
during calibration on the weight-8 realisation bug described in section 3.

  * closed-form objective vs brute force over all 91770 errors: 400 random isotropic
    prefixes of ranks 1..6, 0 mismatches
  * weight counting cross-checked independently
  * rank and commutation computed by independent routines
  * monotonicity |U_{r+1}| <= |U_r|: 100 random chains, 0 violations
  * canonical form is GL(3,2)-invariant: all 168 group elements agree
  * a deliberately different multiset is NOT merged
  * the coverage identity (COV) vs brute force on random isotropic triples, per class
  * every realised Stage-1 class is re-checked: three generators of weight exactly 8,
    pairwise commuting, rank 3, and |U_3| from brute force equal to the enumerator's value

-----------------------------------------------------------------------------------------
8.  HOW TO RUN
-----------------------------------------------------------------------------------------
    cd "...\stabilizer_14q\staged_3w8_3w10_3w10_2w10"
    .\staged_search.exe --selftest 400        correctness suite      (seconds)
    .\staged_search.exe --analyze             search-space accounting
    .\staged_search.exe --stage1              EXHAUSTIVE Stage 1     (0.26 s)
    .\staged_search.exe --calibrate2 8 3000000    Stage-2 measurement on 8 classes
    .\staged_search.exe --calibrate2 1057 2       apply (F2) to all 1057 classes

    --degenerate      score |U \ S| and switch off (F1)
    --out DIR         output directory

Stage 1 writes stage1/classes.txt (all 1057 with objective, min weight, raw multiplicity
and column-type counts) and stage1/representatives.txt (explicit weight-8 generators).

No long search is launched by any of these; --calibrate2 is bounded by its node cap.
