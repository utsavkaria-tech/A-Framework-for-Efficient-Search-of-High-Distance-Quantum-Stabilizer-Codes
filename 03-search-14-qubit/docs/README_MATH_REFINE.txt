=========================================================================================
stab14_refine -- targeted refinement around the known 3+8 solution (Stage-1 orbit 0006)
=========================================================================================
This is a SEPARATE program.  stab14_2stage.cpp, stab14_2stage.exe, orbit_comparison*.txt,
stage1_orbit_0001..0045/ and overnight_orbit0006/ are not read for input and never written
to.  Everything this program produces goes under the directory given by --out (default
"refine").

INPUT (hard-coded, override by editing G3_STR / H8_STR at the top of the source)
  Stage 1, frozen:  g1 = IIIIXXXXIIZZXX
                    g2 = IIIIXXIIXXXXZZ
                    g3 = IIIIIIZZZZZZZZ
  Stage 2 incumbent: the eight weight-10 generators with M_8 = 11181 / 11186.
The program re-derives U_3 from g1..g3, checks |U_3| = 11186, and re-verifies the incumbent
from scratch before searching:  weights 10, [h,g] = 0, [h,h] = 0, rank 11, and a
brute-force recount over all 91770 errors that must equal the algebraic score.

-----------------------------------------------------------------------------------------
1.  WHAT MAKES AN EXACT k-REPLACEMENT SEARCH AFFORDABLE
-----------------------------------------------------------------------------------------
In V = C(S3)/S3 (dim 22) the eight generators span an 8-dimensional totally isotropic D and

      covered(D) = (1/128) * SUM_{v in D\{0}} cov1(v),   cov1(v) = #{E in U_3 : <v,E> = 1}.

The objective is a SUM OVER THE 255 NON-ZERO ELEMENTS OF THE SUBSPACE.  Let K be the span of
the generators we keep and

      f_K(w) = SUM_{u in K} cov1(w ^ u)        (the coset sum at w).

Then for any replacements a_1..a_k,   S(D') = S(K) + SUM over the 2^k - 1 non-zero
combinations c of f_K(c).   f_K is the XOR-convolution of cov1 with the indicator of K, and
because XOR by a fixed b is an involution one basis vector is absorbed by a single in-place
pass -- F(w) and F(w^b) both become F(w) + F(w^b).  So dim K passes tabulate f_K over all
2^22 classes in about 20 ms, and thereafter each candidate costs 2^k - 1 array lookups
instead of 2^(8-k) symplectic sums.  That is the whole trick, and it changes the arithmetic
from "sample the neighbourhood" to "sweep it".

  k = 1:  S = S(K) + f(a)                      -> exhaustive over ~2.0M classes, 0.5 s
  k = 2:  S = S(K) + f(a) + f(b) + f(a^b)      -> exhaustive over all C(8,2) = 28 drops and
          all ~33K^2/2 pairs each, 48 s.  Sorting the candidates by f turns the bound
          S <= S(K) + f(a) + f(b) + max f into a hard break, so the swept region collapses
          as soon as the incumbent is good.
  k = 3:  fix a and put g(w) = f(w) + f(w^a); then
          S = S(K) + f(a) + g(b) + g(c) + g(b^c)
          -- literally the same three-term pair form, so ONE extra tabulation pass reuses the
          k = 2 routine unchanged.  All 56 drop-triples are swept; the first replacement runs
          over the top-M classes by coset gain and the remaining two are exhaustive.

Independence is enforced exactly (a new pair (a,b) is admissible iff a^b is not in K), and
isotropy is a single popcount, so no constraint is ever relaxed.

-----------------------------------------------------------------------------------------
2.  WHY THE SEARCH IS ANCHORED AND BROADENS, RATHER THAN RESTARTING
-----------------------------------------------------------------------------------------
The first design used ruin-and-recreate kicks.  Measured behaviour: a uniformly random
weight-10 class has cov1 ~ 5593 against ~5790 for a good one, so a random replacement moved
the residual from 5 to 14 and the local search settled back only around 10-17 -- the kicks
were leaving the good basin and never returning.  Two changes followed.

  * The kick now rebuilds each dropped generator greedily-randomly, uniformly among the
    top-R classes by coset gain against what is already kept, so it lands in a DIFFERENT
    basin of COMPARABLE quality.  It is off by default (--kick).
  * The default strategy is ANCHORED ITERATIVE BROADENING: every sweep restarts from the
    verified incumbent and the only thing that grows is the breadth M of the first
    replacement in the 3-opt phase.  One cycle runs an untargeted sweep and then one sweep
    aimed at each residual error in turn; when the cycle fails, M doubles.  M = |pool| would
    be exhaustive 3-opt, so this is a monotone ladder towards exhaustiveness instead of a
    random walk away from a certified-good point.

RESIDUAL TARGETING is exact, not heuristic: a residual error E survives iff every generator
commutes with its class, so killing it REQUIRES some new generator v with <v, class(E)> = 1.
Restricting the first replacement to that half of the candidate pool is therefore a sound
way to aim the search at a specific residual error rather than hoping one falls out.

-----------------------------------------------------------------------------------------
3.  OBJECTIVE AND TIE HANDLING
-----------------------------------------------------------------------------------------
Primary: fewest residual errors.  Then, lexicographically, fewest weight-4, weight-3,
weight-2, weight-1 residuals.  The span sum fixes the residual COUNT but not its weight
profile, so equal-score candidates are still enumerated -- a tie in the count can be a
strict improvement in the profile.  The comparison is deliberately NOT a total order: two
solutions with the same count and the same profile compare equal, so a mere relabelling
never counts as an improvement, never triggers a NEW BEST report, and cannot cycle.

Nothing is proxied.  The score above is the exact union |C(h1) u ... u C(h8)| n U_3, and
every candidate that would become the incumbent is re-verified by rebuilding all 91770
errors and recomputing the union directly before it is allowed to replace best_known/.

-----------------------------------------------------------------------------------------
4.  RESULT SO FAR
-----------------------------------------------------------------------------------------
The supplied incumbent verifies exactly:
      M_8 = 11181 / 11186, residual 5 (weights 0/0/2/3), M_final = 91765 / 91770,
      weight-4: 81078 detected, 3 undetected of 81081 (total recomputed, not hard-coded),
      residual span dimension 5, NOT a subgroup, all 14 qubits involved.

  phase 1  exhaustive 1-opt  -- no improvement  (0.5 s)
  phase 2  exhaustive 2-opt  -- no improvement  (48 s)

THE INCUMBENT IS THEREFORE CERTIFIED OPTIMAL AGAINST EVERY REPLACEMENT OF ONE OR TWO OF ITS
EIGHT GENERATORS, exhaustively over all 2 037 794 weight-10 classes of C(S3)/S3 -- not over
a sample.  This is a local optimality certificate only: it says nothing about 3-or-more
simultaneous replacements, and NO impossibility claim is made about residual 4 or 0.

-----------------------------------------------------------------------------------------
5.  OUTPUT LAYOUT
-----------------------------------------------------------------------------------------
  <out>/best_known/          stage1.txt  stage2_best.txt  remaining_errors.txt
                             final_matrix.txt  summary.txt
                             written from the supplied incumbent BEFORE any search, and
                             replaced only by a fully re-verified improvement
  <out>/improved_NNNN/       one directory per improvement, same five files plus
                             new_best_report.txt, written before best_known/ is updated
  <out>/FULL_COVERAGE.txt    only if residual reaches 0
  <out>/refinement_report.txt  the final report

-----------------------------------------------------------------------------------------
6.  COMMAND LINE
-----------------------------------------------------------------------------------------
  --seed S           RNG seed (only matters with --kick)
  --time-limit S     wall-clock budget, seconds (default 3600)
  --threads N        default: all cores
  --top-m N          starting breadth of the 3-opt first replacement (default 48)
  --deterministic    no kicks; the broadening ladder only
  --kick             diversify with ruin-and-recreate instead of broadening M
  --no-3opt          stop after the 1-opt / 2-opt certification (about 50 s)
  --out DIR          output root (default refine)
  --verbose

Cost model for planning a long run: one 3-opt sweep costs roughly M * 56 * 1.2 s, and a full
cycle covers 1 untargeted + 5 targeted sweeps, so a cycle at breadth M is about M * 400 s.
M = 8 -> ~55 min, M = 16 -> ~1.8 h, M = 32 -> ~3.6 h.  Set --time-limit to the hours you
have; the ladder simply stops wherever it has got to, and best_known/ is always current.

-----------------------------------------------------------------------------------------
7.  MEASURED COST MODEL (replaces the estimate in section 6)
-----------------------------------------------------------------------------------------
Timings from a real M = 16 cycle (12 threads):
      untargeted sweep            1350 s
      the five targeted sweeps    56 + 19 + 153 + 140 + 70 s  =  438 s
      full cycle                  ~1790 s  ~=  M * 112 s
The targeted sweeps are an order of magnitude cheaper than the untargeted one, because
restricting the first replacement to classes that anticommute with a chosen residual halves
the pool AND raises the score floor, so the sorted-order break fires much earlier.  Ladder
cost from the incumbent, cumulative:
      M = 16  ~0.5 h      M = 32  ~1 h       M = 64  ~2 h
      M = 128 ~4 h        M = 256 ~8 h
so a 9-hour budget reaches roughly M = 256, i.e. the first replacement is swept over the
256 highest-coset-gain classes for each of the 56 drop-triples, with the other two
replacements exhaustive throughout.

-----------------------------------------------------------------------------------------
8.  RESULT OF THE FIRST REFINEMENT RUN
-----------------------------------------------------------------------------------------
      phase 1  exhaustive 1-opt                       no improvement    0.5 s
      phase 2  exhaustive 2-opt (all 28 drop-pairs)   no improvement     48 s
      phase 3  full 3-opt cycle at breadth M = 16
               (untargeted + one sweep per residual)  no improvement   ~30 min
               plus a partial M = 32 untargeted sweep

      Best Stage-2 residual: 5      Improvement: 0      Status: BEST KNOWN

The incumbent M_8 = 11181 / 11186 (residual 0/0/2/3, M_final = 91765 / 91770) therefore
stands, and it is certified optimal against all 1- and 2-generator replacements.  Reaching
residual 4 requires changing at least three generators simultaneously; that region has been
swept only to breadth M = 16 so far, which is why the overnight ladder continues from there.
No impossibility claim is made.
