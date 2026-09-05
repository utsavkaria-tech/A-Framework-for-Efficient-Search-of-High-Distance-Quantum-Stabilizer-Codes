=========================================================================================
stab14_kopt -- exhaustive 3-opt and targeted 4-opt around the best-known 3+8 solutions
=========================================================================================
SELF-CONTAINED EXPERIMENT.  Everything is inside targeted_3opt4opt/.  stab14_2stage.cpp,
stab14_2stage.exe, stage1_orbit_0001..0045/, overnight_orbit0006/, overnight_refine/,
stab14_refine.cpp and every comparison/result file outside this directory are neither
modified nor written to.  The seven orbit incumbents are compiled into the source; they
were READ once from stage1_orbit_*/summary.txt (orbit 0006 uses the newer 11181 solution
supplied by the user, not the 11179 one in its directory).

-----------------------------------------------------------------------------------------
1.  THE REFORMULATION THAT MAKES EXHAUSTIVE 3-OPT POSSIBLE
-----------------------------------------------------------------------------------------
With S3 frozen the eight Stage-2 generators span an 8-dimensional totally isotropic
D <= V = C(S3)/S3 (dim 22), and  M_8(D) = (1/128) SUM_{v in D\{0}} cov1(v).

A k-opt move keeps K = span of the 8-k retained generators and asks for an 8-dimensional
isotropic D' with K <= D' <= K^perp.  Those D' are in bijection with the k-dimensional
TOTALLY ISOTROPIC SUBSPACES of the quotient

        Q = K^perp / K ,     dim Q = 22 - 2(8-k) = 6 + 2k ,

which carries a non-degenerate symplectic form.  So

        k = 3  ->  dim Q = 12,  |Q| =  4096   (phi table 16 KB, L1-resident)
        k = 4  ->  dim Q = 14,  |Q| = 16384   (64 KB)

rather than the 2^22 = 4.2M-class, 16 MB tables the previous refiner probed.  f_K is
constant on K-cosets, so it descends to phi(q) = SUM_{u in K} cov1(lift(q)^u) on Q, and

        S(D') = S(K) + SUM over the 2^k - 1 non-zero q of D'/K of phi(q).

Fixing q1 and putting psi(w) = phi(w) + phi(w^q1) gives the three-term pair form
        S = S(K) + phi(q1) + psi(q2) + psi(q3) + psi(q2^q3),
so ONE table rebuild per q1 (4096 adds) serves the whole inner double loop.  For k = 4 a
second pass psi2(w) = psi1(w) + psi1(w^q2) reduces the bottom two levels to the identical
pair problem.  Nothing is proxied: this IS the exact union coverage.

CANONICALISATION.  q1 is forced to be the MINIMUM of the seven non-zero span elements via
P(x) = (x > q1) && ((q1^x) > q1), applied to q2, q3 and q2^q3.  That removes most of the
168 bases per subspace for three integer comparisons and does not conflict with the
psi-sorted pruning order.

PRUNING.  Candidates are sorted by psi descending, and since psi(q2^q3) <= max psi the
partial sums give a hard break.  Measured: one drop-triple of orbit 0006 leaves only 104
leaf evaluations out of ~10^9 candidate subspaces, because a strong incumbent makes
base + 2*psi + max psi fall below the threshold almost immediately.

THE 4-OPT TAU FILTER (the fix that made 4-opt practical).  A qualifying pair needs
psi2(a) + psi2(b) + psi2(a^b) > need - base, and psi2(b), psi2(a^b) <= p2max, so BOTH
members must satisfy psi2 > tau = need - base - 2*p2max.  Filtering on tau before sorting
shrinks the inner list from several thousand to a handful.  Sorting the UNFILTERED list was
the real bottleneck: M = 1024 went from "still running after 600 s" to 200 s.

-----------------------------------------------------------------------------------------
2.  WHY THE 3-OPT RESULT IS A REAL CERTIFICATE
-----------------------------------------------------------------------------------------
For each of the C(8,3) = 56 drop-triples EVERY 3-dimensional isotropic subspace of Q is
enumerated -- no breadth limit, so the requested M = 256 is strictly subsumed.  The 3-opt
neighbourhood also CONTAINS the 1-opt and 2-opt neighbourhoods (dropping three generators
and re-adding two of them is a 1-opt move), which is an independent consistency check
against the earlier exhaustive refiner: both say the orbit-0006 incumbent is optimal there.

Weight-10 realisability is checked exactly, never assumed: a q is usable iff its K-coset
contains a class reachable by a weight-10 Pauli of C(S3), and one such representative is
recorded per q.  Only about half the q turn out to be usable, and section 7 shows that is
a parity invariant rather than a defect, so nothing legal is skipped.

POSITIVE CONTROL.  "No improvement found" is only meaningful if the engine can find one.
--degrade N replaces N generators of the incumbent by random valid weight-10 classes and
then runs the search.  Measured on orbit 0006:
        degrade 1 -> M_8 11155 (residual 31)  --> recovered to 11181 in 1.5 s
        degrade 2 -> M_8 11149 (residual 37)  --> recovered to 11181 in 1.5 s
        degrade 3 -> M_8 11125 (residual 61)  --> recovered to 11181 in 1.6 s
The engine climbs back to exactly 11181 every time, so 11181 is a genuine attractor and
the null results below are not an artefact of a broken sweep.

-----------------------------------------------------------------------------------------
3.  MEASURED COST MODEL
-----------------------------------------------------------------------------------------
  per-orbit setup (2^25 Gray sweep + 2^22 Walsh transform)      ~20 s
  exhaustive 3-opt, all 56 drop-triples                         ~0.7 s per orbit
  4-opt, 70 drop-quadruples, breadth M1 = M2 = 1024             200 s   (measured)
  4-opt scales as m1 * m2; unlimited is m1 ~ 16000, m2 ~ 8000
        => fully exhaustive 4-opt  ~ 9.0e9 inner iterations at 2.7 us  ~ 6.7 h
  intermediate points: M = 2048 ~ 13 min, M = 4096 ~ 53 min, M = 8192 ~ 3.6 h

-----------------------------------------------------------------------------------------
4.  RESULTS
-----------------------------------------------------------------------------------------
PHASE 1 -- orbit 0006, exhaustive 3-opt:      11181 -> 11181, NO IMPROVEMENT
PHASE 2 -- orbits 0002/0012/0023/0031/0036/0040, exhaustive 3-opt:
                                              11179 -> 11179 on every one, NO IMPROVEMENT

  Orbit    Starting M8    Best M8    Improvement    Residual    W4 residual
  ---------------------------------------------------------------------------
  0006     11181          11181      0              5           3
  0002     11179          11179      0              7           0
  0012     11179          11179      0              7           0
  0023     11179          11179      0              7           0
  0031     11179          11179      0              7           0
  0036     11179          11179      0              7           0
  0040     11179          11179      0              7           0

  Global best after 3-opt: orbit 0006, M_8 = 11181 / 11186, M_final = 91765 / 91770.

PHASE 3 -- fully exhaustive 4-opt on orbit 0006: see final_report.txt / phase3_exhaustive.log.

-----------------------------------------------------------------------------------------
5.  WHAT IS PROVEN AND WHAT IS NOT
-----------------------------------------------------------------------------------------
PROVEN
  Stage-1 optimum  M_3 = 80584 / 91770, |U_3| = 11186          (earlier work, exhaustive)
  Orbit 0006 incumbent is 1-opt optimal                        (earlier exhaustive refiner)
  Orbit 0006 incumbent is 2-opt optimal                        (earlier exhaustive refiner)
  All seven top orbits are 3-OPT OPTIMAL: for each, every replacement of any three of its
  eight generators -- i.e. every 3-dimensional isotropic subspace of every K^perp/K over
  all 56 drop-triples -- was enumerated and none improves the incumbent.

NOT PROVEN
  Nothing here bounds M_8 from above.  11181 is BEST KNOWN, not optimal.  Larger
  simultaneous replacements (5-opt and beyond) are untouched, and no statement is made
  about whether 11182..11186 is reachable.  A residual of 0 would require the 11-generator
  stabilizer to have no weight-<=4 Pauli in its normaliser, i.e. a pure [[14,3,5]] code
  containing the frozen S3; nothing in this work settles that either way.

-----------------------------------------------------------------------------------------
6.  COMMAND LINE
-----------------------------------------------------------------------------------------
  --phase 3opt|4opt|calibrate      --pipeline            (phase 1 -> 2 -> 3)
  --orbit 0006 | --all-top-orbits
  --M N | --M1 N --M2 N            4-opt breadth; 0 = unlimited (exhaustive)
  --degrade N --seed S             positive control: worsen the incumbent, then search
  --threads N  --time-limit S  --verbose  --out DIR  --checkpoint N  --resume
  --phase calibrate                times one drop-triple and one 4-opt inner sweep and
                                   extrapolates, so a long run can be costed before it starts

OUTPUT LAYOUT
  <out>/orbit<NNNN>_3opt/best_known/     stage1.txt stage2_best.txt remaining_errors.txt
                                         final_matrix.txt summary.txt
  <out>/orbit<NNNN>_3opt/improved_NNNN/  one per improvement, written before best_known
  <out>/final_4opt/best_known/           4-opt result
  <out>/checkpoint.txt                   phase, orbit, drop index, incumbent, residual list
  <out>/phase12_table.txt, <out>/final_report.txt

-----------------------------------------------------------------------------------------
7.  A PARITY INVARIANT ON THE QUOTIENT (checked, not assumed)
-----------------------------------------------------------------------------------------
Only about half the quotient elements are weight-10 usable -- 2079 of 4095 for k = 3 and
8255 of 16383 for k = 4 -- where naive independence would predict 99.999%.  diag_usable.cpp
settles why by histogramming how many A10 classes each K-coset actually contains:

  k = 3, |K| = 32:   0 members: 2016 cosets |  26:6  28:24  29:128  30:336  31:746  32:839
  k = 4, |K| = 16:   0 members: 8128 cosets |  13:90  14:575  15:2391  16:5199

The distribution is bimodal with NOTHING in between: a K-coset holds either no weight-10
class at all or nearly all of its members.  That is the signature of a parity invariant --
the Pauli weight parity of an element of C(S3) is a linear function of its class (all three
Stage-1 generators have even weight 8), so only one parity class can contain weight-10
(even) representatives.  Cosets in the wrong parity class contribute 0; cosets in the right
one contribute 13-16 of 16 (the shortfall being members whose weight is 8 or 12 rather
than 10).

Consequence: the ~50% is a genuine constraint, not a lost search space.  Those q cannot
serve as a generator under any choice of representative, so sweeping the usable half IS the
complete set of legal k-opt moves and the exhaustiveness claim is unaffected -- it just
makes the search twice as cheap as the worst case.

-----------------------------------------------------------------------------------------
8.  PHASE 3 RESULT -- FULLY EXHAUSTIVE 4-OPT
-----------------------------------------------------------------------------------------
Run: stab14_kopt.exe --phase 4opt --orbit 0006 --M 0   (M = 0 means unlimited breadth)

  70 / 70 drop-quadruples swept, each with M1 = 8255 (the complete weight-10-usable set of
  its quotient) and unlimited M2.  Total hits: 0.  Wall clock 9762.3 s = 2.71 h on 12
  threads -- 2.5x faster than the 6.7 h predicted from the M = 1024 scaling, because the tau
  filter bites harder as the breadth grows.

  M_8 = 11181 -> 11181.  NO IMPROVEMENT.

  Note: the final_report.txt written by that run mislabels this as "breadth-limited ... NOT
  exhaustive"; that was a labelling bug (M = 0 IS unlimited), fixed in the source and
  documented in exhaustive_4opt/CORRECTION.txt.  The numbers are unaffected.

-----------------------------------------------------------------------------------------
9.  FINAL STATUS
-----------------------------------------------------------------------------------------
PROVEN
  M_3 = 80584 / 91770, |U_3| = 11186                       Stage-1 global optimum
  orbit 0006 incumbent is 1-opt optimal                    (earlier exhaustive refiner)
  orbit 0006 incumbent is 2-opt optimal                    (earlier exhaustive refiner)
  all seven top orbits are 3-opt optimal                   (this program, exhaustive)
  orbit 0006 incumbent is 4-opt optimal                    (this program, exhaustive)

  So the best-known solution cannot be improved by simultaneously replacing any one, two,
  three or four of its eight generators.  Escaping 11181 requires changing at least FIVE
  of the eight at once.

BEST KNOWN, NOT PROVEN
  M_8 = 11181 / 11186, M_final = 91765 / 91770 = 99.99455160 %
  residual 5, weights 0/0/2/3, span dimension 5, not a subgroup, all 14 qubits involved
  weight-4: 81078 of 81081 detected, 3 undetected

  Nothing here bounds M_8 from above.  5-opt and beyond are untouched.  Whether
  11182..11186 is reachable is open, and a residual of 0 would require a pure [[14,3,5]]
  code containing this S3 -- which nothing in this work settles either way.
