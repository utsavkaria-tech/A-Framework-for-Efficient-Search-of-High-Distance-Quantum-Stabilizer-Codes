=========================================================================================
stab14_6opt -- 6-opt search around the orbit-0006 incumbent (M_8 = 11181 / 11186)
=========================================================================================
SELF-CONTAINED.  Everything this experiment produces lives in targeted_6opt/.  Nothing in
targeted_3opt4opt/, overnight_refine/, overnight_orbit0006/, stage1_orbit_*/ or any other
existing file or binary is read for input or written to.  The incumbent is compiled into
the source and re-verified from scratch at every startup.

-----------------------------------------------------------------------------------------
0.  THE HEADLINE YOU NEED BEFORE CHOOSING A RUN
-----------------------------------------------------------------------------------------
A genuinely EXHAUSTIVE 6-opt search is NOT computationally possible.  Measured, not guessed:

    6-dimensional isotropic subspaces of Q per retained pair : 4.887e17
    x C(8,2) = 28 retained pairs                             : 1.368e19 subspaces
    4-dimensional prefixes a nested DFS must walk            : 1.024e17
    measured cost per prefix (this machine, 1 thread)        : 831.6 us
    => exhaustive 6-opt on 12 threads                        : 7.1e12 s = 225,000 YEARS

So --M 0 exists and is correct, but it will never finish.  The program is therefore built
to run a BREADTH-LIMITED 6-opt search whose scope is stated explicitly in every report and
never described as exhaustive.  See section 6 for the runnable configurations.

-----------------------------------------------------------------------------------------
1.  THE NEIGHBOURHOOD, EXACTLY
-----------------------------------------------------------------------------------------
A 6-opt move retains exactly two of the eight weight-10 generators and replaces the other
six.  With K = span of the two retained classes (dim 2) the admissible results are the
8-dimensional isotropic D' with K <= D' <= K^perp, i.e. the 6-dimensional totally isotropic
subspaces of

        Q = K^perp / K ,    dim Q = 22 - 2*2 = 18 ,    |Q| = 262144 ,

which carries a non-degenerate symplectic form.  All 28 retained pairs are swept.

Score telescoping (exact -- this IS the union coverage, not a proxy).  With
phi(q) = SUM_{u in K} cov1(lift(q)^u) and psi_{d+1}(w) = psi_d(w) + psi_d(w ^ q_{d+1}),

    S(D') = S(K) + phi(q1) + psi1(q2) + psi2(q3) + psi3(q4)
                 + psi4(q5) + psi4(q6) + psi4(q5^q6),

so the last two replacements are again the three-term pair problem and are ALWAYS swept
exhaustively.  M_8 = S/128, and every accepted improvement is re-verified by rebuilding all
91770 errors and recomputing the union directly.

-----------------------------------------------------------------------------------------
2.  HOW THE 4-OPTIMALITY RESULT WAS (AND WAS NOT) USED
-----------------------------------------------------------------------------------------
The incumbent is exhaustively 1-, 2-, 3- and 4-opt optimal.  What that licenses:

  VALID.  The 4-opt sweep enumerated every 8-dimensional isotropic D' containing span(T)
  for every 4-element subset T of the incumbent's basis {h1..h8}.  So any 6-opt candidate
  D' that happens to contain the span of FOUR OR MORE of the incumbent's basis generators
  was already examined and is known to be non-improving, and may be skipped.

  WHY IT IS NOT IMPLEMENTED AS A PRUNE.  A 6-opt candidate already contains the two
  retained generators; for the rule to fire it would have to contain two MORE incumbent
  generators as well.  Among the 4.887e17 candidate subspaces per retained pair that is a
  vanishing fraction, and testing for it costs more per node than it saves.  It is recorded
  here as mathematically valid but computationally worthless.

  INVALID, AND DELIBERATELY NOT USED.  "Every 4-subset of this 6-replacement is
  non-improving, therefore the 6-replacement is non-improving."  False: six generators can
  have collective synergy that no four of them show.  Nothing in this program uses that
  inference, and no branch is discarded on the strength of a smaller replacement's score.

  A SUBTLETY WORTH RECORDING.  Even the valid rule must say "contains the span of four
  BASIS generators", not "meets the incumbent in dimension >= 4".  The 4-opt sweep covered
  D' containing span(T) for coordinate subsets T only; a D' meeting D in a 4-dimensional
  subspace that is not spanned by four of the chosen basis vectors was never examined.

  CONCLUSION: 4-optimality contributes nothing usable to the 6-opt bound.  The pruning that
  actually works is the psi-sorted bound and the tau filter of section 3.

-----------------------------------------------------------------------------------------
3.  PRUNING, AND WHY IT IS WEAKER HERE THAN AT 4-OPT
-----------------------------------------------------------------------------------------
tau filter: a qualifying pair needs psi4(a) + psi4(b) + psi4(a^b) > need - base4, and both
psi4(b) and psi4(a^b) are at most p4max, so BOTH members must satisfy
psi4 > tau = need - base4 - 2*p4max.  Candidates are then sorted by psi4 descending and the
partial sums give a hard break.  Both are exact; nothing that could beat the incumbent is
discarded.

Why it bites less than in the 4-opt program: there K had dimension 4-6, so base already
accounted for most of the 255-element span sum and the threshold was tight.  Here K has
dimension 2 and at depth 4 only 15 of the 63 quotient elements are fixed, so base4 sits far
below `need` and most candidates survive tau.  Measured consequence: one outer block at
M2=M3=M4=64 scores 2.0e9 leaves.  This is the honest reason 6-opt is expensive -- not an
implementation defect.

-----------------------------------------------------------------------------------------
4.  SYMMETRY -- WHAT IS USED AND WHAT IS DELIBERATELY NOT
-----------------------------------------------------------------------------------------
USED: the objective depends only on the SUBSPACE D', so the six replacements are chosen as
a subspace of Q rather than as an ordered tuple of Pauli operators, and the weight-10
representative of each quotient element is picked once (the lexicographically smallest
weight-10 class in the K-coset).  Different Pauli lifts of the same class are never
enumerated separately.

USED: only about 50% of Q is weight-10 usable, and that is a genuine parity invariant, not
a defect -- the Pauli weight parity of an element of C(S3) is a linear function of its
class (all three Stage-1 generators have even weight 8), so only one parity class can hold
weight-10 representatives.  This was verified by histogram in the 3-opt/4-opt project
(diag_usable), where K-cosets held either 0 or nearly all of the weight-10 classes and
nothing in between.  Elements in the wrong parity class cannot serve as a generator under
ANY choice of representative, so skipping them removes no legal move.

NOT USED: a canonical ordering on the six replacement directions (for instance forcing q1
to be the minimum of the span, as the 3-opt code does).  It would remove redundant bases of
the same subspace, but it conflicts with the psi-sorted order that the pruning bound needs,
and at a finite breadth the ordering would change WHICH subspaces get visited.  Visiting a
subspace more than once is wasteful; silently excluding one would be wrong.  Redundancy was
chosen over risk, and is documented rather than hidden.

-----------------------------------------------------------------------------------------
5.  CALIBRATION (measured on this machine, 12 threads)
-----------------------------------------------------------------------------------------
    per-prefix cost                        831.6 us (single thread)
    one outer block at M2=M3=M4=32          26.7 s      32768 prefixes
    one outer block at M2=M3=M4=64         218.0 s     262144 prefixes
    cost scales as M1 * M2 * M3 * M4

    M (all four)   wall clock (12 threads)   fraction of the exhaustive prefix space
    ------------   -----------------------   ---------------------------------------
      32            0h 33m                    2.9e-10
      48            2h 52m                    1.5e-09
      64            9h 03m                    4.6e-09
      96           45h 50m                    2.3e-08
     128          144h (6.0 days)             7.3e-08
       0          225,000 years               1  (not runnable)

The ETA printed during the run is recomputed from measured throughput, not from this table.
Note the calibration uses the q1 with the largest phi, which has the least-filtered inner
lists, so the real run is likely somewhat FASTER than the table suggests.

-----------------------------------------------------------------------------------------
6.  HOW TO RUN IT
-----------------------------------------------------------------------------------------
Recommended overnight run (about 9 hours):

    cd "C:\Users\utsav\OneDrive\ドキュメント\Sidon Set Search\stabilizer_14q\targeted_6opt"
    .\stab14_6opt.exe --M 64

Shorter (about 3 hours):   --M 48
Longer  (about 6 days):    --M 128
Resume after any interruption:   add --resume

Other flags: --M1/--M2/--M3/--M4 to shape the breadth per level (cost is the product),
--threads N, --time-limit SECONDS, --drop-pair N (0..27) to run one retained pair,
--calibrate, --selftest, --out DIR.

Checkpointing: checkpoint.txt is rewritten after every completed outer block and records
the breadth, the retained-pair index, the number of CONTIGUOUSLY completed outer blocks
within it, all counters, the best solution and its residual list.  --resume skips exactly
the completed blocks -- it never silently repeats or skips an unverified region.  Power
loss, sleep or a kill costs at most one outer block.

Logging: every line goes to 6opt_search.log (flushed immediately) and to stdout.  A full
status block -- elapsed, position, percentage, prefixes, leaves, pruned, best M_8,
improvement, residual count and list, throughput, estimated total and remaining time,
timestamp, checkpoint location -- is written every 30 minutes and at every drop-pair
boundary.  No interactive terminal is needed.

-----------------------------------------------------------------------------------------
7.  WHAT MAY AND MAY NOT BE CLAIMED AFTERWARDS
-----------------------------------------------------------------------------------------
If an improvement is found it is real: weights, both commutation blocks, rank 11, the
rebuilt target set, the explicit test of all 11186 targets and an independent recount of
the union must all pass before it is accepted, and each improvement is written to its own
directory under best_known/ (never overwritten).  If M_8 reaches 11186 the program also
re-tests all 91770 weight-<=4 errors against the full 11-generator stabilizer and certifies
that N(L) contains no non-zero Pauli of weight <= 4, i.e. PURE DISTANCE >= 5 -- a [[14,3,5]]
construction -- writing FULL_COVERAGE_[[14,3,5]].txt.

If no improvement is found at a finite breadth, the correct statement is only:

    "The configured breadth-limited 6-opt search found no improvement; M_8 = 11181 stands
     as best known."

It is NOT 6-optimality, and it is certainly not global optimality.  At --M 64 the run
covers roughly 4.6e-09 of the 6-opt neighbourhood.  The program prints exactly this
distinction in final_report.txt and refuses to print the word "exhausted" unless every M
was 0 and every outer block completed.
