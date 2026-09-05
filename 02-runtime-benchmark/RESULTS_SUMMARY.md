# Runtime benchmark results — [[5,1,3]]_2 stabilizer search

Standard workload `--work 5`, NumPy backend, `--target span`,
Python 3.13.7 / NumPy 2.5.2. Candidate space: 4 x 10 binary matrix,
S_0 = 2^40 = 1,099,511,627,776.

## Final results

| Constraint | Runtime (s) | T / T_none | S_R / S_0 |
|---|---|---|---|
| None (unconstrained) | 13851.05 | 1 | 1 |
| Universal stabilizer | 342.99 | 0.02476280 | 0.01434920 |
| Weight spectrum W = {4} | 361.91 | 0.02612860 | 0.02446920 |
| Cyclicity | 0.000227 | 1.63887e-8 | 9.31323e-10 |

Speedups vs unconstrained: universal 40.4x, weight-4 38.3x,
cyclicity 6.10e7x.

## Why T/T_none exceeds S_R/S_0

S_R/S_0 is the idealised candidate-count ratio. The measured runtime ratio
must be LARGER, because the constrained search also pays to evaluate the
constraint. The two differ by exactly one factor: the per-candidate cost
c = T / (candidates processed),

    (T / T_none) / (S_R / S_0)  =  c_constrained / c_none

so the last column below IS the deviation from theory — nothing else
contributes, and it is >= 1 for every constraint.

| Constraint | candidates processed | ns per candidate | c / c_none |
|---|---|---|---|
| None | 1,099,511,627,776 | 12.60 | 1 |
| Universal | 15,777,115,200 | 21.74 | 1.726 |
| Weight W={4} | 26,904,200,625 | 13.45 | 1.068 |
| Cyclicity | 1,024 | 221.7 | 17.60 |

The ordering follows the cost of the constraint test:

- **Weight W={4}** is row-local. A whole block of 1,048,576 candidates is
  accepted or rejected on two table lookups, so the test is nearly free and
  c/c_none sits just above 1.
- **Universal** must resolve 33,131,941,920 candidates at the pair-stage
  test to accept 15,777,115,200 — 2.10 tested per accepted — plus extraction
  of the accepted set. Hence 1.73.
- **Cyclicity** processes only 1024 candidates, far too few to amortise
  NumPy's fixed per-call overhead, so its per-candidate cost is dominated by
  that overhead rather than by the constraint.

### Keeping per-candidate cost comparable

For the identity above to be meaningful the engines must do the same work
per candidate. Three asymmetries were found and removed; all three were
implementation artefacts, and none changed any computed result (every
checksum is unchanged):

1. **Tile size.** The unconstrained engine originally processed
   1,048,576-element arrays while the constrained engines processed far
   smaller, cache-resident ones. The identical work-5 operation sequence
   measured 21.45 ns/candidate on 1,048,576 elements versus 19.05 ns on
   164,025 — an 11% difference from array size alone, which had pushed the
   weight-4 ratio BELOW S_R/S_0. All engines now process in a common tile
   of `_PROC_TILE = 32768` candidates.
2. **Operand form.** The weight-4 engine materialised full (405,405) operand
   arrays where the unconstrained engine used broadcast (n,1)/(1,n)
   operands, costing extra memory reads per candidate. Both now use the
   broadcast form.
3. **Reduction width.** The unconstrained engine summed its base processing
   unit in uint64 while the constrained engines used uint32, worth about 2%.
   All engines now use the same width.

## What each S_R is

- **Universal** — rank(H) = 4 over GF(2) AND H_X H_Z^T = H_Z H_X^T (mod 2).
  S_R = prod_{i=0..3} (2^(10-i) - 2^i) = 1023*510*252*120 = 15,777,115,200.
- **Weight W={4}** — every row has Pauli weight 4, i.e. the row-weight
  vector is [4,4,4,4]. 405 row values have weight 4, so
  S_R = 405^4 = 26,904,200,625.
- **Cyclicity** — row_{i+1} = cyclic_shift(row_i, 1), the shift acting on
  qubit indices independently in the X and Z blocks. The matrix is
  determined by its first row, so S_R = 2^10 = 1024.

No constraint is combined with any other. None imposes distance, code
validity, canonical forms, or duplicate removal.

## Methodology

**Rows 1-3** enumerate the full 2^40 space and test each candidate,
processing only those that pass. The full search was NOT executed; each
figure comes from bounded samples over contiguous deterministic regions of
the same lexicographic enumeration, extrapolated by the exact ratio of
candidate counts (`sample_runtime * 2^40 / sample_candidate_count`).
Sample sizes: unconstrained 1/1024 of the enumeration, universal and
weight-4 1/8 each.

**Row 4 (cyclicity)** uses constructive generation: every possible first row
is enumerated and the remaining rows are generated from it by repeated
shifts, with no cyclicity test afterwards. Every matrix produced is cyclic by
construction and each is produced exactly once, so the enumerated set IS S_R.
This is a complete measurement of 1024 candidates, not an extrapolation.

**Row 4 is therefore methodologically different from rows 1-3** and that must
be stated wherever the number is used. Under the same traversal-and-test
methodology as the others, cyclicity measures 0.235 s: the traversal must
still visit 2^20 blocks and reject 261,888 of every 262,144 on the
row-0/row-1 prefix, so its runtime is floored at Theta(2^20) while
S_R = 2^10. Generating from the first row removes that floor.

### Measurement conditions

The final numbers are minima over interleaved repetitions: the four
benchmarks were run in rotation so all saw the same machine conditions, and
the minimum over repetitions was taken. This matters because the machine had
intermittent background load (OneDrive sync on the working directory) during
the session — observed spreads reached 21% on the unconstrained sample. For
publication the absolute runtimes should be re-taken on a quiet machine with
the working directory outside OneDrive; the ratios were stable throughout.

### Weight-4 sampling correction

Acceptance under W={4} requires row 0 to have weight 4, and weight-4 values
are not uniformly distributed over contiguous ranges of the row alphabet
(a calibration sample on row 0 in [0,8) accepted zero candidates). The
sampled region [128,256) holds 48 weight-4 values against a global density of
405/1024 = 0.39551, and runtime tracks that count almost exactly
(a 60-value region ran 1.258x longer, vs a count ratio of 1.25). The estimate
is therefore scaled by (global density / region density) = 1.0546875. Raw
uncorrected estimate: 343.15 s.

The correction is independently validated: it maps both sampled retention
fractions (0.02320045 from a 48-value region and 0.02900057 from a 60-value
region) onto 0.02446920, the exact theoretical value.

## Validation invariants

- Enumerated candidate count equals S_0 exactly in every traversal run.
- Retention matches the combinatorial S_R for all three constraints.
- Weight-4 accepts four identical weight-4 rows (rank 1) and an
  anticommuting weight-4 pair, confirming no rank or commutation test.
- Cyclicity accepts the all-zero matrix (rank 0) and the XZIII shift chain
  (which fails commutation), confirming the same.
- Constructive generation reproduces the filter version exactly: same 1024
  matrices, same checksum (40,960), same 15 span-mode target occurrences.
- Target-code occurrences are self-consistent: each nonzero span value of
  row 0 admits 14*12*8 = 1344 completions, and 15 * 1344 = 20,160.
- All engines agree with independent brute force on examined, retained and
  checksum at work = 1, 2, 5 across three toy specifications.
- Tiling changed no checksum, confirming it altered only cost, not results.

## Programs

| File | Experiment |
|---|---|
| `stabilizer_search_benchmark.py` | A — no constraints |
| `stabilizer_search_benchmark_universal.py` | B — universal stabilizer constraints |
| `stabilizer_search_benchmark_weight4.py` | C — weight spectrum W={4} |
| `stabilizer_search_benchmark_cyclic.py` | D — cyclicity |
| `stabilizer_search_benchmark_unconstrained_ORIGINAL.py` | pristine copy of A |
| `stabilizer_search_benchmark_BACKUP.py`, `..._universal_BACKUP.py` | pre-modification copies |

Reproduce any row, e.g.:

    py stabilizer_search_benchmark_weight4.py --n 5 --k 1 --d 3 \
       --constraint weight --weights 4 --backend numpy --target span \
       --work 5 --r0-start 128 --r0-count 128 --repetitions 2

    py stabilizer_search_benchmark_cyclic.py --n 5 --k 1 --d 3 \
       --constraint cyclic --mode generate --backend numpy --target span \
       --work 5 --repetitions 200
