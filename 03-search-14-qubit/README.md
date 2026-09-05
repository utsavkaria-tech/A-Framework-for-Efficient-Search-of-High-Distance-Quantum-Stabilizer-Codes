# 03 — The 14-qubit search: applying the framework to `(n, k, q) = (14, 3, 2)`

Section 5 of the paper.

The framework is only useful if it makes a real search tractable. This section applies it
to the `[[14,3,5]]_2` parameters — a set for which no code is known — and produces a
`[[14,3]]_2` stabilizer code that detects **91,765 of the 91,770** Pauli errors of weight
`<= 4`, leaving **5 undetected**.

At physical error probability `p = 0.01` that reaches **99.999999933%** and **99.9733%**
of the detection and correction performance of the theoretically optimal `[[14,3,5]]_2`
code.

---

## The idea in one paragraph

Because semi-perfect codes have finite weight spectra, the search decomposes into `|W|`
sequential stages. **Stage 1** finds the inequivalent classes of the first `f_1` generators
of weight `w_1` that maximise detection of errors of weight `w < d_max`. **Each subsequent
stage** freezes what came before and maximises detection over only the errors still
undetected. For `(14,3,2)` the estimated spectrum is `{8:3, 10:8}`, so there are two
stages: three weight-8 generators, then eight weight-10 ones.

The scoring identity is what makes it fast. The objective depends only on the **weight
enumerator of the stabilizer group**, via a Walsh–Hadamard transform: for a stabilizer
group `L` of rank `r`,

```
    #{E : wt(E) = w, s(E) = 0} = (1/2^r) SUM_{h in L} Ew[wt(h)][w]
```

so scoring a candidate costs 15 table lookups (~70 ns) instead of 91,770 symplectic
products (~950 us) — a measured **13,000x** speed-up per node. Combined with local-Clifford
and permutation symmetry reduction (662M nodes -> 45M, a further 14.6x) and admissible
branch-and-bound pruning, Stage 1 becomes an exhaustive proof rather than a heuristic.

Full derivations: [`docs/README_MATH.txt`](docs/README_MATH.txt) and its companions.

---

## Results

| Stage | Generators | Detected | Status |
|---|---|---|---|
| Stage 1 | 3 x weight-8 | 80,584 / 91,770 (87.811%) | **Proven global optimum** |
| Stage 2 (best of 45 orbits) | + 8 x weight-10 | 91,763 / 91,770 (99.99237%) | Best known |
| **After refinement (orbit 0006)** | + 8 x weight-10 | **91,765 / 91,770 (99.99455%)** | **Best known — the published code** |

Stage 1 solves exhaustively in a fraction of a second, yielding 148 distinct optimal column
multisets which reduce to **45 inequivalent orbits**. A complete Stage-2 sweep over all 45
takes about 19 minutes. Seven orbits reach `M_8 = 11,179`; a deeper pass on the one with
residual weight profile `1/1/2/3` (orbit 0006) improved it to `0/0/2/3` — 5 residual errors.

### The 5 residual undetected errors

```
weight 3:  IIIXIIIIIIIIZZ = X_4  Z_13 Z_14
           IIIIIIXZIIIXII = X_7  Z_8  X_12
weight 4:  YIIIIZIIIYIXII = Y_1  Z_6  Y_10 X_12
           IXIIIIZYXIIIII = X_2  Z_7  Y_8  X_9
           IIZIXIIIIIXXII = Z_3  X_5  X_11 X_12
```

(1-indexed, as in the paper; the result files use 0-indexed qubit labels.)

---

## Contents

### `src/` — the search programs

| File | Program | Role |
|---|---|---|
| `stab14_2stage.cpp` | `stab14_2stage` | **The main pipeline.** Stage 1 exhaustive optimum + complete orbit enumeration + Stage 2 per orbit |
| `stab14_refine.cpp` | `stab14_refine` | **Produces the published code.** Targeted 1/2/3-opt refinement of a frozen incumbent |
| `stab14.cpp` | `stab14` | Earlier single-stage search: 4 weight-8 generators, one shot |
| `stab14_ext.cpp` | `stab14_ext` | Extends the 4-generator solution with 7 weight-10 (or `{10,12}`) generators |
| `stab14_6opt.cpp` | `stab14_6opt` | Breadth-limited 6-opt neighbourhood search around the incumbent |
| `stab14_kopt.cpp` | `stab14_kopt` | General `k`-opt search with control runs |
| `staged_search.cpp` | `staged_search` | Four-stage variant: `3xw8 + 3xw10 + 3xw10 + 2xw10` |
| `reverse_8w10_3w8.cpp` | `reverse_8w10_3w8` | Stage order reversed: eight weight-10 first, then three weight-8 |
| `hitset14.cpp` | `hitset14` | Residual-column hitting-set formulation |
| `diag_usable.cpp` | `diag_usable` | Diagnostic for usable-class counting |

### `build/`

`Makefile` and `build_gcc.sh` for GCC/Clang (`-std=c++20 -O3 -march=native -flto -pthread`);
`build_*.bat` for MSVC (`/std:c++20 /O2 /Oi /Ot /GL /arch:AVX2 /LTCG`). The batch files call
`vcvars64.bat` from a hard-coded Visual Studio path — edit it for your install.

### `docs/` — the mathematics

One derivation per program, written before the code and kept in sync with it.

| File | Covers |
|---|---|
| `README_MATH.txt` | The core: column view, local-Clifford invariants, the Walsh-transform objective, symmetry removal, both branch-and-bound bounds, complexity |
| `README_MATH_2STAGE.txt` | The two-stage decomposition and the identity both stages are built on |
| `README_MATH_ORBITS.txt` | What the orbit experiment does and — explicitly — what it does *not* claim |
| `README_MATH_REFINE.txt` | Why an exact `k`-replacement search is affordable in the quotient space |
| `README_MATH_EXT.txt` | The 20-dimensional symplectic quotient the extension search lives in |
| `README_MATH_6OPT.txt`, `README_MATH_KOPT.txt`, `README_MATH_STAGED.txt`, `README_MATH_REVERSE.txt`, `README_MATH_HITSET.txt` | The corresponding alternative formulations |

### `results/`

| Directory | Content |
|---|---|
| `00-single-stage-4-generator/` | The earlier 4-generator line: `best_stabilizer_14q.txt` (proven optimum, 86,340/91,770, all 15 non-identity group elements of weight exactly 8, perfectly flat syndrome distribution) and its weight-10 extensions. **Not the paper's construction** — kept because it motivates the two-stage split |
| `01-stage1/` | `stage1_best.txt` — the proven Stage-1 optimum, with full verification block. `stage1_undetected_errors.txt` — the 11,186 errors it misses |
| `02-stage2/` | `stage2_best.txt`, `stage2_remaining_errors.txt` (the 7-residual solution), `final_11_generator_code.txt`, `final_summary.txt` with the identity checks `M_final == M_3 + M_8` |
| `03-orbit-sweep/` | `orbit_comparison.txt` — all 45 orbits ranked by `M_8`. `deep_pass_refinement.txt` — the 6x-budget second pass over the 7 best. `stage1_orbit_0001..0045/` — per-orbit Stage 1, Stage 2, final matrix, residual errors and summary |
| `04-final-code/` | **The published code.** `best_known/final_matrix.txt` (Pauli form + `H = [X\|Z]`), `best_known/remaining_errors.txt`, `best_known/summary.txt`, and `refinement_report.txt` with the 1/2-opt optimality certificate |
| `05-further-searches/` | Everything tried that did not improve on 5: seed-robustness sweep, 3/4-opt and 6-opt neighbourhoods, the reversed stage order, the four-stage variant, and the hitting-set formulation |

---

## Proven vs best known

Every result file states its own status in its header. Two statuses are used, and the
distinction is load-bearing:

**`PROVEN GLOBAL OPTIMUM`** — Stage 1 only. The branch-and-bound bounds are admissible
(each relaxes an exact minimisation), so pruning never removes a strictly better solution
and a completed search *is* an optimality proof. Stage 1 also carries an independent
brute-force recount that must equal the algebraic score before the result is written.

**`BEST KNOWN — OPTIMUM NOT PROVEN`** — Stage 2 and everything after it, including the
published code.

The strongest claim made about the final code, from `04-final-code/refinement_report.txt`:

> the incumbent is optimal against every 1-generator and every 2-generator replacement,
> exhaustively over all 2,037,794 weight-10 classes. No impossibility claim is made.

Reaching a residual below 5 therefore requires changing at least **three** generators at
once. That is what `05-further-searches/` sweeps — 3-opt, 4-opt, a breadth-limited 6-opt
(1,792 outer blocks, 470M prefixes, 1.18e12 leaves scored over 10h 52m), 16 independent
seeds, a reversed stage order and a four-stage decomposition. None improved on 5.

There is also a provable upper bound on Stage 2 recorded in `stage2_best.txt`:
`M_8 <= min(|U_3|, 255 * max_v cov1(v) / 128) = 11,186`, i.e. the 7-residual solution is
within 7 of the trivial ceiling.

---

## Verification built into the programs

Nothing is reported that has not been independently recomputed. Each run re-derives the
target error set from the generators and refuses to proceed unless its size and weight
distribution match; checks generator weights, all pairwise commutations and the rank; and
recounts coverage three independent ways — bitset popcount, brute force over all 91,770
errors, and the algebraic Walsh identity — which must all agree. `stab14_refine` re-verifies
its hard-coded incumbent from scratch before searching. `--selftest N` runs the whole check
suite on small instances.

---

## Running

```bash
make -f build/Makefile
./stab14_2stage --mode hybrid          # Stage 1 + 45-orbit Stage 2 sweep (~19 min)
./stab14_refine --time-limit 2400      # -> the published 5-residual code (~40 min)
```

Common flags: `--threads N`, `--seed S`, `--time-limit SEC`, `--deterministic`,
`--selftest N`, `--out DIR`. `stab14_2stage` adds `--stage1-index N` to run one orbit,
`--stage2-time-limit`, `--beam-size`, `--restarts` and `--save-all`. `stab14_refine` adds
`--top-m`, `--no-3opt` and `--kick`.

`stab14_refine` writes only under `--out` and never reads or modifies the Stage-1/Stage-2
artifacts. Its starting incumbent is hard-coded as `G3_STR` / `H8_STR` at the top of the
source.

See [`../docs/REPRODUCING.md`](../docs/REPRODUCING.md) §3 for detail.
