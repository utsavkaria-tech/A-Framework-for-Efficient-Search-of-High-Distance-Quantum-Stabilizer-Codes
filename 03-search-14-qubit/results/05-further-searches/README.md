# Further searches — everything that did not beat 5 residual errors

The published code leaves 5 of the 91,770 weight-`<=4` errors undetected. Reaching 4 would
require changing at least **three** of the eight weight-10 generators at once, because the
incumbent is already certified optimal against every 1- and 2-generator replacement,
exhaustively over all 2,037,794 weight-10 classes.

These are the attempts to do that. All of them are negative results, and they are here
because they are the evidence for calling 5 a robust best-known value rather than merely
the first thing that worked.

| Directory | What was tried | Outcome |
|---|---|---|
| `seed-robustness-orbit0006/` | The Stage-2 heuristic re-run on orbit 0006 across 16 Fibonacci seeds (1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597) | `M_8 = 11,179` never exceeded by the two-stage construction |
| `overnight-refine-harness/` | Long-running refinement harness with a doubling breadth ladder; includes the 1-opt/2-opt certification run | No improvement. The certification is what licenses the "at least three generators" claim |
| `kopt-3opt-4opt/` | 3-opt and exhaustive 4-opt neighbourhoods, with four independent control runs (`ctrl_1`, `ctrl_2`, `ctrl_3`, `ctrl_re`) and a 4-opt probe | No improvement |
| `6opt/` | Breadth-limited 6-opt: 1,792 outer blocks, 470,082,124 prefixes examined, 1.18e12 leaves scored, 10h 52m | No improvement. **Not exhaustive** — with breadth `M1..M4 = 64` only part of the 6-opt neighbourhood was visited, so no optimality claim |
| `reverse-order-8w10-then-3w8/` | Stage order reversed: eight weight-10 generators first, then three weight-8 | Did not improve on the forward order |
| `staged-4-stage/` | Four-stage decomposition `3xw8 + 3xw10 + 3xw10 + 2xw10` instead of two stages | Calibration and self-test only; the search did not complete |
| `residual-column-hitting-set/` | The residual errors recast as a column hitting-set problem, searched to `r = 6` | Calibration, self-test and checkpoints; no improving solution |

Each directory keeps its own `final_report.txt` / `CALIBRATION.txt` / `SELFTEST.txt` and
run log. Where a run was aborted or produced no solution, the empty result directories have
been dropped rather than committed.

The corresponding mathematical derivations are in
[`../../docs/`](../../docs) — `README_MATH_6OPT.txt`, `README_MATH_KOPT.txt`,
`README_MATH_STAGED.txt`, `README_MATH_REVERSE.txt`, `README_MATH_HITSET.txt`.
