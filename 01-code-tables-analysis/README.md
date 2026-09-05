# 01 — Code-tables analysis: the finite-size correction factors `C` and `C'`

Sections 4.1 and 4.2 of the paper, and Section 2 of the supplementary material.

The asymptotic search-space retention factor `R` is a limit. Real codes are finite, and
the gap is the **finite-size correction factor**

```
    C  = H_{q^2}(delta_max) - R          delta_max = max stabilizer row weight / n
    C' = H_{q^2}(delta_avg) - R          delta_avg = mean stabilizer row weight / n
    R  = log_q(X) / (2n),   X = sum_{w in W} binom(n,w) (q^2 - 1)^w
```

This section measures `C` and `C'` empirically across every quantum code in Grassl's
tables that publishes a usable stabilizer matrix, and establishes their scaling:
`C = O(log n / n)` for general codes, `C' = O(1/n)` for semi-perfect ones.

---

## Dataset at a glance

| | `q=2` | `q=3` | `q=4` | `q=5` | `q=7` | `q=8` | total |
|---|---|---|---|---|---|---|---|
| `n` range | 1–256 | 1–100 | 1–100 | 1–100 | 1–100 | 1–100 | |
| examined | 33,152 | 5,150 | 5,150 | 5,150 | 5,150 | 5,150 | **58,902** |
| processed | 17,882 | 4,511 | 3,707 | 2,891 | 2,971 | 2,813 | **34,775** |
| skipped | 15,270 | 639 | 1,443 | 2,259 | 2,179 | 2,337 | **24,127** |
| semi-perfect | 384 | 249 | 236 | 245 | 353 | 442 | **1,909** |

Skip reasons: 23,361 had no parity-check matrix, 756 no stabilizer rows, 9 a matrix
belonging to a different code, 1 malformed. Every skipped set and its reason is logged in
`logs/skipped_codes.csv` — nothing is silently dropped.

---

## Contents

### `src/`

| File | What it does |
|---|---|
| `qecc_pipeline.py` | The whole data pipeline: **fetch** (codetables.de -> SQLite cache, rate-limited and resumable), **parse** (bounds + parity-check matrices out of cached HTML), **classify** (semi-perfect vs normal, using the complete upper-bound grid for that `q`), **export** (one `.xlsx` per `q`, 9 plots per `q`, skip log, run summary) |
| `analysis_scaling.py` | Power-law and `log(n)/n` envelope fits for `C`; nested linear predictor models from `(n)`, `(n,k)`, `(n,q)`, `(n,k,q)` |
| `analysis_semiperfect.py` | The same restricted to semi-perfect codes, plus the exact identity `D = C' - C = H(delta_avg) - H(delta_max)` (in which `R` cancels completely) verified to machine precision |
| `analysis_Cprime_semiperfect.py` | Focused envelope fit of `C'` on semi-perfect codes: two estimators of the same `(a ln n + b)/n` model, plus a power law for comparison, flagged unusable where its linear-space `R^2` goes negative |
| `correction_factor_estimator.py` | Cross-validated KNN / RF / Extra Trees / GP regressors predicting `C'` from `(n,k)`. Feeds the first column of Table 4.2.3 |

### `data/`

| File | Rows | Content |
|---|---|---|
| `all_codes_combined.csv` | 34,775 | Every processed code, all `q`. Columns: `q, n, k, Code Parameters, Type of Code, k/n, Weight Spectrum, Average Weight, Max Weight, R, C, C'` |
| `QECC_codes_q{2,3,4,5,7,8}.xlsx` | per `q` | The same data as filterable workbooks, one per field size |

`Type of Code` is either `Semi-Perfect Quantum Code` or `Normal Code`. `Weight Spectrum`
is a dict `{weight: frequency}`.

### `results/fits/`

Every fitted coefficient in the paper, as CSV.

| File | Content |
|---|---|
| `power_law_fits.csv` | **Table 4.1.4.** `alpha, beta, R2_log` for `C_max(n) ~ alpha n^beta`, per `q`, for both the envelope (`max_k`) and the per-`n` mean |
| `envelope_fits.csv` | Theory (`log n / n`) vs power-law fits side by side, with RMSE, bounding scale and coverage |
| `envelope_coverage.csv` | **§4.1.5.** Fraction of codes under the parameter-free bound (99.0–99.7%), violation counts, and the largest violating `n` |
| `q_collapse.csv` | **§4.1.4.** The `alpha ln(q^2)` collapse across field sizes |
| `predictor_R2.csv` | **§4.1.8.** `R^2` per feature set, per `q`, split at the entropy peak |
| `Cprime_sp_envelope_fits.csv` | **Table 4.2.1.** `b_q` is column `b_0`; `R^2` is column `R2_0_on_M` |
| `Cprime_sp_envelope_points.csv` | The `(q, n, M_q(n))` points those fits were made on |
| `Cprime_sp_D_summary.csv`, `Cprime_sp_D_points.csv`, `Cprime_sp_D_vs_n_bands.csv` | **§4.2.2.** The `C' - C` comparison: signs, magnitudes, Spearman correlation, per-length bands |
| `sp_envelope_fits.csv`, `sp_power_law_fits.csv`, `sp_predictor_R2.csv`, `sp_structure.csv`, `sp_codes.csv` | Semi-perfect-only versions, plus the full per-code table with `delta_avg`, `delta_max`, `H_gap`, `past_peak`, `uniform` |
| `_semi_perfect_working.csv` | Intermediate working table kept for auditability |

### `results/workbooks/` and `results/reports/`

The same three analyses as Excel workbooks, and as long-form generated Markdown reports
(`SCALING_ANALYSIS.md`, `SEMIPERFECT_ANALYSIS.md`, `CPRIME_SEMIPERFECT_RESULTS.md`).
`PIPELINE_NOTES.md` documents the pipeline's conventions; `run_summary.txt` / `.csv` is the
run ledger reproduced in the table above.

### `figures/`

- `by-q/q{2,3,4,5,7,8}/` — 9 scatter plots per field size (the `01_`–`09_` series). Paper
  Figures 4.1.2 and 4.2.1 are the `q2` files; supplementary Figures 2.1 and 2.2 are the
  `03_` and `06_` series across all `q`.
- `analysis/` — the publication figures: the `nC` vs `log n` collapse, the `C'` envelope in
  collapse coordinates, `C' - C` vs `n` (PNG **and** PDF), and the predictor `R^2` bars.

### `logs/`

`run.log` (full run log) and `skipped_codes.csv` (every skipped parameter set with its
reason).

---

## Running it

See [`../docs/REPRODUCING.md`](../docs/REPRODUCING.md) §1 for the full recipe. Short version:

```bash
cd src
python qecc_pipeline.py                   # builds codetables_cache.sqlite, then everything
python analysis_scaling.py
python analysis_semiperfect.py
python analysis_Cprime_semiperfect.py
```

---

## Two things worth knowing

**The SQLite cache is not in this repository.** `codetables_cache.sqlite` is 167 MB, above
GitHub's per-file limit. It is a pure download cache — no original results live in it — and
the fetch phase rebuilds it. Every figure, table and CSV derived from it *is* committed.

**`correction_factor_estimator.py` has a hard-coded absolute path.** It reads
`all_codes_combined.csv` from a Windows path near the top of the file. Change it to
`../data/all_codes_combined.csv` before running the script anywhere else. The `q` value is
also set by editing a variable rather than a flag.

---

## Conventions

- **Row weight** is the **symplectic** weight: the number of coordinates `j` in `1..n` with
  `(X_j, Z_j) != (0,0)`. This keeps `delta = W_i / n` in `[0,1]` so `H_{q^2}(delta)` is
  defined. `--weight-mode raw` counts non-zero symbols across all `2n` printed columns instead.
- **`X`** is summed over the **distinct** weights in the spectrum; the frequencies `f_i`
  enter only through the average weight. `--x-sum-mode rows` weights each term by `f_i`.
- `H_{q^2}(x) = x log_{q^2}(q^2-1) - x log_{q^2}(x) - (1-x) log_{q^2}(1-x)`.
- Nothing is invented: a page without a usable parity-check matrix is skipped and the
  reason recorded.
