# Paper -> file map

Every numbered table and figure in the main paper and the supplementary material, and the
file in this repository that produced or contains it. Paths are relative to the repository
root.

---

## Main paper

### Section 4.1 — Asymptotic behaviour of the correction factor `C`

| Item | Content | Produced by | Data / figure file |
|---|---|---|---|
| §4.1.1 | 58,902 sets examined, 34,775 usable, per-`q` breakdown | `01-code-tables-analysis/src/qecc_pipeline.py` | `01-code-tables-analysis/results/reports/run_summary.txt` and `run_summary.csv` |
| Figure 4.1.2 (a) | `C` vs block length `n`, `q = 2`, 17,882 points | `analysis_scaling.py` | `01-code-tables-analysis/figures/by-q/q2/01_C_vs_n_all_gradient_q2.png` |
| Figure 4.1.2 (b) | `C` vs encoded qubits `k`, `q = 2` | `analysis_scaling.py` | `01-code-tables-analysis/figures/by-q/q2/02_C_vs_k_all_q2.png` |
| **Table 4.1.4** | Power-law fits `C_max(n) ~ alpha n^beta`, columns `alpha, beta, R^2_log` | `analysis_scaling.py` | `01-code-tables-analysis/results/fits/power_law_fits.csv` — rows with `series=C`, `target=max_k C(n,k)` (`alpha`, `beta`, `R2_log`) |
| §4.1.4 | Collapse `alpha ln(q^2) = 0.774 +/- 0.060`, `beta = -0.752 +/- 0.023` | `analysis_scaling.py` | `01-code-tables-analysis/results/fits/q_collapse.csv` (`alpha_x_lnq2` column) |
| §4.1.5 | Parameter-free bound satisfied by 99.0–99.7% of codes; no violation beyond `n = 118` | `analysis_scaling.py` | `01-code-tables-analysis/results/fits/envelope_coverage.csv` (`coverage_all`, `max_n_violation`) |
| §4.1.7 | 100% of codes with `C < 0` satisfy `delta_max > delta*` | `analysis_semiperfect.py` | `01-code-tables-analysis/results/reports/SEMIPERFECT_ANALYSIS.md` |
| §4.1.8 | Predictability from `(n, k, q)`: `R^2 ~ 0.948` below the entropy peak, `~ 0.122` above | `analysis_scaling.py` | `01-code-tables-analysis/results/fits/predictor_R2.csv`; figure `figures/analysis/predictor_R2.png` |

### Section 4.2 — Semi-perfect quantum codes and `C'`

| Item | Content | Produced by | Data / figure file |
|---|---|---|---|
| §4.2.1 | All 1,909 semi-perfect codes satisfy `lim C' = 0` | `analysis_Cprime_semiperfect.py` | `01-code-tables-analysis/results/reports/CPRIME_SEMIPERFECT_RESULTS.md` |
| Figure 4.2.1 (a) | `C'` vs `n`, semi-perfect, `q = 2`, 384 points | `qecc_pipeline.py` | `01-code-tables-analysis/figures/by-q/q2/04_Cprime_vs_n_semiperfect_gradient_q2.png` |
| Figure 4.2.1 (b) | `C'` vs `k`, semi-perfect, `q = 2` | `qecc_pipeline.py` | `01-code-tables-analysis/figures/by-q/q2/05_Cprime_vs_k_semiperfect_q2.png` |
| **Table 4.2.1** | Envelope fits `C'_max(n) ~ b_q / n`, columns `b_q, R^2` | `analysis_Cprime_semiperfect.py` | `01-code-tables-analysis/results/fits/Cprime_sp_envelope_fits.csv` — columns `b_0` (= `b_q`) and `R2_0_on_M` (= `R^2`) |
| §4.2.2 | `C' - C = 0` for 1,289 codes (67.52%), `> 0` for 344 (18.01%), `< 0` for 276 (14.45%); pooled Spearman `rho = -0.43`; range `-0.0919` to `0.0280` | `analysis_Cprime_semiperfect.py` | `01-code-tables-analysis/results/fits/Cprime_sp_D_summary.csv`, `Cprime_sp_D_points.csv`, `Cprime_sp_D_vs_n_bands.csv` |
| **Table 4.2.3** | Estimated weight spectra, `R`, `C'`, `w_avg` and runtime speedups for 4 unresolved semi-perfect codes | `correction_factor_estimator.py` (ML estimate of `C'`), then the manual procedure in Supplementary §1.8.3 | `01-code-tables-analysis/results/reports/CPRIME_SEMIPERFECT_RESULTS.md`; see the note below |

> **Table 4.2.3 is partly manual by design.** The `C'` estimates come from the
> cross-validated regressors in `correction_factor_estimator.py`; the weight spectra
> *without* frequencies were estimated by inspecting neighbouring codes, and the
> frequencies `f_i` were then selected by hand to match the derived `w_avg`. The
> supplementary material (§1.8.3) states the procedure explicitly. Nothing in this step is
> automated end-to-end, so the table is reproduced by following that procedure, not by
> running a single script.

### Section 4.3 — Runtime analysis

| Item | Content | Produced by | Data file |
|---|---|---|---|
| **Table 4.3** | Retention factor, runtime, experimental speedup, theoretical speedup and efficiency ratio `eta` for the `[[5,1,3]]` code under 4 constraint settings | `02-runtime-benchmark/src/*.py` | `02-runtime-benchmark/RESULTS_SUMMARY.md` (final-results table, plus per-candidate cost breakdown, sampling methodology and validation invariants) |

### Section 5 — Application to the `[[14,3,5]]_2` parameters

| Item | Content | Produced by | Data file |
|---|---|---|---|
| §5.2, Stage 1 | Exhaustive Stage-1 optimum, 3 weight-8 generators, 80,584 / 91,770 detected, **proven optimal** | `stab14_2stage` | `03-search-14-qubit/results/01-stage1/stage1_best.txt` |
| §5.2 | 45 inequivalent Stage-1 orbits; complete Stage-2 sweep; per-orbit `M_8`, residual and weight-4 coverage | `stab14_2stage --mode hybrid` | `03-search-14-qubit/results/03-orbit-sweep/orbit_comparison.txt` |
| §5.2 | Deep second pass over the 7 best orbits; residual weight profiles 3/3/1/0, 1/3/3/0, 1/1/2/3 | `stab14_2stage` | `03-search-14-qubit/results/03-orbit-sweep/deep_pass_refinement.txt` |
| §5.2 | 7-residual result, `M_8 = 11,179`, `M_final = 91,763` | `stab14_2stage` | `03-search-14-qubit/results/02-stage2/stage2_best.txt`, `final_summary.txt` |
| **§5.2.1** | Best-known generators `g1..g3` (weight 8) and `h1..h8` (weight 10) | `stab14_refine` | `03-search-14-qubit/results/04-final-code/best_known/final_matrix.txt` |
| **§5.2.2** | The 5 residual undetected errors (2 of weight 3, 3 of weight 4) | `stab14_refine` | `03-search-14-qubit/results/04-final-code/best_known/remaining_errors.txt` |
| §5.3 | 91,765 / 91,770 detected (99.99455%) | `stab14_refine` | `03-search-14-qubit/results/04-final-code/best_known/summary.txt` |
| **§5.3.1** | Canonical generator matrix | derived from the above | `03-search-14-qubit/results/04-final-code/best_known/final_matrix.txt` (`H = [X|Z]` block form) |
| §5.2 | 1/2-opt certification: incumbent optimal against every 1- and 2-generator replacement over all 2,037,794 weight-10 classes | `stab14_refine` | `03-search-14-qubit/results/04-final-code/refinement_report.txt` |

---

## Supplementary material

### §1 — Retention-factor derivations

Sections 1.1–1.8 are analytical derivations with no computational artifact. The
constraint definitions they derive are implemented in
`02-runtime-benchmark/src/` (universal, weight-spectrum, cyclicity) and used as search
constraints in `03-search-14-qubit/src/`.

### §2 — Additional plots

| Figure | Content | File |
|---|---|---|
| Figure 2.1 | `C` vs `k/n`, all codes, `q = 2,3,4,5,7,8` | `01-code-tables-analysis/figures/by-q/q{2,3,4,5,7,8}/03_C_vs_k_over_n_all_q*.png` |
| Figure 2.2 | `C'` vs `k/n`, semi-perfect codes, all `q` | `figures/by-q/q*/06_Cprime_vs_k_over_n_semiperfect_q*.png` |
| Figure 2.3 | `nC` vs `log n` collapse, codes below the entropy peak | `figures/analysis/collapse_nC_vs_logn.png` |
| Figure 2.4 | Upper envelope of `C'` in collapse coordinates | `figures/analysis/pub_Cprime_envelope.png` |
| Figure 2.5.1 | `C' - C` vs block length `n` | `figures/analysis/pub_Dprime_vs_n.png` (and `.pdf`) |
| Figure 2.5.2 | Sign of `C' - C` vs `delta_max` | `figures/analysis/sp_D_vs_deltamax.png` |
| Figure 2.6.1 | `R^2` for predicting `C'` from feature sets, semi-perfect codes | `figures/analysis/sp_predictor_R2.png` |
| Figure 2.6.2 | Explanatory power of feature sets for `C` | `figures/analysis/predictor_R2.png` |

Supporting per-`q` plots not printed in the paper: `07_H_delta_max_vs_R_q*.png`,
`08_H_delta_avg_vs_R_q*.png`, `09_Cprime_vs_n_semiperfect_vs_normal_q*.png`, plus
`figures/analysis/{fit_C_loglog,fit_C_prime_loglog,theory_vs_C,sp_Cprime_vs_n,sp_D_vs_n}.png`.

---

## Notation differences between paper and code

| Paper | Code / CSV | Note |
|---|---|---|
| `R` (retention factor) | `rate_R`, `R` | In the supplementary §1.1 only, `R` denotes the **code rate**; the main paper writes the rate as `r`. The CSVs use `rate_R` for the retention factor and `k/n` for the rate. |
| `delta_max`, `delta_avg` | `delta_max`, `delta_avg` | `max` / `mean` stabilizer row weight divided by `n` |
| `C` | `C` | `H_{q^2}(delta_max) - R` |
| `C'` | `C_prime`, `C'` | `H_{q^2}(delta_avg) - R` |
| `C' - C` | `D` | exact identity `D = H_{q^2}(delta_avg) - H_{q^2}(delta_max)`; `R` cancels |
| `b_q` (Table 4.2.1) | `b_0` | in `Cprime_sp_envelope_fits.csv` |
| `R^2` (Table 4.2.1) | `R2_0_on_M` | measured on `M_q(n) = max_k C'(n,k)` |
