# C' on semi-perfect quantum codes

All 1909 semi-perfect codes over q in {2,3,4,5,7,8}. Non-semi-perfect codes are excluded everywhere in this document. Every number was recomputed from the cache by `analysis_Cprime_semiperfect.py`.


---

## A. Results for the main paper


### A1. The exact identity (Section 3 of the brief)

```
C' - C = [H_q2(delta_avg) - R] - [H_q2(delta_max) - R]
       =  H_q2(delta_avg) - H_q2(delta_max)
```
R cancels identically. Verified numerically over all 1909 semi-perfect codes:

**max |(C' - C) - (H_q2(delta_avg) - H_q2(delta_max))| = 0.000e+00**

The comparison is between two independently stored quantities, so this is a genuine check of the pipeline's arithmetic, not a tautology.

| sign of C' - C | count | fraction | equivalent condition |
|---|---|---|---|
| = 0 | 1289 | 67.5% | delta_avg = delta_max: every stabilizer row has the same weight |
| > 0 | 344 | 18.0% | delta_max > delta* = 1 - 1/q^2 |
| < 0 | 276 | 14.5% | both ratios on the rising branch of H_q2 |

Each of the three characterisations holds on 100% of the codes concerned. The mechanism: delta_avg <= delta_max always, and H_q2 increases on [0, delta*] then decreases on [delta*, 1]. If both ratios sit on the rising branch the larger one scores higher, so C' - C <= 0. A positive value requires the falling branch, hence requires delta_max > delta*. Equality requires the two ratios to coincide.


### A2. Envelope of C' - the requested table (Section 1 of the brief)

Empirical upper envelope M_q(n) = max{ C'(n,k) : semi-perfect code of length n }, fitted in collapse coordinates (Model A):

```
    n * M_q(n)  =  a_q * ln n + b_q          (ordinary least squares)
```
| q | a_q | b_q | R^2 | N_points |
|---|---|---|---|---|
| 2 | +0.1319 | +0.3494 | 0.1280 | 243 |
| 3 | +0.1339 | +0.0149 | 0.2065 | 93 |
| 4 | -0.0018 | +0.2649 | 0.0001 | 90 |
| 5 | -0.1023 | +0.6152 | 0.3623 | 95 |
| 7 | -0.1299 | +0.7307 | 0.3346 | 77 |
| 8 | -0.0735 | +0.5692 | 0.1233 | 81 |

N_points is the number of distinct lengths n with M_q(n) > 0. R^2 is measured on n*M_q(n), the quantity the regression actually minimises.


### A3. Verdict: the logarithmic term is NOT supported

The fitted a_q are +0.132, +0.134, -0.002, -0.102, -0.130, -0.074. They are small and **change sign across q**, which is not the behaviour of a real universal logarithmic correction. Testing the ln n term against the reduced model n*M = b (i.e. a pure inverse law M = b/n) by nested F test:

| q | R^2 of M = b/n (on M) | R^2 of (a ln n + b)/n (on M) | F | p | ln n term |
|---|---|---|---|---|---|
| 2 | 0.8353 | 0.8553 | 35.39 | 0.000 | helps |
| 3 | 0.7455 | 0.6200 | 23.68 | 0.000 | significant but WORSE on M |
| 4 | 0.7935 | 0.7992 | 0.01 | 0.937 | not significant |
| 5 | 0.7601 | 0.7509 | 52.83 | 0.000 | significant but WORSE on M |
| 7 | 0.8430 | 0.0148 | 37.72 | 0.000 | significant but WORSE on M |
| 8 | 0.8459 | 0.1925 | 11.11 | 0.001 | significant but WORSE on M |

The ln n term improves the description of M at only 1 of the six field sizes. At q = 4 it is not significant at all (p = 0.94), and at q = 7 and q = 8 it is statistically significant yet makes the fit to M dramatically **worse** (R^2 falls from 0.84 to 0.01 and from 0.85 to 0.19) - the regression is buying accuracy in n*M at the cost of the quantity of interest. Combined with the sign flip, the conclusion is that the data does not resolve a logarithmic correction.


### A4. Recommended form: the reduced inverse law

```
    C'_max(n)  ~  b_q / n            b_q = mean over lengths of n*M_q(n)
```
| q | b_q | R^2 (on M) | scatter of n*M (CV) | N_points | n range |
|---|---|---|---|---|---|
| 2 | 0.9634 | 0.8353 | 31% | 243 | 5-256 |
| 3 | 0.5161 | 0.7455 | 43% | 93 | 4-100 |
| 4 | 0.2584 | 0.7935 | 62% | 90 | 4-100 |
| 5 | 0.2326 | 0.7601 | 56% | 95 | 4-100 |
| 7 | 0.2628 | 0.8430 | 68% | 77 | 4-100 |
| 8 | 0.3046 | 0.8459 | 52% | 81 | 4-99 |

One parameter per field size, R^2 = 0.75-0.85 on M. Two caveats that belong in the text, not the footnotes:

1. **The R^2 flatters the model.** It is measured on M, whose variance is dominated by the 1/n trend that any inverse-type model reproduces. In collapse coordinates the envelope still scatters by 31-68% (coefficient of variation of n*M_q(n)), with visible branch structure - see `pub_Cprime_envelope.png`. b_q/n captures the trend of the envelope, not its spread.
2. **This is a finite-range empirical description, not an asymptotic claim.** It covers 2 <= n <= 256, and for q != 2 only n <= 100, on a sparse set of lengths (N_points = 77-243). Nothing here constrains n -> infinity.

b_q is not monotone in q (0.963, 0.516, 0.258, 0.233, 0.263, 0.305), so no clean law in q is claimed for it.


---

## B. Supplementary results


### B1. Envelope fits - both estimators and the power law

Model A fits n*M = a ln n + b (equal weight per length). Model B fits M = (a ln n + b)/n directly (small n dominate). Same model, different loss; the paper should quote one and name it.

| q | 0: b | 0: R^2 on M | A: a | A: b | A: R^2 on nM | A: R^2 on M | B: a | B: b | B: R^2 on M | pow alpha | pow beta | pow R^2 lin | pow usable | N |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 0.9634 | 0.8353 | +0.1319 | +0.3494 | 0.1280 | 0.8553 | +0.0547 | +0.6642 | 0.8981 | 0.4351 | -0.8414 | 0.8446 | yes | 243 |
| 3 | 0.5161 | 0.7455 | +0.1339 | +0.0149 | 0.2065 | 0.6200 | -0.0020 | +0.4190 | 0.8387 | 0.1508 | -0.7018 | 0.6551 | yes | 93 |
| 4 | 0.2584 | 0.7935 | -0.0018 | +0.2649 | 0.0001 | 0.7992 | -0.0470 | +0.4347 | 0.8727 | 0.1415 | -0.9101 | 0.5493 | yes | 90 |
| 5 | 0.2326 | 0.7601 | -0.1023 | +0.6152 | 0.3623 | 0.7509 | +0.0077 | +0.3168 | 0.8844 | 0.7659 | -1.3677 | 0.7466 | yes | 95 |
| 7 | 0.2628 | 0.8430 | -0.1299 | +0.7307 | 0.3346 | 0.0148 | +0.0669 | +0.1567 | 0.9127 | 8.5870 | -2.1435 | -16.2940 | NO | 77 |
| 8 | 0.3046 | 0.8459 | -0.0735 | +0.5692 | 0.1233 | 0.1925 | +0.0760 | +0.1207 | 0.9375 | 6.0120 | -1.9782 | -16.7255 | NO | 81 |

Model B reports a high R^2 (0.84-0.94) but that figure is inflated: it is measured on M, whose variance is dominated by the shared 1/n trend that every candidate model reproduces. The reduced one-parameter model M = b/n already attains R^2 = 0.75-0.85 on the same target, so Model B's apparent quality is almost entirely the 1/n factor and not evidence for the logarithm. Quoting Model B's R^2 as support for the (a ln n + b)/n form would be a mistake.


The power law is **unusable** at q = 7, 8: its linear-space R^2 is negative, i.e. the fitted curve is a worse description of M_q(n) than a horizontal line, despite a superficially healthy log-log R^2. It should not be quoted for those field sizes.


### B2. Lengths excluded from the envelope fit

| q | lengths with a semi-perfect code | usable (M>0) | dropped (M=0) |
|---|---|---|---|
| 2 | 249 | 243 | 6 |
| 3 | 99 | 93 | 6 |
| 4 | 99 | 90 | 9 |
| 5 | 99 | 95 | 4 |
| 7 | 99 | 77 | 22 |
| 8 | 99 | 81 | 18 |

A dropped length is one where every semi-perfect code has C' exactly zero. Across the whole set 623 of 1909 semi-perfect codes have C' = 0, and all of them have delta_max = 1 (full-support stabilizer rows), for which R = H_q2(1) = log_q2(q^2-1) identically.


### B3. C' - C statistics per q

| q | N | zero | >0 | <0 | min | max | mean abs (all) | median abs (non-zero) | n range |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 384 | 148 (39%) | 184 (48%) | 52 (14%) | -0.09194 | +0.02797 | 2.46e-03 | 7.44e-04 | 2-256 |
| 3 | 249 | 126 (51%) | 51 (20%) | 72 (29%) | -0.01421 | +0.01410 | 1.28e-03 | 1.16e-03 | 2-100 |
| 4 | 236 | 130 (55%) | 49 (21%) | 57 (24%) | -0.01755 | +0.00692 | 8.95e-04 | 1.32e-03 | 2-100 |
| 5 | 245 | 172 (70%) | 47 (19%) | 26 (11%) | -0.00606 | +0.00405 | 3.37e-04 | 7.01e-04 | 2-100 |
| 7 | 353 | 305 (86%) | 13 (4%) | 35 (10%) | -0.00298 | +0.00067 | 1.23e-04 | 8.14e-04 | 2-100 |
| 8 | 442 | 408 (92%) | 0 (0%) | 34 (8%) | -0.00200 | +0.00000 | 7.32e-05 | 8.85e-04 | 2-100 |

### B4. How |C' - C| varies with n (non-zero codes only)

| q | band | codes | median abs | max abs |
|---|---|---|---|---|
| 2 | 2-10 | 1 | 4.52e-03 | 4.52e-03 |
| 2 | 11-25 | 15 | 1.31e-02 | 6.54e-02 |
| 2 | 26-50 | 25 | 4.65e-03 | 9.19e-02 |
| 2 | 51-100 | 43 | 7.24e-04 | 2.15e-02 |
| 2 | 101-175 | 71 | 7.54e-04 | 6.21e-03 |
| 2 | 176-256 | 81 | 4.75e-04 | 4.93e-03 |
| 3 | 11-25 | 15 | 2.89e-03 | 1.42e-02 |
| 3 | 26-50 | 40 | 4.28e-03 | 1.41e-02 |
| 3 | 51-100 | 68 | 5.82e-04 | 3.19e-03 |
| 4 | 11-25 | 9 | 2.20e-03 | 1.75e-02 |
| 4 | 26-50 | 26 | 1.51e-03 | 6.92e-03 |
| 4 | 51-100 | 71 | 1.04e-03 | 3.75e-03 |
| 5 | 26-50 | 23 | 6.50e-04 | 6.06e-03 |
| 5 | 51-100 | 50 | 7.24e-04 | 4.05e-03 |
| 7 | 51-100 | 48 | 8.14e-04 | 2.98e-03 |
| 8 | 51-100 | 34 | 8.85e-04 | 2.00e-03 |

Spearman rank correlation of |C' - C| with n, non-zero codes:

| q | non-zero codes | Spearman rho |
|---|---|---|
| 2 | 236 | -0.444 |
| 3 | 123 | -0.695 |
| 4 | 106 | -0.193 |
| 5 | 73 | -0.276 |
| 7 | 48 | -0.783 |
| 8 | 34 | -0.638 |
| pooled | 620 | -0.433 |

### B5. Structural dependence on q

| q | frac C'-C = 0 | max abs C'-C |
|---|---|---|
| 2 | 38.5% | 9.19e-02 |
| 3 | 50.6% | 1.42e-02 |
| 4 | 55.1% | 1.75e-02 |
| 5 | 70.2% | 6.06e-03 |
| 7 | 86.4% | 2.98e-03 |
| 8 | 92.3% | 2.00e-03 |

The uniform-weight fraction rises monotonically with q while the largest attainable |C' - C| falls by more than an order of magnitude from q = 2 to q = 8. Over larger fields, C' and C coincide for most semi-perfect codes.


### B6. Files

* `analysis/plots/pub_Dprime_vs_n.png`
* `analysis/plots/pub_Dprime_vs_n.pdf`
* `analysis/plots/pub_Cprime_envelope.png`
* `analysis/Cprime_sp_envelope_fits.csv`
* `analysis/Cprime_sp_envelope_points.csv`
* `analysis/Cprime_sp_D_points.csv`
* `analysis/Cprime_sp_D_summary.csv`
* `analysis/Cprime_sp_D_vs_n_bands.csv`
