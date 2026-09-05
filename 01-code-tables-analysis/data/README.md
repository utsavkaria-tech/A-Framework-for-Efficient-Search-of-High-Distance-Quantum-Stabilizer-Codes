# Data dictionary

Source: **Markus Grassl's code tables**, <https://www.codetables.de/QECC/>. Every row
corresponds to one `(q, n, k)` parameter set that published a usable stabilizer
(parity-check) matrix. Sets without one are not here — they are logged with a reason in
[`../logs/skipped_codes.csv`](../logs/skipped_codes.csv).

## `all_codes_combined.csv` — 34,775 rows, all field sizes

| Column | Meaning |
|---|---|
| `q` | Field size: 2, 3, 4, 5, 7 or 8 |
| `n` | Block length (physical qudits) |
| `k` | Encoded qudits |
| `Code Parameters` | The parameter string, e.g. `[[5,1,3]]_2` |
| `Type of Code` | `Semi-Perfect Quantum Code` or `Normal Code` |
| `k/n` | Code rate |
| `Weight Spectrum` | `{weight: frequency}` over the stabilizer generators |
| `Average Weight` | Mean stabilizer row weight, `w_avg` |
| `Max Weight` | Maximum stabilizer row weight, `w_max` |
| `R` | Search-space retention factor, `log_q(X) / (2n)` |
| `C` | Finite-size correction, `H_{q^2}(delta_max) - R`, with `delta_max = w_max / n` |
| `C'` | Average-weight variant, `H_{q^2}(delta_avg) - R`, with `delta_avg = w_avg / n` |

## `QECC_codes_q{2,3,4,5,7,8}.xlsx`

The same data split by field size, with AutoFilter enabled. One workbook per `q`.

## Conventions

- **Row weight is symplectic**: the number of coordinates `j` in `1..n` with
  `(X_j, Z_j) != (0,0)`. This keeps `delta = W_i / n` in `[0,1]` so the `q^2`-ary entropy
  `H_{q^2}(delta)` is defined.
- `X = sum_{w in W} binom(n, w) (q^2 - 1)^w`, summed over the **distinct** weights in the
  spectrum. Frequencies enter only through `Average Weight`.
- `H_{q^2}(x) = x log_{q^2}(q^2-1) - x log_{q^2}(x) - (1-x) log_{q^2}(1-x)`.
- **Semi-perfect** classification uses the *complete* published upper-bound grid for that
  `q`: a code is semi-perfect when no `[[n',k,d]]_q` exists with `n' < n`, no `[[n,k,d']]_q`
  with `d' > d`, and no `[[n,k',d]]_q` with `k' > k`.

A richer per-code table for the semi-perfect subset — including `delta_avg`, `delta_max`,
`H_gap`, `past_peak` and `uniform` flags — is in
[`../results/fits/sp_codes.csv`](../results/fits/sp_codes.csv).
