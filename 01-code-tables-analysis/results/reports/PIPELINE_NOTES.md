# QECC Pipeline — codetables.de scraper, analyser and visualiser

`qecc_pipeline.py` automatically retrieves every available quantum stabilizer
code from Markus Grassl's [code tables](https://www.codetables.de/QECC/),
computes the requested weight-spectrum / entropy quantities, classifies each
code against the published upper bounds, and produces one Excel workbook and
nine scatter plots for each of the six field sizes *q* = 2, 3, 4, 5, 7, 8.

No parity-check matrix is ever typed in by hand, and nothing is invented: if a
page has no usable stabilizer matrix the parameter set is skipped and the
reason is recorded in the skip log.

---

## 1. Quick start

```bash
python qecc_pipeline.py
```

That processes `q = 2, 3, 4, 5, 7, 8` over the full published ranges. The
first run downloads ~58,900 pages (about 3 hours at the default 6
requests/second); every later run reads them straight from the local cache in
about a minute.

Useful options:

```bash
python qecc_pipeline.py --q 2 --n-max 40        # a smaller slice
python qecc_pipeline.py --skip-fetch            # work offline from the cache
python qecc_pipeline.py --skip-fetch --force-reparse   # recompute all metrics
python qecc_pipeline.py --rps 3 --workers 3     # gentler on the server
```

| Option | Meaning | Default |
|---|---|---|
| `--q` | field sizes to process (2, 3, 4, 5, 7, 8) | all six |
| `--n-min` / `--n-max` | restrict the length range | 1 … table maximum |
| `--rps` | target requests per second | 6 |
| `--workers` | concurrent HTTP workers | 6 |
| `--max-attempts` | retries per URL before giving up | 5 |
| `--weight-mode` | `symplectic` or `raw` (see §4) | `symplectic` |
| `--x-sum-mode` | `distinct` or `rows` (see §4) | `distinct` |
| `--avg-denominator` | `rows` or `n-k` (see §4) | `rows` |
| `--skip-fetch` | never touch the network | off |
| `--force-reparse` | discard stored metrics and recompute | off |
| `--out-dir` | where cache, Excel, plots and logs go | script directory |

---

## 2. Data source and ranges

Pages are requested from

```
https://www.codetables.de/QECC/QECC.php?q=<q^2>&n=<n>&k=<k>
```

so `q=4` is the binary (qubit) table, `q=9` is `q = 3`, and so on. The valid
ranges were determined by probing the server (anything outside returns
`wrong input`):

| q | `q^2` in URL | n | k |
|---|---|---|---|
| 2 | 4 | 1 … 256 | 0 … n |
| 3 | 9 | 1 … 100 | 0 … n |
| 4 | 16 | 1 … 100 | 0 … n |
| 5 | 25 | 1 … 100 | 0 … n |
| 7 | 49 | 1 … 100 | 0 … n |
| 8 | 64 | 1 … 100 | 0 … n |

That is **58,902 parameter sets** in total. `q^2 = 81` and `q^2 = 121` return
an empty document rather than a table, so these six field sizes are the
complete set the server offers.

---

## 3. How a page is read

A typical page contains a bounds table and one construction block:

```
Bounds on [[5,1]]_2      lower bound: 3    upper bound: 3

Construction of a [[5,1,3]] quantum code:
[1]:  [[5, 1, 3]] quantum code over GF(2^2)
     Construction from a stored generator matrix

    stabilizer matrix:

      [1 0 1 0 1|0 0 1 1 0]
      ...
```

* **Bounds** — `lower bound` and `upper bound` are stored for *every* page,
  including pages that are later skipped. The complete upper-bound grid is
  what makes the classification in §5 exact.
* **Multiple constructions** — the body is split at each
  `Construction of a [[n,k,d]] quantum code:` header, so a page carrying more
  than one parity-check matrix yields one candidate per header. Candidates are
  tried in **descending order of d** and the first that validates is used;
  lower-distance matrices are ignored. The `[1]:`, `[2]:` … lines are
  *intermediate* codes in a construction chain, not separate matrices, and are
  never mistaken for one.
* **Validation** — a candidate is accepted only if it has exactly the expected
  number of rows and exactly `n` symbols on each side of the `|`. Anything else
  is logged and the next candidate is tried.
* **Row count depends on the field.** The stabilizer of an `[[n,k]]_q` code is
  an *additive* subgroup of GF(q²)ⁿ of size `q^(n-k)`, and codetables.de prints
  a basis over the **prime** field GF(p). For `q = p^m` that is `m·(n − k)`
  rows with entries in GF(q):

  | q | m | rows printed | entries |
  |---|---|---|---|
  | 2 | 1 | `n − k` | `0 1` |
  | 3 | 1 | `n − k` | `0 1 2` |
  | 4 | 2 | `2(n − k)` | `0 1 a a^2` |
  | 5 | 1 | `n − k` | `0 … 4` |
  | 7 | 1 | `n − k` | `0 … 6` |
  | 8 | 3 | `3(n − k)` | `0 1 w … w^6` |

  Measured across every published matrix the ratio `rows/(n−k)` is exactly
  1.000 for the prime q (2, 3, 5, 7), exactly 2.000 for q = 4 and exactly
  3.000 for q = 8, with only a handful of outliers in total (see below).
  Note the primitive element is printed as `a` for GF(4) but `w` for GF(8);
  the parser never needs to know, because it only ever compares a symbol
  against the literal `0`.

### Skip reasons

| Reason | Meaning |
|---|---|
| `no parity-check matrix` | `Missing Construction`, a Gilbert–Varshamov existence entry, or a construction with no matrix printed |
| `no stabilizer rows` | `k = n`, so the stabilizer matrix is empty and `n − k = 0` would divide by zero |
| `parity-check matrix belongs to a different code` | the row count implies a different `k` than the page claims |
| `malformed parity-check matrix` | a row has the wrong number of symbols |
| `unparseable page` | no bounds could be recovered |
| `page unavailable` | the URL never downloaded successfully |

Every skipped parameter set is written to `logs/skipped_codes.csv` with its
reason and detail.

---

## 4. Quantities computed

For the selected matrix of a code `[[n,k,d]]_q`, with `q2 = q²`:

| Symbol | Definition |
|---|---|
| `W` | weight spectrum `{weight: frequency}` over the stabilizer rows |
| `Avg Weight` | `Σ Wᵢ·fᵢ / (n − k)` |
| `Max Weight` | `max Wᵢ` |
| `X` | `Σᵢ C(n, Wᵢ)·(q²−1)^{Wᵢ}` |
| `R` | `log_q(X) / (2n)` |
| `δ_avg`, `δ_max` | `Avg/n`, `Max/n` |
| `H_{q²}(x)` | `x·log_{q²}(q²−1) − x·log_{q²}(x) − (1−x)·log_{q²}(1−x)` |
| `C` | `H_{q²}(δ_max) − R` |
| `C'` | `H_{q²}(δ_avg) − R` |

`X` is computed in exact integer arithmetic and its logarithm is taken with a
big-integer-safe routine, so lengths up to n = 256 do not overflow.
`H_{q²}` uses the correct limits `H(0) = 0` and `H(1) = log_{q²}(q²−1)`.

### Three interpretation choices

All three are exposed as switches rather than baked in.

**Row weight (`--weight-mode`, default `symplectic`).** A stabilizer row is a
length-`n` vector over GF(q²) printed as a symplectic pair `[X-part | Z-part]`
of `n + n` symbols. The weight counted is the number of coordinates
`j ∈ 1…n` with `(X_j, Z_j) ≠ (0,0)` — the number of non-zero entries of the
length-`n` vector. This is what keeps `W ≤ n`, so `δ ∈ [0,1]` and
`H_{q²}(δ)` is defined; counting non-zero symbols across all `2n` printed
positions instead can give `W > n` and an undefined entropy. (For `[[5,1,3]]`
the symplectic weights are `4,4,4,4` whereas the raw counts are `5,5,4,4`.)
Use `--weight-mode raw` to count all `2n` symbols anyway.

**Summation range of `X` (`--x-sum-mode`, default `distinct`).** The `Avg
Weight` formula names the frequencies `fᵢ` explicitly while the `X` formula
does not, so `X` is summed over the **distinct** weights present in the
spectrum. Use `--x-sum-mode rows` to weight each term by its frequency.

**Denominator of `Avg Weight` (`--avg-denominator`, default `rows`).** The
specified denominator is `n − k`, which is exactly the number of stabilizer
rows for q = 2, 3, 5 — so for three of the four field sizes the two readings
are **identical**. They differ only for q = 4, where the matrix has `2(n − k)`
rows: dividing by `n − k` there would report twice the mean row weight, push
`δ_avg` above 1, and leave `H_{q²}(δ_avg)` — and hence `C'` — undefined for
most codes. The default therefore divides by the actual row count, which is
the mean weight of a stabilizer generator in every case. Use
`--avg-denominator n-k` to force the literal formula.

---

## 5. Classification

A code is a **Semi-Perfect Quantum Code** when, according to the *upper bounds*
published in the table for the same `q`, no code with strictly better
parameters **can exist** — not merely that none has been found:

1. `d` cannot be improved: `d = ub(n, k)`;
2. `k` cannot be improved: `ub(n, k') < d` for every `k' > k`;
3. `n` cannot be reduced: `ub(n', k) < d` for every `n' < n`.

Anything else is a **Normal Code**. All three tests run against the full
upper-bound grid for that `q`, which is why the fetch phase downloads every
`(n, k)` page even when the code itself will be skipped for lacking a matrix.
If any upper bound is missing the run logs a warning naming the number of
affected parameter sets.

Sanity checks on the binary table: `[[5,1,3]]`, `[[8,3,3]]`, `[[11,1,5]]`,
`[[12,0,6]]` are semi-perfect, while `[[7,1,3]]` and `[[9,1,3]]` are normal
because the shorter `[[5,1,3]]` can exist.

---

## 6. Outputs

```
QECC_codes_q2.xlsx   QECC_codes_q3.xlsx   QECC_codes_q4.xlsx
QECC_codes_q5.xlsx   QECC_codes_q7.xlsx   QECC_codes_q8.xlsx
all_codes_combined.csv
run_summary.txt / run_summary.csv
codetables_cache.sqlite
plots/q<q>/01…09_*.png
logs/run.log
logs/skipped_codes.csv
```

Each workbook has a single `Codes` sheet with AutoFilter enabled over the
header row and the first row frozen, so **Type of Code** (and every other
column) is interactively filterable. Columns, in order:

`q`, `n`, `k`, `Code Parameters`, `Type of Code`, `k/n`, `Weight Spectrum`,
`Average Weight`, `Max Weight`, `R`, `C`, `C'`.

The nine plots per `q` are all **discrete scatter plots** — no smoothing, no
interpolation, no trend lines. Every marker carries a thin black outline
(0.3 pt), and each panel gets a dotted reference line: `C = 0` or `C' = 0` on
the panels whose vertical axis is `C`/`C'`, and `H_{q²}(δ) = R` on panels 7
and 8 — which marks the same locus, since `C = 0` exactly when
`H_{q²}(δ) = R`.

| # | Plot | Notes |
|---|---|---|
| 1 | `C` vs `n`, all codes | colour gradient on `k/n` |
| 2 | `C` vs `k`, all codes | |
| 3 | `C` vs `k/n`, all codes | |
| 4 | `C'` vs `n`, semi-perfect only | colour gradient on `k/n`; enlarged markers |
| 5 | `C'` vs `k`, semi-perfect only | enlarged markers |
| 6 | `C'` vs `k/n`, semi-perfect only | enlarged markers |
| 7 | `H_{q²}(δ_max)` vs `R` | dotted `H_{q²}(δ) = R` reference line |
| 8 | `H_{q²}(δ_avg)` vs `R` | dotted `H_{q²}(δ) = R` reference line |
| 9 | `C'` vs `n`, all codes | semi-perfect vibrant, normal muted, with legend, no gradient |

---

## 7. Robustness, caching and resumability

* **One SQLite file** (`codetables_cache.sqlite`) holds the gzipped HTML of
  every page plus the parsed bounds, code metrics and skip log. A single file
  avoids scattering ~58,900 small files across a synced folder.
* **No URL is fetched twice.** Before each run the set of already-cached
  `(n, k)` keys is read back and subtracted from the work list, so an
  interrupted run resumes exactly where it stopped.
* **Rate limiting** is an evenly-spaced global limiter shared by all workers,
  with automatic back-off on HTTP 429/503 and a gradual return to the target
  rate.
* **Retries** use exponential back-off with jitter (5 attempts by default) for
  connection errors, timeouts and 5xx responses. A URL that still fails is
  recorded as an error and reported, never silently dropped.
* **Periodic saves** — pages are committed in batches of 240 during fetching
  and metrics in batches of 500 during parsing, so an interrupted run loses at
  most one batch.
* **Per-code error isolation** — a malformed or unparseable entry is logged and
  skipped; it never aborts the run.
* **Offline re-analysis** — once the cache is populated, `--skip-fetch
  --force-reparse` recomputes every quantity, or re-runs with a different
  `--weight-mode` / `--x-sum-mode`, without touching the network.

## 8. Requirements

Python 3.9+ with `requests`, `beautifulsoup4`-free parsing (regex only),
`pandas`, `openpyxl`, `matplotlib`, `seaborn`, `tqdm`, `numpy`.

```bash
pip install requests pandas openpyxl matplotlib seaborn tqdm numpy
```
