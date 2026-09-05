# A Framework for Efficient Search of High-Distance Quantum Stabilizer Codes

Code, data and results accompanying the paper **"A Framework for Efficient Search of
High-Distance Quantum Stabilizer Codes"** by **Utsav Karia** (Independent Researcher),
together with its supplementary material.

---

## What this is about

Systematic searches for high-distance quantum stabilizer codes become computationally
prohibitive as the code length grows. Known code constructions are far from arbitrary —
they show strong structural regularities — but there was no general framework that says
*how much* such a regularity actually shrinks the search space.

This work supplies one. Stabilizer generators are represented as **q-Galois matrices**, and
for every structural constraint we define a **search-space retention factor**

```
    R = lim (n -> inf)  log_{q^2} S_R / log_{q^2} S_O
```

the asymptotic fraction of the *logarithmic* search space that survives the constraint
(`S_O` = unrestricted space, `S_R` = reduced space). Retention factors are derived for
weight-spectrum, qudit-participation, overlap, support, Pauli-distribution, local
non-commutation, graph-theoretic and symmetry constraints, and for combinations of them.

Finite-length codes deviate from the asymptotic value by a **finite-size correction factor**
`C` (and its average-weight variant `C'`). This repository contains the empirical study of
`C` and `C'` over the published quantum-code tables, the runtime experiments that validate
the predicted savings, and the full 14-qubit code search that applies the framework.

---

## Headline results

| Result | Value | Where |
|---|---|---|
| Code parameter sets examined | 58,902 across `q = 2,3,4,5,7,8` | [`01-code-tables-analysis`](01-code-tables-analysis) |
| Sets with usable stabilizer matrix data | 34,775 | same |
| Semi-perfect quantum codes analysed | 1,909 | same |
| Finite-size correction scaling | `C = O(log n / n)`, `C' = O(1/n)` | same |
| Runtime speedups for 4 unresolved semi-perfect codes | `10^8.3608` to `10^193.9298` | same |
| Measured speedup, `[[5,1,3]]_2` search | 40.4x (universal), 38.3x (weight-4), 6.10e7x (cyclic) | [`02-runtime-benchmark`](02-runtime-benchmark) |
| Best-known `[[14,3]]_2` code | detects **91,765 / 91,770** weight-<=4 errors (99.99455%) | [`03-search-14-qubit`](03-search-14-qubit) |
| Residual undetected errors | **5** (weight profile 0/0/2/3) | same |

At physical error probability `p = 0.01` the constructed `(14,3)` code reaches
**99.999999933%** and **99.9733%** of the theoretical detection and correction performance
of the optimal `[[14,3,5]]_2` code.

---

## Repository layout

```
.
├── paper/                        The paper and supplementary material (PDF)
│
├── 01-code-tables-analysis/      Sections 4.1 and 4.2 - empirical study of C and C'
│   ├── src/                      Scraper/pipeline + three analysis scripts + ML estimator
│   ├── data/                     Per-q workbooks and the combined 34,775-code CSV
│   ├── results/fits/             All fitted coefficients as CSV (the paper's tables)
│   ├── results/workbooks/        The same analyses as Excel workbooks
│   ├── results/reports/          Long-form generated reports (Markdown)
│   ├── figures/                  Every scatter plot and publication figure
│   └── logs/                     Full run log and the skipped-code log
│
├── 02-runtime-benchmark/         Section 4.3 - runtime validation on [[5,1,3]]_2
│   ├── src/                      Four independent benchmark engines
│   └── RESULTS_SUMMARY.md        Measured runtimes, methodology, validation invariants
│
├── 03-search-14-qubit/           Section 5 - applying the framework to (14,3,2)
│   ├── src/                      C++17/20 search programs
│   ├── build/                    Makefile, GCC script, MSVC batch files
│   ├── docs/                     Per-program mathematical derivations (README_MATH*)
│   └── results/                  Stage 1, Stage 2, 45-orbit sweep, final code
│
└── docs/
    ├── PAPER_MAP.md              Every table and figure in the paper -> the file it came from
    └── REPRODUCING.md            Step-by-step reproduction instructions
```

---

## Quick start

```bash
git clone https://github.com/<your-username>/high-distance-stabilizer-code-search.git
cd high-distance-stabilizer-code-search
pip install -r requirements.txt
```

**Reproduce the correction-factor analysis** (needs the SQLite cache; see below):

```bash
cd 01-code-tables-analysis/src
python qecc_pipeline.py                     # fetch + parse + classify + export
python analysis_scaling.py                  # Table 4.1.4, Figures 4.1.2, 2.1, 2.3
python analysis_semiperfect.py              # Section 4.2, Figures 2.2, 2.5.x
python analysis_Cprime_semiperfect.py       # Table 4.2.1, Figure 2.4
```

**Reproduce a runtime row** (Table 4.3):

```bash
cd 02-runtime-benchmark/src
python stabilizer_search_benchmark_universal.py --n 5 --k 1 --d 3 \
    --constraint universal --backend numpy --target span --work 5 --repetitions 2
```

**Reproduce the 14-qubit search** (Section 5):

```bash
cd 03-search-14-qubit
make -f build/Makefile              # or: sh build/build_gcc.sh   (MSVC: build/*.bat)
./stab14_2stage --mode hybrid       # Stage 1 optimum + 45-orbit Stage 2 sweep
./stab14_refine --time-limit 2400   # targeted refinement -> the final 5-residual code
```

Full instructions, expected runtimes and expected outputs are in
[`docs/REPRODUCING.md`](docs/REPRODUCING.md).

---

## Data provenance and the one file that is not here

All quantum-code parameters and stabilizer matrices are taken from **Markus Grassl's
code tables**, <https://www.codetables.de/QECC/>. The pipeline downloads each
`(q, n, k)` page once into a local SQLite cache and does everything else offline.

That cache, `codetables_cache.sqlite`, is **167 MB and is deliberately not committed** —
it exceeds GitHub's 100 MB per-file limit. It is a pure download cache and contains no
original results. Regenerate it by running the fetch phase:

```bash
cd 01-code-tables-analysis/src
python qecc_pipeline.py --q 2 3 4 5 7 8 --n-min 1 --n-max 256 --rps 6
```

The fetch is rate-limited (6 requests/second by default), resumable, and never downloads
the same URL twice. Everything downstream of it — the workbooks, CSVs, fits and figures in
this repository — is committed, so **no result here requires re-downloading anything.**

---

## Requirements

| Component | Needs |
|---|---|
| `01-code-tables-analysis` | Python 3.10+, `numpy pandas scipy matplotlib openpyxl requests beautifulsoup4 seaborn scikit-learn` |
| `02-runtime-benchmark` | Python 3.10+ (3.13 used), `numpy` (2.5.2 used). Pure-Python fallback backend available |
| `03-search-14-qubit` | A C++20 compiler (GCC/Clang with `-O3 -march=native`, or MSVC with `/O2 /arch:AVX2`), pthreads |

`pip install -r requirements.txt` covers both Python sections.

---

## A note on what is proven and what is not

The repository is careful to distinguish the two, and so is the paper:

- **Proven optimal.** The Stage-1 result (three weight-8 generators detecting 80,584 of
  91,770 errors) is a completed branch-and-bound search with admissible bounds — an
  optimality proof, not a heuristic best.
- **Best known.** Every Stage-2 result, and therefore the final 11-generator code, is
  *best known*. The incumbent is certified optimal against every 1-generator and every
  2-generator replacement, exhaustively over all 2,037,794 weight-10 classes, but no
  impossibility claim is made beyond that.

Each result file states its own status in its header. See
[`03-search-14-qubit/README.md`](03-search-14-qubit/README.md).

---

## Citation

```bibtex
@misc{karia_stabilizer_search,
  author = {Karia, Utsav},
  title  = {A Framework for Efficient Search of High-Distance Quantum Stabilizer Codes},
  year   = {2026},
  note   = {Code and data: https://github.com/<your-username>/high-distance-stabilizer-code-search}
}
```

Please also cite the source of the underlying code parameters:
M. Grassl, *Bounds on the minimum distance of linear codes and quantum codes*,
<http://www.codetables.de>.

## License

MIT — see [LICENSE](LICENSE). The derived datasets are released under the same terms;
the upstream code-table data remains subject to codetables.de's own terms of use.
