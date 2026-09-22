# Paper

| File                | Content                                                                              |
| ------------------- | ------------------------------------------------------------------------------------ |
| `main.pdf`          | A Framework for Efficient Search of High-Distance Quantum Stabilizer Codes           |
| `supplementary.pdf` | Supplementary Material — retention-factor derivations (§1) and additional plots (§2) |
| `synopsis.pdf`      | Research Synopsis — summary of the project and its key findings                      |


For a table-by-table and figure-by-figure map from the paper to the files in this
repository, see [`../docs/PAPER_MAP.md`](../docs/PAPER_MAP.md).

## Structure of the paper

| Section | Content | Repository section |
|---|---|---|
| 2 | Notation and definitions: `q`-Galois strings and matrices, local non-commutation count, qudit participation, `p`-wise overlap, the stabilizer overlap graph, the search-space retention factor | — |
| 3 | Structural assumptions and their retention factors: weight spectrum, qudit participation, overlap, SOG constraints, non-contained supports, Pauli distribution, local non-commutation, automorphism groups, mixed constraints | derivations in supplementary §1 |
| 4.1–4.2 | Empirical study of the finite-size corrections `C` and `C'` | [`01-code-tables-analysis`](../01-code-tables-analysis) |
| 4.3 | Runtime analysis on `[[5,1,3]]` | [`02-runtime-benchmark`](../02-runtime-benchmark) |
| 5 | Application to the `[[14,3,5]]_2` parameters | [`03-search-14-qubit`](../03-search-14-qubit) |
