#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stabilizer_search_benchmark_weight4.py
======================================

EXPERIMENT (C): CONSTANT-ROW-WEIGHT RUNTIME BENCHMARK for [[5,1,3]]_2.

Built from stabilizer_search_benchmark_universal.py so that the enumeration,
matrix representation, processing workload (--work), timing methodology,
target recognition and reporting are identical.  The only difference is which
structural constraint is imposed.

THE CONSTRAINT (the ONLY thing imposed)
----------------------------------------
    wt(row_i) in W    for every row i

supplied explicitly with --weights.  For this experiment W = {4}, i.e. the
row-weight vector must be exactly [4,4,4,4].  Row weight is the Pauli weight
by default (number of qubits on which the row is not the identity); use
--weight-mode binary for the plain popcount of the 2n-bit row.

Explicitly NOT imposed: rank(H) = n-k, H_X H_Z^T = H_Z H_X^T (mod 2),
distance = 3, any [[5,1,3]] validity test, cyclicity, automorphism
conditions, canonical forms, or duplicate-generating-set removal.

RETENTION
---------
For [[5,1,3]]_2 with W={4} there are 405 row values of Pauli weight 4
(5 choose 4 positions times 3 non-identity Paulis to the power 4 = 5*81),
so S_R = 405^4 = 26,904,200,625 and S_R/S_0 = 0.0244692... The search never
uses this figure; --dry-run prints it only so a run can be validated.

BOUNDED SAMPLES
---------------
--r0-start / --r0-count restrict the enumeration to a contiguous
deterministic region of row 0, and the report extrapolates the full-search
runtime by the exact ratio of candidate counts.  Nothing is randomised.
NOTE: acceptance under this constraint depends on the weight of row 0, which
is NOT uniformly distributed over contiguous ranges of the row alphabet, so
use several well-separated regions and read the spread.

WHAT THIS PROGRAM IS
--------------------
It measures the wall-clock cost of *traversing* a candidate matrix space, so
that an experimentally observed runtime ratio

        T_R / T_0

can be compared against the theoretical search-space retention prediction

        S_R / S_0 = S_0^(R - 1),      R = log(S_R) / log(S_0).

It is NOT a code-search algorithm, and it deliberately does not try to be
clever.  There is no sampling, no heuristic, no counting formula, no early
exit, no symmetry/canonical-form reduction and no duplicate elimination.
Every candidate matrix in the relevant space is genuinely visited and
processed.

CANDIDATE SPACE
---------------
For [[n, k, d]]_2 the parity-check matrix has

        rows = n - k
        cols = 2 * n

and each entry is an independent bit.  The unconstrained space therefore has

        S_0 = 2 ** (rows * cols) = (2 ** cols) ** rows

elements.  For [[5,1,3]]_2 that is a 4 x 10 binary matrix and

        S_0 = 2 ** 40 = 1,099,511,627,776.

NO validity conditions of any kind are imposed in the baseline: no
commutation / symplectic-orthogonality check, no rank or independence check,
no distance check, no logical-qubit check, no syndrome check, no canonical
form, no duplicate removal.  The baseline is the raw 2^(rows*cols) space.

ROW ENCODING (fixed, explicit)
------------------------------
A row is stored as one integer.  Column order is

        [ x_0 x_1 ... x_{n-1} | z_0 z_1 ... z_{n-1} ]

and column c occupies bit (cols - 1 - c) of the integer, so the integer
printed in binary reads left-to-right exactly like the printed matrix row.
Enumeration order over full matrices is plain lexicographic on
(row_0, row_1, ..., row_{rows-1}) with row_0 varying slowest.

SECTION MAP (see section banners further down)
----------------------------------------------
  [S1] Code specification and row/Pauli encoding helpers
  [S2] GF(2) linear algebra helpers (rank, span, SYMPLECTIC FORM)
  [S3] TARGET-CODE IDENTIFICATION
  [S4] STRUCTURAL CONSTRAINT EVALUATION (modular, one class per constraint)
       -- includes UniversalStabilizerConstraint: rank + commutation
  [S5] Search plan construction (all pre-computation, OUTSIDE the timer)
  [S6] CANDIDATE ENUMERATION + per-candidate processing  <-- the timed work
       -- includes search_numpy_universal(), the constrained engine
  [S7] RUNTIME MEASUREMENT (perf_counter wrapper, repetitions, statistics)
  [S8] Reporting
  [S9] Command-line interface / self-test

PER-CANDIDATE PROCESSING WORKLOAD (identical to experiment (A))
---------------------------------------------------------------
For every candidate that reaches the processing stage the engine computes the
binary weights of all of its rows, forms their sum (a matrix-level quantity),
and folds that into a running checksum.  This is the identical operation used
by the unconstrained baseline -- same arithmetic, same accumulator, same
order -- which is what makes T_universal / T_unconstrained meaningful.

Checksum semantics differ between the two experiments BY DESIGN, because
rejected candidates are skipped rather than processed:

  (A) unconstrained: every candidate is processed, so the checksum equals
      rows * cols/2 * S_0 = 4 * 5 * 2**40 = 21,990,232,555,520.

  (B) this program: only ACCEPTED candidates are processed, so the checksum
      covers the accepted subspace only.  It is not comparable to (A).

What IS comparable, and what verifies that the enumeration itself is
unchanged, is "Total candidate matrices considered", which must still be
exactly 2^40, together with the accepted count, which must match the
combinatorial reference above.

TARGET-CODE IDENTIFICATION
--------------------------
The default target is the five-qubit perfect code [[5,1,3]] with the cyclic
stabilizer generators

        XZZXI, IXZZX, XIXZZ, ZXIXZ

A candidate is reported as the target when (mode "span") its four rows lie in
the GF(2) row space of those generators AND have full rank, i.e. the matrix
generates exactly the same stabilizer group.  Because an ordered basis of a
4-dimensional GF(2) space can be chosen in

        (2^4-1)(2^4-2)(2^4-4)(2^4-8) = 15*14*12*8 = 20160

ways, a complete unconstrained traversal must report exactly 20160
occurrences.  That is a useful correctness check on the run.

Mode "exact" instead matches only the literal generator matrix above (exactly
1 occurrence in a full traversal).

The target test NEVER alters the enumeration and NEVER terminates the search.
It is gated by an O(1) table lookup so its cost is negligible.

BACKENDS
--------
Both backends enumerate the identical space in the identical order, run the
identical per-candidate processing, and produce the identical checksum.

  numpy  (default)  Traversal is blocked over the last two rows: the outer
                    rows are iterated in Python and the (a x b) grid of the
                    final two rows is processed as one vectorised block.
                    Every candidate is still individually weighted, summed
                    and accumulated -- nothing is hoisted out of the traversal
                    except the per-row-value popcount lookup table, which the
                    pure backend uses too.

  pure              Plain nested enumeration via itertools.product with a
                    Python-level loop body.  Same semantics, far slower;
                    useful for small n and for cross-validating the numpy
                    backend.

Use the SAME backend for the baseline and for every constrained run.

SCALE WARNING
-------------
The enumeration still covers S_0 = 2**40 ~ 1.1e12 candidates, but only about
1.43% of them reach the processing stage, so this run is expected to be much
shorter than the unconstrained baseline.  Because it is short enough to
repeat, use --repetitions 3 (or more) so mean / median / stdev are available.
Use --dry-run first to see the exact space size, block count and the
combinatorial reference for the accepted count.

Running this same file with --constraint none reproduces experiment (A)
exactly: that path is byte-identical to the original program and is the
fair-comparison control if you want both numbers from one binary.

DO-NOT list honoured by this file: no plots, no Excel, no per-matrix output,
no extra analysis, no early termination, no sampling, no randomness.
"""

from __future__ import annotations

import argparse
import statistics
import sys
import time
from dataclasses import dataclass, field
from itertools import product
from typing import Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

try:
    import numpy as _np
except ImportError:  # numpy is optional; the "pure" backend needs nothing
    _np = None


# =====================================================================
# [S1] Code specification and row / Pauli encoding helpers
# =====================================================================

@dataclass(frozen=True)
class CodeSpec:
    """Parameters of the code whose parity-check matrix space is searched."""

    n: int
    k: int
    d: Optional[int] = None

    @property
    def rows(self) -> int:
        """Number of stabilizer generators = number of matrix rows = n - k."""
        return self.n - self.k

    @property
    def cols(self) -> int:
        """Number of matrix columns = 2n (X block followed by Z block)."""
        return 2 * self.n

    @property
    def row_alphabet(self) -> int:
        """Number of distinct binary rows = 2 ** (2n)."""
        return 1 << self.cols

    @property
    def space_size(self) -> int:
        """S_0 = unconstrained candidate-matrix count = 2 ** (rows * cols)."""
        return self.row_alphabet ** self.rows

    @property
    def label(self) -> str:
        d = self.d if self.d is not None else "?"
        return f"[[{self.n},{self.k},{d}]]_2"

    def validate(self) -> None:
        if self.n <= 0:
            raise ValueError("n must be positive")
        if not (0 <= self.k < self.n):
            raise ValueError("require 0 <= k < n")
        if self.rows <= 0:
            raise ValueError("n - k must be positive")


def pauli_string_to_row(pauli: str, n: int) -> int:
    """Encode an n-qubit Pauli string ('XZZXI') as one binary symplectic row."""
    if len(pauli) != n:
        raise ValueError(f"Pauli string {pauli!r} does not have length n={n}")
    cols = 2 * n
    row = 0
    for j, p in enumerate(pauli.upper()):
        if p not in "IXYZ":
            raise ValueError(f"bad Pauli letter {p!r} in {pauli!r}")
        if p in ("X", "Y"):
            row |= 1 << (cols - 1 - j)          # X block, qubit j
        if p in ("Z", "Y"):
            row |= 1 << (cols - 1 - (n + j))    # Z block, qubit j
    return row


def combination_subsets(rows: int) -> Tuple[Tuple[int, ...], ...]:
    """
    Row subsets used by the extra per-candidate processing units (--work).

    Only subsets containing BOTH vectorised (inner) rows are used.  A subset
    missing an inner row would make the combination depend on fewer than all
    of the enumerated rows, and its weight could then be computed once per
    row value instead of once per candidate -- which would stop being a
    per-candidate workload.  Requiring both inner rows guarantees every unit
    is genuine per-candidate matrix processing in every engine.

    For rows = 4 this returns the 4 subsets {2,3}, {0,2,3}, {1,2,3}, {0,1,2,3}.

    Identical to the function of the same name in the unconstrained program,
    so the two benchmarks apply the same workload.
    """
    if rows < 2:
        raise ValueError("need at least 2 rows")
    inner = (rows - 2, rows - 1)
    out: List[Tuple[int, ...]] = []
    for c in range(1 << rows):
        if all((c >> i) & 1 for i in inner):
            out.append(tuple(i for i in range(rows) if (c >> i) & 1))
    return tuple(out)


def build_popcount_table(alphabet: int) -> List[int]:
    """Binary weight (number of 1 bits) of every possible row value."""
    return [v.bit_count() for v in range(alphabet)]


def build_pauli_weight_table(alphabet: int, n: int) -> List[int]:
    """
    Pauli weight of every possible row value: the number of qubits on which
    the row is not the identity, i.e. popcount(x_part OR z_part).
    """
    mask = (1 << n) - 1
    table = []
    for v in range(alphabet):
        x = (v >> n) & mask
        z = v & mask
        table.append((x | z).bit_count())
    return table


def cyclic_shift_row(row: int, n: int, shift: int) -> int:
    """
    Cyclically shift a row by `shift` qubit positions, applying the same
    rotation independently to the X block and to the Z block.
    """
    mask = (1 << n) - 1
    x = (row >> n) & mask
    z = row & mask
    s = shift % n
    x = ((x >> s) | (x << (n - s))) & mask
    z = ((z >> s) | (z << (n - s))) & mask
    return (x << n) | z


def row_to_pauli_string(row: int, n: int) -> str:
    """Inverse of pauli_string_to_row (used only by --self-test)."""
    mask = (1 << n) - 1
    x = (row >> n) & mask
    z = row & mask
    out = []
    for j in range(n):
        xb = (x >> (n - 1 - j)) & 1
        zb = (z >> (n - 1 - j)) & 1
        out.append("IXZY"[xb + 2 * zb] if (xb, zb) != (1, 1) else "Y")
    return "".join(out)


# =====================================================================
# [S2] GF(2) linear algebra helpers
#      Used ONLY by target-code identification.  They are never applied as
#      a search constraint and never influence the enumeration.
# =====================================================================

def gf2_rank(vectors: Sequence[int]) -> int:
    """Rank over GF(2) of integer-encoded row vectors."""
    pivots: Dict[int, int] = {}
    rank = 0
    for v in vectors:
        cur = v
        while cur:
            hb = cur.bit_length() - 1
            piv = pivots.get(hb)
            if piv is None:
                pivots[hb] = cur
                rank += 1
                break
            cur ^= piv
    return rank


def gf2_span(generators: Sequence[int]) -> frozenset:
    """All GF(2) linear combinations of the given integer-encoded vectors."""
    span = {0}
    for g in generators:
        span |= {s ^ g for s in span}
    return frozenset(span)


# --- symplectic form: the GF(2) core of the commutation condition ---------
#
# For rows a = (a_X | a_Z) and b = (b_X | b_Z) the (i,j) entry of
# H_X H_Z^T - H_Z H_X^T is
#
#       <a, b> = a_X . b_Z + a_Z . b_X   (mod 2)
#
# so "H_X H_Z^T = H_Z H_X^T (mod 2)" is exactly "<r_i, r_j> = 0 for every
# pair of rows".  Over GF(2) a dot product mod 2 is the parity of a popcount
# of an AND, so the whole test is three machine instructions per pair -- no
# matrix multiplication is ever performed.

def symplectic_product(a: int, b: int, n: int) -> int:
    """<a,b> over GF(2): 0 iff the two Pauli rows commute."""
    mask = (1 << n) - 1
    ax, az = (a >> n) & mask, a & mask
    bx, bz = (b >> n) & mask, b & mask
    return ((ax & bz).bit_count() + (az & bx).bit_count()) & 1


def rows_pairwise_commute(rows: Sequence[int], n: int) -> bool:
    """H_X H_Z^T == H_Z H_X^T (mod 2), evaluated pairwise."""
    m = len(rows)
    for i in range(m):
        ai = rows[i]
        for j in range(i + 1, m):
            if symplectic_product(ai, rows[j], n):
                return False
    return True


def build_commutation_table(spec: "CodeSpec"):
    """
    COMM[a][b] is True iff row values a and b commute.

    This is a table over ROW VALUES (2^2n x 2^2n), not over candidate
    matrices -- exactly like the popcount table the baseline already builds.
    It is constructed once, OUTSIDE the measured interval, and its build time
    is reported separately.
    """
    if _np is None:
        raise RuntimeError("the commutation table requires numpy")
    np = _np
    n = spec.n
    size = spec.row_alphabet
    if size > (1 << 12):
        raise RuntimeError(
            f"commutation table would need {size * size / 2**20:.0f} MiB for "
            f"n={n}; use --backend pure for this code size"
        )
    mask = (1 << n) - 1
    vals = np.arange(size, dtype=np.int64)
    X = ((vals >> n) & mask)          # X half of every row value
    Z = (vals & mask)                 # Z half of every row value
    par = np.array([v.bit_count() & 1 for v in range(1 << n)], dtype=np.uint8)
    t1 = par[(X[:, None] & Z[None, :])]        # parity(a_X . b_Z)
    t2 = par[(Z[:, None] & X[None, :])]        # parity(a_Z . b_X)
    return (t1 ^ t2) == 0


def expected_universal_count(n: int, m: int) -> int:
    """
    Combinatorial reference for the number of m-row candidates satisfying
    BOTH universal constraints, i.e. ordered m-tuples of independent,
    pairwise-commuting rows in a 2n-dimensional symplectic GF(2) space:

        prod_{i=0}^{m-1} ( 2^(2n-i) - 2^i )

    Step i picks a vector in the perp of the i already-chosen rows
    (2^(2n-i) elements) that is outside their span (2^i elements).

    VALIDATION ONLY.  The search never calls this; it is printed by
    --dry-run so a completed run can be checked against it.
    """
    total = 1
    for i in range(m):
        total *= (1 << (2 * n - i)) - (1 << i)
    return total


# =====================================================================
# [S3] TARGET-CODE IDENTIFICATION
#      Recognises a candidate that represents the target code.  It records
#      the hit and NEVER stops the search (see run_search: no break/return
#      is ever executed from inside the traversal).
# =====================================================================

DEFAULT_TARGET_GENERATORS: Dict[Tuple[int, int], Tuple[str, ...]] = {
    # The five-qubit perfect code, cyclic generator set.
    (5, 1): ("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"),
}


@dataclass
class TargetProbe:
    """
    Target-code recogniser.

    mode "span"  : the candidate's rows span exactly the stabilizer group of
                   the target generators (same code, any generating set /
                   ordering).  Expected hit count in a full traversal is the
                   number of ordered bases of the group's GF(2) space.
    mode "exact" : the candidate is literally the target generator matrix, in
                   the given row order.  Expected hit count is 1.
    """

    mode: str
    generators: Tuple[int, ...]
    dim: int
    span: frozenset
    # O(1) gate: values that a candidate row is allowed to take at all.
    row_flag: Tuple[bool, ...]
    # Per-position gate used by the blocked (numpy) traversal.
    position_values: Tuple[frozenset, ...]

    def test(self, rows: Sequence[int]) -> bool:
        if self.mode == "exact":
            return tuple(rows) == self.generators
        # mode == "span"
        for r in rows:
            if r not in self.span:
                return False
        return gf2_rank(rows) == self.dim

    def expected_hits_full_space(self) -> int:
        """Combinatorial invariant for validating a completed full run."""
        if self.mode == "exact":
            return 1
        m = self.dim
        total = 1
        for i in range(len(self.generators)):
            total *= (1 << m) - (1 << i)
        return total


def build_target_probe(spec: CodeSpec, mode: str,
                       pauli_strings: Optional[Sequence[str]]) -> Optional[TargetProbe]:
    if mode == "off":
        return None
    if pauli_strings is None:
        pauli_strings = DEFAULT_TARGET_GENERATORS.get((spec.n, spec.k))
    if pauli_strings is None:
        raise ValueError(
            f"no default target generators known for n={spec.n}, k={spec.k}; "
            f"pass --target-generators or use --target off"
        )
    if len(pauli_strings) != spec.rows:
        raise ValueError(
            f"target needs {spec.rows} generators, got {len(pauli_strings)}"
        )
    gens = tuple(pauli_string_to_row(p, spec.n) for p in pauli_strings)
    span = gf2_span(gens)
    dim = gf2_rank(gens)

    if mode == "span":
        allowed_per_position = tuple(span for _ in range(spec.rows))
        flag = tuple(v in span for v in range(spec.row_alphabet))
    elif mode == "exact":
        allowed_per_position = tuple(frozenset((g,)) for g in gens)
        flag = tuple(v == gens[0] for v in range(spec.row_alphabet))
    else:
        raise ValueError(f"unknown target mode {mode!r}")

    return TargetProbe(mode=mode, generators=gens, dim=dim, span=span,
                       row_flag=flag, position_values=allowed_per_position)


# =====================================================================
# [S4] STRUCTURAL CONSTRAINT EVALUATION
#      One class per constraint.  Constraints are NEVER combined implicitly;
#      a combination happens only when several names are passed explicitly.
#      Adding a new constraint means adding one class plus one registry entry
#      -- the enumeration engine below does not change.
# =====================================================================

class Constraint:
    """Interface for a structural constraint on candidate matrices."""

    name = "abstract"

    def describe(self) -> str:
        raise NotImplementedError

    @property
    def is_trivial(self) -> bool:
        """True for the baseline (no constraint at all)."""
        return False

    def accepts(self, rows: Sequence[int]) -> bool:
        """Matrix-level predicate, used in --mode filter."""
        return True

    def allowed_row_values(self, spec: CodeSpec) -> Optional[List[List[int]]]:
        """
        Row-local form of the constraint: for each row position, the list of
        row values that satisfy it.  Returned only when the constraint really
        is row-local; otherwise None.  Used by --mode generate.
        """
        return None

    def generate_matrices(self, spec: CodeSpec) -> Optional[Iterator[Tuple[int, ...]]]:
        """
        Matrix-level constructive enumeration of exactly the retained space,
        for constraints that are not row-local.  None if unavailable.
        """
        return None

    def retained_size(self, spec: CodeSpec) -> Optional[int]:
        """S_R, the size of the retained subspace, if it is known exactly."""
        return None


class NoConstraint(Constraint):
    """(A) Baseline: the raw candidate space, nothing imposed."""

    name = "none"

    def describe(self) -> str:
        return "None (unconstrained baseline)"

    @property
    def is_trivial(self) -> bool:
        return True

    def allowed_row_values(self, spec: CodeSpec) -> Optional[List[List[int]]]:
        full = list(range(spec.row_alphabet))
        return [full for _ in range(spec.rows)]

    def retained_size(self, spec: CodeSpec) -> Optional[int]:
        return spec.space_size


class WeightSpectrumConstraint(Constraint):
    """
    (C) Weight-spectrum constraint.

    A candidate is retained iff EVERY row weight lies in the supplied
    spectrum W.  The spectrum is never assumed -- it must be given with
    --weights.

    weight_mode:
      "pauli"  (default) row weight = number of qubits on which the row's
               Pauli is not the identity = popcount(x | z).  This is the
               reading under which W={4} is the correct spectrum for the
               [[5,1,3]] code.
      "binary" row weight = popcount of the whole 2n-bit row.
    """

    name = "weight"

    def __init__(self, weights: Iterable[int], weight_mode: str = "pauli"):
        self.weights = frozenset(int(w) for w in weights)
        if not self.weights:
            raise ValueError("weight spectrum must not be empty")
        if weight_mode not in ("pauli", "binary"):
            raise ValueError("weight-mode must be 'pauli' or 'binary'")
        self.weight_mode = weight_mode
        self._table: Optional[List[int]] = None
        self._spec: Optional[CodeSpec] = None

    def describe(self) -> str:
        w = "{" + ",".join(str(x) for x in sorted(self.weights)) + "}"
        return f"Weight spectrum W={w} ({self.weight_mode} row weight)"

    def _weight_table(self, spec: CodeSpec) -> List[int]:
        if self._table is None or self._spec != spec:
            self._table = (build_pauli_weight_table(spec.row_alphabet, spec.n)
                           if self.weight_mode == "pauli"
                           else build_popcount_table(spec.row_alphabet))
            self._spec = spec
        return self._table

    def bind(self, spec: CodeSpec) -> None:
        """Pre-compute the row table (done outside the timed region)."""
        self._weight_table(spec)

    def accepts(self, rows: Sequence[int]) -> bool:
        tbl = self._table
        wts = self.weights
        for r in rows:
            if tbl[r] not in wts:
                return False
        return True

    def allowed_row_values(self, spec: CodeSpec) -> Optional[List[List[int]]]:
        tbl = self._weight_table(spec)
        ok = [v for v in range(spec.row_alphabet) if tbl[v] in self.weights]
        return [ok for _ in range(spec.rows)]

    def retained_size(self, spec: CodeSpec) -> Optional[int]:
        tbl = self._weight_table(spec)
        m = sum(1 for v in range(spec.row_alphabet) if tbl[v] in self.weights)
        return m ** spec.rows


class CyclicConstraint(Constraint):
    """
    (D) Cyclicity / automorphism constraint.

    Imposed condition, exactly as implemented (nothing more general, and no
    additional automorphism test):

        row_{i+1} == cyclic_shift(row_i, shift)     for i = 0 .. rows-2

    where cyclic_shift rotates the qubit indices by `shift` positions,
    applying the same rotation independently to the X block and the Z block.

    NOTE: confirm this is the exact cyclic condition you want before using it
    for a published measurement; it is a single well-defined rule and is easy
    to swap for a different one in this class alone.
    """

    name = "cyclic"

    def __init__(self, shift: int = 1):
        self.shift = int(shift)
        self._n: Optional[int] = None
        self._shift_table: Optional[List[int]] = None

    def describe(self) -> str:
        return (f"Cyclicity: row_(i+1) = qubit-cyclic-shift(row_i, {self.shift}) "
                f"applied to the X and Z blocks independently")

    def bind(self, spec: CodeSpec) -> None:
        self._n = spec.n
        self._shift_table = [cyclic_shift_row(v, spec.n, self.shift)
                             for v in range(spec.row_alphabet)]

    def accepts(self, rows: Sequence[int]) -> bool:
        tbl = self._shift_table
        for i in range(len(rows) - 1):
            if rows[i + 1] != tbl[rows[i]]:
                return False
        return True

    def generate_matrices(self, spec: CodeSpec) -> Optional[Iterator[Tuple[int, ...]]]:
        tbl = self._shift_table
        if tbl is None:
            self.bind(spec)
            tbl = self._shift_table
        rows = spec.rows

        def _gen() -> Iterator[Tuple[int, ...]]:
            for first in range(spec.row_alphabet):
                cand = [first]
                cur = first
                for _ in range(rows - 1):
                    cur = tbl[cur]
                    cand.append(cur)
                yield tuple(cand)

        return _gen()

    def retained_size(self, spec: CodeSpec) -> Optional[int]:
        return spec.row_alphabet


class UniversalStabilizerConstraint(Constraint):
    """
    (B) Universal stabilizer constraints.  EXACTLY two conditions, nothing
    else:

        (1) rank(H) = n - k          over GF(2)
        (2) H_X H_Z^T = H_Z H_X^T    (mod 2)

    Condition (2) is evaluated as the pairwise symplectic form
    <r_i, r_j> = 0 (see [S2]); it is mathematically identical to the matrix
    identity but costs a few bit operations per pair instead of a matrix
    multiplication.

    Explicitly NOT imposed: distance, [[5,1,3]] validity, weight spectrum,
    cyclicity, automorphisms, canonical form, duplicate-generating-set
    removal, or anything else.

    Not row-local (both conditions couple rows), so there is no
    allowed_row_values() form; the dedicated engine search_numpy_universal()
    evaluates it hierarchically over the blocked traversal.
    """

    name = "universal"

    def __init__(self) -> None:
        self._n: Optional[int] = None
        self._rows: Optional[int] = None
        self._comm = None            # numpy bool table, built by bind()
        self.table_build_seconds: float = 0.0

    def describe(self) -> str:
        m = self._rows if self._rows is not None else "n-k"
        return (f"Universal stabilizer constraints: rank(H)={m} over GF(2) "
                f"AND H_X H_Z^T = H_Z H_X^T (mod 2)")

    def bind(self, spec: CodeSpec) -> None:
        """Pre-computation, performed OUTSIDE the measured interval."""
        self._n = spec.n
        self._rows = spec.rows
        if _np is not None:
            t0 = time.perf_counter()
            try:
                self._comm = build_commutation_table(spec)
            except RuntimeError:
                self._comm = None    # too large; pure backend still works
            self.table_build_seconds = time.perf_counter() - t0

    # ---- condition (2) then condition (1); both must hold ----------------
    def accepts(self, rows: Sequence[int]) -> bool:
        n = self._n
        if not rows_pairwise_commute(rows, n):      # H_X H_Z^T = H_Z H_X^T
            return False
        return gf2_rank(rows) == len(rows)          # rank(H) = n - k

    def allowed_row_values(self, spec: CodeSpec) -> Optional[List[List[int]]]:
        return None                                  # not row-local

    def retained_size(self, spec: CodeSpec) -> Optional[int]:
        # Validation reference only; never used to shortcut the search.
        return expected_universal_count(spec.n, spec.rows)


class CombinedConstraint(Constraint):
    """(E) Explicit conjunction of several constraints. Never implicit."""

    name = "combined"

    def __init__(self, parts: Sequence[Constraint]):
        self.parts = list(parts)

    def describe(self) -> str:
        return " AND ".join(p.describe() for p in self.parts)

    def bind(self, spec: CodeSpec) -> None:
        for p in self.parts:
            binder = getattr(p, "bind", None)
            if binder is not None:
                binder(spec)

    def accepts(self, rows: Sequence[int]) -> bool:
        for p in self.parts:
            if not p.accepts(rows):
                return False
        return True

    def allowed_row_values(self, spec: CodeSpec) -> Optional[List[List[int]]]:
        per_part = [p.allowed_row_values(spec) for p in self.parts]
        if any(x is None for x in per_part):
            return None
        out: List[List[int]] = []
        for i in range(spec.rows):
            acc = set(per_part[0][i])
            for other in per_part[1:]:
                acc &= set(other[i])
            out.append(sorted(acc))
        return out

    def retained_size(self, spec: CodeSpec) -> Optional[int]:
        allowed = self.allowed_row_values(spec)
        if allowed is None:
            return None
        total = 1
        for vals in allowed:
            total *= len(vals)
        return total


def build_constraint(names: Sequence[str], args: argparse.Namespace) -> Constraint:
    """Constraint registry.  Add new constraints here and in [S4] only."""
    parts: List[Constraint] = []
    for raw in names:
        nm = raw.strip().lower()
        if nm in ("none", "baseline", ""):
            parts.append(NoConstraint())
        elif nm in ("weight", "weight_spectrum", "spectrum"):
            if not args.weights:
                raise ValueError(
                    "--constraint weight requires --weights (e.g. --weights 4); "
                    "no spectrum is ever assumed"
                )
            wts = [int(x) for x in args.weights.split(",") if x.strip() != ""]
            parts.append(WeightSpectrumConstraint(wts, args.weight_mode))
        elif nm in ("cyclic", "cyclicity"):
            parts.append(CyclicConstraint(args.cyclic_shift))
        elif nm == "universal":
            parts.append(UniversalStabilizerConstraint())
        else:
            raise ValueError(f"unknown constraint {raw!r}")

    real = [p for p in parts if not p.is_trivial]
    if not real:
        return NoConstraint()
    if len(real) == 1:
        return real[0]
    return CombinedConstraint(real)


# =====================================================================
# [S5] Search plan construction
#      Everything here is pre-computation and happens strictly OUTSIDE the
#      measured interval: lookup tables, the enumeration alphabets and the
#      target probe.  It is rebuilt identically before every repetition so
#      that no repetition can benefit from work done by a previous one, and
#      it is identical in structure for the baseline and constrained runs.
# =====================================================================

@dataclass
class SearchPlan:
    spec: CodeSpec
    constraint: Constraint
    mode: str                       # "filter" or "generate"
    backend: str                    # "numpy" or "pure"
    row_values: Optional[List[List[int]]]      # per-row-position alphabets
    matrix_iter_factory: Optional[Callable[[], Iterator[Tuple[int, ...]]]]
    pop_table: List[int]            # per-candidate processing workload table
    constraint_row_ok: Optional[List[bool]]    # row-local constraint mask
    check_constraint: bool          # evaluate the constraint per candidate?
    target: Optional[TargetProbe]
    announce_all: bool
    progress: bool
    constraint_eval: str = "hierarchical"   # or "per-candidate"
    work: int = 1                   # per-candidate processing units (--work)


@dataclass
class SearchResult:
    examined: int = 0               # candidate matrices enumerated
    retained: int = 0               # candidates satisfying the constraint
    target_hits: int = 0            # target-code occurrences
    first_target_index: Optional[int] = None
    checksum: int = 0               # per-candidate processing accumulator
    elapsed: float = 0.0
    # Constrained-run diagnostics (0 for the unconstrained baseline path).
    blocks_total: int = 0
    blocks_short_circuited: int = 0
    # Candidates whose constraint status was resolved by the pair-stage test
    # rather than by a cheaper prefix condition. This is the quantity that
    # explains why the runtime ratio exceeds S_R / S_0.
    pair_stage_evaluated: int = 0

    @property
    def rejected(self) -> int:
        return self.examined - self.retained


def build_plan(spec: CodeSpec, constraint: Constraint, mode: str, backend: str,
               target: Optional[TargetProbe], announce_all: bool,
               progress: bool,
               constraint_eval: str = "hierarchical",
               work: int = 1) -> SearchPlan:
    binder = getattr(constraint, "bind", None)
    if binder is not None:
        binder(spec)

    pop_table = build_popcount_table(spec.row_alphabet)
    full = list(range(spec.row_alphabet))

    row_values: Optional[List[List[int]]] = None
    matrix_iter_factory: Optional[Callable[[], Iterator[Tuple[int, ...]]]] = None
    check_constraint = False
    constraint_row_ok: Optional[List[bool]] = None

    if mode == "filter" or constraint.is_trivial:
        # ---- traverse the FULL space, evaluating the constraint per candidate
        row_values = [full for _ in range(spec.rows)]
        check_constraint = not constraint.is_trivial
        if check_constraint:
            local = constraint.allowed_row_values(spec)
            # The fast row-mask path is only valid when the constraint imposes
            # the SAME row-local condition at every row position.  If a future
            # constraint is position-dependent, drop back to accepts() rather
            # than silently applying position 0's mask everywhere.
            if local is not None and all(list(v) == list(local[0]) for v in local):
                ok = [False] * spec.row_alphabet
                for v in local[0]:
                    ok[v] = True
                constraint_row_ok = ok
    else:
        # ---- "generate": enumerate exactly the retained subspace
        local = constraint.allowed_row_values(spec)
        if local is not None:
            row_values = [list(vals) for vals in local]
        else:
            gen = constraint.generate_matrices(spec)
            if gen is None:
                raise ValueError(
                    f"constraint {constraint.name!r} has no constructive form; "
                    f"use --mode filter"
                )
            matrix_iter_factory = lambda: constraint.generate_matrices(spec)

    return SearchPlan(
        spec=spec, constraint=constraint, mode=mode, backend=backend,
        row_values=row_values, matrix_iter_factory=matrix_iter_factory,
        pop_table=pop_table, constraint_row_ok=constraint_row_ok,
        check_constraint=check_constraint, target=target,
        announce_all=announce_all, progress=progress,
        constraint_eval=constraint_eval, work=max(1, int(work)),
    )


def plan_enumeration_size(plan: SearchPlan) -> Optional[int]:
    """How many candidates the plan will enumerate (a counter, not a result)."""
    if plan.row_values is not None:
        total = 1
        for vals in plan.row_values:
            total *= len(vals)
        return total
    return plan.constraint.retained_size(plan.spec)


# =====================================================================
# [S6] CANDIDATE ENUMERATION  (this is the timed work)
#
#      Three engines with identical semantics.  In each of them the loop body
#      performs, for EVERY enumerated candidate and in this order:
#         (1) per-candidate processing  -- row weights -> matrix weight sum
#                                          -> checksum fold
#         (2) structural constraint evaluation (only when one is imposed)
#         (3) target-code identification (O(1)-gated, never terminates)
#      There is no break, no return and no early exit inside any traversal.
# =====================================================================

_PROGRESS_STEPS = 200

# Candidates are processed in tiles of this many elements in EVERY engine.
# Without it the unconstrained engine works on 1,048,576-element arrays while
# a constrained engine works on far smaller, cache-resident ones, so the same
# arithmetic costs a different amount per candidate and T/T_none stops being
# a clean measure of the search-space reduction. Equal tiles remove that
# confound: the only thing left to distinguish the engines is how many
# candidates they process and what the constraint test costs.
_PROC_TILE = 32768


def _progress_line(done: int, total: Optional[int], t0: float) -> None:
    el = time.perf_counter() - t0
    if total:
        pct = 100.0 * done / total
        sys.stderr.write(f"\r  [progress] {done:,}/{total:,} blocks "
                         f"({pct:6.2f}%)  {el:9.1f} s elapsed")
    else:
        sys.stderr.write(f"\r  [progress] {done:,} blocks  {el:9.1f} s elapsed")
    sys.stderr.flush()


def search_pure(plan: SearchPlan) -> SearchResult:
    """Plain Python traversal: one loop iteration per candidate matrix."""
    res = SearchResult()
    POP = plan.pop_table
    check = plan.check_constraint
    accepts = plan.constraint.accepts
    row_ok = plan.constraint_row_ok
    probe = plan.target
    tflag = probe.row_flag if probe is not None else None
    ttest = probe.test if probe is not None else None
    announce_all = plan.announce_all

    if plan.matrix_iter_factory is not None:
        candidates: Iterable[Tuple[int, ...]] = plan.matrix_iter_factory()
    else:
        candidates = product(*plan.row_values)

    examined = 0
    retained = 0
    hits = 0
    checksum = 0
    first_hit: Optional[int] = None

    work = plan.work
    subsets = combination_subsets(plan.spec.rows) if work > 1 else ()
    nsub = len(subsets)

    for cand in candidates:
        examined += 1

        # ---------- (1) structural constraint evaluation, FIRST ----------
        # A rejected candidate is skipped immediately and never reaches the
        # candidate-processing workload below.  With no constraint imposed
        # (`check` is False) this branch vanishes and the loop body is the
        # unconstrained baseline unchanged.
        if check:
            if row_ok is not None:
                ok = True
                for r in cand:
                    if not row_ok[r]:
                        ok = False
                        break
            else:
                ok = accepts(cand)
            if not ok:
                continue

        # ---------- (2) per-candidate processing (accepted only) ----------
        w = 0
        for r in cand:
            w += POP[r]
        checksum += w
        # extra processing units: weights of row-space combinations
        for t in range(1, work):
            v = 0
            for i in subsets[(t - 1) % nsub]:
                v ^= cand[i]
            checksum += POP[v]
        retained += 1

        # ---------- (3) target-code identification (no early exit) ----------
        if tflag is not None and tflag[cand[0]] and ttest(cand):
            hits += 1
            if first_hit is None:
                first_hit = examined - 1
                print("Code Found")
            elif announce_all:
                print("Code Found")

    res.examined = examined
    res.retained = retained
    res.target_hits = hits
    res.checksum = checksum
    res.first_target_index = first_hit
    return res


def search_numpy(plan: SearchPlan) -> SearchResult:
    """
    Vectorised traversal, identical space and identical order.

    The last two row positions form a block of size (na x nb); the remaining
    row positions are iterated in Python.  Every candidate in every block is
    individually weighted and folded into the checksum -- the block is a
    memory layout, not a shortcut.  The only thing shared across candidates
    is the per-row-value popcount table, which the pure backend uses too.
    """
    if _np is None:
        raise RuntimeError("numpy backend requested but numpy is not installed")
    if plan.matrix_iter_factory is not None:
        raise RuntimeError("numpy backend needs a product-form enumeration")

    spec = plan.spec
    np = _np
    res = SearchResult()

    row_values = plan.row_values
    assert row_values is not None
    if spec.rows < 2:
        return _search_numpy_single_row(plan)

    outer_values = row_values[:-2]
    A = np.asarray(row_values[-2], dtype=np.int64)
    B = np.asarray(row_values[-1], dtype=np.int64)
    na, nb = int(A.size), int(B.size)
    block = na * nb
    if block == 0:
        return res

    # Accumulator width: max possible per-candidate weight sum is rows*cols.
    max_sum = spec.rows * spec.cols
    acc_dtype = (np.uint8 if max_sum < 256
                 else np.uint16 if max_sum < 65536 else np.uint32)

    POP_np = np.asarray(plan.pop_table, dtype=np.int64)
    POP = plan.pop_table
    wA = POP_np[A].astype(acc_dtype).reshape(na, 1)   # weights of row rows-2
    wB = POP_np[B].astype(acc_dtype).reshape(1, nb)   # weights of row rows-1
    trows = max(1, min(na, _PROC_TILE // max(1, nb)))
    wbuf = np.empty((trows, nb), dtype=acc_dtype)

    # Reduction accumulator for every processing unit. Defined here, not
    # inside the --work branch, so the base unit uses the SAME width as the
    # constrained engines do -- otherwise the baseline would pay for a
    # uint64 reduction that the constrained runs avoid, which is worth about
    # 2% and would bias T_R / T_0.
    red_dtype = (np.uint32 if na * nb * spec.cols < (1 << 32) else np.uint64)

    # ---- extra per-candidate processing units (--work) ------------------
    work = plan.work
    if work > 1:
        subsets = combination_subsets(spec.rows)
        nsub = len(subsets)
        val_dtype = np.int16 if spec.cols <= 15 else np.int32
        POP_val = np.asarray(plan.pop_table, dtype=np.uint8)
        Acol_v = A.astype(val_dtype).reshape(na, 1)
        Brow_v = B.astype(val_dtype).reshape(1, nb)
        vbuf = np.empty((trows, nb), dtype=val_dtype)
        obuf = np.empty((trows, nb), dtype=np.uint8)

    # Constraint masks for the two blocked row positions (row-local only).
    check = plan.check_constraint
    row_ok = plan.constraint_row_ok
    if check and row_ok is None:
        raise RuntimeError(
            "numpy backend can only filter row-local constraints; "
            "use --backend pure for this constraint"
        )
    if check:
        ok_np = np.asarray(row_ok, dtype=bool)
        mA = ok_np[A].reshape(na, 1)
        mB = ok_np[B].reshape(1, nb)
        mbuf = np.empty((na, nb), dtype=bool)

    # Target-probe pre-computation: positions inside A / B that could
    # possibly participate in a target matrix.  Sorted so that hits are found
    # in the same lexicographic order the traversal visits them.
    probe = plan.target
    if probe is not None:
        pos_a = probe.position_values[spec.rows - 2]
        pos_b = probe.position_values[spec.rows - 1]
        idx_a = [i for i in range(na) if int(A[i]) in pos_a]
        idx_b = [i for i in range(nb) if int(B[i]) in pos_b]
        val_a = [int(A[i]) for i in idx_a]
        val_b = [int(B[i]) for i in idx_b]
        outer_sets = probe.position_values[:spec.rows - 2]
        ttest = probe.test
    announce_all = plan.announce_all

    outer_total = 1
    for vals in outer_values:
        outer_total *= len(vals)
    prog_every = max(1, outer_total // _PROGRESS_STEPS)

    examined = 0
    retained = 0
    hits = 0
    checksum = 0
    first_hit: Optional[int] = None
    blocks_done = 0
    t0 = time.perf_counter()

    outer_iter: Iterable[Tuple[int, ...]] = (product(*outer_values)
                                             if outer_values else ((),))

    for outer in outer_iter:
        # ---------- (1) per-candidate processing, whole block ----------
        w_outer = 0
        for r in outer:
            w_outer += POP[r]
        # Per-block XOR constants for the extra units, hoisted out of the
        # tile loop so tiling adds no arithmetic.
        if work > 1:
            consts = []
            for t in range(1, work):
                sel = subsets[(t - 1) % nsub]
                const = 0
                for i in sel:
                    if i < len(outer):
                        const ^= outer[i]
                consts.append(const)

        # Candidate processing, in fixed-size tiles (see _PROC_TILE).
        for i0 in range(0, na, trows):
            i1 = i0 + trows
            if i1 > na:
                i1 = na
            h = i1 - i0
            wt = wbuf[:h]
            np.add(wA[i0:i1], wB, out=wt)        # weight of the two block rows
            if w_outer:
                wt += acc_dtype(w_outer)         # + weight of the outer rows
            checksum += int(wt.sum(dtype=red_dtype))

            # extra processing units: weights of row-space combinations
            if work > 1:
                vt = vbuf[:h]
                ot = obuf[:h]
                ac = Acol_v[i0:i1]
                for const in consts:
                    np.bitwise_xor(ac, Brow_v, out=vt)
                    if const:
                        vt ^= val_dtype(const)
                    np.take(POP_val, vt, out=ot)
                    checksum += int(ot.sum(dtype=red_dtype))
        examined += block

        # ---------- (2) structural constraint evaluation ----------
        if check:
            np.logical_and(mA, mB, out=mbuf)     # evaluated for every candidate
            cnt = int(mbuf.sum())
            outer_ok = True
            for r in outer:
                if not row_ok[r]:
                    outer_ok = False
                    break
            retained += cnt if outer_ok else 0
        else:
            retained += block

        # ---------- (3) target-code identification (no early exit) ----------
        if probe is not None:
            outer_ok_t = True
            for i, r in enumerate(outer):
                if r not in outer_sets[i]:
                    outer_ok_t = False
                    break
            if outer_ok_t and idx_a and idx_b:
                base = examined - block
                for pa, ra in zip(idx_a, val_a):
                    row_base = base + pa * nb
                    for pb, rb in zip(idx_b, val_b):
                        cand = outer + (ra, rb)
                        if ttest(cand):
                            hits += 1
                            if first_hit is None:
                                first_hit = row_base + pb
                                print("Code Found")
                            elif announce_all:
                                print("Code Found")

        blocks_done += 1
        if plan.progress and blocks_done % prog_every == 0:
            _progress_line(blocks_done, outer_total, t0)

    if plan.progress:
        sys.stderr.write("\n")
        sys.stderr.flush()

    res.examined = examined
    res.retained = retained
    res.target_hits = hits
    res.checksum = checksum
    res.first_target_index = first_hit
    return res


def _search_numpy_single_row(plan: SearchPlan) -> SearchResult:
    """Degenerate case rows == 1, kept so small specs work in both backends."""
    np = _np
    spec = plan.spec
    res = SearchResult()
    A = np.asarray(plan.row_values[0], dtype=np.int64)
    POP_np = np.asarray(plan.pop_table, dtype=np.int64)
    w = POP_np[A]
    res.checksum = int(w.sum())
    res.examined = int(A.size)
    if plan.check_constraint:
        ok_np = np.asarray(plan.constraint_row_ok, dtype=bool)
        res.retained = int(ok_np[A].sum())
    else:
        res.retained = res.examined
    probe = plan.target
    if probe is not None:
        for i in range(int(A.size)):
            cand = (int(A[i]),)
            if probe.row_flag[cand[0]] and probe.test(cand):
                res.target_hits += 1
                if res.first_target_index is None:
                    res.first_target_index = i
                    print("Code Found")
                elif plan.announce_all:
                    print("Code Found")
    return res


def search_numpy_weight(plan: SearchPlan) -> SearchResult:
    """
    CONSTRAINED ENGINE for the constant-row-weight constraint (EXPERIMENT C).

    The ONLY condition imposed is

        wt(row_i) in W        for every row i of the candidate matrix

    with W supplied explicitly by --weights (W={4} for [[5,1,3]]_2, i.e. the
    row-weight vector must be [4,4,4,4]).  Nothing else is imposed: no
    rank(H) = n-k, no H_X H_Z^T = H_Z H_X^T, no distance, no cyclicity, no
    automorphism condition, no [[5,1,3]] validity test of any kind.

    Enumeration order, matrix representation, processing workload, timing
    methodology and target recognition are identical to the other programs
    in this set, so the runtimes are directly comparable.

    The condition is ROW-LOCAL, so it decomposes as:

      Stage 1  every outer row must satisfy the weight condition.  If one
               does not, every candidate sharing that prefix violates the
               constraint and the whole block is rejected.
      Stage 2  each inner row must lie in S = { v : wt(v) in W }.
      Stage 3  there is NO pair condition.  The accepted set of a surviving
               block is therefore exactly the rectangle S x S, so the
               processing workload runs on it directly: no acceptance mask
               and no extraction are needed, because no candidate inside the
               rectangle is rejected.

    S is derived once from the row alphabet, outside the measured interval,
    exactly like the popcount table.  Rejected candidates never reach the
    processing workload.
    """
    if _np is None:
        raise RuntimeError("numpy backend requested but numpy is not installed")
    np = _np
    spec = plan.spec
    res = SearchResult()

    if plan.row_values is None or spec.rows < 2:
        raise RuntimeError("weight engine needs rows >= 2 and product form")
    alphabet = spec.row_alphabet
    if (len(plan.row_values[-2]) != alphabet
            or len(plan.row_values[-1]) != alphabet):
        raise RuntimeError("weight engine expects the full row alphabet")
    row_ok = plan.constraint_row_ok
    if row_ok is None:
        raise RuntimeError("weight engine needs a row-local constraint mask")

    outer_values = plan.row_values[:-2]
    p = spec.rows - 2
    block = alphabet * alphabet

    # S = row values satisfying the weight condition (405 values for W={4}).
    S = np.asarray([v for v in range(alphabet) if row_ok[v]], dtype=np.int64)
    S_list = [int(v) for v in S]
    ns = int(S.size)

    max_sum = spec.rows * spec.cols
    acc_dtype = (np.uint8 if max_sum < 256
                 else np.uint16 if max_sum < 65536 else np.uint32)
    red_dtype = (np.uint32 if max(1, ns * ns) * max_sum < (1 << 32)
                 else np.uint64)

    POP = plan.pop_table
    POP_np = np.asarray(POP, dtype=acc_dtype)
    POPa = POP_np[S].reshape(ns, 1) if ns else None
    POPb = POP_np[S].reshape(1, ns) if ns else None
    wtrows = max(1, min(ns, _PROC_TILE // max(1, ns))) if ns else 1
    wbuf = np.empty((wtrows, ns), dtype=acc_dtype) if ns else None

    # ---- extra per-candidate processing units (--work) ------------------
    work = plan.work
    subsets = combination_subsets(spec.rows) if work > 1 else ()
    nsub = len(subsets)
    val_dtype = np.int16 if spec.cols <= 15 else np.int32
    POP_val = np.asarray(plan.pop_table, dtype=np.uint8)
    if work > 1 and ns:
        # Same broadcast form the unconstrained engine uses, so the tile does
        # the same reads and writes per candidate there as it does here.
        av = S.reshape(ns, 1).astype(val_dtype)
        bv = S.reshape(1, ns).astype(val_dtype)
        vbuf = np.empty((wtrows, ns), dtype=val_dtype)
        obuf = np.empty((wtrows, ns), dtype=np.uint8)

    probe = plan.target
    if probe is not None:
        tgt_a = [v for v in S_list
                 if v in probe.position_values[spec.rows - 2]]
        tgt_b = [v for v in S_list
                 if v in probe.position_values[spec.rows - 1]]
    announce_all = plan.announce_all

    outer_total = 1
    for vals in outer_values:
        outer_total *= len(vals)
    prog_every = max(1, outer_total // _PROGRESS_STEPS)

    examined = 0
    retained = 0
    hits = 0
    checksum = 0
    first_hit: Optional[int] = None
    blocks_done = 0
    short_circuited = 0
    t0 = time.perf_counter()
    outer_iter: Iterable[Tuple[int, ...]] = (product(*outer_values)
                                             if outer_values else ((),))

    for outer in outer_iter:
        examined += block                  # traversal covers every candidate
        blocks_done += 1
        base = examined - block

        # ---------- Stage 1: weight condition on the outer rows -----------
        ok = True
        for r in outer:
            if not row_ok[r]:
                ok = False
                break
        if not ok:
            short_circuited += 1
            if plan.progress and blocks_done % prog_every == 0:
                _progress_line(blocks_done, outer_total, t0)
            continue

        # ---------- Stages 2 and 3: inner rows must lie in S --------------
        n_acc = ns * ns
        retained += n_acc
        if n_acc == 0:
            if plan.progress and blocks_done % prog_every == 0:
                _progress_line(blocks_done, outer_total, t0)
            continue

        # ---------- candidate processing: ACCEPTED candidates only --------
        w_outer = 0
        for r in outer:
            w_outer += POP[r]
        if work > 1:
            consts = []
            for t in range(1, work):
                sel = subsets[(t - 1) % nsub]
                const = 0
                for i in sel:
                    if i < p:
                        const ^= outer[i]
                consts.append(const)

        # Candidate processing, in the SAME fixed-size tiles the
        # unconstrained engine uses (see _PROC_TILE).
        for i0 in range(0, ns, wtrows):
            i1 = i0 + wtrows
            if i1 > ns:
                i1 = ns
            h = i1 - i0
            wt = wbuf[:h]
            np.add(POPa[i0:i1], POPb, out=wt)
            if w_outer:
                wt += acc_dtype(w_outer)
            checksum += int(wt.sum(dtype=red_dtype))
            if work > 1:
                vt = vbuf[:h]
                ot = obuf[:h]
                ac = av[i0:i1]
                for const in consts:
                    np.bitwise_xor(ac, bv, out=vt)
                    if const:
                        vt ^= val_dtype(const)
                    np.take(POP_val, vt, out=ot)
                    checksum += int(ot.sum(dtype=red_dtype))

        # ---------- target-code identification (no early exit) ------------
        if probe is not None:
            ok_t = True
            for i, r in enumerate(outer):
                if r not in probe.position_values[i]:
                    ok_t = False
                    break
            if ok_t:
                for ra in tgt_a:                  # ascending == lex order
                    row_base = base + ra * alphabet
                    for rb in tgt_b:
                        cand = outer + (ra, rb)
                        if probe.test(cand):
                            hits += 1
                            if first_hit is None:
                                first_hit = row_base + rb
                                print("Code Found")
                            elif announce_all:
                                print("Code Found")

        if plan.progress and blocks_done % prog_every == 0:
            _progress_line(blocks_done, outer_total, t0)

    if plan.progress:
        sys.stderr.write("\n")
        sys.stderr.flush()

    res.examined = examined
    res.retained = retained
    res.target_hits = hits
    res.checksum = checksum
    res.first_target_index = first_hit
    res.blocks_total = blocks_done
    res.blocks_short_circuited = short_circuited
    res.pair_stage_evaluated = 0        # row-local: there is no pair stage
    return res


def run_search(plan: SearchPlan) -> SearchResult:
    """Dispatch to the selected engine. Contains no timing and no I/O setup."""
    if plan.backend == "numpy":
        if (isinstance(plan.constraint, WeightSpectrumConstraint)
                and plan.check_constraint):
            return search_numpy_weight(plan)
        return search_numpy(plan)
    return search_pure(plan)


# =====================================================================
# [S7] RUNTIME MEASUREMENT
#      The measured interval starts immediately before enumeration and ends
#      immediately after it.  Plan construction, table building, argument
#      parsing, imports and all printing of results sit outside it.
# =====================================================================

def measure_runtime(plan: SearchPlan) -> SearchResult:
    t_start = time.perf_counter()          # <-- timer starts here
    result = run_search(plan)
    t_end = time.perf_counter()            # <-- timer stops here
    result.elapsed = t_end - t_start
    return result


def run_benchmark(spec: CodeSpec, constraint: Constraint, mode: str, backend: str,
                  target: Optional[TargetProbe], repetitions: int,
                  announce_all: bool, progress: bool,
                  constraint_eval: str = "hierarchical",
                  work: int = 1, r0_start: int = 0,
                  r0_count: Optional[int] = None) -> List[SearchResult]:
    results: List[SearchResult] = []
    for rep in range(repetitions):
        # The plan (and every lookup table in it, including the commutation
        # table) is rebuilt from scratch for each repetition, outside the
        # timed region, so that no repetition can reuse anything computed by
        # a previous one.
        plan = build_plan(spec, constraint, mode, backend, target,
                          announce_all, progress, constraint_eval, work)
        if r0_count is not None and plan.row_values is not None:
            # BOUNDED SAMPLE: keep a contiguous slice of the row-0 alphabet.
            # Every other row still ranges over its full alphabet, so the
            # sampled region is an exact contiguous sub-range of the same
            # enumeration, in the same order.
            plan.row_values = ([plan.row_values[0][r0_start:r0_start + r0_count]]
                               + list(plan.row_values[1:]))
        if repetitions > 1:
            print(f"-- repetition {rep + 1}/{repetitions}")
        results.append(measure_runtime(plan))
    return results


# =====================================================================
# [S8] Reporting  (aggregate only: never a single candidate matrix)
# =====================================================================

def report(spec: CodeSpec, constraint: Constraint, mode: str, backend: str,
           target: Optional[TargetProbe], results: List[SearchResult],
           constraint_eval: str = "hierarchical", work: int = 1) -> None:
    r0 = results[0]
    print()
    print(f"Code: {spec.label}")
    print(f"Per-candidate processing units (--work): {work}")
    print(f"Matrix shape: {spec.rows} x {spec.cols} binary "
          f"(entries independent, 2 values each)")
    print(f"Constraint: {constraint.describe()}")
    print(f"Constraint application mode: {mode}")
    if not constraint.is_trivial:
        print(f"Constraint evaluation strategy: {constraint_eval}")
    print(f"Backend: {backend}")
    print(f"Unconstrained search-space size S_0: {spec.space_size:,}")
    print(f"Total candidate matrices considered: {r0.examined:,}")
    print(f"Candidates satisfying constraint: {r0.retained:,}")
    print(f"Candidates rejected by constraint: {r0.rejected:,}")
    if r0.examined:
        frac = r0.retained / r0.examined
        print(f"Fraction satisfying constraint: {frac:.12g}")
    if r0.blocks_total:
        print(f"Traversal blocks: {r0.blocks_total:,} "
              f"(prefix-rejected: {r0.blocks_short_circuited:,})")
    if r0.pair_stage_evaluated:
        print(f"Candidates resolved by the pair-stage test: "
              f"{r0.pair_stage_evaluated:,}")
        if r0.examined:
            print(f"  as a fraction of S_0: "
                  f"{r0.pair_stage_evaluated / r0.examined:.12g}")
        if r0.retained:
            print(f"  per accepted candidate: "
                  f"{r0.pair_stage_evaluated / r0.retained:.6g}")
    build_s = getattr(constraint, "table_build_seconds", 0.0)
    if build_s:
        print(f"Constraint table build time (EXCLUDED from runtime): "
              f"{build_s:.6f} seconds")
    if target is not None:
        print(f"Target code: {target.mode} match against "
              f"{[row_to_pauli_string(g, spec.n) for g in target.generators]}")
        print(f"Target code occurrences: {r0.target_hits:,}")
        if r0.first_target_index is not None:
            print(f"First target occurrence at enumeration index: "
                  f"{r0.first_target_index:,}")
        else:
            print("First target occurrence at enumeration index: (none)")
    else:
        print("Target code occurrences: (target check disabled)")
    if constraint.is_trivial:
        print(f"Processing checksum: {r0.checksum:,}")
    else:
        print(f"Processing checksum (ACCEPTED candidates only, not "
              f"comparable to the unconstrained run): {r0.checksum:,}")

    times = [r.elapsed for r in results]
    for i, t in enumerate(times, 1):
        rate = r0.examined / t if t > 0 else float("inf")
        print(f"Runtime[{i}]: {t:.6f} seconds "
              f"({rate:,.0f} candidates/second)")
    if len(times) > 1:
        print(f"Mean runtime: {statistics.mean(times):.6f} seconds")
        print(f"Median runtime: {statistics.median(times):.6f} seconds")
        print(f"Stdev runtime: {statistics.stdev(times):.6f} seconds")
    mean_t = statistics.mean(times)
    if mean_t > 0:
        print(f"Candidate-processing rate: {r0.examined / mean_t:,.0f} "
              f"candidates/second (enumerated)")
        if not constraint.is_trivial:
            print(f"Accepted-candidate rate: {r0.retained / mean_t:,.0f} "
                  f"candidates/second (processed)")

    # Consistency notes (facts about this run, not analysis of the results).
    consistent = all(r.examined == r0.examined and r.retained == r0.retained
                     and r.target_hits == r0.target_hits
                     and r.checksum == r0.checksum for r in results)
    print(f"Repetitions identical: {'yes' if consistent else 'NO - INVESTIGATE'}")


def dry_run_report(spec: CodeSpec, constraint: Constraint, mode: str,
                   backend: str, target: Optional[TargetProbe]) -> None:
    plan = build_plan(spec, constraint, mode, backend, target, False, False)
    enum_size = plan_enumeration_size(plan)
    print(f"Code: {spec.label}")
    print(f"Matrix shape: {spec.rows} x {spec.cols} binary")
    print(f"Constraint: {constraint.describe()}")
    print(f"Constraint application mode: {mode}")
    print(f"Backend: {backend}")
    print(f"Unconstrained search-space size S_0: {spec.space_size:,}")
    sr = constraint.retained_size(spec)
    if sr is not None:
        print(f"Retained subspace size S_R (combinatorial reference for "
              f"VALIDATION only; the search never uses it): {sr:,}")
        print(f"Expected fraction S_R / S_0: {sr / spec.space_size:.12g}")
    print(f"Candidates this run will enumerate: "
          f"{enum_size:,}" if enum_size is not None else "unknown")
    if backend == "numpy" and spec.rows >= 2 and plan.row_values is not None:
        na = len(plan.row_values[-2])
        nb = len(plan.row_values[-1])
        outer = 1
        for vals in plan.row_values[:-2]:
            outer *= len(vals)
        print(f"Vectorised blocks: {outer:,} blocks of {na * nb:,} candidates")
    if target is not None:
        print(f"Target generators: "
              f"{[row_to_pauli_string(g, spec.n) for g in target.generators]}")
        print(f"Target occurrences expected in a FULL unconstrained traversal: "
              f"{target.expected_hits_full_space():,}")
    print("(dry run: nothing was enumerated, no timing performed)")


def sample_report(spec: CodeSpec, results: List[SearchResult],
                  r0_start: int, r0_count: Optional[int], work: int) -> None:
    """Bounded-sample accounting and the full-search extrapolation."""
    if r0_count is None or r0_count >= spec.row_alphabet:
        return
    r0 = results[0]
    full = spec.space_size
    factor = full / r0.examined if r0.examined else 0.0
    print()
    print("--- BOUNDED SAMPLE ---")
    print(f"Sampled region: row 0 in [{r0_start}, {r0_start + r0_count}) "
          f"of [0, {spec.row_alphabet})  (contiguous, deterministic)")
    print(f"Sample fraction of the enumeration: {r0_count}/"
          f"{spec.row_alphabet} = {r0_count / spec.row_alphabet:.10g}")
    print(f"Full search-space size S_0: {full:,}")
    print(f"Sample candidate count: {r0.examined:,}")
    print(f"Sample candidates satisfying constraint: {r0.retained:,}")
    if r0.examined:
        print(f"Sample retention fraction: {r0.retained / r0.examined:.12g}")
    print(f"Extrapolation factor (S_0 / sample): {factor:.10g}")
    times = [r.elapsed for r in results]
    ests = [t * factor for t in times]
    for i, (t, e) in enumerate(zip(times, ests), 1):
        print(f"  sample runtime[{i}]: {t:.6f} s  ->  estimated full "
              f"runtime: {e:.3f} s ({e / 60:.2f} min, {e / 3600:.4f} h)")
    mean_e = statistics.mean(ests)
    print(f"Estimated full-search runtime (mean): {mean_e:.3f} s "
          f"({mean_e / 60:.2f} min, {mean_e / 3600:.4f} h)")
    if len(ests) > 1:
        print(f"Estimated full-search runtime (median): "
              f"{statistics.median(ests):.3f} s")
        print(f"Estimated full-search runtime (stdev): "
              f"{statistics.stdev(ests):.3f} s")
    print(f"Workload units used (--work): {work}")
    print("NOTE: the full 2**{} search was NOT executed; the number above is "
          "an extrapolation from the bounded sample.".format(
              spec.rows * spec.cols))


# =====================================================================
# [S9] CLI and self-test
# =====================================================================

def self_test() -> int:
    """
    Cheap static/behavioural checks that do NOT run the big search.
    Verifies the encoders, the GF(2) helpers, the target probe and that the
    two backends agree on a tiny specification.
    """
    failures = 0

    def check(label: str, got, want) -> None:
        nonlocal failures
        ok = got == want
        print(f"  [{'ok ' if ok else 'FAIL'}] {label}: got {got!r}, want {want!r}")
        if not ok:
            failures += 1

    print("Encoding:")
    check("XZZXI -> int", pauli_string_to_row("XZZXI", 5), 0b1001001100)
    check("IXZZX -> int", pauli_string_to_row("IXZZX", 5), 0b0100100110)
    check("XIXZZ -> int", pauli_string_to_row("XIXZZ", 5), 0b1010000011)
    check("ZXIXZ -> int", pauli_string_to_row("ZXIXZ", 5), 0b0101010001)
    check("round trip", row_to_pauli_string(pauli_string_to_row("XZZXI", 5), 5),
          "XZZXI")
    check("Y round trip", row_to_pauli_string(pauli_string_to_row("YIYZX", 5), 5),
          "YIYZX")

    print("Weights:")
    pw = build_pauli_weight_table(1 << 10, 5)
    check("Pauli weight of XZZXI", pw[pauli_string_to_row("XZZXI", 5)], 4)
    check("binary weight of XZZXI",
          build_popcount_table(1 << 10)[pauli_string_to_row("XZZXI", 5)], 4)

    print("GF(2) / target probe:")
    spec513 = CodeSpec(5, 1, 3)
    probe = build_target_probe(spec513, "span", None)
    check("stabilizer group size", len(probe.span), 16)
    check("stabilizer group dimension", probe.dim, 4)
    check("expected hits in full space", probe.expected_hits_full_space(), 20160)
    check("canonical generators recognised", probe.test(probe.generators), True)
    g = probe.generators
    check("reordered generators still recognised",
          probe.test((g[1], g[0], g[3], g[2])), True)
    check("products of generators still recognised",
          probe.test((g[0], g[0] ^ g[1], g[2], g[3] ^ g[0])), True)
    check("dependent set rejected", probe.test((g[0], g[1], g[2], g[0] ^ g[1])),
          False)
    check("off-code row rejected", probe.test((g[0], g[1], g[2], 0b1111111111)),
          False)
    probe_exact = build_target_probe(spec513, "exact", None)
    check("exact mode matches generator matrix", probe_exact.test(g), True)
    check("exact mode rejects reorder", probe_exact.test((g[1], g[0], g[2], g[3])),
          False)

    print("Cyclic shift:")
    check("shift(XZZXI,1) == IXZZX",
          cyclic_shift_row(pauli_string_to_row("XZZXI", 5), 5, 1),
          pauli_string_to_row("IXZZX", 5))

    print("Backend agreement on a tiny spec [[3,1]] (4096 candidates):")
    tiny = CodeSpec(3, 1)
    tiny_target = build_target_probe(tiny, "exact", ("XZI", "IXZ"))
    outs = {}
    for backend in (("pure", "numpy") if _np is not None else ("pure",)):
        plan = build_plan(tiny, NoConstraint(), "filter", backend,
                          tiny_target, False, False)
        r = run_search(plan)
        outs[backend] = (r.examined, r.retained, r.checksum, r.target_hits)
        print(f"  {backend}: examined={r.examined} retained={r.retained} "
              f"checksum={r.checksum} hits={r.target_hits}")
    check("tiny examined", outs["pure"][0], tiny.space_size)
    # every bit is 1 in half the space: checksum = rows * cols/2 * S_0
    check("tiny checksum invariant", outs["pure"][2],
          tiny.rows * (tiny.cols // 2) * tiny.space_size)
    check("tiny exact-target hits", outs["pure"][3], 1)
    if "numpy" in outs:
        check("numpy == pure", outs["numpy"], outs["pure"])

    print("Weight-spectrum sanity on [[3,1]] with W={2} (filter vs generate):")
    wc = WeightSpectrumConstraint([2], "pauli")
    wc.bind(tiny)
    p_filter = build_plan(tiny, wc, "filter", "pure", None, False, False)
    r_filter = run_search(p_filter)
    p_gen = build_plan(tiny, wc, "generate", "pure", None, False, False)
    r_gen = run_search(p_gen)
    check("filter retained == generate examined", r_filter.retained,
          r_gen.examined)
    check("generate examined == S_R", r_gen.examined, wc.retained_size(tiny))

    # =================================================================
    # UNIVERSAL STABILIZER CONSTRAINT
    # =================================================================
    print("Symplectic form (condition 2 in bit form):")
    g513 = probe.generators
    check("XZZXI vs IXZZX commute", symplectic_product(g513[0], g513[1], 5), 0)
    check("[[5,1,3]] generators pairwise commute",
          rows_pairwise_commute(g513, 5), True)
    x0 = pauli_string_to_row("XIIII", 5)
    z0 = pauli_string_to_row("ZIIII", 5)
    check("XIIII vs ZIIII anticommute", symplectic_product(x0, z0, 5), 1)
    check("XIIII vs IZIII commute", symplectic_product(
        x0, pauli_string_to_row("IZIII", 5), 5), 0)
    check("every row commutes with itself",
          all(symplectic_product(v, v, 5) == 0 for v in range(1024)), True)
    check("symplectic form is symmetric",
          all(symplectic_product(a, b, 5) == symplectic_product(b, a, 5)
              for a in range(0, 1024, 37) for b in range(0, 1024, 53)), True)

    print("Bit form agrees with the literal matrix identity H_X H_Z^T = H_Z H_X^T:")

    def literal_commutation(rows, n):
        """Direct GF(2) evaluation of H_X H_Z^T - H_Z H_X^T, entry by entry."""
        mask = (1 << n) - 1
        m = len(rows)
        HX = [(r >> n) & mask for r in rows]
        HZ = [r & mask for r in rows]
        for i in range(m):
            for j in range(m):
                lhs = (HX[i] & HZ[j]).bit_count() & 1      # (H_X H_Z^T)_ij
                rhs = (HZ[i] & HX[j]).bit_count() & 1      # (H_Z H_X^T)_ij
                if lhs != rhs:
                    return False
        return True

    # Deterministic and exhaustive on small specs -- no sampling anywhere.
    for (nn, mm) in ((2, 2), (2, 3), (2, 4), (3, 2)):
        alpha = 1 << (2 * nn)
        agree = all(rows_pairwise_commute(rs, nn) == literal_commutation(rs, nn)
                    for rs in product(range(alpha), repeat=mm))
        check(f"bit test == matrix identity, all {alpha ** mm} "
              f"{mm}x{2 * nn} matrices", agree, True)
    # And on a deterministic grid of five-qubit rows.
    grid = list(range(0, 1024, 87))
    agree5 = all(rows_pairwise_commute(rs, 5) == literal_commutation(rs, 5)
                 for rs in product(grid, repeat=4))
    check(f"bit test == matrix identity, {len(grid) ** 4} five-qubit cases",
          agree5, True)

    print("Universal constraint predicate on [[5,1,3]]:")
    uni = UniversalStabilizerConstraint()
    uni.bind(spec513)
    check("canonical generator matrix accepted", uni.accepts(g513), True)
    check("reordered generators accepted",
          uni.accepts((g513[2], g513[0], g513[3], g513[1])), True)
    check("dependent set rejected (rank 3)",
          uni.accepts((g513[0], g513[1], g513[2], g513[0] ^ g513[1])), False)
    check("zero row rejected (rank < 4)",
          uni.accepts((0, g513[1], g513[2], g513[3])), False)
    check("anticommuting pair rejected",
          uni.accepts((x0, z0, g513[2], g513[3])), False)
    check("expected count formula, n=5 m=4", expected_universal_count(5, 4),
          1023 * 510 * 252 * 120)
    check("expected count value", expected_universal_count(5, 4), 15777115200)

    print("Formula vs brute force on toy specs:")
    for (nn, kk) in ((2, 0), (3, 1), (3, 0)):
        ts = CodeSpec(nn, kk)
        u = UniversalStabilizerConstraint()
        u.bind(ts)
        brute = sum(1 for cand in product(range(ts.row_alphabet),
                                          repeat=ts.rows)
                    if u.accepts(cand))
        check(f"brute force [[{nn},{kk}]] ({ts.space_size} candidates)",
              brute, expected_universal_count(nn, ts.rows))

    print("All 20160 target matrices satisfy the universal constraint:")
    span_sorted = sorted(probe.span)
    both = sum(1 for c in product(span_sorted, repeat=4)
               if probe.test(c) and uni.accepts(c))
    check("target hits that also pass the universal constraint", both, 20160)

    if _np is not None:
        print("Universal engine: numpy vs pure, both evaluation strategies:")
        for (nn, kk) in ((2, 0), (3, 1), (3, 0)):
            ts = CodeSpec(nn, kk)
            u = UniversalStabilizerConstraint()
            ref = None
            outs = {}
            for be, ev in (("pure", "hierarchical"),
                           ("numpy", "hierarchical"),
                           ("numpy", "per-candidate")):
                pl = build_plan(ts, u, "filter", be, None, False, False, ev)
                r = run_search(pl)
                outs[(be, ev)] = (r.examined, r.retained, r.checksum)
                if ref is None:
                    ref = outs[(be, ev)]
            check(f"[[{nn},{kk}]] examined == S_0", ref[0], ts.space_size)
            check(f"[[{nn},{kk}]] retained == formula", ref[1],
                  expected_universal_count(nn, ts.rows))
            check(f"[[{nn},{kk}]] numpy hierarchical == pure",
                  outs[("numpy", "hierarchical")], ref)
            check(f"[[{nn},{kk}]] numpy per-candidate == pure",
                  outs[("numpy", "per-candidate")], ref)

        print("Universal engine with 2 outer rows (same structure as [[5,1,3]]):")
        ts = CodeSpec(4, 0)
        u = UniversalStabilizerConstraint()
        pl = build_plan(ts, u, "filter", "numpy", None, False, False)
        r = run_search(pl)
        check("[[4,0]] examined == S_0", r.examined, ts.space_size)
        check("[[4,0]] retained == formula", r.retained,
              expected_universal_count(4, 4))
        check("[[4,0]] blocks == alphabet^2", r.blocks_total,
              ts.row_alphabet ** 2)

        print("Universal engine preserves target recognition:")
        # The toy target must itself satisfy the universal constraint,
        # otherwise it is correctly rejected and never reported. XXI/ZZI
        # commute and are independent; XZI/IXZ anticommute and would not.
        ts = CodeSpec(3, 1)
        u = UniversalStabilizerConstraint()
        u.bind(ts)
        toy = (pauli_string_to_row("XXI", 3), pauli_string_to_row("ZZI", 3))
        check("toy target commutes", rows_pairwise_commute(toy, 3), True)
        check("toy target passes the universal constraint", u.accepts(toy), True)
        check("XZI/IXZ anticommute, so they are correctly rejected",
              u.accepts((pauli_string_to_row("XZI", 3),
                         pauli_string_to_row("IXZ", 3))), False)
        tp = build_target_probe(ts, "exact", ("XXI", "ZZI"))
        pl = build_plan(ts, u, "filter", "numpy", tp, False, False)
        r = run_search(pl)
        check("[[3,1]] exact target still found once", r.target_hits, 1)
        pl2 = build_plan(ts, u, "filter", "pure", tp, False, False)
        r2 = run_search(pl2)
        check("[[3,1]] pure agrees on target hits", r2.target_hits, r.target_hits)
        check("[[3,1]] pure agrees on target index",
              r2.first_target_index, r.first_target_index)
        check("[[3,1]] target index is the lexicographic position",
              r.first_target_index, toy[0] * ts.row_alphabet + toy[1])

        # Span-mode recognition through the constrained engine.
        tp_span = build_target_probe(ts, "span", ("XXI", "ZZI"))
        pl3 = build_plan(ts, u, "filter", "numpy", tp_span, False, False)
        r3 = run_search(pl3)
        pl4 = build_plan(ts, u, "filter", "pure", tp_span, False, False)
        r4 = run_search(pl4)
        check("[[3,1]] span-mode hits == ordered bases of a 2-dim space",
              r3.target_hits, (4 - 1) * (4 - 2))
        check("[[3,1]] span-mode numpy == pure", r3.target_hits, r4.target_hits)

    print()
    print("SELF-TEST FAILURES:" if failures else "SELF-TEST: all checks passed.",
          failures if failures else "")
    return 1 if failures else 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Exhaustive stabilizer parity-check-matrix search "
                    "runtime benchmark (no sampling, no early exit).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--n", type=int, default=5, help="number of physical qubits")
    p.add_argument("--k", type=int, default=1, help="number of logical qubits")
    p.add_argument("--d", type=int, default=3, help="distance, label only")
    p.add_argument("--constraint", default="weight",
                   help="comma-separated constraint names: universal "
                        "(default here: rank(H)=n-k AND H_X H_Z^T = H_Z H_X^T), "
                        "none (reproduces the unconstrained baseline), weight, "
                        "cyclic. Constraints are combined ONLY when explicitly "
                        "listed together.")
    p.add_argument("--constraint-eval", default="hierarchical",
                   choices=("hierarchical", "per-candidate"),
                   help="hierarchical (default): reject a candidate with the "
                        "cheapest sub-condition it violates, short-circuiting "
                        "on shared row prefixes. per-candidate: materialise "
                        "the full acceptance mask for every block so the "
                        "constraint cost stays Theta(S_0). These measure "
                        "different things; report which you used.")
    p.add_argument("--weights", default="",
                   help="weight spectrum for --constraint weight, "
                        "e.g. '4' or '2,4'. Never assumed.")
    p.add_argument("--weight-mode", default="pauli", choices=("pauli", "binary"),
                   help="row weight definition (default: pauli weight)")
    p.add_argument("--cyclic-shift", type=int, default=1,
                   help="shift for --constraint cyclic (default 1)")
    p.add_argument("--mode", default="filter", choices=("filter", "generate"),
                   help="filter: traverse the full space and test every "
                        "candidate. generate: enumerate only the retained "
                        "subspace (this is the mode whose runtime should track "
                        "S_R/S_0).")
    p.add_argument("--backend", default="auto", choices=("auto", "numpy", "pure"),
                   help="enumeration engine; use the SAME one for baseline and "
                        "constrained runs")
    p.add_argument("--repetitions", type=int, default=1,
                   help="number of independent timed repetitions")
    p.add_argument("--r0-start", type=int, default=0,
                   help="BOUNDED SAMPLE: first value of row 0 to enumerate. "
                        "The sample is a contiguous deterministic region of "
                        "the same lexicographic enumeration the full search "
                        "uses; nothing is randomised and the candidate "
                        "distribution inside the region is unchanged.")
    p.add_argument("--r0-count", type=int, default=None,
                   help="BOUNDED SAMPLE: how many values of row 0 to "
                        "enumerate (default: all of them = the full search). "
                        "Sample candidate count = r0-count * alphabet**"
                        "(rows-1); the report extrapolates the full-search "
                        "runtime by the exact ratio of candidate counts.")
    p.add_argument("--work", type=int, default=1,
                   help="per-candidate processing units. 1 (default) is the "
                        "original workload -- the sum of the matrix's row "
                        "weights. Each extra unit computes the weight of one "
                        "element of the candidate's row space. Raising it "
                        "makes per-candidate processing dominate the fixed "
                        "constraint-evaluation cost, driving the observed "
                        "runtime ratio toward S_R/S_0. MUST match the value "
                        "used for the unconstrained baseline run.")
    p.add_argument("--target", default="span", choices=("span", "exact", "off"),
                   help="target-code recognition mode (default: span)")
    p.add_argument("--target-generators", default="",
                   help="comma-separated Pauli strings overriding the default "
                        "target generators")
    p.add_argument("--announce-all", action="store_true",
                   help="print 'Code Found' for EVERY occurrence instead of "
                        "only the first (adds I/O inside the timed region)")
    p.add_argument("--no-progress", action="store_true",
                   help="suppress the stderr progress line")
    p.add_argument("--dry-run", action="store_true",
                   help="print the configuration and space sizes, enumerate "
                        "nothing")
    p.add_argument("--self-test", action="store_true",
                   help="run cheap correctness checks (does not run the search)")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)

    if args.self_test:
        return self_test()

    spec = CodeSpec(args.n, args.k, args.d)
    spec.validate()

    backend = args.backend
    if backend == "auto":
        backend = "numpy" if _np is not None else "pure"
    if backend == "numpy" and _np is None:
        print("numpy not available; falling back to the pure backend",
              file=sys.stderr)
        backend = "pure"

    gens = ([s.strip() for s in args.target_generators.split(",") if s.strip()]
            if args.target_generators else None)
    target = build_target_probe(spec, args.target, gens)

    constraint = build_constraint(args.constraint.split(","), args)

    if args.dry_run:
        dry_run_report(spec, constraint, args.mode, backend, target)
        return 0

    universal = isinstance(constraint, WeightSpectrumConstraint)
    if universal and args.mode != "filter":
        print("This benchmark enumerates the full candidate space and tests "
              "each candidate; use --mode filter.", file=sys.stderr)
        return 2

    # The numpy backend cannot vectorise an arbitrary non-row-local
    # constraint; fall back rather than silently changing what is measured.
    # The universal constraint is the exception: search_numpy_universal()
    # handles it natively, so it keeps the numpy backend.
    if backend == "numpy" and not constraint.is_trivial and not universal:
        binder = getattr(constraint, "bind", None)
        if binder is not None:
            binder(spec)
        local = constraint.allowed_row_values(spec)
        probe_plan_ok = (local is not None
                         and all(list(v) == list(local[0]) for v in local))
        if not probe_plan_ok:
            print(f"NOTE: constraint {constraint.name!r} is not row-local; "
                  f"switching to the pure backend. Re-run the BASELINE with "
                  f"--backend pure so the comparison stays fair.",
                  file=sys.stderr)
            backend = "pure"

    print(f"Starting exhaustive traversal: {spec.label}, "
          f"{spec.rows}x{spec.cols} binary matrix, backend={backend}, "
          f"mode={args.mode}", file=sys.stderr)
    print(f"S_0 = {spec.space_size:,} candidate matrices. "
          f"No early termination.", file=sys.stderr)
    if universal:
        print(f"Constraint: {constraint.describe()}", file=sys.stderr)

    results = run_benchmark(spec, constraint, args.mode, backend, target,
                            max(1, args.repetitions), args.announce_all,
                            not args.no_progress, args.constraint_eval,
                            max(1, args.work),
                            args.r0_start, args.r0_count)
    report(spec, constraint, args.mode, backend, target, results,
           args.constraint_eval, max(1, args.work))
    sample_report(spec, results, args.r0_start, args.r0_count, max(1, args.work))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
