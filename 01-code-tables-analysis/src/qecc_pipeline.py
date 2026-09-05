#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qecc_pipeline.py
================

Automated retrieval, analysis, export and visualisation of quantum stabilizer
codes from Markus Grassl's code tables at

    https://www.codetables.de/QECC/QECC.php?q=<q^2>&n=<n>&k=<k>

The pipeline runs in four phases:

  1. FETCH    - download every (q, n, k) page in range into a local SQLite
                cache.  Rate limited, retrying, resumable, never downloads the
                same URL twice.
  2. PARSE    - parse bounds + stabilizer (parity check) matrices out of the
                cached HTML, compute every requested quantity, and store the
                results back into SQLite.  Re-runnable offline.
  3. CLASSIFY - using the *complete* upper-bound grid for a given q, label each
                code "Semi-Perfect Quantum Code" or "Normal Code".
  4. EXPORT   - one .xlsx per q (AutoFilter enabled) + 9 scatter plots per q
                + skip log + run summary.

Nothing is invented.  If a page has no usable parity-check matrix the code is
skipped and the reason recorded in the skip log.

--------------------------------------------------------------------------
Conventions / interpretation notes (see README notes printed by --explain)
--------------------------------------------------------------------------
* Row weight  W_i is the SYMPLECTIC weight of stabilizer row i: the number of
  coordinates j in 1..n for which (X_j, Z_j) != (0, 0).  A stabilizer matrix
  row is a vector of length n over GF(q^2) written as an n+n binary/GF(q)
  symplectic pair, so this is the number of non-zero entries of that length-n
  vector.  This keeps W_i <= n so that delta = W_i/n lies in [0,1] and the
  q^2-ary entropy H_{q^2}(delta) is defined.  Use --weight-mode raw to count
  non-zero entries across all 2n printed symbols instead.
* X = sum_i C(n, W_i) * (q^2 - 1)^{W_i} is summed over the DISTINCT weights
  appearing in the weight spectrum (the frequencies f_i appear only in the
  average-weight formula).  Use --x-sum-mode rows to weight each term by f_i.
* R  = log_q(X) / (2n)            (logarithm base q, as specified)
* H_{q^2}(x) = x log_{q^2}(q^2-1) - x log_{q^2}(x) - (1-x) log_{q^2}(1-x)
* C  = H_{q^2}(delta_max) - R,    C' = H_{q^2}(delta_avg) - R

Author: generated for the Code Tables analysis project.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html as html_mod
import json
import logging
import math
import os
import random
import re
import sqlite3
import sys
import threading
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# --------------------------------------------------------------------------
# Third party
# --------------------------------------------------------------------------
try:
    import requests
    from requests.adapters import HTTPAdapter
except ImportError as exc:  # pragma: no cover
    sys.exit("Missing dependency: %s  (pip install requests)" % exc)

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover
    def tqdm(iterable=None, **kwargs):  # type: ignore
        return iterable if iterable is not None else None

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

try:
    import seaborn as sns
except ImportError:  # pragma: no cover
    sns = None


# ==========================================================================
# Configuration
# ==========================================================================

BASE_URL = "https://www.codetables.de/QECC/QECC.php"

#: Maximum length n published by codetables.de for each prime power q.
#: (verified by probing: q=2 -> n<=256, every other q -> n<=100; n+1 gives
#: "wrong input".  q^2 = 81 and 121 are not served at all, so q = 2, 3, 4,
#: 5, 7, 8 is the complete set of field sizes the table offers.)
DEFAULT_N_MAX: Dict[int, int] = {2: 256, 3: 100, 4: 100, 5: 100,
                                 7: 100, 8: 100}

SEMI_PERFECT = "Semi-Perfect Quantum Code"
NORMAL = "Normal Code"

EXCEL_COLUMNS = [
    "q",
    "n",
    "k",
    "Code Parameters",
    "Type of Code",
    "k/n",
    "Weight Spectrum",
    "Average Weight",
    "Max Weight",
    "R",
    "C",
    "C'",
]


@dataclass
class Config:
    out_dir: str
    q_values: Tuple[int, ...] = (2, 3, 4, 5, 7, 8)
    n_max: Dict[int, int] = field(default_factory=lambda: dict(DEFAULT_N_MAX))
    n_min: int = 1
    workers: int = 6
    rps: float = 6.0
    timeout: float = 30.0
    max_attempts: int = 5
    weight_mode: str = "symplectic"      # symplectic | raw
    x_sum_mode: str = "distinct"         # distinct | rows
    avg_denominator: str = "rows"        # rows | n-k
    skip_fetch: bool = False
    force_reparse: bool = False
    dpi: int = 160

    @property
    def db_path(self) -> str:
        return os.path.join(self.out_dir, "codetables_cache.sqlite")

    @property
    def plot_dir(self) -> str:
        return os.path.join(self.out_dir, "plots")

    @property
    def log_dir(self) -> str:
        return os.path.join(self.out_dir, "logs")


log = logging.getLogger("qecc")


# ==========================================================================
# Small utilities
# ==========================================================================

def ln_bigint(value: int) -> float:
    """Natural logarithm of an arbitrarily large positive integer."""
    if value <= 0:
        raise ValueError("ln of non-positive integer")
    bits = value.bit_length()
    if bits <= 1000:
        return math.log(value)
    shift = bits - 900
    return math.log(value >> shift) + shift * math.log(2.0)


def entropy_q2(x: float, q2: int) -> float:
    """q^2-ary entropy function H_{q^2}(x), with the correct limits at 0 and 1."""
    if x <= 0.0:
        return 0.0
    log_q2 = math.log(q2)
    base_term = math.log(q2 - 1) / log_q2
    if x >= 1.0:
        return base_term
    return (
        x * base_term
        - x * math.log(x) / log_q2
        - (1.0 - x) * math.log(1.0 - x) / log_q2
    )


def format_spectrum(spectrum: Dict[int, int]) -> str:
    """'{3: 2, 4: 3}' with weights in ascending order."""
    return "{" + ", ".join("%d: %d" % (w, spectrum[w]) for w in sorted(spectrum)) + "}"


def chunks(seq: Sequence, size: int) -> Iterable[Sequence]:
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


# ==========================================================================
# Rate limiting
# ==========================================================================

class RateLimiter:
    """Simple thread-safe evenly-spaced rate limiter with adaptive back-off."""

    def __init__(self, rps: float) -> None:
        self._base_interval = 1.0 / max(rps, 0.05)
        self._interval = self._base_interval
        self._next_slot = time.monotonic()
        self._lock = threading.Lock()

    def acquire(self) -> None:
        with self._lock:
            now = time.monotonic()
            if self._next_slot < now:
                self._next_slot = now
            wait = self._next_slot - now
            self._next_slot += self._interval
        if wait > 0:
            time.sleep(wait)

    def penalise(self, factor: float = 2.0, cap: float = 4.0) -> None:
        """Server pushed back - slow down (bounded)."""
        with self._lock:
            self._interval = min(self._interval * factor, self._base_interval * cap)
            log.warning("Rate limit backed off to %.2f req/s", 1.0 / self._interval)

    def relax(self) -> None:
        with self._lock:
            if self._interval > self._base_interval:
                self._interval = max(self._base_interval, self._interval * 0.9)


# ==========================================================================
# Persistent store (cache + results + skip log), single SQLite file
# ==========================================================================

class Store:
    """SQLite-backed page cache and result store.  Resumable across runs."""

    def __init__(self, path: str) -> None:
        self.path = path
        self.conn = sqlite3.connect(path, check_same_thread=False)
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA synchronous=NORMAL")
        self._lock = threading.Lock()
        self._create_schema()

    def _create_schema(self) -> None:
        cur = self.conn.cursor()
        cur.execute(
            """CREATE TABLE IF NOT EXISTS pages (
                   q INTEGER NOT NULL,
                   n INTEGER NOT NULL,
                   k INTEGER NOT NULL,
                   url TEXT NOT NULL,
                   fetch_status TEXT NOT NULL,   -- ok | invalid | error
                   http_status INTEGER,
                   body BLOB,                    -- gzipped utf-8 html
                   note TEXT,
                   fetched_at TEXT NOT NULL,
                   PRIMARY KEY (q, n, k))"""
        )
        cur.execute(
            """CREATE TABLE IF NOT EXISTS bounds (
                   q INTEGER NOT NULL,
                   n INTEGER NOT NULL,
                   k INTEGER NOT NULL,
                   lower_bound INTEGER,
                   upper_bound INTEGER,
                   has_matrix INTEGER NOT NULL DEFAULT 0,
                   PRIMARY KEY (q, n, k))"""
        )
        cur.execute(
            """CREATE TABLE IF NOT EXISTS codes (
                   q INTEGER NOT NULL,
                   n INTEGER NOT NULL,
                   k INTEGER NOT NULL,
                   d INTEGER NOT NULL,
                   lower_bound INTEGER,
                   upper_bound INTEGER,
                   label TEXT NOT NULL,
                   n_rows INTEGER NOT NULL,
                   spectrum TEXT NOT NULL,
                   avg_weight REAL NOT NULL,
                   max_weight INTEGER NOT NULL,
                   x_value TEXT NOT NULL,
                   rate_R REAL NOT NULL,
                   delta_avg REAL NOT NULL,
                   delta_max REAL NOT NULL,
                   H_delta_avg REAL NOT NULL,
                   H_delta_max REAL NOT NULL,
                   C REAL NOT NULL,
                   C_prime REAL NOT NULL,
                   n_constructions INTEGER NOT NULL,
                   matrix_label TEXT,
                   PRIMARY KEY (q, n, k))"""
        )
        cur.execute(
            """CREATE TABLE IF NOT EXISTS skips (
                   q INTEGER NOT NULL,
                   n INTEGER NOT NULL,
                   k INTEGER NOT NULL,
                   reason TEXT NOT NULL,
                   detail TEXT,
                   PRIMARY KEY (q, n, k))"""
        )
        self.conn.commit()

    # ---------------- page cache ----------------

    def cached_keys(self, q: int) -> set:
        cur = self.conn.execute(
            "SELECT n, k FROM pages WHERE q=? AND fetch_status IN ('ok','invalid')", (q,)
        )
        return set(cur.fetchall())

    def store_pages(self, rows: List[tuple]) -> None:
        """rows: (q, n, k, url, fetch_status, http_status, body, note, ts)"""
        with self._lock:
            self.conn.executemany(
                "INSERT OR REPLACE INTO pages "
                "(q,n,k,url,fetch_status,http_status,body,note,fetched_at) "
                "VALUES (?,?,?,?,?,?,?,?,?)",
                rows,
            )
            self.conn.commit()

    def iter_pages(self, q: int):
        cur = self.conn.execute(
            "SELECT n, k, fetch_status, body FROM pages WHERE q=? ORDER BY n, k", (q,)
        )
        while True:
            batch = cur.fetchmany(256)
            if not batch:
                return
            for n, k, status, body in batch:
                text = gzip.decompress(body).decode("utf-8", "replace") if body else ""
                yield n, k, status, text

    def count_pages(self, q: int) -> int:
        return self.conn.execute(
            "SELECT COUNT(*) FROM pages WHERE q=?", (q,)
        ).fetchone()[0]

    # ---------------- results ----------------

    def clear_results(self, q: int) -> None:
        with self._lock:
            self.conn.execute("DELETE FROM bounds WHERE q=?", (q,))
            self.conn.execute("DELETE FROM codes WHERE q=?", (q,))
            self.conn.execute("DELETE FROM skips WHERE q=?", (q,))
            self.conn.commit()

    def has_results(self, q: int) -> bool:
        return self.conn.execute(
            "SELECT COUNT(*) FROM bounds WHERE q=?", (q,)
        ).fetchone()[0] > 0

    def store_bounds(self, rows: List[tuple]) -> None:
        with self._lock:
            self.conn.executemany(
                "INSERT OR REPLACE INTO bounds (q,n,k,lower_bound,upper_bound,has_matrix)"
                " VALUES (?,?,?,?,?,?)", rows)
            self.conn.commit()

    def store_codes(self, rows: List[tuple]) -> None:
        with self._lock:
            self.conn.executemany(
                "INSERT OR REPLACE INTO codes VALUES (" + ",".join(["?"] * 21) + ")", rows)
            self.conn.commit()

    def store_skips(self, rows: List[tuple]) -> None:
        with self._lock:
            self.conn.executemany(
                "INSERT OR REPLACE INTO skips (q,n,k,reason,detail) VALUES (?,?,?,?,?)",
                rows)
            self.conn.commit()

    def upper_bound_grid(self, q: int) -> Dict[Tuple[int, int], int]:
        cur = self.conn.execute(
            "SELECT n, k, upper_bound FROM bounds WHERE q=? AND upper_bound IS NOT NULL",
            (q,))
        return {(n, k): ub for n, k, ub in cur.fetchall()}

    def codes_dataframe(self, q: int) -> pd.DataFrame:
        return pd.read_sql_query(
            "SELECT * FROM codes WHERE q=? ORDER BY n, k", self.conn, params=(q,))

    def skips_dataframe(self, q: Optional[int] = None) -> pd.DataFrame:
        if q is None:
            return pd.read_sql_query(
                "SELECT * FROM skips ORDER BY q, n, k", self.conn)
        return pd.read_sql_query(
            "SELECT * FROM skips WHERE q=? ORDER BY n, k", self.conn, params=(q,))

    def close(self) -> None:
        self.conn.close()


# ==========================================================================
# Phase 1 - fetching
# ==========================================================================

_thread_local = threading.local()


def _session(cfg: Config) -> requests.Session:
    sess = getattr(_thread_local, "session", None)
    if sess is None:
        sess = requests.Session()
        # We run our own retry / back-off loop, so the adapter must not retry.
        adapter = HTTPAdapter(pool_connections=cfg.workers + 2,
                              pool_maxsize=cfg.workers + 2,
                              max_retries=0)
        sess.mount("https://", adapter)
        sess.mount("http://", adapter)
        sess.headers.update({
            "User-Agent": ("QECC-research-bot/1.0 (academic study of quantum code "
                           "weight spectra; contact via codetables.de)"),
            "Accept": "text/html,application/xhtml+xml",
            "Accept-Encoding": "gzip, deflate",
            "Connection": "keep-alive",
        })
        _thread_local.session = sess
    return sess


WRONG_INPUT_RE = re.compile(r"wrong\s+input", re.I)


def fetch_one(cfg: Config, limiter: RateLimiter, q: int, n: int, k: int) -> tuple:
    """Download a single page.  Returns a row for the `pages` table."""
    params = {"q": q * q, "n": n, "k": k}
    url = "%s?q=%d&n=%d&k=%d" % (BASE_URL, q * q, n, k)
    last_error = ""

    for attempt in range(1, cfg.max_attempts + 1):
        limiter.acquire()
        try:
            resp = _session(cfg).get(BASE_URL, params=params, timeout=cfg.timeout)
        except Exception as exc:                       # noqa: BLE001 - network layer
            last_error = "%s: %s" % (type(exc).__name__, exc)
            log.debug("q=%d n=%d k=%d attempt %d failed: %s", q, n, k, attempt, last_error)
            time.sleep(min(2.0 ** attempt, 30.0) * (0.5 + random.random()))
            continue

        if resp.status_code in (429, 503):
            limiter.penalise()
            last_error = "HTTP %d" % resp.status_code
            time.sleep(min(2.0 ** attempt, 30.0))
            continue
        if resp.status_code >= 500:
            last_error = "HTTP %d" % resp.status_code
            time.sleep(min(2.0 ** attempt, 30.0) * (0.5 + random.random()))
            continue
        if resp.status_code != 200:
            return (q, n, k, url, "error", resp.status_code, None,
                    "HTTP %d" % resp.status_code, _now())

        text = resp.text
        limiter.relax()
        if WRONG_INPUT_RE.search(text) and "Bounds on [[" not in text:
            return (q, n, k, url, "invalid", 200, None, "wrong input (out of range)",
                    _now())
        blob = gzip.compress(text.encode("utf-8"), 6)
        return (q, n, k, url, "ok", 200, blob, None, _now())

    return (q, n, k, url, "error", None, None,
            "giving up after %d attempts: %s" % (cfg.max_attempts, last_error), _now())


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")


def fetch_phase(cfg: Config, store: Store) -> Dict[str, int]:
    """Make sure every (q, n, k) page in range is present in the cache."""
    stats = {"requested": 0, "already_cached": 0, "downloaded": 0, "errors": 0}
    limiter = RateLimiter(cfg.rps)

    for q in cfg.q_values:
        n_max = cfg.n_max[q]
        all_keys = [(n, k) for n in range(cfg.n_min, n_max + 1) for k in range(0, n + 1)]
        cached = store.cached_keys(q)
        todo = [(n, k) for (n, k) in all_keys if (n, k) not in cached]
        stats["requested"] += len(all_keys)
        stats["already_cached"] += len(all_keys) - len(todo)

        if not todo:
            log.info("q=%d: all %d pages already cached", q, len(all_keys))
            continue

        log.info("q=%d: %d pages in range, %d cached, downloading %d",
                 q, len(all_keys), len(all_keys) - len(todo), len(todo))

        bar = tqdm(total=len(todo), desc="fetch q=%d" % q, unit="pg",
                   ascii=True, dynamic_ncols=False, mininterval=5.0)
        with ThreadPoolExecutor(max_workers=cfg.workers) as pool:
            for batch in chunks(todo, 240):
                rows = list(pool.map(
                    lambda nk: fetch_one(cfg, limiter, q, nk[0], nk[1]), batch))
                store.store_pages(rows)
                for r in rows:
                    if r[4] == "error":
                        stats["errors"] += 1
                        log.warning("fetch error q=%d n=%d k=%d: %s", q, r[1], r[2], r[7])
                    else:
                        stats["downloaded"] += 1
                if bar is not None:
                    bar.update(len(batch))
        if bar is not None:
            bar.close()

    return stats


# ==========================================================================
# Phase 2 - parsing
# ==========================================================================

BOUNDS_HDR_RE = re.compile(
    r"Bounds\s+on\s+\[\[\s*(\d+)\s*,\s*(\d+)\s*\]\]", re.I)
LOWER_RE = re.compile(r"lower\s*bound:\s*</TD>\s*<TD[^>]*>(.*?)</TD>", re.I | re.S)
UPPER_RE = re.compile(r"upper\s*bound:\s*</TD>\s*<TD[^>]*>(.*?)</TD>", re.I | re.S)
PRE_RE = re.compile(r"<PRE>(.*?)</PRE>", re.I | re.S)
CONSTRUCTION_HDR_RE = re.compile(
    r"Construction\s+of\s+an?\s+\[\[\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\]\]", re.I)
MATRIX_LABEL_RE = re.compile(
    r"^[ \t]*((?:stabilizer|generator|parity[ \-]?check)\s+matrix)\s*:\s*$", re.I)
ROW_RE = re.compile(r"^\s*\[([^\[\]|]*)\|([^\[\]|]*)\]\s*$")
INT_RE = re.compile(r"-?\d+")
MISSING_CONSTRUCTION_RE = re.compile(r"Missing\s+Construction", re.I)


@dataclass
class MatrixCandidate:
    d: int
    label: str
    rows: List[Tuple[List[str], List[str]]]


@dataclass
class ParsedPage:
    lower_bound: Optional[int] = None
    upper_bound: Optional[int] = None
    candidates: List[MatrixCandidate] = field(default_factory=list)
    missing_construction: bool = False
    has_pre: bool = False
    note: str = ""


def _first_int(text: str) -> Optional[int]:
    plain = re.sub(r"<[^>]*>", " ", text)
    m = INT_RE.search(plain)
    return int(m.group(0)) if m else None


def _extract_matrix_blocks(text: str) -> List[Tuple[str, List[Tuple[List[str], List[str]]]]]:
    """Find every '<something> matrix:' block and its bracketed rows."""
    lines = text.splitlines()
    blocks: List[Tuple[str, List[Tuple[List[str], List[str]]]]] = []
    i = 0
    while i < len(lines):
        m = MATRIX_LABEL_RE.match(lines[i])
        if not m:
            i += 1
            continue
        label = " ".join(m.group(1).split()).lower()
        j = i + 1
        rows: List[Tuple[List[str], List[str]]] = []
        # allow blank lines before the first row
        while j < len(lines) and not lines[j].strip() and not rows:
            j += 1
        while j < len(lines):
            line = lines[j]
            rm = ROW_RE.match(line)
            if rm:
                left = rm.group(1).split()
                right = rm.group(2).split()
                rows.append((left, right))
                j += 1
                continue
            break
        if rows:
            blocks.append((label, rows))
        i = max(j, i + 1)
    return blocks


def parse_page(text: str, n: int, k: int) -> ParsedPage:
    """Parse one QECC page into bounds + all matrix candidates."""
    out = ParsedPage()

    hdr = BOUNDS_HDR_RE.search(text)
    if hdr is None:
        out.note = "no bounds header"
        return out
    if int(hdr.group(1)) != n or int(hdr.group(2)) != k:
        out.note = "bounds header [[%s,%s]] does not match requested [[%d,%d]]" % (
            hdr.group(1), hdr.group(2), n, k)
        return out

    lm = LOWER_RE.search(text)
    um = UPPER_RE.search(text)
    out.lower_bound = _first_int(lm.group(1)) if lm else None
    out.upper_bound = _first_int(um.group(1)) if um else None

    pre = PRE_RE.search(text)
    if pre is None:
        out.note = "no <PRE> construction block"
        return out
    out.has_pre = True
    body = html_mod.unescape(pre.group(1))

    if MISSING_CONSTRUCTION_RE.search(body):
        out.missing_construction = True

    # Split the PRE body into construction segments, one per
    # "Construction of a [[n,k,d]] quantum code:" header.
    marks = list(CONSTRUCTION_HDR_RE.finditer(body))
    segments: List[Tuple[Optional[int], str]] = []
    if not marks:
        segments.append((None, body))
    else:
        if marks[0].start() > 0:
            segments.append((None, body[:marks[0].start()]))
        for idx, mk in enumerate(marks):
            end = marks[idx + 1].start() if idx + 1 < len(marks) else len(body)
            segments.append((int(mk.group(3)), body[mk.start():end]))

    for seg_d, seg_text in segments:
        for label, rows in _extract_matrix_blocks(seg_text):
            d = seg_d if seg_d is not None else out.lower_bound
            if d is None:
                continue
            out.candidates.append(MatrixCandidate(d=d, label=label, rows=rows))

    return out


def prime_power(q: int) -> Tuple[int, int]:
    """Write q = p^m and return (p, m)."""
    for p in (2, 3, 5, 7, 11, 13):
        if q % p == 0:
            m, x = 0, q
            while x % p == 0:
                x //= p
                m += 1
            if x != 1:
                raise ValueError("q=%d is not a prime power" % q)
            return p, m
    raise ValueError("q=%d is not a supported prime power" % q)


def expected_row_count(q: int, n: int, k: int) -> int:
    """Number of rows codetables.de prints for a [[n,k]]_q stabilizer matrix.

    The stabilizer of an [[n,k]]_q code is an additive subgroup of GF(q^2)^n of
    size q^(n-k).  codetables.de prints a basis over the PRIME field GF(p), so
    for q = p^m the matrix has m*(n-k) rows whose entries lie in GF(q).
    For q = 2, 3, 5 (m = 1) that is the familiar n-k rows; for q = 4 = 2^2 it
    is 2(n-k) rows written over GF(4) = {0, 1, a, a^2}.
    """
    _, m = prime_power(q)
    return m * (n - k)


def validate_candidate(cand: MatrixCandidate, n: int, k: int,
                       expected_rows: int) -> Optional[str]:
    """Return None if the matrix is usable, else a human readable problem."""
    if expected_rows <= 0:
        return "k = n, stabilizer matrix has no rows"
    if len(cand.rows) != expected_rows:
        # The printed matrix belongs to a different code than the one this page
        # is supposed to be about (codetables.de occasionally shows the matrix
        # of the last code in the construction chain, e.g. [[27,15]] displays a
        # [[27,14,4]] matrix).  Never substitute it - say so and skip.
        return ("row count %d != expected %d, i.e. the printed matrix is for "
                "[[%d,%s]] rather than the requested [[%d,%d]] "
                "(different code)"
                % (len(cand.rows), expected_rows, n,
                   ("%g" % (n - len(cand.rows) * (n - k) / float(expected_rows))),
                   n, k))
    for idx, (left, right) in enumerate(cand.rows):
        if len(left) != n or len(right) != n:
            return ("row %d has %d|%d symbols, expected %d|%d"
                    % (idx, len(left), len(right), n, n))
    return None


def row_weights(cand: MatrixCandidate, n: int, mode: str) -> List[int]:
    """Weight of each stabilizer row."""
    weights = []
    for left, right in cand.rows:
        if mode == "raw":
            w = sum(1 for tok in left if tok != "0") + \
                sum(1 for tok in right if tok != "0")
        else:  # symplectic: non-zero coordinates of the length-n GF(q^2) vector
            w = sum(1 for j in range(n) if left[j] != "0" or right[j] != "0")
        weights.append(w)
    return weights


def compute_metrics(q: int, n: int, k: int, d: int, weights: List[int],
                    x_sum_mode: str, avg_denominator: str = "rows") -> dict:
    """All requested derived quantities for one code."""
    q2 = q * q
    spectrum = dict(Counter(weights))
    n_rows = len(weights)

    # Avg Weight = sum_i W_i f_i / (n - k).  The stabilizer matrix has n-k rows
    # whenever q is prime, so the two denominators agree for q = 2, 3, 5.  For
    # q = 4 the matrix is an m(n-k)-row GF(p)-basis (see expected_row_count),
    # and dividing by n-k would double the mean and push delta_avg past 1,
    # leaving H_{q^2} undefined - so the row count is the right denominator.
    denom = float(n_rows) if avg_denominator == "rows" else float(n - k)
    avg_weight = sum(w * f for w, f in spectrum.items()) / denom
    max_weight = max(spectrum)

    if x_sum_mode == "rows":
        x_value = sum(math.comb(n, w) * (q2 - 1) ** w * f for w, f in spectrum.items())
    else:  # distinct weights
        x_value = sum(math.comb(n, w) * (q2 - 1) ** w for w in spectrum)

    rate_r = ln_bigint(x_value) / (math.log(q) * 2.0 * n) if x_value > 0 else 0.0

    delta_avg = avg_weight / n
    delta_max = max_weight / float(n)
    h_avg = entropy_q2(delta_avg, q2)
    h_max = entropy_q2(delta_max, q2)

    return {
        "spectrum": spectrum,
        "n_rows": n_rows,
        "avg_weight": avg_weight,
        "max_weight": max_weight,
        "x_value": x_value,
        "rate_R": rate_r,
        "delta_avg": delta_avg,
        "delta_max": delta_max,
        "H_delta_avg": h_avg,
        "H_delta_max": h_max,
        "C": h_max - rate_r,
        "C_prime": h_avg - rate_r,
    }


def parse_phase(cfg: Config, store: Store) -> Dict[int, Dict[str, int]]:
    """Parse every cached page for every q and store bounds / codes / skips."""
    per_q: Dict[int, Dict[str, int]] = {}

    for q in cfg.q_values:
        if store.has_results(q) and not cfg.force_reparse:
            log.info("q=%d: results already present, skipping parse "
                     "(use --force-reparse to redo)", q)
            per_q[q] = _summarise_existing(store, q)
            continue

        store.clear_results(q)
        total_pages = store.count_pages(q)
        stats = Counter()
        bounds_rows: List[tuple] = []
        code_rows: List[tuple] = []
        skip_rows: List[tuple] = []

        bar = tqdm(total=total_pages, desc="parse q=%d" % q, unit="pg",
                   ascii=True, dynamic_ncols=False, mininterval=5.0)

        def flush():
            if bounds_rows:
                store.store_bounds(bounds_rows)
                bounds_rows.clear()
            if code_rows:
                store.store_codes(code_rows)
                code_rows.clear()
            if skip_rows:
                store.store_skips(skip_rows)
                skip_rows.clear()

        for n, k, status, text in store.iter_pages(q):
            stats["examined"] += 1
            if bar is not None:
                bar.update(1)

            if status != "ok":
                reason = ("page not downloaded (%s)" % status)
                skip_rows.append((q, n, k, "page unavailable", reason))
                stats["skipped"] += 1
                stats["skip_page_unavailable"] += 1
                continue

            page = parse_page(text, n, k)
            bounds_rows.append((q, n, k, page.lower_bound, page.upper_bound, 0))

            if page.lower_bound is None and page.upper_bound is None:
                skip_rows.append((q, n, k, "unparseable page",
                                  page.note or "no bounds found"))
                stats["skipped"] += 1
                stats["skip_unparseable"] += 1
                continue

            exp_rows = expected_row_count(q, n, k)
            if exp_rows <= 0:
                skip_rows.append((q, n, k, "no stabilizer rows",
                                  "k = n, the stabilizer matrix is empty"))
                stats["skipped"] += 1
                stats["skip_no_rows"] += 1
                continue

            if not page.candidates:
                reason = ("construction listed as 'Missing Construction'"
                          if page.missing_construction else
                          ("no parity-check matrix printed"
                           if page.has_pre else (page.note or "no construction block")))
                skip_rows.append((q, n, k, "no parity-check matrix", reason))
                stats["skipped"] += 1
                stats["skip_no_matrix"] += 1
                continue

            # Highest distance first; fall back to the next candidate if a
            # matrix turns out to be malformed.
            ordered = sorted(page.candidates, key=lambda c: (-c.d, c.label != "stabilizer matrix"))
            chosen: Optional[MatrixCandidate] = None
            problems: List[str] = []
            for cand in ordered:
                problem = validate_candidate(cand, n, k, exp_rows)
                if problem is None:
                    chosen = cand
                    break
                problems.append("d=%d: %s" % (cand.d, problem))

            if chosen is None:
                if problems and all("different code" in p for p in problems):
                    reason = "parity-check matrix belongs to a different code"
                    stats["skip_wrong_code"] += 1
                else:
                    reason = "malformed parity-check matrix"
                    stats["skip_malformed"] += 1
                skip_rows.append((q, n, k, reason, "; ".join(problems)[:500]))
                stats["skipped"] += 1
                continue

            if problems:
                stats["fallback_used"] += 1
                log.debug("q=%d n=%d k=%d used fallback construction: %s",
                          q, n, k, "; ".join(problems))

            try:
                weights = row_weights(chosen, n, cfg.weight_mode)
                metrics = compute_metrics(q, n, k, chosen.d, weights,
                                          cfg.x_sum_mode, cfg.avg_denominator)
            except Exception as exc:                    # noqa: BLE001
                skip_rows.append((q, n, k, "computation failed",
                                  "%s: %s" % (type(exc).__name__, exc)))
                stats["skipped"] += 1
                stats["skip_computation"] += 1
                continue

            bounds_rows[-1] = (q, n, k, page.lower_bound, page.upper_bound, 1)
            label = "[[%d,%d,%d]]_%d" % (n, k, chosen.d, q)
            code_rows.append((
                q, n, k, chosen.d, page.lower_bound, page.upper_bound, label,
                metrics["n_rows"], format_spectrum(metrics["spectrum"]),
                metrics["avg_weight"], metrics["max_weight"], str(metrics["x_value"]),
                metrics["rate_R"], metrics["delta_avg"], metrics["delta_max"],
                metrics["H_delta_avg"], metrics["H_delta_max"],
                metrics["C"], metrics["C_prime"],
                len(page.candidates), chosen.label,
            ))
            stats["processed"] += 1
            if len(page.candidates) > 1:
                stats["multi_construction"] += 1

            if len(code_rows) + len(skip_rows) >= 500:
                flush()

        flush()
        if bar is not None:
            bar.close()
        per_q[q] = dict(stats)
        log.info("q=%d parsed: examined=%d processed=%d skipped=%d",
                 q, stats["examined"], stats["processed"], stats["skipped"])

    return per_q


def _summarise_existing(store: Store, q: int) -> Dict[str, int]:
    processed = store.conn.execute(
        "SELECT COUNT(*) FROM codes WHERE q=?", (q,)).fetchone()[0]
    skipped = store.conn.execute(
        "SELECT COUNT(*) FROM skips WHERE q=?", (q,)).fetchone()[0]
    return {"examined": processed + skipped, "processed": processed, "skipped": skipped}


# ==========================================================================
# Phase 3 - classification
# ==========================================================================

def classify_codes(df: pd.DataFrame, ub_grid: Dict[Tuple[int, int], int],
                   n_max: int) -> pd.Series:
    """Semi-Perfect iff no code CAN exist (per the table's upper bounds) with
    strictly better parameters for the same q:

        (a) d cannot be improved : d == ub(n, k)
        (b) k cannot be improved : ub(n, k') < d  for every k' > k
        (c) n cannot be reduced  : ub(n', k) < d  for every n' < n (n' >= k)
    """
    # Pre-compute, for each n, the best upper bound achievable with k' > k.
    best_above_k: Dict[int, List[int]] = {}
    for (n, k), ub in ub_grid.items():
        arr = best_above_k.get(n)
        if arr is None:
            arr = [-1] * (n + 2)
            best_above_k[n] = arr
        if k < len(arr):
            arr[k] = max(arr[k], ub)
    # suffix maxima over k
    for n, arr in best_above_k.items():
        for k in range(len(arr) - 2, -1, -1):
            arr[k] = max(arr[k], arr[k + 1])

    # Pre-compute, for each k, the best upper bound achievable with n' < n.
    best_below_n: Dict[int, List[int]] = {}
    for k in {kk for (_, kk) in ub_grid}:
        arr = [-1] * (n_max + 2)
        running = -1
        for n in range(0, n_max + 1):
            arr[n] = running                       # best over n' < n
            ub = ub_grid.get((n, k))
            if ub is not None:
                running = max(running, ub)
        best_below_n[k] = arr

    labels = []
    for n, k, d in zip(df["n"], df["k"], df["d"]):
        ub_here = ub_grid.get((n, k))
        if ub_here is None or d < ub_here:
            labels.append(NORMAL)
            continue
        arr_k = best_above_k.get(n)
        if arr_k is not None and k + 1 < len(arr_k) and arr_k[k + 1] >= d:
            labels.append(NORMAL)
            continue
        arr_n = best_below_n.get(k)
        if arr_n is not None and n < len(arr_n) and arr_n[n] >= d:
            labels.append(NORMAL)
            continue
        labels.append(SEMI_PERFECT)
    return pd.Series(labels, index=df.index)


# ==========================================================================
# Phase 4 - export
# ==========================================================================

def build_export_frame(df: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame({
        "q": df["q"].astype(int),
        "n": df["n"].astype(int),
        "k": df["k"].astype(int),
        "Code Parameters": df["label"],
        "Type of Code": df["code_type"],
        "k/n": df["k"] / df["n"],
        "Weight Spectrum": df["spectrum"],
        "Average Weight": df["avg_weight"],
        "Max Weight": df["max_weight"].astype(int),
        "R": df["rate_R"],
        "C": df["C"],
        "C'": df["C_prime"],
    })
    return out[EXCEL_COLUMNS]


def write_excel(df_export: pd.DataFrame, path: str) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df_export.to_excel(writer, sheet_name="Codes", index=False)
        ws = writer.sheets["Codes"]

        # AutoFilter over the whole table so 'Type of Code' is filterable.
        ws.auto_filter.ref = ws.dimensions
        ws.freeze_panes = "A2"

        widths = {"q": 6, "n": 7, "k": 7, "Code Parameters": 20,
                  "Type of Code": 26, "k/n": 11, "Weight Spectrum": 44,
                  "Average Weight": 16, "Max Weight": 12,
                  "R": 13, "C": 13, "C'": 13}
        num_fmt = {"k/n": "0.000000", "Average Weight": "0.000000",
                   "R": "0.000000", "C": "0.000000", "C'": "0.000000"}

        from openpyxl.styles import Alignment, Font, PatternFill
        header_font = Font(bold=True, color="FFFFFF")
        header_fill = PatternFill("solid", fgColor="2F4F6F")
        for col_idx, name in enumerate(df_export.columns, start=1):
            letter = ws.cell(row=1, column=col_idx).column_letter
            ws.column_dimensions[letter].width = widths.get(name, 14)
            cell = ws.cell(row=1, column=col_idx)
            cell.font = header_font
            cell.fill = header_fill
            cell.alignment = Alignment(horizontal="center", vertical="center")
            fmt = num_fmt.get(name)
            if fmt:
                for row_idx in range(2, ws.max_row + 1):
                    ws.cell(row=row_idx, column=col_idx).number_format = fmt
    log.info("wrote %s (%d rows)", path, len(df_export))


# ==========================================================================
# Phase 4 - plots
# ==========================================================================

#: Semi-perfect-only panels hold far fewer points than the all-codes panels,
#: so their markers get an extra bump on top of the density-based sizing.
SEMI_SIZE_SCALE = 1.6

#: Thin black outline on every marker, and the dotted reference-line style.
MARKER_EDGE_COLOR = "black"
MARKER_EDGE_WIDTH = 0.3
REFLINE_COLOR = "#444444"
REFLINE_WIDTH = 1.7

VIBRANT = "#E8412C"
MUTED = "#CBD5DF"
GRADIENT_CMAP = "viridis"


def _marker_size(n_points: int) -> float:
    if n_points <= 0:
        return 20.0
    return float(max(4.0, min(28.0, 2200.0 / math.sqrt(n_points))))


def _finish(fig, ax, title: str, xlabel: str, ylabel: str, path: str, dpi: int) -> None:
    ax.set_title(title, fontsize=13, pad=12)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.grid(True, alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _add_identity_line(ax) -> None:
    """Dotted y = x reference, drawn without letting it expand the view."""
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    lo, hi = min(xlim[0], ylim[0]), max(xlim[1], ylim[1])
    ax.plot([lo, hi], [lo, hi], linestyle=":", color=REFLINE_COLOR,
            linewidth=REFLINE_WIDTH, zorder=2, label="$H_{q^2}(\\delta) = \\mathcal{R}$")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def _add_zero_line(ax, label: str) -> None:
    """Dotted horizontal reference at C = 0 (resp. C' = 0)."""
    ax.axhline(0.0, linestyle=":", color=REFLINE_COLOR, linewidth=REFLINE_WIDTH,
               zorder=2, label=label)


def _legend_if_labelled(ax, **kwargs) -> None:
    """Draw a legend only when something on the axes carries a label."""
    handles, _ = ax.get_legend_handles_labels()
    if handles:
        leg = ax.legend(loc="best", frameon=True, fontsize=10, **kwargs)
        leg.get_frame().set_alpha(0.95)


def _scatter_gradient(x, y, c, title, xlabel, ylabel, path, dpi, cbar_label="k/n",
                      size_scale=1.0, identity_line=False, zero_line=None):
    fig, ax = plt.subplots(figsize=(9, 5.6))
    sc = ax.scatter(x, y, c=c, cmap=GRADIENT_CMAP,
                    s=_marker_size(len(x)) * size_scale,
                    edgecolors=MARKER_EDGE_COLOR, linewidths=MARKER_EDGE_WIDTH,
                    alpha=0.85, rasterized=len(x) > 8000, zorder=3)
    if identity_line:
        _add_identity_line(ax)
    if zero_line:
        _add_zero_line(ax, zero_line)
    _legend_if_labelled(ax)
    cbar = fig.colorbar(sc, ax=ax, pad=0.015)
    cbar.set_label(cbar_label, fontsize=11)
    _finish(fig, ax, title, xlabel, ylabel, path, dpi)


def _scatter_plain(x, y, title, xlabel, ylabel, path, dpi, color="#2F6FA8",
                   integer_x=False, size_scale=1.0, zero_line=None):
    fig, ax = plt.subplots(figsize=(9, 5.6))
    ax.scatter(x, y, s=_marker_size(len(x)) * size_scale, color=color,
               edgecolors=MARKER_EDGE_COLOR, linewidths=MARKER_EDGE_WIDTH,
               alpha=0.7, rasterized=len(x) > 8000, zorder=3)
    if zero_line:
        _add_zero_line(ax, zero_line)
    _legend_if_labelled(ax)
    if integer_x:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    _finish(fig, ax, title, xlabel, ylabel, path, dpi)


def make_plots(q: int, df: pd.DataFrame, out_dir: str, dpi: int) -> List[str]:
    """All nine requested discrete scatter plots for a single q."""
    os.makedirs(out_dir, exist_ok=True)
    paths: List[str] = []
    if df.empty:
        log.warning("q=%d: no codes, no plots produced", q)
        return paths

    if sns is not None:
        sns.set_theme(style="whitegrid", context="notebook")
    plt.rcParams["axes.facecolor"] = "#FFFFFF"
    plt.rcParams["figure.facecolor"] = "#FFFFFF"

    kn = df["k"] / df["n"]
    semi = df[df["code_type"] == SEMI_PERFECT]
    kn_semi = semi["k"] / semi["n"] if not semi.empty else pd.Series(dtype=float)
    suffix = "q%d" % q
    tag = "  (q = %d)" % q

    def p(name: str) -> str:
        path = os.path.join(out_dir, "%s_%s.png" % (name, suffix))
        paths.append(path)
        return path

    # 1. C vs n, all codes, gradient by k/n
    _scatter_gradient(df["n"], df["C"], kn,
                      "C vs n - all codes" + tag, "n", "C = $H_{q^2}(\\delta_{max}) - \\mathcal{R}$",
                      p("01_C_vs_n_all_gradient"), dpi, zero_line="$C = 0$")

    # 2. C vs k, all codes
    _scatter_plain(df["k"], df["C"],
                   "C vs k - all codes" + tag, "k", "C = $H_{q^2}(\\delta_{max}) - \\mathcal{R}$",
                   p("02_C_vs_k_all"), dpi, integer_x=True,
                   zero_line="$C = 0$")

    # 3. C vs k/n, all codes
    _scatter_plain(kn, df["C"],
                   "C vs k/n - all codes" + tag, "k/n",
                   "C = $H_{q^2}(\\delta_{max}) - \\mathcal{R}$",
                   p("03_C_vs_k_over_n_all"), dpi, zero_line="$C = 0$")

    # 4-6. Semi-perfect only
    if semi.empty:
        log.warning("q=%d: no semi-perfect codes, plots 4-6 skipped", q)
    else:
        _scatter_gradient(semi["n"], semi["C_prime"], kn_semi,
                          "C' vs n - semi-perfect quantum codes" + tag, "n",
                          "C' = $H_{q^2}(\\delta_{avg}) - \\mathcal{R}$",
                          p("04_Cprime_vs_n_semiperfect_gradient"), dpi,
                          size_scale=SEMI_SIZE_SCALE, zero_line="$C' = 0$")
        _scatter_plain(semi["k"], semi["C_prime"],
                       "C' vs k - semi-perfect quantum codes" + tag, "k",
                       "C' = $H_{q^2}(\\delta_{avg}) - \\mathcal{R}$",
                       p("05_Cprime_vs_k_semiperfect"), dpi, color=VIBRANT,
                       integer_x=True, size_scale=SEMI_SIZE_SCALE,
                       zero_line="$C' = 0$")
        _scatter_plain(kn_semi, semi["C_prime"],
                       "C' vs k/n - semi-perfect quantum codes" + tag, "k/n",
                       "C' = $H_{q^2}(\\delta_{avg}) - \\mathcal{R}$",
                       p("06_Cprime_vs_k_over_n_semiperfect"), dpi, color=VIBRANT,
                       size_scale=SEMI_SIZE_SCALE, zero_line="$C' = 0$")

    # 7. H(delta_max) vs R
    _scatter_gradient(df["rate_R"], df["H_delta_max"], kn,
                      "$H_{q^2}(\\delta_{max})$ vs $\\mathcal{R}$ - all codes" + tag,
                      "$\\mathcal{R} = \\log_q(X)/(2n)$", "$H_{q^2}(\\delta_{max})$",
                      p("07_H_delta_max_vs_R"), dpi, identity_line=True)

    # 8. H(delta_avg) vs R
    _scatter_gradient(df["rate_R"], df["H_delta_avg"], kn,
                      "$H_{q^2}(\\delta_{avg})$ vs $\\mathcal{R}$ - all codes" + tag,
                      "$\\mathcal{R} = \\log_q(X)/(2n)$", "$H_{q^2}(\\delta_{avg})$",
                      p("08_H_delta_avg_vs_R"), dpi, identity_line=True)

    # 9. C' vs n, all codes, categorical colouring + legend
    normal = df[df["code_type"] != SEMI_PERFECT]
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    size = _marker_size(len(df))
    if not normal.empty:
        ax.scatter(normal["n"], normal["C_prime"], s=size, color=MUTED,
                   edgecolors=MARKER_EDGE_COLOR, linewidths=MARKER_EDGE_WIDTH,
                   alpha=0.75, label="Normal Code (%d)" % len(normal),
                   rasterized=len(normal) > 8000, zorder=2)
    if not semi.empty:
        ax.scatter(semi["n"], semi["C_prime"], s=size * 1.6, color=VIBRANT,
                   edgecolors=MARKER_EDGE_COLOR, linewidths=MARKER_EDGE_WIDTH,
                   alpha=0.95,
                   label="Semi-Perfect Quantum Code (%d)" % len(semi),
                   rasterized=len(semi) > 8000, zorder=3)
    _add_zero_line(ax, "$C' = 0$")
    _legend_if_labelled(ax, markerscale=2.2)
    _finish(fig, ax, "C' vs n - semi-perfect vs normal codes" + tag, "n",
            "C' = $H_{q^2}(\\delta_{avg}) - \\mathcal{R}$",
            p("09_Cprime_vs_n_semiperfect_vs_normal"), dpi)

    log.info("q=%d: wrote %d plots to %s", q, len(paths), out_dir)
    return paths


# ==========================================================================
# Orchestration
# ==========================================================================

def setup_logging(cfg: Config, verbose: bool) -> None:
    os.makedirs(cfg.log_dir, exist_ok=True)
    log.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s",
                            "%H:%M:%S")

    fh = logging.FileHandler(os.path.join(cfg.log_dir, "run.log"), encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    log.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.DEBUG if verbose else logging.INFO)
    sh.setFormatter(fmt)
    log.addHandler(sh)

    logging.getLogger("urllib3").setLevel(logging.WARNING)


def run(cfg: Config, verbose: bool) -> None:
    os.makedirs(cfg.out_dir, exist_ok=True)
    os.makedirs(cfg.plot_dir, exist_ok=True)
    setup_logging(cfg, verbose)

    started = time.time()
    log.info("=" * 74)
    log.info("QECC pipeline  |  q = %s  |  weight-mode=%s  x-sum-mode=%s  avg-denom=%s",
             ", ".join(map(str, cfg.q_values)), cfg.weight_mode, cfg.x_sum_mode,
             cfg.avg_denominator)
    log.info("output directory: %s", cfg.out_dir)
    log.info("=" * 74)

    store = Store(cfg.db_path)
    try:
        # ---- phase 1 ----
        if cfg.skip_fetch:
            fetch_stats = {"requested": 0, "already_cached": 0,
                           "downloaded": 0, "errors": 0}
            log.info("fetch phase skipped (--skip-fetch)")
        else:
            fetch_stats = fetch_phase(cfg, store)
            log.info("fetch done: %d in range, %d already cached, %d downloaded, "
                     "%d errors", fetch_stats["requested"],
                     fetch_stats["already_cached"], fetch_stats["downloaded"],
                     fetch_stats["errors"])

        # ---- phase 2 ----
        parse_stats = parse_phase(cfg, store)

        # ---- phases 3 + 4 ----
        summary_rows = []
        combined = []
        for q in cfg.q_values:
            df = store.codes_dataframe(q)
            ub_grid = store.upper_bound_grid(q)
            n_max = cfg.n_max[q]

            missing_bounds = sum(
                1 for n in range(cfg.n_min, n_max + 1) for k in range(0, n + 1)
                if (n, k) not in ub_grid)
            if missing_bounds:
                log.warning("q=%d: upper bound unknown for %d parameter sets; "
                            "classification of affected codes may be optimistic",
                            q, missing_bounds)

            if df.empty:
                log.warning("q=%d: no processed codes", q)
                df["code_type"] = []
            else:
                df["code_type"] = classify_codes(df, ub_grid, n_max)

            n_semi = int((df["code_type"] == SEMI_PERFECT).sum()) if not df.empty else 0
            st = parse_stats.get(q, {})
            log.info("q=%d: %d codes processed, %d semi-perfect, %d normal",
                     q, len(df), n_semi, len(df) - n_semi)

            export = build_export_frame(df) if not df.empty else \
                pd.DataFrame(columns=EXCEL_COLUMNS)
            xlsx = os.path.join(cfg.out_dir, "QECC_codes_q%d.xlsx" % q)
            write_excel(export, xlsx)

            make_plots(q, df, os.path.join(cfg.plot_dir, "q%d" % q), cfg.dpi)

            if not df.empty:
                combined.append(export)

            summary_rows.append({
                "q": q,
                "n_range": "%d..%d" % (cfg.n_min, n_max),
                "parameter_sets_examined": st.get("examined", len(df)),
                "codes_processed": len(df),
                "codes_skipped": st.get("skipped", 0),
                "semi_perfect": n_semi,
                "normal": len(df) - n_semi,
                "skip_no_matrix": st.get("skip_no_matrix", 0),
                "skip_no_rows": st.get("skip_no_rows", 0),
                "skip_malformed": st.get("skip_malformed", 0),
                "skip_wrong_code": st.get("skip_wrong_code", 0),
                "skip_unparseable": st.get("skip_unparseable", 0),
                "skip_page_unavailable": st.get("skip_page_unavailable", 0),
                "pages_with_multiple_constructions": st.get("multi_construction", 0),
                "excel": os.path.basename(xlsx),
            })

        # ---- skip log + summary ----
        skips = store.skips_dataframe()
        skip_csv = os.path.join(cfg.log_dir, "skipped_codes.csv")
        skips.to_csv(skip_csv, index=False, encoding="utf-8")
        log.info("skip log written: %s (%d entries)", skip_csv, len(skips))

        if combined:
            all_df = pd.concat(combined, ignore_index=True)
            all_csv = os.path.join(cfg.out_dir, "all_codes_combined.csv")
            all_df.to_csv(all_csv, index=False, encoding="utf-8")
            log.info("combined dataset written: %s (%d rows)", all_csv, len(all_df))

        summary = pd.DataFrame(summary_rows)
        summary_path = os.path.join(cfg.out_dir, "run_summary.csv")
        summary.to_csv(summary_path, index=False, encoding="utf-8")

        elapsed = time.time() - started
        report = _render_summary(cfg, summary, skips, fetch_stats, elapsed)
        with open(os.path.join(cfg.out_dir, "run_summary.txt"), "w",
                  encoding="utf-8") as fh:
            fh.write(report)
        print(report)
        log.info("finished in %.1f s", elapsed)
    finally:
        store.close()


def _render_summary(cfg: Config, summary: pd.DataFrame, skips: pd.DataFrame,
                    fetch_stats: dict, elapsed: float) -> str:
    lines = []
    add = lines.append
    add("=" * 78)
    add("QECC PIPELINE - RUN SUMMARY")
    add("=" * 78)
    add("output directory : %s" % cfg.out_dir)
    add("q values         : %s" % ", ".join(map(str, cfg.q_values)))
    add("weight mode      : %s   X sum mode: %s" % (cfg.weight_mode, cfg.x_sum_mode))
    add("elapsed          : %.1f s" % elapsed)
    add("")
    add("-- fetching ------------------------------------------------------------")
    add("pages in range        : %d" % fetch_stats.get("requested", 0))
    add("already cached        : %d" % fetch_stats.get("already_cached", 0))
    add("downloaded this run   : %d" % fetch_stats.get("downloaded", 0))
    add("download failures     : %d" % fetch_stats.get("errors", 0))
    add("")
    add("-- per q ---------------------------------------------------------------")
    if not summary.empty:
        header = ("  q   n range   examined  processed  skipped  semi-perfect  normal")
        add(header)
        add("  " + "-" * (len(header) - 2))
        for _, r in summary.iterrows():
            add("  %-3d %-9s %9d %10d %8d %13d %7d" % (
                r["q"], r["n_range"], r["parameter_sets_examined"],
                r["codes_processed"], r["codes_skipped"], r["semi_perfect"],
                r["normal"]))
        add("  " + "-" * (len(header) - 2))
        add("  TOTAL           %9d %10d %8d %13d %7d" % (
            summary["parameter_sets_examined"].sum(),
            summary["codes_processed"].sum(),
            summary["codes_skipped"].sum(),
            summary["semi_perfect"].sum(),
            summary["normal"].sum()))
    add("")
    add("-- skip reasons --------------------------------------------------------")
    if skips.empty:
        add("  (none)")
    else:
        for reason, count in skips["reason"].value_counts().items():
            add("  %-34s %8d" % (reason, count))
    add("")
    add("-- outputs -------------------------------------------------------------")
    for _, r in summary.iterrows():
        add("  %s" % r["excel"])
    add("  plots/q<q>/*.png            (9 scatter plots per q)")
    add("  logs/skipped_codes.csv      (every skipped parameter set + reason)")
    add("  logs/run.log                (full run log)")
    add("  run_summary.csv / .txt")
    add("=" * 78)
    return "\n".join(lines)


# ==========================================================================
# CLI
# ==========================================================================

def parse_args(argv: Optional[List[str]] = None) -> Tuple[Config, bool]:
    ap = argparse.ArgumentParser(
        description="Scrape, analyse, export and plot quantum codes from "
                    "codetables.de",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--out-dir", default=os.path.dirname(os.path.abspath(__file__)),
                    help="directory for cache, Excel files, plots and logs")
    ap.add_argument("--q", type=int, nargs="+", default=[2, 3, 4, 5, 7, 8],
                    choices=[2, 3, 4, 5, 7, 8],
                    help="prime powers q to process")
    ap.add_argument("--n-min", type=int, default=1)
    ap.add_argument("--n-max", type=int, default=None,
                    help="cap n (default: the table maximum, 256 for q=2 else 100)")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--rps", type=float, default=6.0,
                    help="target HTTP requests per second (politeness limit)")
    ap.add_argument("--timeout", type=float, default=30.0)
    ap.add_argument("--max-attempts", type=int, default=5)
    ap.add_argument("--weight-mode", choices=["symplectic", "raw"],
                    default="symplectic")
    ap.add_argument("--x-sum-mode", choices=["distinct", "rows"], default="distinct")
    ap.add_argument("--avg-denominator", choices=["rows", "n-k"], default="rows",
                    help="denominator of Avg Weight; identical for q=2,3,5")
    ap.add_argument("--skip-fetch", action="store_true",
                    help="work entirely from the local cache")
    ap.add_argument("--force-reparse", action="store_true",
                    help="recompute all metrics even if results already exist")
    ap.add_argument("--dpi", type=int, default=160)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    n_max = dict(DEFAULT_N_MAX)
    if args.n_max is not None:
        for q in n_max:
            n_max[q] = min(n_max[q], args.n_max)

    cfg = Config(
        out_dir=args.out_dir,
        q_values=tuple(sorted(set(args.q))),
        n_max=n_max,
        n_min=args.n_min,
        workers=args.workers,
        rps=args.rps,
        timeout=args.timeout,
        max_attempts=args.max_attempts,
        weight_mode=args.weight_mode,
        x_sum_mode=args.x_sum_mode,
        avg_denominator=args.avg_denominator,
        skip_fetch=args.skip_fetch,
        force_reparse=args.force_reparse,
        dpi=args.dpi,
    )
    return cfg, args.verbose


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")  # type: ignore
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")  # type: ignore
    except Exception:                                   # noqa: BLE001
        pass
    cfg, verbose = parse_args()
    run(cfg, verbose)


if __name__ == "__main__":
    main()
