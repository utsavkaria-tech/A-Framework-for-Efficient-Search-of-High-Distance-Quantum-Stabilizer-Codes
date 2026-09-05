// ===========================================================================================
//  reverse_8w10_3w8.cpp
//
//  REVERSE two-stage optimisation of a 14-qubit stabilizer code.
//
//      STAGE 1 :  8 independent, mutually commuting generators, each of Pauli weight 10.
//      STAGE 2 :  3 further generators of Pauli weight 8, commuting with everything,
//                 chosen to detect the errors Stage 1 leaves undetected.
//
//  The earlier experiments in this project did it the other way round (3 x weight-8 first,
//  then 8 x weight-10).  Nothing from those experiments is read or written here.
//
// -------------------------------------------------------------------------------------------
//  1.  REPRESENTATION
//
//  A Pauli modulo phase on 14 qubits is a pair (x|z) in F_2^14 x F_2^14, packed into one
//  uint32 as  x | (z << 16).  The Pauli (qubit) weight is  wt = popcount(x | z), so
//  wt(I)=0 and wt(X)=wt(Y)=wt(Z)=1.  This is NOT the binary Hamming weight of (x|z).
//
//  Two Paulis commute iff the symplectic form  <a,b> = x_a.z_b + z_a.x_b  vanishes mod 2.
//  The syndrome bit of error e under generator g is exactly <g,e>, so
//
//        e is DETECTED  <=>  e does not commute with the whole stabilizer group.
//
// -------------------------------------------------------------------------------------------
//  2.  THE EXACT OBJECTIVE, IN CLOSED FORM  (this is the engine of the whole program)
//
//  Let L <= F_2^28 be the F_2-span of the chosen generators, dim L = r.  The number of
//  weight-w Paulis that commute with every element of L is, by Fourier inversion over L,
//
//        #{e : wt(e)=w, s(e)=0}  =  (1/2^r) SUM_{h in L} SUM_{wt(e)=w} (-1)^{<h,e>} .
//
//  The inner sum factorises over the 14 qubits.  On a qubit where h acts trivially all three
//  non-identity Paulis commute with it and contribute +1 each, giving a factor 3; on a qubit
//  where h acts non-trivially exactly one of the three commutes, giving 1-2 = -1.  Hence with
//  m = wt(h) the generating function is  (1 - y)^m (1 + 3y)^{14-m}  and
//
//        EW[m][w]  =  SUM_k C(m,k) (-1)^k C(14-m, w-k) 3^{w-k} .
//
//  Summing the objective weight window w = 1..4 gives a single table
//
//        P[m] = SUM_{w=1..4} EW[m][w] ,
//
//  and the number of UNDETECTED errors of weight 1..4 is exactly
//
//        U(L)  =  (1/2^r) SUM_{h in L} P(wt(h)) ,        M = TOTAL - U.
//
//  So the objective of an r-generator code is a sum of 2^r table lookups instead of
//  91770 * r symplectic products.  It depends only on the SUBSPACE L -- not on the basis --
//  which is what makes the whole k-opt reformulation of section 4 legitimate.
//  P[] is computed at run time and checked against brute force; nothing is hard-coded.
//
// -------------------------------------------------------------------------------------------
//  3.  COLUMN TYPES  (the mandated qubit-structured formulation)
//
//  Write the Stage-1 check matrix as H = [X|Z] with X,Z in F_2^{8x14}.  Qubit j contributes
//  the pair of columns (x_j, z_j) in F_2^8 x F_2^8 and we set
//
//        W_j  =  span{ x_j , z_j }  <=  F_2^8 ,      dim W_j <= 2 .
//
//  For h = sum_i a_i g_i the Pauli of h on qubit j is (a.x_j, a.z_j), which is the identity
//  iff a is orthogonal to W_j.  Therefore
//
//        wt(h_a)  =  14 - #{ j : a _|_ W_j } ,
//
//  i.e. the entire weight function -- hence the entire objective -- is determined by the
//  MULTISET { W_1, ..., W_14 } of subspaces of F_2^8.  Consequences used below:
//
//    * a local Clifford on qubit j replaces (x_j,z_j) by (ax_j+bz_j, cx_j+dz_j) with
//      ad+bc = 1; this fixes W_j, so it preserves every weight, the objective, and (because
//      the qubit-j contribution to <g_i,g_k> is multiplied by det = 1) all commutators;
//    * permuting qubits permutes the multiset, so only the multiset matters;
//    * commutation of the whole family is  SUM_j plucker(W_j) = 0  with
//      plucker(W) = x z^T + z x^T (symmetric, zero diagonal, basis-independent, zero when
//      dim W <= 1);
//    * generator i has weight 10  <=>  exactly 4 columns satisfy W_j <= H_i := {v : v_i = 0};
//    * rank H = 8  <=>  W_1 + ... + W_14 = F_2^8  <=>  wt(h_a) > 0 for every a != 0.
//
//  Raw column types: 2^16 = 65536.  Reduced types (subspaces of F_2^8 of dim <= 2):
//  1 + 255 + 10795 = 11051.  Both counts are recomputed by the program.  The reduction is
//  exact but, unlike the r=3 and r=4 cases of the earlier experiments (15 and 51 types),
//  the resulting multiset space is astronomically large -- see section 5 and the README.
//
// -------------------------------------------------------------------------------------------
//  4.  k-OPT AS A SUBSPACE SEARCH IN A QUOTIENT
//
//  Keep t = 8-k of the Stage-1 generators, K = their span (dim t, isotropic).  Every valid
//  completion is an 8-dimensional isotropic D' with K <= D' <= K^perp, i.e. exactly a
//  k-dimensional totally isotropic subspace of the symplectic space
//
//        Q = K^perp / K ,      dim Q = 28 - 2t = 12 + 2k .
//
//  Writing phi(q) = SUM_{u in K} P(wt(lift(q) + u)), the objective telescopes exactly:
//
//        S(D')  =  SUM_{v in span(q_1..q_k)} phi(v) ,        U = S / 256 .
//
//  k=1 : dim Q = 14 (16384 classes)  -- swept exhaustively, cost ~ms.
//  k=2 : dim Q = 16 (65536 classes)  -- swept exhaustively over all 28 retained pairs.
//
//  The weight-10 requirement is applied to the SUBSPACE, not to individual replacements:
//  D' is admissible iff its weight-10 elements have rank 8, which is the exact condition for
//  "D' has a basis of eight weight-10 generators".  This is a strictly larger (and correct)
//  neighbourhood than demanding that each replaced class carry a weight-10 lift.
//
// -------------------------------------------------------------------------------------------
//  5.  WHAT IS PROVEN AND WHAT IS NOT
//
//  Stage 1 CANNOT be exhausted.  The number of 8-dimensional isotropic subspaces of F_2^28 is
//  prod_{i=0}^{7}(2^{28-i} - 2^i) / |GL(8,2)|; even after dividing by the full 14! * 6^14
//  local-Clifford-and-permutation group an astronomical number remains.  The column multiset
//  formulation is no better: C(11064,14) multisets.  Both counts are computed and printed.
//  Every Stage-1 number this program prints is therefore BEST KNOWN, never OPTIMUM.  A
//  rigorous LP lower bound on the residual |U_8| is computed separately so that the gap
//  between best-known and provable is explicit.
//
//  Stage 2 IS exhausted.  With Stage 1 frozen, Q2 = L8^perp/L8 has dim 12 (4096 classes) and
//  the program enumerates every 3-dimensional isotropic subspace spanned by weight-8-usable
//  classes.  Its output is a PROVEN OPTIMUM conditional on the frozen Stage-1 code.
// ===========================================================================================

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cstdarg>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>
#include <random>
#include <atomic>
#include <mutex>
#include <thread>
#include <chrono>
#include <bit>
#include <filesystem>
#include <ctime>

constexpr int NQ   = 14;          // physical qubits
constexpr int R1   = 8;           // Stage-1 generators
constexpr int R2   = 3;           // Stage-2 generators
constexpr int RTOT = R1 + R2;     // 11
constexpr int MAXW = 4;           // error weight window
constexpr int W1   = 10;          // Stage-1 generator weight
constexpr int W2   = 8;           // Stage-2 generator weight
constexpr uint32_t QM   = (1u << NQ) - 1u;
constexpr uint32_t NONE = 0xFFFFFFFFu;
constexpr int S1SZ = 1 << R1;     // 256 elements of the Stage-1 group

static inline uint32_t pk(uint32_t x, uint32_t z) { return x | (z << 16); }
static inline uint32_t px(uint32_t p) { return p & QM; }
static inline uint32_t pz(uint32_t p) { return (p >> 16) & QM; }
static inline int pwt(uint32_t p) { return std::popcount((p | (p >> 16)) & QM); }
static inline int par(unsigned a) { return std::popcount(a) & 1; }
static inline int symp(uint32_t a, uint32_t b) { return par((px(a) & pz(b)) ^ (pz(a) & px(b))); }

static const char PCH[2][2] = { {'I','Z'}, {'X','Y'} };
static std::string pstr(uint32_t p) {
    std::string s; s.reserve(NQ);
    for (int j = 0; j < NQ; ++j) s += PCH[(px(p) >> j) & 1][(pz(p) >> j) & 1];
    return s;
}
static std::string bstr(uint32_t v) {
    std::string s; for (int j = 0; j < NQ; ++j) s += char('0' + ((v >> j) & 1)); return s;
}
static std::string supp(uint32_t E) {
    std::string s; bool first = true;
    for (int j = 0; j < NQ; ++j) {
        int x = (px(E) >> j) & 1, z = (pz(E) >> j) & 1;
        if (!x && !z) continue;
        if (!first) s += " ";
        s += PCH[x][z]; s += "_" + std::to_string(j); first = false;
    }
    return s.empty() ? std::string("I") : s;
}
static uint32_t parse_pauli(const char* s) {
    uint32_t x = 0, z = 0;
    for (int j = 0; j < NQ; ++j) {
        switch (s[j]) {
            case 'I': break;
            case 'X': x |= 1u << j; break;
            case 'Z': z |= 1u << j; break;
            case 'Y': x |= 1u << j; z |= 1u << j; break;
            default: fprintf(stderr, "bad Pauli char\n"); exit(1);
        }
    }
    return pk(x, z);
}

// ---------------------------------------------------------------- timing / text buffers
// Start time, initialised on first use.  A function-local static is initialised exactly once
// and thread-safely; a namespace-scope time_point assigned at the top of main() was reading
// back as the clock epoch here, which made every absolute "elapsed" figure the machine
// uptime instead of the run time (the deltas were always right, which is what hid it).
static long long now_ns() {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(
               std::chrono::steady_clock::now().time_since_epoch()).count();
}
static std::atomic<long long> T0_NS{0};
static void set_start_time() { T0_NS.store(now_ns(), std::memory_order_relaxed); }
// Self-healing: if the origin was never set, the first caller sets it, so "elapsed" is never
// the clock epoch (which on this platform is boot time and would report the machine uptime).
static double elapsed() {
    long long t0 = T0_NS.load(std::memory_order_relaxed);
    if (t0 == 0) {
        long long n = now_ns(), expect = 0;
        t0 = T0_NS.compare_exchange_strong(expect, n) ? n : expect;
    }
    return double(now_ns() - t0) * 1e-9;
}
static std::string RB;
static void ap(const char* fmt, ...) {
    char buf[16384]; va_list a; va_start(a, fmt); vsnprintf(buf, sizeof(buf), fmt, a); va_end(a); RB += buf;
}
static void dump(const std::string& path, bool echo = false) {
    if (echo) { fputs(RB.c_str(), stdout); fflush(stdout); }
    FILE* f = fopen(path.c_str(), "w"); if (f) { fputs(RB.c_str(), f); fclose(f); }
}
static FILE* LOGF = nullptr;
static std::mutex LOGMTX;
static void logline(const char* fmt, ...) {
    char buf[16384]; va_list a; va_start(a, fmt); vsnprintf(buf, sizeof(buf), fmt, a); va_end(a);
    std::lock_guard<std::mutex> lk(LOGMTX);
    fputs(buf, stdout); fflush(stdout);
    if (LOGF) { fputs(buf, LOGF); fflush(LOGF); }
}
static std::string now_stamp() {
    std::time_t t = std::time(nullptr);
    char b[64]; std::strftime(b, sizeof(b), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
    return std::string(b);
}
static std::string hms(double s) {
    if (s < 0) s = 0;
    long long x = (long long)s;
    char b[64]; snprintf(b, sizeof(b), "%lldh %02lldm %02llds", x / 3600, (x / 60) % 60, x % 60);
    return std::string(b);
}

// ---------------------------------------------------------------- the closed-form tables
// EW[m][w] for the FULL weight range w = 0..14, not just the objective window w = 0..4:
// the MacWilliams/duality certificate of section 7 needs every K_w.  (Sizing this [MAXW+1]
// and filling the tail later overruns the array and silently corrupts neighbouring globals.)
static int64_t EWT[NQ + 1][NQ + 1];     // EW[m][w]
static int64_t PTAB[NQ + 1];            // P[m] = sum_{w=1..4} EW[m][w]
static int64_t BINOM[NQ + 1][NQ + 1];
static void build_tables() {
    for (int n = 0; n <= NQ; ++n) {
        for (int k = 0; k <= NQ; ++k) BINOM[n][k] = 0;
        BINOM[n][0] = 1;
        for (int k = 1; k <= n; ++k) BINOM[n][k] = BINOM[n-1][k-1] + BINOM[n-1][k];
    }
    for (int m = 0; m <= NQ; ++m) {
        for (int w = 0; w <= NQ; ++w) {
            int64_t s = 0;
            for (int k = 0; k <= w; ++k) {
                if (k > m) break;
                int r = w - k; if (r > NQ - m) continue;
                int64_t p3 = 1; for (int i = 0; i < r; ++i) p3 *= 3;
                int64_t term = BINOM[m][k] * BINOM[NQ - m][r] * p3;
                s += (k & 1) ? -term : term;
            }
            EWT[m][w] = s;
        }
        PTAB[m] = 0; for (int w = 1; w <= MAXW; ++w) PTAB[m] += EWT[m][w];
    }
}

// ---------------------------------------------------------------- error enumeration
template <class F> static void for_each_error(F&& f) {
    static const uint32_t PX3[3] = {1,1,0}, PZ3[3] = {0,1,1};
    int idx[MAXW];
    for (int w = 1; w <= MAXW; ++w) {
        for (int i = 0; i < w; ++i) idx[i] = i;
        for (;;) {
            int npat = 1; for (int i = 0; i < w; ++i) npat *= 3;
            for (int p = 0; p < npat; ++p) {
                uint32_t x = 0, z = 0; int q = p;
                for (int i = 0; i < w; ++i) { int a = q % 3; q /= 3;
                    x |= PX3[a] << idx[i]; z |= PZ3[a] << idx[i]; }
                f(pk(x, z), w);
            }
            int i = w - 1; while (i >= 0 && idx[i] == NQ - w + i) --i;
            if (i < 0) break;
            ++idx[i]; for (int k = i + 1; k < w; ++k) idx[k] = idx[k-1] + 1;
        }
    }
}
static std::vector<uint32_t> ALLERR;      // all 91770, built once
static std::vector<uint8_t>  ALLW;
static int64_t TOT_ERR = 0, TOT_W4 = 0, TOT_BY_W[MAXW + 1];
static void build_errors() {
    ALLERR.clear(); ALLW.clear(); TOT_ERR = 0; TOT_W4 = 0;
    for (int w = 0; w <= MAXW; ++w) TOT_BY_W[w] = 0;
    for_each_error([&](uint32_t E, int w) {
        ALLERR.push_back(E); ALLW.push_back((uint8_t)w);
        ++TOT_ERR; ++TOT_BY_W[w]; if (w == MAXW) ++TOT_W4;
    });
}

// ===========================================================================================
//  LINEAR ALGEBRA OVER F_2, THE STAGE-1 SOLUTION OBJECT, AND ITS EXACT SCORE
// ===========================================================================================
// row-echelon rank over F_2, pivots indexed by highest set bit
static int rank_of(const uint32_t* v, int n) {
    uint32_t piv[32] = {0}; int r = 0;
    for (int i = 0; i < n; ++i) {
        uint32_t x = v[i];
        while (x) {
            int b = 31 - std::countl_zero(x);
            if (!piv[b]) { piv[b] = x; ++r; break; }
            x ^= piv[b];
        }
    }
    return r;
}
// greedily pick an independent subset of size <= want, writing the ORIGINAL vectors to out
static int extract_independent(const uint32_t* v, int n, uint32_t* out, int want) {
    uint32_t piv[32] = {0}; int r = 0;
    for (int i = 0; i < n && r < want; ++i) {
        uint32_t x = v[i];
        while (x) {
            int b = 31 - std::countl_zero(x);
            if (!piv[b]) { piv[b] = x; out[r++] = v[i]; break; }
            x ^= piv[b];
        }
    }
    return r;
}
// is x in the span of the r independent vectors v?
static bool in_span(const uint32_t* v, int r, uint32_t x) {
    uint32_t piv[32] = {0};
    for (int i = 0; i < r; ++i) {
        uint32_t y = v[i];
        while (y) { int b = 31 - std::countl_zero(y);
            if (!piv[b]) { piv[b] = y; break; } y ^= piv[b]; }
    }
    while (x) { int b = 31 - std::countl_zero(x); if (!piv[b]) return false; x ^= piv[b]; }
    return true;
}

struct S1 {
    uint32_t g[R1];
    int64_t  sumP = 0;          // SUM_{h in span} P(wt h)      -- minimise
    int64_t  U8   = 0;          // sumP / 256  -- undetected weight-<=4 errors
    int32_t  A[NQ + 1];         // weight enumerator of the span (A[0] = 1)
};

// enumerate the 2^r span of g by Gray code; fills elems (optional), sumP and A
static void span_stats(const uint32_t* g, int r, int64_t& sumP, int32_t* A, uint32_t* elems) {
    const int n = 1 << r;
    if (A) for (int m = 0; m <= NQ; ++m) A[m] = 0;
    uint32_t cur = 0; sumP = PTAB[0]; if (A) A[0] = 1; if (elems) elems[0] = 0;
    for (int i = 1; i < n; ++i) {
        cur ^= g[std::countr_zero((unsigned)i)];
        int m = pwt(cur); sumP += PTAB[m]; if (A) ++A[m];
        if (elems) elems[i ^ (i >> 1)] = cur;
    }
}
static void fill_s1(S1& S) {
    span_stats(S.g, R1, S.sumP, S.A, nullptr);
    S.U8 = S.sumP / S1SZ;
}
static bool is_isotropic(const uint32_t* g, int r) {
    for (int i = 0; i < r; ++i) for (int j = i + 1; j < r; ++j) if (symp(g[i], g[j])) return false;
    return true;
}
// exact admissibility of an 8-dimensional subspace: its weight-10 elements must span it
static bool w10_spans(const uint32_t* elems, int n, uint32_t* basis_out) {
    uint32_t cand[512]; int nc = 0;
    for (int i = 1; i < n; ++i) if (pwt(elems[i]) == W1 && nc < 512) cand[nc++] = elems[i];
    uint32_t tmp[R1];
    int r = extract_independent(cand, nc, tmp, R1);
    if (r != R1) return false;
    if (basis_out) for (int i = 0; i < R1; ++i) basis_out[i] = tmp[i];
    return true;
}
static bool s1_valid(const S1& S, std::string* why) {
    for (int i = 0; i < R1; ++i) if (pwt(S.g[i]) != W1) {
        if (why) *why = "generator weight != 10"; return false; }
    if (!is_isotropic(S.g, R1)) { if (why) *why = "generators do not commute"; return false; }
    if (rank_of(S.g, R1) != R1) { if (why) *why = "generators dependent"; return false; }
    return true;
}

// ===========================================================================================
//  THE QUOTIENT  Q = K^perp / K   (used for both the Stage-1 k-opt and the Stage-2 sweep)
// ===========================================================================================
struct QSpace {
    int t = 0, h = 0, dim = 0, size = 0;
    uint32_t kb[R1 + 1];                  // basis of K
    uint32_t UB[NQ], WB[NQ];              // hyperbolic pairs
    uint32_t basisv[2 * NQ];              // UB[0..h-1] then WB[0..h-1]; class bit = index
    std::vector<uint32_t> lift;           // class -> a Pauli representative
    std::vector<int64_t>  phi;            // class -> SUM_{u in K} P(wt(lift ^ u))
    std::vector<uint32_t> repW;           // class -> smallest representative of weight wantW
    int64_t sK = 0;
    uint32_t lowmask = 0;
};
static inline int qsymp(const QSpace& Q, uint32_t a, uint32_t b) {
    uint32_t m = ((b >> Q.h) & Q.lowmask) | ((b & Q.lowmask) << Q.h);
    return std::popcount(a & m) & 1;
}
static uint32_t qcls(const QSpace& Q, uint32_t v) {
    uint32_t c = 0;
    for (int i = 0; i < Q.h; ++i) {
        c |= uint32_t(symp(Q.WB[i], v)) << i;
        c |= uint32_t(symp(Q.UB[i], v)) << (i + Q.h);
    }
    return c;
}
// keep[0..t-1] must be independent and isotropic.  withPhi builds the phi table (needs P[]).
// wantW > 0 also records, for every class, the lexicographically smallest lift of that weight.
static bool build_qspace(const uint32_t* keep, int t, QSpace& Q, bool withPhi, int wantW) {
    Q.t = t;
    for (int i = 0; i < t; ++i) Q.kb[i] = keep[i];
    std::vector<uint32_t> b;
    for (int j = 0; j < NQ; ++j) { b.push_back(pk(1u << j, 0)); b.push_back(pk(0, 1u << j)); }
    for (int gI = 0; gI < t; ++gI) {
        int piv = -1;
        for (size_t i = 0; i < b.size(); ++i) if (symp(keep[gI], b[i])) { piv = (int)i; break; }
        if (piv < 0) return false;                       // generator was zero / dependent
        uint32_t pv = b[piv]; b.erase(b.begin() + piv);
        for (auto& q : b) if (symp(keep[gI], q)) q ^= pv;
    }
    if ((int)b.size() != 2 * NQ - t) return false;
    int np = 0;
    while (true) {
        int pi = -1, qi = -1;
        for (size_t i = 0; i < b.size() && pi < 0; ++i)
            for (size_t j = i + 1; j < b.size(); ++j)
                if (symp(b[i], b[j])) { pi = (int)i; qi = (int)j; break; }
        if (pi < 0) break;
        uint32_t u = b[pi], w = b[qi];
        b.erase(b.begin() + qi); b.erase(b.begin() + pi);
        for (auto& r : b) { if (symp(r, w)) r ^= u; if (symp(r, u)) r ^= w; }
        Q.UB[np] = u; Q.WB[np] = w; ++np;
    }
    Q.h = np; Q.dim = 2 * np; Q.size = 1 << Q.dim; Q.lowmask = (1u << Q.h) - 1u;
    if (Q.dim != 2 * NQ - 2 * t) return false;
    if ((int)b.size() != t) return false;
    // the radical must be exactly K
    {
        std::vector<uint32_t> kel(1 << t, 0);
        uint32_t cur = 0;
        for (int i = 1; i < (1 << t); ++i) { cur ^= keep[std::countr_zero((unsigned)i)];
            kel[i ^ (i >> 1)] = cur; }
        for (uint32_t r : b) if (std::find(kel.begin(), kel.end(), r) == kel.end()) return false;
    }
    for (int i = 0; i < Q.h; ++i) { Q.basisv[i] = Q.UB[i]; Q.basisv[Q.h + i] = Q.WB[i]; }
    Q.lift.assign(Q.size, 0);
    { uint32_t cur = 0;
      for (int i = 1; i < Q.size; ++i) { cur ^= Q.basisv[std::countr_zero((unsigned)i)];
          Q.lift[i ^ (i >> 1)] = cur; } }
    if (withPhi || wantW > 0) {
        if (withPhi) Q.phi.assign(Q.size, 0);
        if (wantW > 0) Q.repW.assign(Q.size, NONE);
        std::vector<uint32_t> gb, gc;
        for (int i = 0; i < t; ++i) { gb.push_back(keep[i]); gc.push_back(0); }
        for (int i = 0; i < Q.dim; ++i) { gb.push_back(Q.basisv[i]); gc.push_back(1u << i); }
        const int nb = (int)gb.size();
        uint32_t cur = 0, ccl = 0;
        if (withPhi) Q.phi[0] += PTAB[0];
        for (uint64_t i = 1; i < (1ull << nb); ++i) {
            int bb = std::countr_zero(i);
            cur ^= gb[bb]; ccl ^= gc[bb];
            int w = pwt(cur);
            if (withPhi) Q.phi[ccl] += PTAB[w];
            if (wantW > 0 && w == wantW && (Q.repW[ccl] == NONE || cur < Q.repW[ccl]))
                Q.repW[ccl] = cur;
        }
        if (withPhi) Q.sK = Q.phi[0];
    }
    return true;
}

// ===========================================================================================
//  COLUMN-TYPE MODULE  (section 3 of the header:  the mandated qubit-structured formulation)
// ===========================================================================================
struct ColumnView {
    uint32_t x[NQ], z[NQ];        // columns of X and Z as elements of F_2^8
    int      dimW[NQ];            // dim span{x_j, z_j}
    uint64_t plucker[NQ];         // C(8,2) = 28 bits
    uint64_t pluckerXor = 0;
    int      zerosPerGen[R1];     // #{ j : W_j <= H_i }
    int      spanRank = 0;        // dim (W_1 + ... + W_14)
};
static void column_view(const uint32_t* g, ColumnView& C) {
    for (int j = 0; j < NQ; ++j) {
        uint32_t xa = 0, za = 0;
        for (int i = 0; i < R1; ++i) {
            if ((px(g[i]) >> j) & 1) xa |= 1u << i;
            if ((pz(g[i]) >> j) & 1) za |= 1u << i;
        }
        C.x[j] = xa; C.z[j] = za;
        C.dimW[j] = (xa == 0 && za == 0) ? 0 : ((xa == za || xa == 0 || za == 0) ? 1 : 2);
        uint64_t pl = 0; int bit = 0;
        for (int a = 0; a < R1; ++a) for (int c = a + 1; c < R1; ++c) {
            int v = (((xa >> a) & 1) & ((za >> c) & 1)) ^ (((xa >> c) & 1) & ((za >> a) & 1));
            if (v) pl |= 1ull << bit;
            ++bit;
        }
        C.plucker[j] = pl;
    }
    C.pluckerXor = 0; for (int j = 0; j < NQ; ++j) C.pluckerXor ^= C.plucker[j];
    for (int i = 0; i < R1; ++i) {
        int c = 0;
        for (int j = 0; j < NQ; ++j) if (!(((C.x[j] >> i) & 1) | ((C.z[j] >> i) & 1))) ++c;
        C.zerosPerGen[i] = c;
    }
    uint32_t all[2 * NQ];
    for (int j = 0; j < NQ; ++j) { all[2*j] = C.x[j]; all[2*j+1] = C.z[j]; }
    C.spanRank = rank_of(all, 2 * NQ);
}
// weight of the element a of the span, computed ONLY from the column multiset
static int weight_from_columns(const ColumnView& C, uint32_t a) {
    int triv = 0;
    for (int j = 0; j < NQ; ++j)
        if (!par(a & C.x[j]) && !par(a & C.z[j])) ++triv;
    return NQ - triv;
}
static int64_t gauss_binom(int n, int k) {          // number of k-subspaces of F_2^n
    long double num = 1.0L, den = 1.0L;
    for (int i = 0; i < k; ++i) {
        num *= (powl(2.0L, n - i) - 1.0L);
        den *= (powl(2.0L, k - i) - 1.0L);
    }
    return (int64_t)llroundl(num / den);
}

// ===========================================================================================
//  SYMMETRY FINGERPRINT OF A STAGE-1 SOLUTION
//
//  Two Stage-1 codes are equivalent when a qubit permutation composed with per-qubit local
//  Cliffords carries one span onto the other.  Both operations preserve, for every element of
//  the span, its Pauli weight and its exact support.  The fingerprint below is therefore a
//  genuine invariant of the equivalence class:
//
//      (a) the weight enumerator A[0..14] of the 256-element span;
//      (b) for every qubit j the vector n_j[m] = #{h in L : wt(h)=m, h acts on j}, with the
//          14 vectors sorted (kills the qubit permutation);
//      (c) for every qubit pair {j,k} the vector n_{jk}[m] = #{h : wt(h)=m, h acts on both},
//          with the 91 vectors sorted.
//
//  DIRECTION OF THE ERROR.  Equal fingerprints do not PROVE equivalence; different
//  fingerprints DO prove inequivalence.  So the reported class count is a lower bound on the
//  number of inequivalent solutions found, and the program guards against over-merging by
//  keeping several distinct members per fingerprint and running Stage 2 on each of them: if
//  two members of one fingerprint class produce different Stage-2 optima or different
//  residual structure they are provably inequivalent and the program says so.
// ===========================================================================================
static std::string fingerprint(const uint32_t* g) {
    uint32_t el[S1SZ]; int64_t sp; int32_t A[NQ + 1];
    span_stats(g, R1, sp, A, el);
    static const int NPAIR = NQ * (NQ - 1) / 2;
    std::vector<std::array<int32_t, NQ + 1>> nq(NQ), npr(NPAIR);
    for (auto& v : nq)  v.fill(0);
    for (auto& v : npr) v.fill(0);
    for (int i = 1; i < S1SZ; ++i) {
        uint32_t h = el[i]; uint32_t sup = (px(h) | pz(h)) & QM; int m = pwt(h);
        for (int j = 0; j < NQ; ++j) if ((sup >> j) & 1) {
            ++nq[j][m];
            int base = j * (2 * NQ - j - 1) / 2 - j - 1;
            for (int k = j + 1; k < NQ; ++k) if ((sup >> k) & 1) ++npr[base + k][m];
        }
    }
    std::sort(nq.begin(), nq.end());
    std::sort(npr.begin(), npr.end());
    std::string s;
    auto push = [&](int32_t v) { char b[4]; memcpy(b, &v, 4); s.append(b, 4); };
    for (int m = 0; m <= NQ; ++m) push(A[m]);
    for (auto& v : nq)  for (int m = 0; m <= NQ; ++m) push(v[m]);
    for (auto& v : npr) for (int m = 0; m <= NQ; ++m) push(v[m]);
    return s;
}

// ===========================================================================================
//  STAGE-1 CONSTRUCTION AND LOCAL SEARCH
// ===========================================================================================
static bool random_start(S1& S, std::mt19937_64& rng, int ncand) {
    std::vector<uint32_t> pb;
    for (int j = 0; j < NQ; ++j) { pb.push_back(pk(1u << j, 0)); pb.push_back(pk(0, 1u << j)); }
    uint32_t el[S1SZ]; el[0] = 0; int nel = 1;
    int64_t cur = PTAB[0];
    for (int m = 0; m < R1; ++m) {
        uint32_t bestx = NONE; int64_t bestv = 0;
        int tries = 0, found = 0;
        while (found < ncand && tries < ncand * 200) {
            ++tries;
            uint64_t mask = rng();
            uint32_t x = 0;
            for (size_t i = 0; i < pb.size(); ++i) if ((mask >> (i & 63)) & 1) x ^= pb[i];
            if (pwt(x) != W1) continue;
            bool dup = false;
            for (int i = 0; i < nel; ++i) if (el[i] == x) { dup = true; break; }
            if (dup) continue;
            ++found;
            int64_t v = 0;
            for (int i = 0; i < nel; ++i) v += PTAB[pwt(x ^ el[i])];
            if (bestx == NONE || v < bestv) { bestv = v; bestx = x; }
        }
        if (bestx == NONE) return false;
        S.g[m] = bestx; cur += bestv;
        for (int i = 0; i < nel; ++i) el[nel + i] = el[i] ^ bestx;
        nel <<= 1;
        int piv = -1;
        for (size_t i = 0; i < pb.size(); ++i) if (symp(bestx, pb[i])) { piv = (int)i; break; }
        if (piv < 0) return false;
        uint32_t pv = pb[piv]; pb.erase(pb.begin() + piv);
        for (auto& q : pb) if (symp(bestx, q)) q ^= pv;
    }
    fill_s1(S);
    return s1_valid(S, nullptr) && S.sumP == cur;
}

// build the 8-dim subspace K + span(lifts) and, if its weight-10 elements span it, write a
// weight-10 basis into out.  This is the exact admissibility test of section 4.
static bool realise(const QSpace& Q, const uint32_t* qcl, int k, uint32_t* out) {
    uint32_t basis[R1];
    for (int i = 0; i < Q.t; ++i) basis[i] = Q.kb[i];
    for (int i = 0; i < k; ++i) basis[Q.t + i] = Q.lift[qcl[i]];
    if (rank_of(basis, R1) != R1) return false;
    uint32_t el[S1SZ]; int64_t sp; span_stats(basis, R1, sp, nullptr, el);
    return w10_spans(el, S1SZ, out);
}

struct Cand { int64_t v; uint32_t c; };

// exhaustive 1-opt: every 8-dim isotropic D' containing seven of the current generators.
// Returns true and overwrites S when a strictly better admissible D' exists.
static bool one_opt(S1& S, uint64_t* nodes) {
    int64_t bestv = S.sumP; uint32_t bg[R1]; bool got = false;
    for (int drop = 0; drop < R1; ++drop) {
        uint32_t keep[R1 - 1]; int t = 0;
        for (int i = 0; i < R1; ++i) if (i != drop) keep[t++] = S.g[i];
        QSpace Q;
        if (!build_qspace(keep, R1 - 1, Q, true, 0)) continue;
        std::vector<Cand> cs;
        for (int c = 1; c < Q.size; ++c) {
            int64_t v = Q.sK + Q.phi[c];
            if (v < bestv) cs.push_back({v, (uint32_t)c});
        }
        if (nodes) *nodes += Q.size;
        std::sort(cs.begin(), cs.end(), [](const Cand& a, const Cand& b) { return a.v < b.v; });
        for (const Cand& cd : cs) {
            if (cd.v >= bestv) break;
            uint32_t one = cd.c, ob[R1];
            if (!realise(Q, &one, 1, ob)) continue;
            bestv = cd.v; for (int i = 0; i < R1; ++i) bg[i] = ob[i]; got = true;
            break;
        }
    }
    if (!got) return false;
    for (int i = 0; i < R1; ++i) S.g[i] = bg[i];
    fill_s1(S);
    return true;
}

// exhaustive 2-opt over all C(8,2) = 28 retained sextuples.  The filter is exact: with
// phimin = min over non-zero classes of phi, a pair can only beat the incumbent if
// phi(c1) < need - sK - 2*phimin and phi(c2) < need - sK - phi(c1) - phimin.
static bool two_opt(S1& S, uint64_t* nodes, uint64_t* pruned) {
    int64_t bestv = S.sumP; uint32_t bg[R1]; bool got = false;
    for (int da = 0; da < R1; ++da) for (int db = da + 1; db < R1; ++db) {
        uint32_t keep[R1 - 2]; int t = 0;
        for (int i = 0; i < R1; ++i) if (i != da && i != db) keep[t++] = S.g[i];
        QSpace Q;
        if (!build_qspace(keep, R1 - 2, Q, true, 0)) continue;
        std::vector<uint32_t> ord;
        ord.reserve(Q.size);
        for (int c = 1; c < Q.size; ++c) ord.push_back((uint32_t)c);
        std::sort(ord.begin(), ord.end(),
                  [&](uint32_t a, uint32_t b) { return Q.phi[a] < Q.phi[b]; });
        // pre[k] = sum of the k smallest phi over non-zero classes.  Any k DISTINCT non-zero
        // classes have phi-sum at least pre[k] -- valid whatever else has been chosen, since
        // excluding classes can only raise the sum.  This is far tighter than k*phimin.
        int64_t pre[4] = {0,0,0,0};
        for (int k = 1; k <= 3 && k <= (int)ord.size(); ++k) pre[k] = pre[k-1] + Q.phi[ord[k-1]];
        for (uint32_t c1 : ord) {
            int64_t rest = bestv - Q.sK - Q.phi[c1];
            if (Q.phi[c1] + pre[2] >= bestv - Q.sK) { if (pruned) ++*pruned; break; }
            for (uint32_t c2 : ord) {
                if (Q.phi[c2] + pre[1] >= rest) break;
                if (c2 <= c1) continue;
                if (qsymp(Q, c1, c2)) continue;
                if (nodes) ++*nodes;
                int64_t v = Q.sK + Q.phi[c1] + Q.phi[c2] + Q.phi[c1 ^ c2];
                if (v >= bestv) continue;
                uint32_t two[2] = { c1, c2 }, ob[R1];
                if (!realise(Q, two, 2, ob)) continue;
                bestv = v; for (int i = 0; i < R1; ++i) bg[i] = ob[i]; got = true;
                rest = bestv - Q.sK - Q.phi[c1];
            }
        }
    }
    if (!got) return false;
    for (int i = 0; i < R1; ++i) S.g[i] = bg[i];
    fill_s1(S);
    return true;
}

// exhaustive 3-opt over all C(8,3) = 56 retained quintuples.  Q has dimension 18 here, so a
// naive triple loop over 262144 classes is hopeless.  The fix is to telescope one level:
// with psi(w) = phi(w) + phi(w ^ c1) the seven non-zero classes of span(c1,c2,c3) regroup as
//
//     S = sK + phi(c1) + psi(c2) + psi(c3) + psi(c2 ^ c3) ,
//
// i.e. once c1 is fixed the remaining problem is exactly the cheap three-term pair problem
// that 2-opt already solves.  Let pre2[k] be the sum of the k smallest psi over the classes
// other than 0 and c1.  If a triple beats the incumbent then psi(c2)+psi(c3)+psi(c2^c3) is
// below need2 = bestv - sK - phi(c1), and since any TWO of those three distinct classes sum
// to at least pre2[2], each one individually satisfies psi(x) < need2 - pre2[2].  So all
// three live in the single candidate list L = { x : psi(x) + pre2[2] < need2 }, which is
// built once per c1 and is normally tiny.  Sorting L by psi ascending then gives exact
// breaks at both remaining levels.  Nothing that could beat the incumbent is discarded.
static bool three_opt(S1& S, uint64_t* nodes, uint64_t* pruned, double deadline, bool* complete,
                      int nthreads, uint64_t* c1done, uint64_t* c1total) {
    int64_t bestv = S.sumP; uint32_t bg[R1]; bool got = false;
    if (complete) *complete = true;
    std::vector<std::array<int,3>> drops;
    for (int da = 0; da < R1; ++da)
      for (int db = da + 1; db < R1; ++db)
        for (int dc = db + 1; dc < R1; ++dc) drops.push_back({da, db, dc});
    std::mutex bmtx;
    std::atomic<size_t> nextDrop{0};
    std::atomic<uint64_t> anodes{0}, apruned{0}, ac1{0};
    std::atomic<bool> incomplete{false};
    auto worker = [&]() {
      for (;;) {
        size_t di = nextDrop.fetch_add(1);
        if (di >= drops.size()) break;
        const int da = drops[di][0], db = drops[di][1], dc = drops[di][2];
        if (deadline > 0 && elapsed() > deadline) { incomplete = true; break; }
        uint32_t keep[R1 - 3]; int t = 0;
        uint32_t basis[R1];
        { std::lock_guard<std::mutex> lk(bmtx);
          for (int i = 0; i < R1; ++i) basis[i] = S.g[i]; }
        for (int i = 0; i < R1; ++i) if (i != da && i != db && i != dc) keep[t++] = basis[i];
        QSpace Q;
        if (!build_qspace(keep, R1 - 3, Q, true, 0)) continue;
        std::vector<uint32_t> ordPhi; ordPhi.reserve(Q.size);
        for (int c = 1; c < Q.size; ++c) ordPhi.push_back((uint32_t)c);
        std::sort(ordPhi.begin(), ordPhi.end(),
                  [&](uint32_t a, uint32_t b) { return Q.phi[a] < Q.phi[b]; });
        int64_t prePhi = 0;                       // sum of the 6 smallest phi
        for (int k = 0; k < 6 && k < (int)ordPhi.size(); ++k) prePhi += Q.phi[ordPhi[k]];
        std::vector<int64_t> psi(Q.size);
        std::vector<uint32_t> L;
        uint64_t seen = 0;
        size_t ci = 0;
        for (uint32_t c1 : ordPhi) {
            ++ci;
            int64_t bnow;
            { std::lock_guard<std::mutex> lk(bmtx); bnow = bestv; }
            // the list is sorted by phi ascending, so this break retires the whole tail
            if (Q.phi[c1] + prePhi >= bnow - Q.sK) {
                ++apruned; ac1 += (uint64_t)(ordPhi.size() - ci + 1); break; }
            ++ac1;
            if ((++seen & 0x3FF) == 0 && deadline > 0 && elapsed() > deadline) {
                incomplete = true; break; }
            const int64_t need2 = bnow - Q.sK - Q.phi[c1];
            for (int c = 0; c < Q.size; ++c) psi[c] = Q.phi[c] + Q.phi[c ^ c1];
            int64_t m1 = INT64_MAX, m2 = INT64_MAX;
            for (int c = 1; c < Q.size; ++c) {
                if ((uint32_t)c == c1) continue;
                int64_t v = psi[c];
                if (v < m1) { m2 = m1; m1 = v; } else if (v < m2) m2 = v;
            }
            if (m2 == INT64_MAX) continue;
            const int64_t pre2 = m1 + m2;
            L.clear();
            for (int c = 1; c < Q.size; ++c) {
                if ((uint32_t)c == c1) continue;
                if (psi[c] + pre2 < need2) L.push_back((uint32_t)c);
            }
            if (L.size() < 2) continue;
            std::sort(L.begin(), L.end(),
                      [&](uint32_t a, uint32_t b) { return psi[a] < psi[b]; });
            for (size_t i2 = 0; i2 < L.size(); ++i2) {
                uint32_t c2 = L[i2];
                if (psi[c2] + pre2 >= need2) break;
                if (qsymp(Q, c1, c2)) continue;
                for (size_t i3 = i2 + 1; i3 < L.size(); ++i3) {
                    uint32_t c3 = L[i3];
                    if (psi[c3] + m1 >= need2 - psi[c2]) break;
                    if (c3 == (c1 ^ c2)) continue;
                    if (qsymp(Q, c1, c3) || qsymp(Q, c2, c3)) continue;
                    ++anodes;
                    int64_t v = Q.sK + Q.phi[c1] + psi[c2] + psi[c3] + psi[c2 ^ c3];
                    {
                        std::lock_guard<std::mutex> lk(bmtx);
                        if (v >= bestv) continue;
                    }
                    uint32_t three[3] = { c1, c2, c3 }, ob[R1];
                    if (!realise(Q, three, 3, ob)) continue;
                    std::lock_guard<std::mutex> lk(bmtx);
                    if (v < bestv) {
                        bestv = v; for (int i = 0; i < R1; ++i) bg[i] = ob[i]; got = true;
                    }
                }
            }
        }
      }
    };
    {
        int nt = std::max(1, nthreads);
        std::vector<std::thread> th;
        for (int i = 0; i < nt; ++i) th.emplace_back(worker);
        for (auto& t : th) t.join();
    }
    if (nodes) *nodes += anodes.load();
    if (pruned) *pruned += apruned.load();
    if (c1done) *c1done = ac1.load();
    if (c1total) *c1total = (uint64_t)drops.size() * (uint64_t)((1u << (2*NQ-2*(R1-3))) - 1);
    if (complete && incomplete.load()) *complete = false;
    if (!got) return false;
    for (int i = 0; i < R1; ++i) S.g[i] = bg[i];
    fill_s1(S);
    return true;
}

// ===========================================================================================
//  STAGE 2 -- EXHAUSTIVE, AND THEREFORE A PROVEN OPTIMUM ONCE STAGE 1 IS FROZEN
//
//  Freeze L8.  A Stage-2 generator must lie in L8^perp, have weight 8, and be independent of
//  L8; the three of them must commute.  Since the symplectic form descends to
//  Q2 = L8^perp / L8 (dim 28 - 16 = 12, only 4096 classes), the admissible Stage-2 additions
//  are exactly the 3-dimensional totally isotropic subspaces of Q2 spanned by classes that
//  carry a weight-8 lift.  That space is small enough to enumerate completely.
//
//  Coverage identity.  For e in U8 the syndrome bit <h,e> depends only on the class of h, so
//  with cov1(q) = #{e in U8 : <q,e> = 1},
//
//        #covered by a d-dimensional D  =  (1 / 2^{d-1}) * SUM_{q in D} cov1(q) ,
//
//  which for d = 3 is a quarter of the sum over the seven non-zero classes.  Derivation:
//  #{e : <q,e> = 0 for all q in D} = (1/2^d) SUM_{q in D} SUM_e (-1)^{<q,e>}
//                                  = (1/2^d) SUM_q (|U8| - 2 cov1(q)) .
//  cov1 itself is obtained by one Walsh-Hadamard transform over the 4096 classes.
//
//  Prune (admissible): with base = cov1(a)+cov1(b)+cov1(a^b), the remaining four terms are
//  each at most COVMAX, so base + 4*COVMAX <= incumbent kills the pair outright.
// ===========================================================================================
struct S2Result {
    uint32_t h[R2] = {0,0,0};
    int64_t  cover = -1;            // errors of U8 detected by the three new generators
    int64_t  residual = -1;         // |U8| - cover
    uint64_t nodes = 0;
    int64_t  nOptimalSubspaces = 0;
    bool     exhaustive = false;
    bool     found = false;
    int      usable = 0, qdim = 0, qsize = 0;
    int      covmax = 0;
};

static void stage2_exhaustive(const S1& S, const std::vector<uint32_t>& U8,
                              S2Result& R, int threads) {
    QSpace Q;
    if (!build_qspace(S.g, R1, Q, false, W2)) { R.found = false; return; }
    R.qdim = Q.dim; R.qsize = Q.size;
    const int h = Q.h; const uint32_t low = Q.lowmask;
    // ---- cov1 by Walsh-Hadamard over the classes
    std::vector<int32_t> f(Q.size, 0);
    for (uint32_t e : U8) {
        uint32_t c = qcls(Q, e);
        f[((c >> h) & low) | ((c & low) << h)] += 1;
    }
    for (int b = 0; b < Q.dim; ++b) {
        int step = 1 << b;
        for (int i = 0; i < Q.size; i += step << 1)
            for (int k = i; k < i + step; ++k) {
                int32_t a = f[k], d = f[k + step]; f[k] = a + d; f[k + step] = a - d; }
    }
    const int32_t N = (int32_t)U8.size();
    std::vector<int32_t> cov1(Q.size);
    int covmax = 0;
    for (int v = 0; v < Q.size; ++v) {
        cov1[v] = (N - f[v]) / 2;
        if (v && cov1[v] > covmax) covmax = cov1[v];
    }
    R.covmax = covmax;
    std::vector<uint32_t> usable;
    for (int c = 1; c < Q.size; ++c) if (Q.repW[c] != NONE) usable.push_back((uint32_t)c);
    R.usable = (int)usable.size();
    if (usable.size() < 3) { R.exhaustive = true; R.found = false; return; }

    std::mutex mtx;
    int64_t bestTot = -1;
    uint32_t bestTrip[3] = {0,0,0};
    std::set<std::array<uint32_t, 7>> optSet;
    std::atomic<uint64_t> nodes{0};
    std::atomic<size_t> next{0};

    auto worker = [&]() {
        std::vector<uint32_t> basis;
        uint64_t localNodes = 0;
        for (;;) {
            size_t ia = next.fetch_add(1);
            if (ia >= usable.size()) break;
            uint32_t a = usable[ia];
            for (size_t ib = ia + 1; ib < usable.size(); ++ib) {
                uint32_t b = usable[ib];
                if (qsymp(Q, a, b)) continue;
                uint32_t ab = a ^ b;
                int64_t base = (int64_t)cov1[a] + cov1[b] + cov1[ab];
                int64_t snap;
                { std::lock_guard<std::mutex> lk(mtx); snap = bestTot; }
                if (base + 4LL * covmax < snap) continue;
                // basis of a^perp AND b^perp inside Q  (dim 10)
                basis.clear();
                for (int i = 0; i < Q.dim; ++i) basis.push_back(1u << i);
                for (uint32_t g : { a, b }) {
                    int piv = -1;
                    for (size_t i = 0; i < basis.size(); ++i)
                        if (qsymp(Q, g, basis[i])) { piv = (int)i; break; }
                    if (piv < 0) continue;
                    uint32_t pv = basis[piv]; basis.erase(basis.begin() + piv);
                    for (auto& q : basis) if (qsymp(Q, g, q)) q ^= pv;
                }
                const size_t nb = basis.size();
                uint32_t cur = 0;
                for (uint64_t i = 1; i < (1ull << nb); ++i) {
                    cur ^= basis[std::countr_zero(i)];
                    uint32_t c = cur;
                    if (c == a || c == b || c == ab) continue;
                    if (Q.repW[c] == NONE) continue;
                    ++localNodes;
                    int64_t tot = base + cov1[c] + cov1[c ^ a] + cov1[c ^ b] + cov1[c ^ ab];
                    if (tot < snap) continue;
                    std::array<uint32_t, 7> key = { a, b, ab, c, (uint32_t)(c^a),
                                                    (uint32_t)(c^b), (uint32_t)(c^ab) };
                    std::sort(key.begin(), key.end());
                    std::lock_guard<std::mutex> lk(mtx);
                    if (tot > bestTot) {
                        bestTot = tot; optSet.clear();
                        bestTrip[0] = a; bestTrip[1] = b; bestTrip[2] = c;
                    }
                    if (tot == bestTot && optSet.size() < 200000) optSet.insert(key);
                    snap = bestTot;
                }
            }
        }
        nodes += localNodes;
    };
    int nt = std::max(1, threads);
    std::vector<std::thread> th;
    for (int i = 0; i < nt; ++i) th.emplace_back(worker);
    for (auto& t : th) t.join();

    R.nodes = nodes.load();
    R.exhaustive = true;
    if (bestTot < 0) { R.found = false; return; }
    R.found = true;
    R.nOptimalSubspaces = (int64_t)optSet.size();
    R.cover = bestTot / 4;
    R.residual = (int64_t)U8.size() - R.cover;
    for (int i = 0; i < R2; ++i) R.h[i] = Q.repW[bestTrip[i]];
}

// ===========================================================================================
//  RESIDUAL ANALYSIS AND FULLY INDEPENDENT VERIFICATION
// ===========================================================================================
struct ResidualInfo {
    int64_t n = 0, by_w[MAXW + 1] = {0,0,0,0,0};
    bool subgroup = false, all_x = false, subgroupTested = false;
    int dim = 0;
    std::vector<uint32_t> gens, list, qubits;
};
static ResidualInfo analyse(const std::vector<uint32_t>& Rv) {
    ResidualInfo r; r.list = Rv; r.n = (int64_t)Rv.size();
    for (uint32_t E : Rv) { int w = pwt(E); if (w <= MAXW) r.by_w[w]++; }
    uint32_t basis[64]; int d = extract_independent(Rv.data(), (int)Rv.size(), basis, 32);
    r.dim = d;
    for (int i = 0; i < d; ++i) r.gens.push_back(basis[i]);
    // The subgroup test enumerates the span, so it is only attempted when that is cheap.
    // A residual of a few hundred errors can have span dimension up to 20; 2^16 is the point
    // beyond which the test costs more than it tells us, and it is reported as "not tested".
    r.subgroupTested = (d <= 16 && !Rv.empty());
    if (r.subgroupTested) {
        std::set<uint32_t> sp; sp.insert(0);
        uint32_t cur = 0;
        for (uint64_t i = 1; i < (1ull << d); ++i) {
            cur ^= basis[std::countr_zero(i)]; sp.insert(cur);
        }
        r.subgroup = (sp.size() == Rv.size() + 1);
        if (r.subgroup) for (uint32_t E : Rv) if (!sp.count(E)) r.subgroup = false;
    }
    r.all_x = true;
    uint32_t qmask = 0;
    for (uint32_t E : Rv) { if (pz(E)) r.all_x = false; qmask |= (px(E) | pz(E)) & QM; }
    if (Rv.empty()) r.all_x = false;
    for (int j = 0; j < NQ; ++j) if ((qmask >> j) & 1) r.qubits.push_back((uint32_t)j);
    return r;
}

struct Verify {
    bool w1_ok = false, w2_ok = false, cg_ok = false, ch_ok = false, rank_ok = false;
    bool u8_ok = false, cover_ok = false, ptab_ok = false;
    int  rank = 0;
    int64_t u8size = 0, det11 = 0, all_total = 0, w4_total = 0, w4_det = 0, w4_und = 0;
    int64_t cover = 0, residual = 0;
    int weights1[R1] = {0}, weights2[R2] = {0};
    std::vector<uint32_t> leftover;
    int64_t det_by_w[MAXW + 1] = {0,0,0,0,0};
};
// Everything here is recomputed from scratch by direct symplectic products; nothing reuses
// the algebraic score.  The two must agree.
static Verify verify_full(const S1& S, const uint32_t* hh, bool haveStage2) {
    Verify V;
    V.w1_ok = true;
    for (int i = 0; i < R1; ++i) { V.weights1[i] = pwt(S.g[i]); if (V.weights1[i] != W1) V.w1_ok = false; }
    if (haveStage2) { V.w2_ok = true;
        for (int a = 0; a < R2; ++a) { V.weights2[a] = pwt(hh[a]); if (V.weights2[a] != W2) V.w2_ok = false; } }
    V.cg_ok = true;
    for (int i = 0; i < R1; ++i) for (int j = i + 1; j < R1; ++j) if (symp(S.g[i], S.g[j])) V.cg_ok = false;
    if (haveStage2) {
        for (int a = 0; a < R2; ++a) for (int i = 0; i < R1; ++i) if (symp(hh[a], S.g[i])) V.cg_ok = false;
        V.ch_ok = true;
        for (int a = 0; a < R2; ++a) for (int b = a + 1; b < R2; ++b) if (symp(hh[a], hh[b])) V.ch_ok = false;
    } else V.ch_ok = true;
    uint32_t all[RTOT]; int n = 0;
    for (int i = 0; i < R1; ++i) all[n++] = S.g[i];
    if (haveStage2) for (int a = 0; a < R2; ++a) all[n++] = hh[a];
    V.rank = rank_of(all, n);
    V.rank_ok = (V.rank == n);
    V.all_total = 0; V.w4_total = 0; V.det11 = 0; V.u8size = 0; V.w4_det = 0;
    for_each_error([&](uint32_t E, int w) {
        ++V.all_total; if (w == MAXW) ++V.w4_total;
        int s1 = 0;
        for (int i = 0; i < R1; ++i) if (symp(S.g[i], E)) { s1 = 1; break; }
        if (!s1) ++V.u8size;
        int s = s1;
        if (!s && haveStage2)
            for (int a = 0; a < R2; ++a) if (symp(hh[a], E)) { s = 1; break; }
        if (s) { ++V.det11; ++V.det_by_w[w]; if (w == MAXW) ++V.w4_det; }
        else V.leftover.push_back(E);
    });
    V.w4_und = V.w4_total - V.w4_det;
    V.u8_ok = (V.u8size == S.U8);
    V.residual = (int64_t)V.leftover.size();
    V.cover = V.u8size - V.residual;
    V.cover_ok = true;
    // independent check of the closed-form identity used everywhere in the search
    {
        int64_t sp; span_stats(S.g, R1, sp, nullptr, nullptr);
        V.ptab_ok = (sp % S1SZ == 0) && (sp / S1SZ == V.u8size);
    }
    return V;
}

// ===========================================================================================
//  COLUMN-MULTISET BRANCH AND BOUND
//
//  This is the same engine that proved the earlier 3-generator and 4-generator weight-8
//  optima, written generically in the rank r and the generator weight w.  A state is a
//  NON-DECREASING sequence of column-type indices, which is exactly a multiset of 14 columns
//  and therefore already quotients out the S_14 qubit permutations.  Column types are
//  subspaces W <= F_2^r of dimension <= 2 (the complete local-Clifford invariant of a column).
//
//    weight       c(e_i) = w for each generator i, where c(chi) = #columns on which chi acts
//    commutation  XOR of the Plucker invariants of the chosen columns = 0
//    rank         c(chi) >= 1 for every non-zero chi
//    objective    minimise SUM_chi P(c(chi)) ; U = that / 2^r
//
//  BOUND (admissible).  With rem columns still to place, every chi already has c(chi) and can
//  only grow, to at most c(chi)+rem.  The generator points chi = e_i are pinned at exactly w.
//  So  P(0) + r*P(w) + SUM_{other chi} min_{0<=d<=rem} P(c(chi)+d)  never exceeds the true
//  optimum of the subtree; pruning on "bound >= incumbent" therefore cannot discard a strictly
//  better solution, and a completed tree is a proof.
//
//  S_r SYMMETRY.  Permuting generator labels permutes F_2^r and hence the column types.  A
//  node is pruned when some permutation makes the sorted image of the partial multiset
//  lexicographically smaller.  Validity: for a canonical complete solution S the k smallest
//  entries of sorted(pi(S)) are entrywise <= sorted(pi(p)) for any k-prefix p, so
//  sorted(pi(p)) < p would force sorted(pi(S)) < S and contradict canonicity.  The test is
//  applied only for the first --sym-depth columns (still valid: it is a subset of the prunes).
//
//  FEASIBILITY.  The type count is 1 + (2^r - 1) + gauss(r,2), i.e. 15 at r=3, 51 at r=4,
//  187 at r=5, 715 at r=6, 2795 at r=7 and 11051 at r=8, and the raw state count is
//  C(types+13, 14).  r <= 5 completes; r = 8 does not, by an enormous margin (printed by
//  --spacesize).  This is why Stage 1 at r=8 is BEST KNOWN and not PROVEN.
// ===========================================================================================
struct BBTypes {
    int r = 0, ntypes = 0, nchi = 0;
    std::vector<uint64_t> contrib;    // bit chi set  <=>  chi is NOT orthogonal to W_t
    std::vector<uint32_t> plu;        // Plucker invariant, C(r,2) bits
    std::vector<int>      dimW;
    std::vector<uint64_t> spanmask;   // the elements of W_t as a bitmask over F_2^r
};
// The per-type "which chi does this column act on" set is a uint64 bitmask over the 2^r
// characters, so the engine is only defined for r <= 6.  That is not a limitation in
// practice: the tree already fails to terminate well before r = 7 (see --spacesize).
static int bb_max_rank() { return 6; }
static int64_t bb_type_count(int r) { return 1 + ((int64_t)1 << r) - 1 + gauss_binom(r, 2); }
static void build_types(int r, BBTypes& T) {
    if (r < 1 || r > bb_max_rank()) {
        fprintf(stderr, "branch and bound is implemented for rank 1..%d only\n", bb_max_rank());
        exit(1);
    }
    T.r = r; T.nchi = 1 << r;
    std::vector<std::pair<uint64_t, std::pair<uint32_t,int>>> tmp;
    std::set<uint64_t> seen;
    auto add = [&](uint32_t x, uint32_t z) {
        uint64_t sm = 1ull;                       // element 0 always in W
        sm |= 1ull << x; sm |= 1ull << z; sm |= 1ull << (x ^ z);
        if (seen.count(sm)) return;
        seen.insert(sm);
        int d = std::popcount(sm) == 1 ? 0 : (std::popcount(sm) == 2 ? 1 : 2);
        uint32_t pl = 0; int bit = 0;
        for (int i = 0; i < r; ++i) for (int j = i + 1; j < r; ++j) {
            int v = (((x >> i) & 1) & ((z >> j) & 1)) ^ (((x >> j) & 1) & ((z >> i) & 1));
            if (v) pl |= 1u << bit;
            ++bit;
        }
        tmp.push_back({ sm, { pl, d } });
    };
    for (uint32_t x = 0; x < (uint32_t)T.nchi; ++x)
        for (uint32_t z = 0; z < (uint32_t)T.nchi; ++z) add(x, z);
    std::sort(tmp.begin(), tmp.end(), [](const auto& a, const auto& b) {
        if (a.second.second != b.second.second) return a.second.second < b.second.second;
        return a.first < b.first; });
    T.ntypes = (int)tmp.size();
    T.contrib.resize(T.ntypes); T.plu.resize(T.ntypes);
    T.dimW.resize(T.ntypes); T.spanmask.resize(T.ntypes);
    for (int t = 0; t < T.ntypes; ++t) {
        T.spanmask[t] = tmp[t].first; T.plu[t] = tmp[t].second.first; T.dimW[t] = tmp[t].second.second;
        uint64_t cm = 0;
        for (int chi = 0; chi < T.nchi; ++chi) {
            bool orth = true;
            for (int v = 0; v < T.nchi; ++v) if ((tmp[t].first >> v) & 1)
                if (par((uint32_t)chi & (uint32_t)v)) { orth = false; break; }
            if (!orth) cm |= 1ull << chi;
        }
        T.contrib[t] = cm;
    }
}
struct BBResult {
    int64_t bestSum = INT64_MAX, U = -1, M = -1, seedIncumbent = INT64_MAX;
    uint64_t nodes = 0, pruned = 0, leaves = 0;
    bool completed = false;
    std::vector<std::array<uint32_t, 16>> gens;   // optimal generator sets (capped)
    int64_t rawOptimal = 0;
    double seconds = 0;
};
static int64_t PMIN_TAB[NQ + 2][NQ + 2];          // PMIN[c][rem] = min_{0<=d<=rem} P[c+d]
static int64_t RATIO_TAB[NQ + 2][NQ + 2];         // 1024 * max_{d<=rem} (P[c]-P[c+d])/d, rounded UP
static void build_pmin() {
    for (int c = 0; c <= NQ + 1; ++c) for (int rem = 0; rem <= NQ + 1; ++rem) {
        int64_t best = INT64_MAX, ratio = 0;
        for (int d = 0; d <= rem; ++d) {
            int m = c + d; if (m > NQ) break;
            best = std::min(best, PTAB[m]);
            if (d > 0) {
                int64_t num = PTAB[c] - PTAB[m];
                if (num > 0) ratio = std::max(ratio, (num * 1024 + d - 1) / d);
            }
        }
        if (best == INT64_MAX) best = PTAB[NQ];
        PMIN_TAB[c][rem] = best;
        RATIO_TAB[c][rem] = ratio;
    }
}
struct BBCtx {
    const BBTypes* T; int r, w, nchi, symDepth, gencap, deltaMax = 1;
    bool keepTies = false;          // explore subtrees that can only equal the incumbent
    std::vector<std::vector<int>> perms;          // type permutations induced by S_r
    int64_t incumbent;
    std::mutex mtx;
    BBResult* R;
};
static void build_typeperms(const BBTypes& T, int r, std::vector<std::vector<int>>& out) {
    out.clear();
    std::vector<int> p(r); std::iota(p.begin(), p.end(), 0);
    std::map<uint64_t,int> idx;
    for (int t = 0; t < T.ntypes; ++t) idx[T.spanmask[t]] = t;
    do {
        std::vector<int> m(T.ntypes);
        for (int t = 0; t < T.ntypes; ++t) {
            uint64_t sm = 0;
            for (int v = 0; v < T.nchi; ++v) if ((T.spanmask[t] >> v) & 1) {
                uint32_t nv = 0;
                for (int i = 0; i < r; ++i) if ((v >> i) & 1) nv |= 1u << p[i];
                sm |= 1ull << nv;
            }
            m[t] = idx[sm];
        }
        out.push_back(m);
    } while (std::next_permutation(p.begin(), p.end()));
}
static bool sym_ok(const BBCtx& C, const int* seq, int k) {
    if (k == 0 || k > C.symDepth) return true;
    std::vector<int> img(k);
    for (const auto& m : C.perms) {
        for (int i = 0; i < k; ++i) img[i] = m[seq[i]];
        std::sort(img.begin(), img.end());
        for (int i = 0; i < k; ++i) {
            if (img[i] < seq[i]) return false;
            if (img[i] > seq[i]) break;
        }
    }
    return true;
}
static double BB_DEADLINE = 0;                 // 0 = no limit
static std::atomic<bool> BB_ABORTED{false};
static std::atomic<int64_t> BB_INC{INT64_MAX}; // incumbent, read lock-free on the hot path
// false: prune subtrees that can only TIE the incumbent (fast; finds the optimum but does not
// enumerate every optimal leaf).  true: keep them, so the optimal-family count is complete.
static bool BB_KEEP_TIES = false;

// A greedy randomised construction at arbitrary rank and generator weight.  Used only to give
// the branch and bound a starting incumbent: any ACHIEVABLE value is a legitimate upper bound
// on the optimum, so seeding can never change which optimum the completed tree reports -- it
// only decides how much of the tree has to be walked.
static bool greedy_rank(uint32_t* g, int r, int w, std::mt19937_64& rng, int ncand) {
    std::vector<uint32_t> pb;
    for (int j = 0; j < NQ; ++j) { pb.push_back(pk(1u << j, 0)); pb.push_back(pk(0, 1u << j)); }
    std::vector<uint32_t> el; el.push_back(0);
    for (int m = 0; m < r; ++m) {
        uint32_t bestx = NONE; int64_t bestv = 0;
        int tries = 0, found = 0;
        while (found < ncand && tries < ncand * 300) {
            ++tries;
            uint64_t mask = rng(); uint32_t x = 0;
            for (size_t i = 0; i < pb.size(); ++i) if ((mask >> (i & 63)) & 1) x ^= pb[i];
            if (pwt(x) != w) continue;
            if (std::find(el.begin(), el.end(), x) != el.end()) continue;
            ++found;
            int64_t v = 0;
            for (uint32_t u : el) v += PTAB[pwt(x ^ u)];
            if (bestx == NONE || v < bestv) { bestv = v; bestx = x; }
        }
        if (bestx == NONE) return false;
        g[m] = bestx;
        size_t n = el.size(); for (size_t i = 0; i < n; ++i) el.push_back(el[i] ^ bestx);
        int piv = -1;
        for (size_t i = 0; i < pb.size(); ++i) if (symp(bestx, pb[i])) { piv = (int)i; break; }
        if (piv < 0) return false;
        uint32_t pv = pb[piv]; pb.erase(pb.begin() + piv);
        for (auto& q : pb) if (symp(bestx, q)) q ^= pv;
    }
    return true;
}
static int64_t bb_seed_incumbent(int r, int w, uint64_t seed, double budget) {
    std::mt19937_64 rng(seed);
    int64_t best = INT64_MAX;
    double t0 = elapsed();
    int iters = 0;
    while (elapsed() - t0 < budget || iters < 50) {
        if (++iters > 200000) break;
        uint32_t g[16];
        if (!greedy_rank(g, r, w, rng, 24)) continue;
        int64_t sp; span_stats(g, r, sp, nullptr, nullptr);
        if (sp < best) best = sp;
        if (elapsed() - t0 > budget && iters >= 50) break;
    }
    return best;
}
static void bb_dfs(BBCtx& C, int depth, int startType, int* seq, int* cnt,
                   uint32_t pluXor, BBResult& L) {
    if (BB_ABORTED.load(std::memory_order_relaxed)) return;
    const BBTypes& T = *C.T;
    const int rem = NQ - depth;
    if (depth == NQ) {
        ++L.leaves;
        if (pluXor) return;
        for (int i = 0; i < C.r; ++i) if (cnt[1 << i] != C.w) return;
        for (int chi = 1; chi < C.nchi; ++chi) if (cnt[chi] == 0) return;
        int64_t s = 0;
        for (int chi = 0; chi < C.nchi; ++chi) s += PTAB[cnt[chi]];
        std::lock_guard<std::mutex> lk(C.mtx);
        if (s < C.R->bestSum) { C.R->bestSum = s; C.R->rawOptimal = 0; C.R->gens.clear();
                                BB_INC.store(s, std::memory_order_relaxed); }
        if (s == C.R->bestSum) {
            ++C.R->rawOptimal;
            if ((int)C.R->gens.size() < C.gencap) {
                std::array<uint32_t,16> g{}; g.fill(0);
                for (int j = 0; j < NQ; ++j) {
                    uint64_t sm = T.spanmask[seq[j]];
                    uint32_t xv = 0, zv = 0; int got = 0;
                    for (int v = 1; v < C.nchi; ++v) if ((sm >> v) & 1) {
                        if (got == 0) { xv = (uint32_t)v; ++got; }
                        else if (got == 1 && (uint32_t)v != xv) { zv = (uint32_t)v; ++got; break; }
                    }
                    if (T.dimW[seq[j]] == 1) zv = 0;
                    for (int i = 0; i < C.r; ++i) {
                        if ((xv >> i) & 1) g[i] |= 1u << j;
                        if ((zv >> i) & 1) g[i] |= 1u << (j + 16);
                    }
                }
                C.R->gens.push_back(g);
            }
        }
        return;
    }
    // ---------------- admissible bound ----------------
    // (a) chi = 0 contributes P(0) and each of the r generator points is pinned at weight w;
    // (b) every remaining column raises SUM_{chi != 0} c(chi) by delta(W) = 2^r - 2^{r-dim W},
    //     at most deltaMax = 3*2^{r-2}, and the generator points already claim
    //     needGen = SUM_i (w - c(e_i)) of that;  what is left is the BUDGET the other
    //     characters can share.  A character starting at count c cannot fall below
    //     PMIN[c][rem], and cannot fall faster than RATIO[c][rem] per unit of budget, so the
    //     total decrease is bounded by both sums.  Taking the smaller keeps the bound valid.
    {
        int64_t needGen = 0;
        for (int i = 0; i < C.r; ++i) {
            int d = C.w - cnt[1 << i];
            if (d < 0 || d > rem) { ++L.pruned; return; }
            needGen += d;
        }
        int64_t budget = (int64_t)C.deltaMax * rem - needGen;
        if (budget < 0) { ++L.pruned; return; }
        int64_t sumP = 0, sumDec = 0, maxRatio = 0, nzero = 0;
        for (int chi = 1; chi < C.nchi; ++chi) {
            if (std::popcount((unsigned)chi) == 1) continue;      // pinned generator points
            int c = cnt[chi];
            if (c == 0) ++nzero;                                  // rank needs c >= 1
            sumP += PTAB[c];
            sumDec += PTAB[c] - PMIN_TAB[c][rem];
            int64_t rr = RATIO_TAB[c][rem];
            if (rr > maxRatio) maxRatio = rr;
        }
        if (nzero > budget) { ++L.pruned; return; }               // cannot reach rank r
        int64_t byBudget = (budget * maxRatio + 1023) / 1024;
        int64_t b = PTAB[0] + (int64_t)C.r * PTAB[C.w] + sumP - std::min(sumDec, byBudget);
        if (b >= BB_INC.load(std::memory_order_relaxed) + (C.keepTies ? 1 : 0)) {
            ++L.pruned; return; }
    }
    for (int t = startType; t < T.ntypes; ++t) {
        uint64_t cm = T.contrib[t];
        bool bad = false;
        for (int i = 0; i < C.r && !bad; ++i) {
            int c = cnt[1 << i] + (int)((cm >> (1 << i)) & 1);
            if (c > C.w || C.w - c > rem - 1) bad = true;
        }
        if (bad) continue;
        seq[depth] = t;
        if (!sym_ok(C, seq, depth + 1)) continue;
        for (int chi = 0; chi < C.nchi; ++chi) cnt[chi] += (int)((cm >> chi) & 1);
        ++L.nodes;
        if (BB_DEADLINE > 0 && (L.nodes & 0xFFFFu) == 0 && elapsed() > BB_DEADLINE)
            BB_ABORTED.store(true, std::memory_order_relaxed);
        bb_dfs(C, depth + 1, t, seq, cnt, pluXor ^ T.plu[t], L);
        for (int chi = 0; chi < C.nchi; ++chi) cnt[chi] -= (int)((cm >> chi) & 1);
    }
}
static void bb_search(int r, int w, int symDepth, int nthreads, int gencap, BBResult& R) {
    double t0 = elapsed();
    BBTypes T; build_types(r, T);
    build_pmin();
    BBCtx C; C.T = &T; C.r = r; C.w = w; C.nchi = 1 << r;
    C.symDepth = symDepth; C.gencap = gencap; C.R = &R;
    C.deltaMax = (r >= 2) ? 3 * (1 << (r - 2)) : 1;
    C.keepTies = BB_KEEP_TIES;
    if (symDepth > 0 && r <= 6) build_typeperms(T, r, C.perms);
    else C.symDepth = 0;
    R.bestSum = INT64_MAX;
    R.seedIncumbent = bb_seed_incumbent(r, w, 0xC0FFEEull + r * 131 + w, 0.35);
    // with tie-pruning the test is "b >= inc", so the seed must be one unit loose or the leaf
    // that achieves the seed value would itself be pruned and the run could return nothing
    BB_INC.store(R.seedIncumbent == INT64_MAX ? INT64_MAX
                 : R.seedIncumbent + (BB_KEEP_TIES ? 0 : 1), std::memory_order_relaxed);
    std::atomic<int> next{0};
    std::vector<BBResult> loc(std::max(1, nthreads));
    auto worker = [&](int id) {
        std::vector<int> seq(NQ, 0), cnt(C.nchi, 0);
        for (;;) {
            int t = next.fetch_add(1);
            if (t >= T.ntypes) break;
            uint64_t cm = T.contrib[t];
            bool bad = false;
            for (int i = 0; i < r && !bad; ++i) {
                int c = (int)((cm >> (1 << i)) & 1);
                if (c > w || w - c > NQ - 1) bad = true;
            }
            if (bad) continue;
            std::fill(cnt.begin(), cnt.end(), 0);
            for (int chi = 0; chi < C.nchi; ++chi) cnt[chi] += (int)((cm >> chi) & 1);
            seq[0] = t;
            if (!sym_ok(C, seq.data(), 1)) continue;
            ++loc[id].nodes;
            bb_dfs(C, 1, t, seq.data(), cnt.data(), T.plu[t], loc[id]);
        }
    };
    int nt = std::max(1, nthreads);
    std::vector<std::thread> th;
    for (int i = 0; i < nt; ++i) th.emplace_back(worker, i);
    for (auto& t : th) t.join();
    for (auto& l : loc) { R.nodes += l.nodes; R.pruned += l.pruned; R.leaves += l.leaves; }
    R.completed = !BB_ABORTED.load();
    if (R.bestSum != INT64_MAX) { R.U = R.bestSum >> r; R.M = TOT_ERR - R.U; }
    R.seconds = elapsed() - t0;
}

// ===========================================================================================
//  A CERTIFIED LOWER BOUND ON |U_8|   (the only rigorous Stage-1 statement available)
//
//  L is an additive self-orthogonal code over the 4-element alphabet {I,X,Y,Z} of length 14
//  with 256 codewords.  Let A_m be its weight distribution and B_w that of its dual L^perp
//  (2^20 words).  Then, with K_w(m) = SUM_k C(m,k)(-1)^k C(14-m,w-k) 3^{w-k},
//
//        B_w = (1/256) SUM_m A_m K_w(m) ,      B_w >= A_w   (because L <= L^perp) ,
//        A_m >= 0 ,   A_0 = 1 ,   SUM_m A_m = 256 ,   A_10 >= 8   (eight weight-10 generators)
//
//  and |U_8| = SUM_{w=1..4} B_w = (1/256) SUM_m A_m P(m).  For ANY multipliers y_w, z_w, s >= 0
//  put   d_m = P(m)/256 - (1/256) SUM_w (y_w + z_w) K_w(m) + z_m - s*[m=10] .  Then
//
//        |U_8|  >=  8 s  +  d_0  +  255 * min_{m>=1} d_m ,
//
//  because SUM_{m>=1} A_m = 255 with A_m >= 0.  EVERY choice of multipliers gives a valid
//  bound -- weak duality -- so the search below cannot produce a wrong answer, only a weak
//  one.  No LP solver is trusted anywhere in this argument.
// ===========================================================================================
struct BoundResult { double value = 0; int64_t intBound = 0; double y[NQ+1]{}, z[NQ+1]{}, s = 0; };
static double bound_eval(const double* y, const double* z, double s) {
    double d[NQ + 1];
    for (int m = 0; m <= NQ; ++m) {
        double acc = (double)PTAB[m];
        for (int w = 0; w <= NQ; ++w) acc -= (y[w] + z[w]) * (double)EWT[m][w];
        d[m] = acc / 256.0 + z[m] - (m == 10 ? s : 0.0);
    }
    double mn = d[1];
    for (int m = 2; m <= NQ; ++m) mn = std::min(mn, d[m]);
    return 8.0 * s + d[0] + 255.0 * mn;
}
static void certified_bound(BoundResult& B, int passes, uint64_t seed) {
    std::mt19937_64 rng(seed);
    double by[NQ+1] = {0}, bz[NQ+1] = {0}, bs = 0;
    double best = bound_eval(by, bz, bs);
    for (int restart = 0; restart < 8; ++restart) {
        double y[NQ+1] = {0}, z[NQ+1] = {0}, s = 0;
        if (restart) {
            for (int w = 0; w <= NQ; ++w) {
                y[w] = (rng() % 1000) / 1000.0 * std::pow(0.25, w);
                z[w] = (rng() % 1000) / 1000.0 * std::pow(0.25, w);
            }
            s = (rng() % 1000) / 10.0;
        }
        double cur = bound_eval(y, z, s);
        for (double step = 64.0; step > 1e-9; step *= 0.5) {
            for (int pass = 0; pass < passes; ++pass) {
                bool moved = false;
                for (int idx = 0; idx <= 2 * NQ + 2; ++idx) {
                    double* v; double scale = 1.0;
                    if (idx <= NQ) { v = &y[idx]; scale = std::pow(0.25, idx); }
                    else if (idx <= 2 * NQ + 1) { v = &z[idx - NQ - 1]; scale = std::pow(0.25, idx - NQ - 1); }
                    else v = &s;
                    for (int sgn = 0; sgn < 2; ++sgn) {
                        double delta = (sgn ? -step : step) * scale;
                        double old = *v, nv = old + delta;
                        if (nv < 0) nv = 0;
                        if (nv == old) continue;
                        *v = nv;
                        double val = bound_eval(y, z, s);
                        if (val > cur + 1e-12) { cur = val; moved = true; break; }
                        *v = old;
                    }
                }
                if (!moved) break;
            }
        }
        if (cur > best) { best = cur;
            for (int w = 0; w <= NQ; ++w) { by[w] = y[w]; bz[w] = z[w]; } bs = s; }
    }
    // final independent re-evaluation of the certificate (never trust the search loop)
    B.value = bound_eval(by, bz, bs);
    for (int w = 0; w <= NQ; ++w) { B.y[w] = by[w]; B.z[w] = bz[w]; }
    B.s = bs;
    B.intBound = (int64_t)std::ceil(B.value - 1e-6);
    if (B.intBound < 0) B.intBound = 0;
}

// ===========================================================================================
//  THE OPTIMAL STAGE-1 FAMILY
// ===========================================================================================
static std::string subspace_key(const uint32_t* g) {          // RREF -> canonical per subspace
    uint32_t b[R1]; int n = 0;
    for (int i = 0; i < R1; ++i) {
        uint32_t x = g[i];
        for (int k = 0; k < n; ++k) { int hb = 31 - std::countl_zero(b[k]);
            if ((x >> hb) & 1) x ^= b[k]; }
        if (x) { b[n++] = x;
            for (int k = 0; k < n - 1; ++k) { int hb = 31 - std::countl_zero(x);
                if ((b[k] >> hb) & 1) b[k] ^= x; }
            std::sort(b, b + n, std::greater<uint32_t>()); }
    }
    std::string s; s.reserve(n * 4);
    for (int i = 0; i < n; ++i) { char q[4]; memcpy(q, &b[i], 4); s.append(q, 4); }
    return s;
}
struct FamilyMember { S1 s; bool twoOpt = false, threeOpt = false; };
struct FamilyClass  { uint64_t raw = 0; std::vector<FamilyMember> members; };

static std::map<std::string, FamilyClass> FAMILY;
static std::set<std::string> SEEN_SUBSPACES;
static int64_t  BEST_SUM = INT64_MAX;
static int64_t  RAW_DISTINCT = 0;
static S1       BEST;
static bool     HAVE_BEST = false;
static std::mutex FMTX;
static int  MEMCAP = 4, CLASSCAP = 512;
static size_t SEENCAP = 2000000;
static std::atomic<uint64_t> RESTARTS{0}, ONEOPT_NODES{0}, TWOOPT_NODES{0}, TWOOPT_PRUNED{0};
static std::atomic<uint64_t> THREEOPT_NODES{0};
static std::atomic<bool> THREEOPT_INCOMPLETE{false};
static std::atomic<uint64_t> IMPROVEMENTS{0};
static std::string OUTDIR = ".";
static double LAST_REPORT = 0, REPORT_EVERY = 1800.0;

// returns 1 if this is a new global best, 2 if it joined the current best family, 0 otherwise
static int offer_s1(const S1& Sin) {
    S1 S = Sin;
    if (!s1_valid(S, nullptr)) return 0;
    fill_s1(S);
    std::lock_guard<std::mutex> lk(FMTX);
    if (S.sumP > BEST_SUM) return 0;
    std::string key = subspace_key(S.g);
    if (S.sumP < BEST_SUM) {
        BEST_SUM = S.sumP; FAMILY.clear(); SEEN_SUBSPACES.clear(); RAW_DISTINCT = 0;
        BEST = S; HAVE_BEST = true; ++IMPROVEMENTS;
        SEEN_SUBSPACES.insert(key); RAW_DISTINCT = 1;
        FamilyClass fc; fc.raw = 1; fc.members.push_back({S, false});
        FAMILY[fingerprint(S.g)] = fc;
        return 1;
    }
    if (SEEN_SUBSPACES.count(key)) return 2;
    if (SEEN_SUBSPACES.size() < SEENCAP) SEEN_SUBSPACES.insert(key);
    ++RAW_DISTINCT;
    std::string fp = fingerprint(S.g);
    auto it = FAMILY.find(fp);
    if (it == FAMILY.end()) {
        if ((int)FAMILY.size() >= CLASSCAP) return 2;
        FamilyClass fc; fc.raw = 1; fc.members.push_back({S, false});
        FAMILY[fp] = fc;
    } else {
        ++it->second.raw;
        if ((int)it->second.members.size() < MEMCAP) it->second.members.push_back({S, false});
    }
    return 2;
}

// ---------------------------------------------------------------- construction with a prefix
static bool random_extend(S1& S, int fixed, std::mt19937_64& rng, int ncand) {
    std::vector<uint32_t> pb;
    for (int j = 0; j < NQ; ++j) { pb.push_back(pk(1u << j, 0)); pb.push_back(pk(0, 1u << j)); }
    uint32_t el[S1SZ]; el[0] = 0; int nel = 1;
    for (int m = 0; m < fixed; ++m) {
        uint32_t x = S.g[m];
        for (int i = 0; i < nel; ++i) el[nel + i] = el[i] ^ x;
        nel <<= 1;
        int piv = -1;
        for (size_t i = 0; i < pb.size(); ++i) if (symp(x, pb[i])) { piv = (int)i; break; }
        if (piv < 0) return false;
        uint32_t pv = pb[piv]; pb.erase(pb.begin() + piv);
        for (auto& q : pb) if (symp(x, q)) q ^= pv;
    }
    for (int m = fixed; m < R1; ++m) {
        uint32_t bestx = NONE; int64_t bestv = 0;
        int tries = 0, found = 0;
        while (found < ncand && tries < ncand * 300) {
            ++tries;
            uint64_t mask = rng();
            uint32_t x = 0;
            for (size_t i = 0; i < pb.size(); ++i) if ((mask >> (i & 63)) & 1) x ^= pb[i];
            if (pwt(x) != W1) continue;
            bool dup = false;
            for (int i = 0; i < nel; ++i) if (el[i] == x) { dup = true; break; }
            if (dup) continue;
            ++found;
            int64_t v = 0;
            for (int i = 0; i < nel; ++i) v += PTAB[pwt(x ^ el[i])];
            if (bestx == NONE || v < bestv) { bestv = v; bestx = x; }
        }
        if (bestx == NONE) return false;
        S.g[m] = bestx;
        for (int i = 0; i < nel; ++i) el[nel + i] = el[i] ^ bestx;
        nel <<= 1;
        int piv = -1;
        for (size_t i = 0; i < pb.size(); ++i) if (symp(bestx, pb[i])) { piv = (int)i; break; }
        if (piv < 0) return false;
        uint32_t pv = pb[piv]; pb.erase(pb.begin() + piv);
        for (auto& q : pb) if (symp(bestx, q)) q ^= pv;
    }
    fill_s1(S);
    return s1_valid(S, nullptr);
}
static void polish1(S1& S) {
    uint64_t nodes = 0;
    while (one_opt(S, &nodes)) {}
    ONEOPT_NODES += nodes;
}

// ---------------------------------------------------------------- progress reporting
static void status_report(bool force, const char* where, int64_t totalRestarts, double tlimit) {
    if (!force && elapsed() - LAST_REPORT < REPORT_EVERY) return;
    LAST_REPORT = elapsed();
    std::lock_guard<std::mutex> lk(FMTX);
    double el = elapsed();
    uint64_t done = RESTARTS.load();
    double rate = done / std::max(1e-9, el);
    logline("\n---------------- STATUS %s ----------------\n", now_stamp().c_str());
    logline("phase                  : %s\n", where);
    logline("elapsed                : %s (%.1f s)\n", hms(el).c_str(), el);
    if (totalRestarts > 0)
        logline("restarts               : %llu / %lld = %.4f %%\n",
                (unsigned long long)done, (long long)totalRestarts,
                100.0 * double(done) / double(totalRestarts));
    else
        logline("restarts               : %llu (no restart budget; time limited)\n",
                (unsigned long long)done);
    logline("1-opt classes swept    : %llu\n", (unsigned long long)ONEOPT_NODES.load());
    logline("2-opt pairs scored     : %llu\n", (unsigned long long)TWOOPT_NODES.load());
    logline("3-opt triples scored   : %llu\n", (unsigned long long)THREEOPT_NODES.load());
    logline("k-opt prefixes pruned  : %llu\n", (unsigned long long)TWOOPT_PRUNED.load());
    if (HAVE_BEST) {
        logline("best |U_8|             : %lld   (M_8 = %lld / %lld)\n",
                (long long)(BEST_SUM / S1SZ), (long long)(TOT_ERR - BEST_SUM / S1SZ),
                (long long)TOT_ERR);
        logline("improvements so far    : %llu\n", (unsigned long long)IMPROVEMENTS.load());
        logline("distinct optimal spans : %lld\n", (long long)RAW_DISTINCT);
        logline("fingerprint classes    : %zu\n", FAMILY.size());
    } else logline("best |U_8|             : (none yet)\n");
    logline("throughput             : %.2f restarts/s\n", rate);
    if (totalRestarts > 0 && rate > 0) {
        double remain = (double(totalRestarts) - double(done)) / rate;
        logline("estimated total run    : %s\n", hms(el + remain).c_str());
        logline("estimated remaining    : %s\n", hms(remain).c_str());
    } else if (tlimit > 0) {
        logline("estimated remaining    : %s (time limit)\n", hms(tlimit - el).c_str());
    }
    logline("checkpoint             : %s/checkpoint.txt\n", OUTDIR.c_str());
    logline("-------------------------------------------------------------\n\n");
}

// ---------------------------------------------------------------- checkpointing
static void save_checkpoint(int64_t restartsDone, const char* phase) {
    std::lock_guard<std::mutex> lk(FMTX);
    RB.clear();
    ap("# reverse_8w10_3w8 checkpoint -- resume with --resume\n");
    ap("version 1\ntimestamp %s\nelapsed %.2f\nphase %s\n",
       now_stamp().c_str(), elapsed(), phase);
    ap("restarts %lld\n", (long long)restartsDone);
    ap("best_sum %lld\n", (long long)(HAVE_BEST ? BEST_SUM : -1));
    ap("raw_distinct %lld\n", (long long)RAW_DISTINCT);
    ap("oneopt_nodes %llu\ntwoopt_nodes %llu\ntwoopt_pruned %llu\nimprovements %llu\n",
       (unsigned long long)ONEOPT_NODES.load(), (unsigned long long)TWOOPT_NODES.load(),
       (unsigned long long)TWOOPT_PRUNED.load(), (unsigned long long)IMPROVEMENTS.load());
    ap("classes %zu\n", FAMILY.size());
    for (auto& kv : FAMILY) {
        ap("class %llu %zu\n", (unsigned long long)kv.second.raw, kv.second.members.size());
        for (auto& m : kv.second.members) {
            ap("member %d", m.twoOpt ? 1 : 0);
            for (int i = 0; i < R1; ++i) ap(" %s", pstr(m.s.g[i]).c_str());
            ap("\n");
        }
    }
    dump(OUTDIR + "/checkpoint.txt");
}
static bool load_checkpoint(int64_t& restartsDone) {
    FILE* f = fopen((OUTDIR + "/checkpoint.txt").c_str(), "r");
    if (!f) return false;
    char line[4096];
    int64_t bs = -1;
    std::vector<S1> loaded; std::vector<int> loadedTwo;
    while (fgets(line, sizeof(line), f)) {
        if (!strncmp(line, "restarts ", 9)) restartsDone = atoll(line + 9);
        else if (!strncmp(line, "best_sum ", 9)) bs = atoll(line + 9);
        else if (!strncmp(line, "oneopt_nodes ", 13)) ONEOPT_NODES = strtoull(line + 13, nullptr, 10);
        else if (!strncmp(line, "twoopt_nodes ", 13)) TWOOPT_NODES = strtoull(line + 13, nullptr, 10);
        else if (!strncmp(line, "twoopt_pruned ", 14)) TWOOPT_PRUNED = strtoull(line + 14, nullptr, 10);
        else if (!strncmp(line, "improvements ", 13)) IMPROVEMENTS = strtoull(line + 13, nullptr, 10);
        else if (!strncmp(line, "member ", 7)) {
            char* p = line + 7; int two = (int)strtol(p, &p, 10);
            S1 S; bool ok = true;
            for (int i = 0; i < R1; ++i) {
                while (*p == ' ') ++p;
                if ((int)strlen(p) < NQ) { ok = false; break; }
                char buf[NQ + 1]; memcpy(buf, p, NQ); buf[NQ] = 0;
                S.g[i] = parse_pauli(buf); p += NQ;
            }
            if (ok) { fill_s1(S); loaded.push_back(S); loadedTwo.push_back(two); }
        }
    }
    fclose(f);
    if (bs < 0) return false;
    BEST_SUM = INT64_MAX;
    for (size_t i = 0; i < loaded.size(); ++i) {
        int r = offer_s1(loaded[i]);
        (void)r;
        if (loadedTwo[i]) {
            std::lock_guard<std::mutex> lk(FMTX);
            auto it = FAMILY.find(fingerprint(loaded[i].g));
            if (it != FAMILY.end()) for (auto& m : it->second.members)
                if (subspace_key(m.s.g) == subspace_key(loaded[i].g)) m.twoOpt = true;
        }
    }
    return HAVE_BEST;
}

// ===========================================================================================
//  REPORTING
// ===========================================================================================
static void build_U8(const S1& S, std::vector<uint32_t>& U8) {
    U8.clear();
    for_each_error([&](uint32_t E, int w) {
        (void)w;
        for (int i = 0; i < R1; ++i) if (symp(S.g[i], E)) return;
        U8.push_back(E);
    });
}
static std::string log10_str(long double log10v) {
    char b[64];
    long double mant = powl(10.0L, log10v - floorl(log10v));
    // print as double: MSVC's CRT does not handle the %Lf length modifier reliably
    snprintf(b, sizeof(b), "%.3fe+%02d", (double)mant, (int)floorl(log10v));
    return std::string(b);
}
static void space_sizes(long double& logIso, long double& logIsoMod, long double& logMultiset,
                        int64_t& ntypes8) {
    logIso = 0;
    for (int i = 0; i < R1; ++i) logIso += log10l(powl(2.0L, 28 - i) - powl(2.0L, i));
    long double logGL = 0;
    for (int i = 0; i < R1; ++i) logGL += log10l(powl(2.0L, R1) - powl(2.0L, i));
    logIso -= logGL;
    long double logGrp = 0;                       // 14! * 6^14
    for (int i = 2; i <= NQ; ++i) logGrp += log10l((long double)i);
    logGrp += NQ * log10l(6.0L);
    logIsoMod = logIso - logGrp;
    ntypes8 = 1 + ((int64_t)1 << R1) - 1 + gauss_binom(R1, 2);
    long double lm = 0;                           // C(ntypes+13, 14)
    for (int i = 0; i < NQ; ++i) lm += log10l((long double)(ntypes8 + 13 - i));
    for (int i = 2; i <= NQ; ++i) lm -= log10l((long double)i);
    logMultiset = lm;
}
static void residual_block(const ResidualInfo& r) {
    ap("Remaining undetected errors: %lld\n\n", (long long)r.n);
    for (int w = 1; w <= MAXW; ++w) {
        ap("Weight %d:\n", w); int n = 0;
        for (uint32_t E : r.list) if (pwt(E) == w) {
            if (n < 4000) ap("  %s = %s\n", pstr(E).c_str(), supp(E).c_str());
            ++n; }
        if (!n) ap("  (none)\n");
        else if (n > 4000) ap("  ... and %d more of this weight\n", n - 4000);
        ap("\n");
    }
    ap("Residual weight distribution: %lld / %lld / %lld / %lld  (weights 1/2/3/4)\n",
       (long long)r.by_w[1], (long long)r.by_w[2], (long long)r.by_w[3], (long long)r.by_w[4]);
    ap("Residual set forms a subgroup (with the identity): %s\n",
       !r.subgroupTested ? "not tested (span dimension too large)" : (r.subgroup ? "YES" : "NO"));
    ap("Residual F_2 span dimension: %d\n", r.dim);
    ap("Generators of the span:\n");
    for (uint32_t g : r.gens) ap("  %s = %s\n", pstr(g).c_str(), supp(g).c_str());
    ap("All residual errors are X-only: %s\n", r.all_x ? "YES" : "NO");
    ap("Physical qubits involved (%zu):", r.qubits.size());
    for (uint32_t q : r.qubits) ap(" %u", q);
    ap("\n");
}
static void write_columns_report(const std::string& path, const S1& S) {
    ColumnView C; column_view(S.g, C);
    RB.clear();
    ap("COLUMN-TYPE VIEW OF THE STAGE-1 CODE  (the qubit-structured formulation)\n\n");
    ap("qubit   x_j (F_2^8)   z_j (F_2^8)   dim W_j\n");
    for (int j = 0; j < NQ; ++j) {
        std::string xs, zs;
        for (int i = 0; i < R1; ++i) { xs += char('0' + ((C.x[j] >> i) & 1));
                                       zs += char('0' + ((C.z[j] >> i) & 1)); }
        ap("  %2d    %s      %s        %d\n", j, xs.c_str(), zs.c_str(), C.dimW[j]);
    }
    int dimcnt[3] = {0,0,0};
    for (int j = 0; j < NQ; ++j) dimcnt[C.dimW[j]]++;
    ap("\ncolumn dimension multiset: dim0 x %d, dim1 x %d, dim2 x %d\n",
       dimcnt[0], dimcnt[1], dimcnt[2]);
    ap("\nCONSTRAINTS EXPRESSED IN THE COLUMN REPRESENTATION\n");
    ap("  weight  : #{ j : W_j <= H_i } must be %d for every generator i\n", NQ - W1);
    ap("            values:");
    for (int i = 0; i < R1; ++i) ap(" %d", C.zerosPerGen[i]);
    bool wok = true; for (int i = 0; i < R1; ++i) if (C.zerosPerGen[i] != NQ - W1) wok = false;
    ap("   -> %s\n", wok ? "OK" : "VIOLATED");
    ap("  commute : XOR of the 28-bit Plucker invariants must vanish -> %s (value %llu)\n",
       C.pluckerXor == 0 ? "OK" : "VIOLATED", (unsigned long long)C.pluckerXor);
    ap("  rank    : W_1 + ... + W_14 must be all of F_2^8 -> dim %d %s\n",
       C.spanRank, C.spanRank == R1 ? "OK" : "VIOLATED");
    ap("\nWEIGHT FUNCTION RECOVERED FROM THE COLUMN MULTISET ALONE\n");
    uint32_t el[S1SZ]; int64_t sp; span_stats(S.g, R1, sp, nullptr, el);
    bool allok = true;
    for (int a = 0; a < S1SZ; ++a) {
        uint32_t h = 0; // element with coefficient vector a
        for (int i = 0; i < R1; ++i) if ((a >> i) & 1) h ^= S.g[i];
        if (weight_from_columns(C, (uint32_t)a) != pwt(h)) allok = false;
    }
    ap("  wt(h_a) = 14 - #{ j : a _|_ W_j } for all 256 a  ->  %s\n", allok ? "OK" : "MISMATCH");
    int64_t sumFromCols = 0;
    for (int a = 0; a < S1SZ; ++a) sumFromCols += PTAB[weight_from_columns(C, (uint32_t)a)];
    ap("  objective from columns  : %lld  (span sum %lld) -> %s\n",
       (long long)sumFromCols, (long long)sp, sumFromCols == sp ? "OK" : "MISMATCH");
    long double li, lim, lms; int64_t nt8;
    space_sizes(li, lim, lms, nt8);
    ap("\nSEARCH-SPACE ACCOUNTING FOR THIS REPRESENTATION\n");
    ap("  raw column types (2^16)                     : %d\n", 1 << (2 * R1));
    ap("  reduced column types (subspaces of F_2^8)   : %lld\n", (long long)nt8);
    ap("  multisets of 14 reduced types               : %s\n", log10_str(lms).c_str());
    dump(path);
}
static void write_class_dir(const std::string& dir, int idx, const S1& S,
                            const std::vector<uint32_t>& U8, const S2Result& R2r,
                            const Verify& V) {
    std::error_code ec; std::filesystem::create_directories(dir, ec);
    RB.clear();
    ap("STAGE-1 REPRESENTATIVE  (class %04d)\n\n", idx);
    for (int i = 0; i < R1; ++i) ap("g%d = %s   (weight %d)\n", i + 1, pstr(S.g[i]).c_str(), pwt(S.g[i]));
    ap("\nM_8 = %lld / %lld     |U_8| = %lld\n",
       (long long)(TOT_ERR - S.U8), (long long)TOT_ERR, (long long)S.U8);
    ap("\nWeight enumerator of the 256-element Stage-1 group:\n");
    for (int m = 0; m <= NQ; ++m) if (S.A[m]) ap("  A[%2d] = %d\n", m, S.A[m]);
    ap("\nStage-1 undetected set by weight: ");
    { int64_t bw[MAXW+1] = {0,0,0,0,0};
      for (uint32_t E : U8) bw[pwt(E)]++;
      for (int w = 1; w <= MAXW; ++w) ap("w%d=%lld ", w, (long long)bw[w]); }
    ap("\n");
    dump(dir + "/stage1_representative.txt");
    RB.clear();
    ap("STAGE-1 UNDETECTED SET U_8   (%zu errors)\n\n", U8.size());
    for (uint32_t E : U8) ap("%s  %s\n", pstr(E).c_str(), supp(E).c_str());
    dump(dir + "/stage1_undetected.txt");
    RB.clear();
    if (R2r.found) {
        ap("STAGE-2 OPTIMUM FOR THIS FROZEN STAGE-1 CODE\n");
        ap("(exhaustive over every 3-dimensional isotropic subspace of Q2 spanned by\n");
        ap(" weight-8-usable classes -- PROVEN OPTIMUM conditional on the frozen Stage 1)\n\n");
        for (int a = 0; a < R2; ++a)
            ap("h%d = %s   (weight %d)\n", a + 1, pstr(R2r.h[a]).c_str(), pwt(R2r.h[a]));
        ap("\n|U_8|                       = %zu\n", U8.size());
        ap("detected by Stage 2         = %lld\n", (long long)R2r.cover);
        ap("remaining                   = %lld\n", (long long)R2r.residual);
        ap("M_final                     = %lld / %lld  (%.8f %%)\n",
           (long long)V.det11, (long long)V.all_total,
           100.0 * double(V.det11) / double(V.all_total));
        ap("\nquotient Q2 dim             = %d  (|Q2| = %d)\n", R2r.qdim, R2r.qsize);
        ap("weight-8-usable classes     = %d\n", R2r.usable);
        ap("max cov1 over classes       = %d\n", R2r.covmax);
        ap("triples scored              = %llu\n", (unsigned long long)R2r.nodes);
        ap("distinct optimal subspaces  = %lld\n", (long long)R2r.nOptimalSubspaces);
    } else ap("STAGE 2: no valid triple exists for this frozen Stage-1 code.\n");
    dump(dir + "/stage2_best.txt");
    RB.clear();
    ap("COMPLETE 11-GENERATOR STABILIZER\n\n");
    for (int i = 0; i < R1; ++i) ap("g%d  = %s   (weight %d)\n", i + 1, pstr(S.g[i]).c_str(), pwt(S.g[i]));
    if (R2r.found) for (int a = 0; a < R2; ++a)
        ap("h%d  = %s   (weight %d)\n", a + 1, pstr(R2r.h[a]).c_str(), pwt(R2r.h[a]));
    ap("\nH = [ X | Z ]   (%d x 28):\n", R2r.found ? RTOT : R1);
    for (int i = 0; i < R1; ++i) ap("  %s | %s\n", bstr(px(S.g[i])).c_str(), bstr(pz(S.g[i])).c_str());
    if (R2r.found) for (int a = 0; a < R2; ++a)
        ap("  %s | %s\n", bstr(px(R2r.h[a])).c_str(), bstr(pz(R2r.h[a])).c_str());
    ap("\nrank = %d\n", V.rank);
    dump(dir + "/final_matrix.txt");
    RB.clear(); residual_block(analyse(V.leftover)); dump(dir + "/remaining_errors.txt");
    RB.clear();
    ap("INDEPENDENT VERIFICATION  (all %lld weight-<=4 errors regenerated from scratch and\n"
       " their syndromes recomputed by direct symplectic products, not by the algebraic score)\n\n",
       (long long)V.all_total);
    ap("  1. eight Stage-1 generators of weight exactly %d : %s  (", W1, V.w1_ok ? "OK" : "FAIL");
    for (int i = 0; i < R1; ++i) ap("%d ", V.weights1[i]); ap(")\n");
    if (R2r.found) {
        ap("  2. three Stage-2 generators of weight exactly %d : %s  (", W2, V.w2_ok ? "OK" : "FAIL");
        for (int a = 0; a < R2; ++a) ap("%d ", V.weights2[a]); ap(")\n");
    } else ap("  2. Stage-2 generators                            : none\n");
    ap("  3. all generator pairs commute                   : %s\n", V.cg_ok && V.ch_ok ? "OK" : "FAIL");
    ap("  4. rank of the stabilizer                        : %d %s\n", V.rank, V.rank_ok ? "OK" : "FAIL");
    ap("  5. |U_8| rebuilt by brute force                  : %lld %s\n",
       (long long)V.u8size, V.u8_ok ? "OK (matches the algebraic score)" : "*** MISMATCH ***");
    ap("  6. closed-form identity  sum P(wt)/256 == |U_8|  : %s\n", V.ptab_ok ? "OK" : "*** MISMATCH ***");
    ap("  7. Stage-2 coverage recomputed independently     : %lld %s\n",
       (long long)V.cover, (R2r.found && V.cover == R2r.cover) ? "OK" : (R2r.found ? "*** MISMATCH ***" : "n/a"));
    ap("  8. final residual                                : %lld\n", (long long)V.residual);
    ap("\nWeight-4 errors total     : %lld\n", (long long)V.w4_total);
    ap("Weight-4 errors detected  : %lld\n", (long long)V.w4_det);
    ap("Weight-4 errors undetected: %lld\n", (long long)V.w4_und);
    if (V.w4_und == 0) ap("CONDITION B SATISFIED: ALL WEIGHT-4 ERRORS DETECTED\n");
    ap("\nDetected by weight: ");
    for (int w = 1; w <= MAXW; ++w) ap("w%d=%lld/%lld  ", w, (long long)V.det_by_w[w], (long long)TOT_BY_W[w]);
    ap("\n");
    if (V.residual == 0) {
        ap("\nCONDITION A SATISFIED: every weight-<=4 Pauli error has a non-zero 11-bit syndrome,\n"
           "so N(L) contains no non-zero Pauli of weight <= 4 and the code has PURE DISTANCE >= 5:\n"
           "a [[14,3,5]] construction.\n");
    }
    dump(dir + "/verification.txt");
    write_columns_report(dir + "/columns.txt", S);
}

// ===========================================================================================
//  SELF-TESTS   (every algebraic shortcut is checked against direct enumeration)
// ===========================================================================================
static bool random_isotropic(uint32_t* g, int r, std::mt19937_64& rng, int wantWeight) {
    std::vector<uint32_t> pb;
    for (int j = 0; j < NQ; ++j) { pb.push_back(pk(1u << j, 0)); pb.push_back(pk(0, 1u << j)); }
    std::vector<uint32_t> el; el.push_back(0);
    for (int m = 0; m < r; ++m) {
        uint32_t x = 0; bool ok = false;
        for (int tries = 0; tries < 200000 && !ok; ++tries) {
            uint64_t mask = rng(); x = 0;
            for (size_t i = 0; i < pb.size(); ++i) if ((mask >> (i & 63)) & 1) x ^= pb[i];
            if (!x) continue;
            if (wantWeight > 0 && pwt(x) != wantWeight) continue;
            if (std::find(el.begin(), el.end(), x) != el.end()) continue;
            ok = true;
        }
        if (!ok) return false;
        g[m] = x;
        size_t n = el.size(); for (size_t i = 0; i < n; ++i) el.push_back(el[i] ^ x);
        int piv = -1;
        for (size_t i = 0; i < pb.size(); ++i) if (symp(x, pb[i])) { piv = (int)i; break; }
        if (piv < 0) return false;
        uint32_t pv = pb[piv]; pb.erase(pb.begin() + piv);
        for (auto& q : pb) if (symp(x, q)) q ^= pv;
    }
    return true;
}
static int64_t brute_undetected(const uint32_t* g, int r) {
    int64_t c = 0;
    for_each_error([&](uint32_t E, int w) { (void)w;
        for (int i = 0; i < r; ++i) if (symp(g[i], E)) return;
        ++c; });
    return c;
}
static bool selftest(int trials, int threads, double bbcap) {
    bool pass = true;
    logline("=========================== SELF TEST ===========================\n");
    int64_t tot = 0;
    for (int w = 1; w <= MAXW; ++w) { int64_t p3 = 1; for (int i = 0; i < w; ++i) p3 *= 3;
        tot += BINOM[NQ][w] * p3; }
    logline("  total weight-<=4 Pauli errors  : enumerated %lld, formula %lld  %s\n",
            (long long)TOT_ERR, (long long)tot, TOT_ERR == tot ? "OK" : "FAIL");
    if (TOT_ERR != tot) pass = false;
    logline("  weight-4 errors                : enumerated %lld, C(14,4)*81 = %lld  %s\n",
            (long long)TOT_BY_W[4], (long long)(BINOM[NQ][4] * 81),
            TOT_BY_W[4] == BINOM[NQ][4] * 81 ? "OK" : "FAIL");
    if (TOT_BY_W[4] != BINOM[NQ][4] * 81) pass = false;
    logline("  P[0] equals the total          : %lld  %s\n", (long long)PTAB[0],
            PTAB[0] == TOT_ERR ? "OK" : "FAIL");
    if (PTAB[0] != TOT_ERR) pass = false;
    { std::string s = "  P[m], m=0..14                  :";
      for (int m = 0; m <= NQ; ++m) { char b[32]; snprintf(b, sizeof(b), " %lld", (long long)PTAB[m]); s += b; }
      logline("%s\n", s.c_str()); }

    std::mt19937_64 rng(12345);
    int bad = 0, done = 0;
    for (int t = 0; t < trials; ++t) {
        int r = 1 + (int)(rng() % R1);
        int ww = (rng() % 2) ? W1 : 0;
        uint32_t g[R1];
        if (!random_isotropic(g, r, rng, ww)) continue;
        ++done;
        int64_t sp; span_stats(g, r, sp, nullptr, nullptr);
        if (sp % (1 << r) != 0) { ++bad; continue; }
        if (sp / (1 << r) != brute_undetected(g, r)) ++bad;
    }
    logline("  closed form vs brute force     : %d random subspaces, %d mismatches  %s\n",
            done, bad, bad == 0 ? "OK" : "FAIL");
    if (bad) pass = false;

    // The duality certificate needs EW[m][w] for the FULL range w = 0..14, not just the
    // objective window.  Check the whole table: the dual distribution must be non-negative
    // integers summing to 2^(28-r), and its w=1..4 part must equal the brute-force count.
    {
        int bad2 = 0, tested2 = 0;
        for (int t = 0; t < 25; ++t) {
            int r = 1 + (int)(rng() % 6);
            uint32_t g[R1];
            if (!random_isotropic(g, r, rng, 0)) continue;
            int32_t A[NQ + 1]; int64_t sp;
            span_stats(g, r, sp, A, nullptr);
            // exact integer arithmetic throughout: the largest partial sum is about 6e10
            bool ok = true; int64_t tot = 0, low = 0;
            const int64_t den = (int64_t)1 << r;
            for (int w = 0; w <= NQ; ++w) {
                int64_t num = 0;
                for (int m = 0; m <= NQ; ++m) num += (int64_t)A[m] * EWT[m][w];
                if (num % den != 0 || num < 0) { ok = false; break; }
                int64_t B = num / den;
                tot += B;
                if (w >= 1 && w <= MAXW) low += B;
            }
            int64_t want = (int64_t)1 << (28 - r);
            int64_t bf = brute_undetected(g, r);
            if (tot != want) ok = false;
            if (low != bf) ok = false;
            ++tested2;
            if (!ok) {
                if (bad2 < 3)
                    logline("      [diag] r=%d tot=%lld want=%lld low=%lld brute=%lld sp=%lld\n",
                            r, (long long)tot, (long long)want, (long long)low,
                            (long long)bf, (long long)sp);
                ++bad2;
            }
        }
        logline("  MacWilliams dual, w = 0..14    : %d subspaces, %d mismatches  %s\n",
                tested2, bad2, bad2 == 0 ? "OK" : "FAIL");
        if (bad2) pass = false;
    }

    // -------- a genuine Stage-1 configuration, then column view and quotient checks
    S1 S; bool got = false;
    for (int t = 0; t < 200 && !got; ++t) got = random_extend(S, 0, rng, 24);
    if (!got) { logline("  could not build a Stage-1 configuration -- FAIL\n"); return false; }
    logline("  built a valid Stage-1 code     : |U_8| = %lld, M_8 = %lld\n",
            (long long)S.U8, (long long)(TOT_ERR - S.U8));
    {
        ColumnView C; column_view(S.g, C);
        bool wok = true, allok = true;
        for (int i = 0; i < R1; ++i) if (C.zerosPerGen[i] != NQ - W1) wok = false;
        int64_t sumc = 0;
        for (int a = 0; a < S1SZ; ++a) {
            uint32_t hh = 0;
            for (int i = 0; i < R1; ++i) if ((a >> i) & 1) hh ^= S.g[i];
            if (weight_from_columns(C, (uint32_t)a) != pwt(hh)) allok = false;
            sumc += PTAB[weight_from_columns(C, (uint32_t)a)];
        }
        logline("  column weight constraint       : %s\n", wok ? "OK" : "FAIL");
        logline("  column Plucker commutation     : %s\n", C.pluckerXor == 0 ? "OK" : "FAIL");
        logline("  column rank constraint         : dim %d %s\n", C.spanRank,
                C.spanRank == R1 ? "OK" : "FAIL");
        logline("  weight function from columns   : %s\n", allok ? "OK" : "FAIL");
        logline("  objective from columns         : %lld vs %lld  %s\n",
                (long long)sumc, (long long)S.sumP, sumc == S.sumP ? "OK" : "FAIL");
        if (!wok || C.pluckerXor || C.spanRank != R1 || !allok || sumc != S.sumP) pass = false;
    }
    {
        int64_t nt = 1 + ((int64_t)1 << 3) - 1 + gauss_binom(3, 2);
        int64_t nf = 1 + ((int64_t)1 << 4) - 1 + gauss_binom(4, 2);
        logline("  reduced column types r=3 / r=4 : %lld / %lld  %s (expect 15 / 51)\n",
                (long long)nt, (long long)nf, (nt == 15 && nf == 51) ? "OK" : "FAIL");
        if (nt != 15 || nf != 51) pass = false;
    }
    {   // quotient consistency for a 7-generator retained set
        uint32_t keep[R1 - 1];
        for (int i = 0; i < R1 - 1; ++i) keep[i] = S.g[i];
        QSpace Q;
        bool ok = build_qspace(keep, R1 - 1, Q, true, 0);
        bool cls_ok = ok, form_ok = ok, phi_ok = ok;
        if (ok) {
            for (int c = 0; c < Q.size; ++c) if (qcls(Q, Q.lift[c]) != (uint32_t)c) { cls_ok = false; break; }
            for (int t = 0; t < 2000; ++t) {
                uint32_t a = (uint32_t)(rng() % Q.size), b = (uint32_t)(rng() % Q.size);
                if (qsymp(Q, a, b) != symp(Q.lift[a], Q.lift[b])) { form_ok = false; break; }
            }
            for (int t = 0; t < 200; ++t) {
                uint32_t c = (uint32_t)(rng() % Q.size);
                int64_t acc = 0; uint32_t cur = Q.lift[c];
                acc += PTAB[pwt(cur)];
                for (int i = 1; i < (1 << (R1 - 1)); ++i) {
                    cur ^= keep[std::countr_zero((unsigned)i)];
                    acc += PTAB[pwt(cur)];
                }
                // cur has walked a Gray cycle; recompute cleanly
                acc = 0;
                for (int i = 0; i < (1 << (R1 - 1)); ++i) {
                    uint32_t u = 0;
                    for (int k = 0; k < R1 - 1; ++k) if ((i >> k) & 1) u ^= keep[k];
                    acc += PTAB[pwt(Q.lift[c] ^ u)];
                }
                if (acc != Q.phi[c]) { phi_ok = false; break; }
            }
        }
        logline("  quotient dim / size            : %d / %d %s\n", Q.dim, Q.size,
                Q.dim == 2 * NQ - 2 * (R1 - 1) ? "OK" : "FAIL");
        logline("  class(lift(c)) == c            : %s\n", cls_ok ? "OK" : "FAIL");
        logline("  descended symplectic form      : %s\n", form_ok ? "OK" : "FAIL");
        logline("  phi table vs direct sum        : %s\n", phi_ok ? "OK" : "FAIL");
        if (!ok || !cls_ok || !form_ok || !phi_ok) pass = false;
    }
    {   // Stage-2 coverage identity
        std::vector<uint32_t> U8; build_U8(S, U8);
        bool sizeok = ((int64_t)U8.size() == S.U8);
        logline("  |U_8| brute force vs algebra   : %zu vs %lld  %s\n", U8.size(),
                (long long)S.U8, sizeok ? "OK" : "FAIL");
        if (!sizeok) pass = false;
        QSpace Q;
        if (build_qspace(S.g, R1, Q, false, W2)) {
            std::vector<int32_t> f(Q.size, 0);
            for (uint32_t e : U8) { uint32_t c = qcls(Q, e);
                f[((c >> Q.h) & Q.lowmask) | ((c & Q.lowmask) << Q.h)] += 1; }
            for (int b = 0; b < Q.dim; ++b) { int step = 1 << b;
                for (int i = 0; i < Q.size; i += step << 1)
                    for (int k = i; k < i + step; ++k) {
                        int32_t a2 = f[k], d2 = f[k + step]; f[k] = a2 + d2; f[k + step] = a2 - d2; } }
            std::vector<int32_t> cov1(Q.size);
            for (int v = 0; v < Q.size; ++v) cov1[v] = ((int32_t)U8.size() - f[v]) / 2;
            int idok = 1, tested = 0;
            for (int t = 0; t < 200; ++t) {
                uint32_t a = (uint32_t)(1 + rng() % (Q.size - 1));
                uint32_t b = (uint32_t)(1 + rng() % (Q.size - 1));
                uint32_t c = (uint32_t)(1 + rng() % (Q.size - 1));
                if (a == b || a == c || b == c || c == (a ^ b)) continue;
                if (qsymp(Q, a, b) || qsymp(Q, a, c) || qsymp(Q, b, c)) continue;
                int64_t alg = (int64_t)cov1[a] + cov1[b] + cov1[c] + cov1[a^b] + cov1[a^c]
                            + cov1[b^c] + cov1[a^b^c];
                if (alg % 4) { idok = 0; break; }
                alg /= 4;
                int64_t bf = 0;
                for (uint32_t e : U8) {
                    if (symp(Q.lift[a], e) || symp(Q.lift[b], e) || symp(Q.lift[c], e)) ++bf; }
                ++tested;
                if (alg != bf) { idok = 0; break; }
            }
            logline("  Stage-2 coverage identity      : %d random 3-spaces  %s\n",
                    tested, idok ? "OK" : "FAIL");
            if (!idok) pass = false;
            int usable = 0;
            for (int c = 1; c < Q.size; ++c) if (Q.repW[c] != NONE) ++usable;
            logline("  weight-8-usable classes in Q2  : %d / %d\n", usable, Q.size - 1);
        } else { logline("  Q2 construction                : FAIL\n"); pass = false; }
    }
    // -------- branch and bound regression against the two optima proven earlier in this project
    {
        const int64_t REF[2] = { 80584, 86340 };   // proven earlier in this project, r=3 and r=4
        for (int r = 3; r <= 4; ++r) {
            BB_ABORTED = false;
            BB_DEADLINE = elapsed() + bbcap;
            BBResult R; bb_search(r, 8, 6, threads, 1, R);
            BB_DEADLINE = 0;
            bool ok = R.completed && R.M == REF[r - 3];
            logline("  B&B r=%d w=8  M = %lld  (previously proven in this project: %lld)  %s\n",
                    r, (long long)R.M, (long long)REF[r - 3],
                    !R.completed ? "TREE CUT SHORT -- inconclusive" : (ok ? "MATCH" : "DIFFERENT"));
            logline("               seed incumbent |U| = %lld, nodes %llu, leaves %llu, "
                    "pruned %llu, %.2f s\n",
                    (long long)(R.seedIncumbent >> r), (unsigned long long)R.nodes,
                    (unsigned long long)R.leaves, (unsigned long long)R.pruned, R.seconds);
            if (!ok) pass = false;
        }
        BB_ABORTED = false;
    }
    logline("=================== SELF TEST %s ===================\n\n", pass ? "PASSED" : "FAILED");
    return pass;
}

// ===========================================================================================
//  DRIVER
// ===========================================================================================
struct Opts {
    int threads = 0, restarts = 0, ncand = 24, symDepth = 6, seedbb = 0;
    int memcap = 4, classcap = 512;
    double tlimit = 0, bbtime = 0, reportEvery = 1800.0;
    bool resume = false, twoopt = true, threeopt = false, forceExact = false, deterministic = false;
    double threeoptTime = 0;
    std::string mode = "hybrid", outdir = ".";
    uint64_t seed = 0x5eed1234ULL;
    int selftest = 0, calibrate = 0, bbr = 0, bbw = 0;
};
static std::vector<std::array<uint32_t,16>> SEEDPOOL;

struct ClassReport {
    int idx = 0; S1 s; int64_t u8 = 0, s2cover = 0, s2res = 0, mfinal = 0;
    int64_t u8byw[MAXW+1] = {0,0,0,0,0}, resbyw[MAXW+1] = {0,0,0,0,0};
    int resDim = 0; bool resSub = false, allW4 = false, s2found = false, verifyOK = false;
    int usable = 0; int64_t nOpt = 0; uint64_t members = 0;
    uint32_t h[R2] = {0,0,0};
};
static std::vector<ClassReport> REPORTS;

static void new_best_block(const S1& S) {
    logline("\n==================================================\n");
    logline("NEW STAGE-1 BEST\n");
    logline("==================================================\n");
    logline("Detected errors: %lld / %lld\n", (long long)(TOT_ERR - S.U8), (long long)TOT_ERR);
    logline("Undetected: %lld\n", (long long)S.U8);
    logline("Detection percentage: %.8f %%\n",
            100.0 * double(TOT_ERR - S.U8) / double(TOT_ERR));
    for (int i = 0; i < R1; ++i) logline("g%d = %s\n", i + 1, pstr(S.g[i]).c_str());
    logline("Restarts: %llu   Elapsed: %.1f s\n",
            (unsigned long long)RESTARTS.load(), elapsed());
    logline("==================================================\n\n");
}

static void phase_multistart(const Opts& O) {
    logline("[phase 2] multi-start construction + exhaustive 1-opt, %d threads\n", O.threads);
    std::atomic<int64_t> issued{0};
    double deadline = O.tlimit > 0 ? elapsed() + O.tlimit : 0;
    auto worker = [&](int id) {
        std::mt19937_64 rng(O.seed + 0x9e3779b97f4a7c15ULL * (uint64_t)(id + 1));
        for (;;) {
            if (O.restarts > 0 && issued.fetch_add(1) >= O.restarts) break;
            if (O.restarts <= 0) issued.fetch_add(1);
            if (deadline > 0 && elapsed() > deadline) break;
            S1 S; bool ok;
            if (!SEEDPOOL.empty() && (rng() % 2)) {
                const auto& sd = SEEDPOOL[rng() % SEEDPOOL.size()];
                int fixed = 0;
                for (int i = 0; i < R1; ++i) if (sd[i]) { S.g[i] = sd[i]; fixed = i + 1; }
                ok = random_extend(S, fixed, rng, O.ncand);
            } else ok = random_extend(S, 0, rng, O.ncand);
            if (!ok) continue;
            polish1(S);
            int r = offer_s1(S);
            RESTARTS.fetch_add(1);
            if (r == 1) new_best_block(S);
            if (id == 0) {
                status_report(false, "multi-start", O.restarts, O.tlimit);
                if (elapsed() - LAST_REPORT < 1.0) save_checkpoint((int64_t)RESTARTS.load(), "multistart");
            }
        }
    };
    std::vector<std::thread> th;
    for (int i = 0; i < O.threads; ++i) th.emplace_back(worker, i);
    for (auto& t : th) t.join();
    status_report(true, "multi-start done", O.restarts, O.tlimit);
    save_checkpoint((int64_t)RESTARTS.load(), "multistart-done");
}

static void phase_twoopt(const Opts& O) {
    logline("[phase 3] exhaustive 2-opt certification of every family member\n");
    for (int round = 0; round < 8; ++round) {
        std::vector<S1> todo;
        {
            std::lock_guard<std::mutex> lk(FMTX);
            for (auto& kv : FAMILY) for (auto& m : kv.second.members)
                if (!m.twoOpt) todo.push_back(m.s);
        }
        if (todo.empty()) break;
        logline("  round %d: %zu members to certify\n", round + 1, todo.size());
        std::atomic<size_t> next{0};
        std::atomic<int> improved{0};
        std::mutex impmtx; std::vector<S1> improvedSols;
        auto worker = [&]() {
            for (;;) {
                size_t i = next.fetch_add(1);
                if (i >= todo.size()) break;
                S1 S = todo[i];
                uint64_t nodes = 0, pruned = 0;
                bool got = two_opt(S, &nodes, &pruned);
                TWOOPT_NODES += nodes; TWOOPT_PRUNED += pruned;
                if (got) {
                    polish1(S);
                    std::lock_guard<std::mutex> lk(impmtx); improvedSols.push_back(S); ++improved;
                } else {
                    std::lock_guard<std::mutex> lk(FMTX);
                    auto it = FAMILY.find(fingerprint(todo[i].g));
                    if (it != FAMILY.end()) for (auto& m : it->second.members)
                        if (subspace_key(m.s.g) == subspace_key(todo[i].g)) m.twoOpt = true;
                }
                status_report(false, "2-opt", 0, 0);
            }
        };
        std::vector<std::thread> th;
        for (int i = 0; i < O.threads; ++i) th.emplace_back(worker);
        for (auto& t : th) t.join();
        for (auto& S : improvedSols) { int r = offer_s1(S); if (r == 1) new_best_block(S); }
        if (improved.load() == 0) break;
        logline("  round %d improved %d members; repeating\n", round + 1, improved.load());
    }
    status_report(true, "2-opt done", 0, 0);
    save_checkpoint((int64_t)RESTARTS.load(), "twoopt-done");
}

// Optional deep certification.  three_opt is internally parallel, so members are taken one at
// a time with every thread on each -- a full 3-opt is roughly an hour and a quarter on 12
// threads, against a fraction of a second for 2-opt.  A member is only marked 3-opt certified
// when its sweep actually COMPLETED; a run cut short by --threeopt-time is reported as such
// and never described as a certification.
static void phase_threeopt(const Opts& O) {
    logline("[phase 3b] exhaustive 3-opt certification (%d threads per member%s)\n",
            O.threads, O.threeoptTime > 0 ? ", time capped" : "");
    for (int round = 0; round < 8; ++round) {
        std::vector<S1> todo;
        { std::lock_guard<std::mutex> lk(FMTX);
          for (auto& kv : FAMILY) for (auto& m : kv.second.members)
              if (!m.threeOpt) todo.push_back(m.s); }
        if (todo.empty()) break;
        int improved = 0;
        for (size_t i = 0; i < todo.size(); ++i) {
            if (O.tlimit > 0 && elapsed() > O.tlimit) {
                logline("  time limit reached; %zu members left uncertified\n", todo.size() - i);
                return;
            }
            S1 S = todo[i];
            uint64_t n3 = 0, p3 = 0, c1d = 0, c1t = 0; bool comp3 = true;
            double t0 = elapsed();
            bool got = three_opt(S, &n3, &p3,
                                 O.threeoptTime > 0 ? elapsed() + O.threeoptTime : 0,
                                 &comp3, O.threads, &c1d, &c1t);
            THREEOPT_NODES += n3;
            if (!comp3) THREEOPT_INCOMPLETE = true;
            logline("  member %zu/%zu: %.1f s, %llu triples, %.3f %% of c1 prefixes, %s%s\n",
                    i + 1, todo.size(), elapsed() - t0, (unsigned long long)n3,
                    c1t ? 100.0 * double(c1d) / double(c1t) : 0.0,
                    comp3 ? "COMPLETE" : "cut short",
                    got ? ", IMPROVED" : "");
            if (got) { polish1(S); int r = offer_s1(S); if (r == 1) new_best_block(S); ++improved; }
            else if (comp3) {
                std::lock_guard<std::mutex> lk(FMTX);
                auto it = FAMILY.find(fingerprint(todo[i].g));
                if (it != FAMILY.end()) for (auto& m : it->second.members)
                    if (subspace_key(m.s.g) == subspace_key(todo[i].g)) m.threeOpt = true;
            } else {
                std::lock_guard<std::mutex> lk(FMTX);      // do not retry an uncertifiable one
                auto it = FAMILY.find(fingerprint(todo[i].g));
                if (it != FAMILY.end()) for (auto& m : it->second.members)
                    if (subspace_key(m.s.g) == subspace_key(todo[i].g)) m.threeOpt = true;
            }
            save_checkpoint((int64_t)RESTARTS.load(), "threeopt");
        }
        if (!improved) break;
        logline("  round %d improved %d members; repeating\n", round + 1, improved);
    }
    status_report(true, "3-opt done", 0, 0);
}

static void phase_stage2(const Opts& O) {
    std::vector<std::pair<std::string, FamilyClass>> classes;
    { std::lock_guard<std::mutex> lk(FMTX);
      for (auto& kv : FAMILY) classes.push_back(kv); }
    logline("[phase 4] exhaustive Stage-2 for %zu fingerprint classes "
            "(up to %d members each)\n", classes.size(), O.memcap);
    int idx = 0;
    for (auto& kv : classes) {
        for (size_t mi = 0; mi < kv.second.members.size(); ++mi) {
            ++idx;
            const S1& S = kv.second.members[mi].s;
            std::vector<uint32_t> U8; build_U8(S, U8);
            S2Result R2r;
            double t0 = elapsed();
            stage2_exhaustive(S, U8, R2r, O.threads);
            Verify V = verify_full(S, R2r.h, R2r.found);
            ClassReport CR;
            CR.idx = idx; CR.s = S; CR.u8 = (int64_t)U8.size(); CR.members = kv.second.raw;
            for (uint32_t E : U8) CR.u8byw[pwt(E)]++;
            CR.s2found = R2r.found;
            if (R2r.found) { CR.s2cover = R2r.cover; CR.s2res = R2r.residual;
                for (int a = 0; a < R2; ++a) CR.h[a] = R2r.h[a]; }
            CR.mfinal = V.det11;
            ResidualInfo ri = analyse(V.leftover);
            for (int w = 1; w <= MAXW; ++w) CR.resbyw[w] = ri.by_w[w];
            CR.resDim = ri.dim; CR.resSub = ri.subgroup;
            CR.allW4 = (V.w4_und == 0);
            CR.usable = R2r.usable; CR.nOpt = R2r.nOptimalSubspaces;
            CR.verifyOK = V.w1_ok && V.cg_ok && V.ch_ok && V.rank_ok && V.u8_ok && V.ptab_ok
                          && (!R2r.found || V.cover == R2r.cover);
            REPORTS.push_back(CR);
            char dn[64]; snprintf(dn, sizeof(dn), "/stage1_orbit_%04d", idx);
            write_class_dir(OUTDIR + dn, idx, S, U8, R2r, V);
            logline("  class %04d: |U_8| = %5lld  Stage-2 detects %5lld  residual %5lld  "
                    "M_final = %lld  (%.1f s, %d usable classes)%s\n",
                    idx, (long long)CR.u8, (long long)CR.s2cover, (long long)CR.s2res,
                    (long long)CR.mfinal, elapsed() - t0, CR.usable,
                    CR.verifyOK ? "" : "  *** VERIFICATION PROBLEM ***");
            if (R2r.found && R2r.residual == 0)
                logline("  *** CONDITION A: FULL COVERAGE -- pure [[14,3,5]] construction ***\n");
            if (CR.allW4)
                logline("  *** CONDITION B: all %lld weight-4 errors detected ***\n",
                        (long long)TOT_BY_W[4]);
            status_report(false, "stage 2", 0, 0);
        }
    }
    save_checkpoint((int64_t)RESTARTS.load(), "stage2-done");
}

// ===========================================================================================
//  FINAL REPORTS
// ===========================================================================================
static void write_reports(const Opts& O, const BoundResult& B, bool timeLimited) {
    const int64_t bestU = HAVE_BEST ? BEST_SUM / S1SZ : -1;
    RB.clear();
    ap("============================================================\n");
    ap("STAGE 1 -- EIGHT WEIGHT-10 GENERATORS\n");
    ap("============================================================\n\n");
    ap("Total errors: %lld\n\n", (long long)TOT_ERR);
    ap("STATUS: BEST KNOWN  (an exhaustive Stage-1 search is not computationally possible;\n");
    ap("        see the search-space accounting below and README_MATH_REVERSE.txt)\n\n");
    ap("Best known:\n");
    ap("M8 = %lld / %lld\n", (long long)(TOT_ERR - bestU), (long long)TOT_ERR);
    ap("Residual = %lld\n\n", (long long)bestU);
    if (HAVE_BEST) {
        std::vector<uint32_t> U8; build_U8(BEST, U8);
        int64_t bw[MAXW+1] = {0,0,0,0,0};
        for (uint32_t E : U8) bw[pwt(E)]++;
        ap("Residual distribution:\n");
        for (int w = 1; w <= MAXW; ++w) ap("weight %d: %lld\n", w, (long long)bw[w]);
        ap("\nBest-known Stage-1 generators:\n");
        for (int i = 0; i < R1; ++i) ap("g%d = %s   (weight %d)\n", i + 1,
                                        pstr(BEST.g[i]).c_str(), pwt(BEST.g[i]));
        ap("\nWeight enumerator of the 256-element group:\n");
        for (int m = 0; m <= NQ; ++m) if (BEST.A[m]) ap("  A[%2d] = %d\n", m, BEST.A[m]);
    }
    ap("\nNumber of distinct optimal Stage-1 subspaces found:\n%lld\n", (long long)RAW_DISTINCT);
    ap("\nNumber of inequivalent optimal Stage-1 classes (fingerprint classes):\n%zu\n",
       FAMILY.size());
    ap("\nCERTIFIED LOWER BOUND on |U_8| over ALL 8-generator weight-10 stabilizer codes:\n");
    ap("  |U_8| >= %lld   (linear-programming/weak-duality certificate, value %.4f)\n",
       (long long)B.intBound, B.value);
    ap("  gap to best known: %lld\n", (long long)(bestU - B.intBound));
    ap("\nSymmetry reduction:\n");
    ap("  qubit permutations S_14           : columns treated as a multiset (B&B) and\n");
    ap("                                      quotiented in the fingerprint (local search)\n");
    ap("  local Cliffords (S_3)^14          : a column is reduced to the subspace\n");
    ap("                                      W_j = span{x_j,z_j} <= F_2^8\n");
    ap("  generator relabelling S_8         : the objective depends only on the SUBSPACE L,\n");
    ap("                                      so the k-opt search never enumerates bases\n");
    ap("  raw column types                  : %d\n", 1 << (2 * R1));
    { long double li, lim, lms; int64_t nt8; space_sizes(li, lim, lms, nt8);
      ap("  reduced column types              : %lld\n", (long long)nt8);
      ap("  multisets of 14 reduced types     : %s\n", log10_str(lms).c_str());
      ap("  8-dim isotropic subspaces of F_2^28: %s\n", log10_str(li).c_str());
      ap("  ... modulo 14! * 6^14             : %s\n", log10_str(lim).c_str()); }
    ap("\nSearch effort:\n");
    ap("  restarts                          : %llu\n", (unsigned long long)RESTARTS.load());
    ap("  1-opt classes swept               : %llu\n", (unsigned long long)ONEOPT_NODES.load());
    ap("  2-opt pairs scored                : %llu\n", (unsigned long long)TWOOPT_NODES.load());
    ap("  elapsed                           : %s\n", hms(elapsed()).c_str());
    ap("  stopped because                   : %s\n",
       timeLimited ? "time limit reached" : "restart budget exhausted");
    dump(OUTDIR + "/stage1_summary.txt", false);

    RB.clear();
    ap("ALL INEQUIVALENT BEST-KNOWN STAGE-1 CLASSES\n");
    ap("(fingerprint = weight enumerator + sorted per-qubit and per-pair weight profiles;\n");
    ap(" different fingerprints PROVE inequivalence, equal fingerprints only suggest\n");
    ap(" equivalence, which is why several members per class are carried through Stage 2)\n\n");
    int ci = 0;
    for (auto& kv : FAMILY) {
        ++ci;
        ap("---- class %04d : %llu distinct spans found, %zu members retained ----\n",
           ci, (unsigned long long)kv.second.raw, kv.second.members.size());
        for (size_t m = 0; m < kv.second.members.size(); ++m) {
            ap("  member %zu%s\n", m + 1, kv.second.members[m].twoOpt ? "  [2-opt certified]" : "");
            for (int i = 0; i < R1; ++i) ap("    g%d = %s\n", i + 1,
                                            pstr(kv.second.members[m].s.g[i]).c_str());
            ap("    |U_8| = %lld\n", (long long)kv.second.members[m].s.U8);
        }
        ap("\n");
    }
    dump(OUTDIR + "/stage1_optimal_family.txt", false);

    RB.clear();
    ap("STAGE-2 COMPARISON ACROSS THE BEST-KNOWN STAGE-1 FAMILY\n");
    ap("(Stage 2 is EXHAUSTIVE for every frozen Stage-1 code, so each Stage-2 optimum below\n");
    ap(" is a PROVEN OPTIMUM conditional on that Stage-1 code.)\n\n");
    ap("Orbit    |U8|   U8 by weight        Stage2 det   Final      Residual  resDim  sub  allW4\n");
    ap("--------------------------------------------------------------------------------------\n");
    for (const auto& r : REPORTS) {
        ap("%04d   %6lld   %3lld/%4lld/%5lld/%5lld  %8lld  %8lld  %8lld  %4d   %-3s  %-3s\n",
           r.idx, (long long)r.u8, (long long)r.u8byw[1], (long long)r.u8byw[2],
           (long long)r.u8byw[3], (long long)r.u8byw[4], (long long)r.s2cover,
           (long long)r.mfinal, (long long)r.s2res, r.resDim,
           r.resSub ? "yes" : "no", r.allW4 ? "YES" : "no");
    }
    ap("\n");
    int64_t bestFinal = -1, bestRes = INT64_MAX; int bestIdx = -1;
    bool anyFull = false, anyW4 = false;
    for (const auto& r : REPORTS) {
        if (r.mfinal > bestFinal) { bestFinal = r.mfinal; bestIdx = r.idx; }
        bestRes = std::min(bestRes, r.s2res);
        if (r.s2found && r.s2res == 0) anyFull = true;
        if (r.allW4) anyW4 = true;
    }
    ap("Stage-1 optimum (BEST KNOWN)          : M_8 = %lld / %lld, |U_8| = %lld\n",
       (long long)(TOT_ERR - bestU), (long long)TOT_ERR, (long long)bestU);
    ap("Best Stage-2 completion               : residual %lld (orbit %04d)\n",
       (long long)(bestRes == INT64_MAX ? -1 : bestRes), bestIdx);
    ap("Best final coverage                   : %lld / %lld  (%.8f %%)\n",
       (long long)bestFinal, (long long)TOT_ERR,
       bestFinal < 0 ? 0.0 : 100.0 * double(bestFinal) / double(TOT_ERR));
    ap("Full coverage found                   : %s\n", anyFull ? "YES" : "NO");
    ap("All weight-4 errors detected          : %s\n", anyW4 ? "YES" : "NO");
    ap("\nCONDITION C -- do different optimal Stage-1 classes behave differently?\n");
    {
        std::set<int64_t> us, cs, rs, ds;
        for (const auto& r : REPORTS) { us.insert(r.u8); cs.insert(r.s2cover);
            rs.insert(r.s2res); ds.insert(r.resDim); }
        ap("  distinct |U_8| values            : %zu\n", us.size());
        ap("  distinct Stage-2 optima          : %zu\n", cs.size());
        ap("  distinct final residual sizes    : %zu\n", rs.size());
        ap("  distinct residual span dimensions: %zu\n", ds.size());
        ap("  => the Stage-1 choice %s the Stage-2 outcome within the best-known family.\n",
           (cs.size() > 1 || rs.size() > 1) ? "DOES change" : "does NOT change");
    }
    dump(OUTDIR + "/stage2_comparison.txt", false);

    RB.clear();
    ap("============================================================\n");
    ap("REVERSE EXPERIMENT  8 x weight-10  ->  3 x weight-8   FINAL SUMMARY\n");
    ap("============================================================\n\n");
    ap("timestamp: %s\nelapsed:   %s\n\n", now_stamp().c_str(), hms(elapsed()).c_str());
    ap("STATUS LABELS\n");
    ap("  Stage 1                : BEST KNOWN     (exhaustive search provably out of reach)\n");
    ap("  Stage 1 lower bound    : PROVEN         (|U_8| >= %lld, weak-duality certificate)\n",
       (long long)B.intBound);
    ap("  Stage 2, each class    : PROVEN OPTIMUM conditional on that frozen Stage-1 code\n");
    ap("  B&B engine             : validated against the r=3 and r=4 optima proven earlier\n");
    ap("  Overall best code      : BEST KNOWN\n\n");
    ap("RESULT LADDER\n");
    ap("  M_8      = %lld / %lld   (|U_8| = %lld)\n",
       (long long)(TOT_ERR - bestU), (long long)TOT_ERR, (long long)bestU);
    if (!REPORTS.empty()) {
        ap("  M_3(2)   = %lld  (best Stage-2 coverage of U_8)\n",
           (long long)(bestRes == INT64_MAX ? -1 : bestU - bestRes));
        ap("  M_final  = %lld / %lld  (%.8f %%)\n", (long long)bestFinal, (long long)TOT_ERR,
           bestFinal < 0 ? 0.0 : 100.0 * double(bestFinal) / double(TOT_ERR));
        ap("  residual = %lld\n", (long long)(bestRes == INT64_MAX ? -1 : bestRes));
    }
    ap("\nCONDITION A (some class admits |U_final| = 0)      : %s\n", anyFull ? "YES" : "NO");
    ap("CONDITION B (all %lld weight-4 errors detected)  : %s\n",
       (long long)TOT_BY_W[4], anyW4 ? "YES" : "NO");
    ap("\nWHAT MAY NOT BE CLAIMED\n");
    { long double li, lim, lms; int64_t nt; space_sizes(li, lim, lms, nt);
      ap("  A negative Condition A here is NOT an impossibility proof.  Stage 2 is exhaustive,\n");
      ap("  but only for the Stage-1 codes that were actually found, and the Stage-1 search is\n");
      ap("  a heuristic over a space of about %s inequivalent 8-dimensional\n", log10_str(lim).c_str());
      ap("  isotropic subspaces.  Nothing here rules out a [[14,3,5]] code built the other way\n");
      ap("  round, or one sitting on a Stage-1 code that was never visited.\n"); }
    dump(OUTDIR + "/final_summary.txt", true);
}

// ===========================================================================================
//  CALIBRATION  --  every runtime estimate the program prints comes from these measurements
// ===========================================================================================
static void calibrate(const Opts& O) {
    logline("========================= CALIBRATION =========================\n");
    logline("threads available: %d (using %d)\n\n", (int)std::thread::hardware_concurrency(), O.threads);
    std::mt19937_64 rng(O.seed);

    // ---- one restart = construction + exhaustive 1-opt to a local optimum
    double t0 = elapsed(); int nres = 0; S1 keep; bool haveKeep = false;
    int64_t sumU = 0, bestU = INT64_MAX;
    while (elapsed() - t0 < 6.0 && nres < 200) {
        S1 S;
        if (!random_extend(S, 0, rng, O.ncand)) continue;
        polish1(S);
        ++nres; sumU += S.U8;
        if (S.U8 < bestU) { bestU = S.U8; keep = S; haveKeep = true; }
    }
    double perRestart = (elapsed() - t0) / std::max(1, nres);
    logline("single-thread restart (construction + 1-opt to a local optimum)\n");
    logline("  restarts measured        : %d\n", nres);
    logline("  seconds per restart      : %.4f\n", perRestart);
    logline("  mean |U_8| reached       : %.1f\n", nres ? double(sumU) / nres : 0.0);
    logline("  best  |U_8| in this batch: %lld  (M_8 = %lld)\n\n",
            (long long)bestU, (long long)(TOT_ERR - bestU));
    if (!haveKeep) { logline("calibration could not build a solution\n"); return; }

    // ---- one 1-opt call (all eight drops)
    { S1 S = keep; uint64_t nodes = 0; double a = elapsed();
      one_opt(S, &nodes); double b = elapsed();
      logline("one 1-opt sweep (8 drops, quotient dim %d, %d classes each)\n",
              2 * NQ - 2 * (R1 - 1), 1 << (2 * NQ - 2 * (R1 - 1)));
      logline("  seconds                  : %.4f\n", b - a);
      logline("  classes swept            : %llu\n\n", (unsigned long long)nodes); }

    // ---- one full 2-opt certification (all 28 retained sextuples)
    { S1 S = keep; uint64_t nodes = 0, pruned = 0; double a = elapsed();
      bool got = two_opt(S, &nodes, &pruned); double b = elapsed();
      logline("one full 2-opt certification (28 drop-pairs, quotient dim %d, %d classes)\n",
              2 * NQ - 2 * (R1 - 2), 1 << (2 * NQ - 2 * (R1 - 2)));
      logline("  seconds (single thread)  : %.3f\n", b - a);
      logline("  pairs scored             : %llu\n", (unsigned long long)nodes);
      logline("  prefixes cut by the filter: %llu\n", (unsigned long long)pruned);
      logline("  improved the incumbent   : %s\n\n", got ? "yes" : "no"); }

    // ---- one full 3-opt certification (all 56 retained quintuples)
    { S1 S = keep; uint64_t nodes = 0, pruned = 0; double a = elapsed();
      bool comp3 = true;
      uint64_t c1d = 0, c1t = 0;
      bool got = three_opt(S, &nodes, &pruned, elapsed() + 120.0, &comp3, O.threads, &c1d, &c1t);
      double b = elapsed();
      logline("one full 3-opt certification (56 drop-triples, quotient dim %d, %d classes, "
              "%d threads)\n", 2 * NQ - 2 * (R1 - 3), 1 << (2 * NQ - 2 * (R1 - 3)), O.threads);
      logline("  seconds                  : %.3f (capped at 120)\n", b - a);
      logline("  triples scored           : %llu\n", (unsigned long long)nodes);
      logline("  drop-triples fully cut   : %llu\n", (unsigned long long)pruned);
      logline("  c1 prefixes retired      : %llu of %llu = %.4f %%\n",
              (unsigned long long)c1d, (unsigned long long)c1t,
              c1t ? 100.0 * double(c1d) / double(c1t) : 0.0);
      logline("  completed                : %s\n", comp3 ? "yes" : "NO (cut short)");
      if (!comp3 && c1d)
          logline("  projected FULL 3-opt     : %s on %d threads\n",
                  hms((b - a) * double(c1t) / double(c1d)).c_str(), O.threads);
      logline("  improved the incumbent   : %s\n", got ? "yes" : "no");
      if (got) { logline("  new |U_8|                : %lld\n", (long long)S.U8); keep = S; }
      logline("\n"); }

    // ---- one exhaustive Stage-2
    { std::vector<uint32_t> U8; build_U8(keep, U8);
      S2Result R; double a = elapsed();
      stage2_exhaustive(keep, U8, R, O.threads); double b = elapsed();
      logline("one exhaustive Stage-2 (frozen Stage-1, quotient dim %d, %d classes)\n",
              R.qdim, R.qsize);
      logline("  |U_8|                    : %zu\n", U8.size());
      logline("  weight-8-usable classes  : %d of %d\n", R.usable, R.qsize - 1);
      logline("  triples scored           : %llu\n", (unsigned long long)R.nodes);
      logline("  seconds (%d threads)      : %.3f\n", O.threads, b - a);
      if (R.found) logline("  Stage-2 optimum          : detects %lld of %zu, residual %lld\n\n",
                           (long long)R.cover, U8.size(), (long long)R.residual);
      else logline("  Stage-2                  : no valid triple\n\n"); }

    // ---- branch and bound scaling
    logline("column-multiset branch and bound (the exhaustive engine), generator weight %d\n", W1);
    logline("  NOTE: these rows are the rank-r problem at weight %d, which is a DIFFERENT\n", W1);
    logline("  problem from the weight-8 regressions in --selftest; the values will differ.\n");
    logline("  rank  types   nodes        leaves       seconds   completed  M at weight %d\n", W1);
    for (int r = 3; r <= 6; ++r) {
        BBTypes T; build_types(r, T);
        BB_ABORTED = false;
        BB_DEADLINE = (r >= 5) ? elapsed() + (O.bbtime > 0 ? O.bbtime : 20.0) : 0;
        BBResult R; bb_search(r, W1, O.symDepth, O.threads, 1, R);
        logline("  %4d  %5d   %-11llu  %-11llu  %7.2f   %-9s  %lld\n", r, T.ntypes,
                (unsigned long long)R.nodes, (unsigned long long)R.leaves, R.seconds,
                R.completed ? "yes" : "NO (cut)", (long long)R.M);
        BB_DEADLINE = 0; BB_ABORTED = false;
    }
    for (int r = 7; r <= 8; ++r) {
        int64_t nt = bb_type_count(r);
        long double lm = 0;
        for (int i = 0; i < NQ; ++i) lm += log10l((long double)(nt + 13 - i));
        for (int i = 2; i <= NQ; ++i) lm -= log10l((long double)i);
        logline("  %4d  %5lld   -- not run: %s states, and the engine stops at rank %d --\n",
                r, (long long)nt, log10_str(lm).c_str(), bb_max_rank());
    }
    logline("\n");

    // ---- projections
    logline("PROJECTED RUNTIMES on %d threads (from the measurements above)\n", O.threads);
    logline("  restarts        wall clock\n");
    const int64_t budgets[] = { 1000, 10000, 100000, 1000000 };
    for (int64_t b : budgets)
        logline("  %-14lld  %s\n", (long long)b, hms(double(b) * perRestart / O.threads).c_str());
    logline("\n  A restart is one construction plus a complete 1-opt descent.  Phase 3 adds one\n");
    logline("  full 2-opt certification per retained family member, and phase 4 one exhaustive\n");
    logline("  Stage-2 per member; both are reported above and are small next to phase 2.\n");
    logline("===============================================================\n");
}

// ===========================================================================================
static void usage() {
    printf(
    "reverse_8w10_3w8 -- reverse two-stage search: 8 x weight-10, then 3 x weight-8\n\n"
    "  --mode heuristic|hybrid|exact   default hybrid (Stage-1 heuristic + exhaustive Stage 2)\n"
    "  --restarts N        Stage-1 multi-start budget (0 = run until --time-limit)\n"
    "  --time-limit S      seconds for the Stage-1 multi-start phase\n"
    "  --threads N         default: all hardware threads\n"
    "  --seed N            RNG seed\n"
    "  --deterministic     fixed seed, single-threaded restart ordering\n"
    "  --candidates N      greedy candidate samples per generator during construction (24)\n"
    "  --members N         members retained per fingerprint class (4)\n"
    "  --classes N         maximum fingerprint classes retained (512)\n"
    "  --no-twoopt         skip the exhaustive 2-opt certification phase\n"
    "  --threeopt          also run the exhaustive 3-opt certification (much slower)\n"
    "  --threeopt-time S   per-member wall-clock cap for 3-opt (0 = none)\n"
    "  --seed-bb R         seed restarts from proven-optimal rank-R weight-10 configurations\n"
    "  --sym-depth N       S_r symmetry-breaking depth in the branch and bound (6)\n"
    "  --bb R W            run ONLY the branch and bound at rank R, generator weight W\n"
    "  --bb-time S         wall-clock cap for a branch-and-bound run\n"
    "  --bb-all            keep tie subtrees so every optimal leaf is enumerated (slower)\n"
    "  --report-every S    status block interval, seconds (1800)\n"
    "  --out DIR           output directory (default: the current directory)\n"
    "  --resume            continue from checkpoint.txt\n"
    "  --selftest [N]      run the correctness suite (N random subspaces, default 400);\n"
    "                      it closes the rank-4 branch-and-bound tree, so allow ~10 minutes\n"
    "  --calibrate         measure the engine and project runtimes, then exit\n"
    "  --spacesize         print the search-space accounting and exit\n"
    "  --force-exact       allow --mode exact at rank 8 even though it cannot finish\n"
    "  --help\n");
}

int main(int argc, char** argv) {
    set_start_time();
    Opts O;
    O.threads = (int)std::thread::hardware_concurrency(); if (O.threads <= 0) O.threads = 4;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&](int d = 1) { return (i + d < argc) ? argv[i + d] : nullptr; };
        if (a == "--help") { usage(); return 0; }
        else if (a == "--mode" && next()) O.mode = argv[++i];
        else if (a == "--restarts" && next()) O.restarts = atoi(argv[++i]);
        else if (a == "--time-limit" && next()) O.tlimit = atof(argv[++i]);
        else if (a == "--threads" && next()) O.threads = std::max(1, atoi(argv[++i]));
        else if (a == "--seed" && next()) O.seed = strtoull(argv[++i], nullptr, 10);
        else if (a == "--deterministic") { O.deterministic = true; O.threads = 1; }
        else if (a == "--candidates" && next()) O.ncand = std::max(1, atoi(argv[++i]));
        else if (a == "--members" && next()) O.memcap = std::max(1, atoi(argv[++i]));
        else if (a == "--classes" && next()) O.classcap = std::max(1, atoi(argv[++i]));
        else if (a == "--no-twoopt") O.twoopt = false;
        else if (a == "--threeopt") O.threeopt = true;
        else if (a == "--threeopt-time" && next()) O.threeoptTime = atof(argv[++i]);
        else if (a == "--seed-bb" && next()) O.seedbb = atoi(argv[++i]);
        else if (a == "--sym-depth" && next()) O.symDepth = atoi(argv[++i]);
        else if (a == "--bb" && next(2)) { O.bbr = atoi(argv[i+1]); O.bbw = atoi(argv[i+2]); i += 2; }
        else if (a == "--bb-time" && next()) O.bbtime = atof(argv[++i]);
        else if (a == "--report-every" && next()) O.reportEvery = atof(argv[++i]);
        else if (a == "--out" && next()) O.outdir = argv[++i];
        else if (a == "--resume") O.resume = true;
        else if (a == "--selftest") { O.selftest = 400;
            if (next() && argv[i+1][0] != '-') O.selftest = atoi(argv[++i]); }
        else if (a == "--calibrate") O.calibrate = 1;
        else if (a == "--spacesize") O.calibrate = 2;
        else if (a == "--force-exact") O.forceExact = true;
        else if (a == "--bb-all") BB_KEEP_TIES = true;
        else { fprintf(stderr, "unknown option: %s\n", a.c_str()); usage(); return 1; }
    }
    OUTDIR = O.outdir; MEMCAP = O.memcap; CLASSCAP = O.classcap; REPORT_EVERY = O.reportEvery;
    std::error_code ec; std::filesystem::create_directories(OUTDIR, ec);
    LOGF = fopen((OUTDIR + "/reverse_search.log").c_str(), O.resume ? "a" : "w");

    build_tables();
    build_errors();
    build_pmin();

    logline("=========================================================================\n");
    logline("REVERSE TWO-STAGE SEARCH   8 x weight-%d  ->  3 x weight-%d\n", W1, W2);
    logline("started %s   threads %d   mode %s\n", now_stamp().c_str(), O.threads, O.mode.c_str());
    logline("=========================================================================\n");
    logline("weight-<=4 Pauli errors enumerated: %lld   (weight 4 alone: %lld)\n",
            (long long)TOT_ERR, (long long)TOT_BY_W[4]);

    if (O.calibrate == 2) {
        long double li, lim, lms; int64_t nt8; space_sizes(li, lim, lms, nt8);
        logline("\nSEARCH-SPACE ACCOUNTING\n");
        logline("  8-dim isotropic subspaces of F_2^28          : %s\n", log10_str(li).c_str());
        logline("  ... modulo qubit permutations and local Cliffords: %s\n", log10_str(lim).c_str());
        logline("  raw column types (2^16)                      : %d\n", 1 << (2 * R1));
        logline("  reduced column types (subspaces of F_2^8)    : %lld\n", (long long)nt8);
        logline("  multisets of 14 reduced types                : %s\n", log10_str(lms).c_str());
        logline("\n  For comparison the same accounting at the ranks that WERE proven earlier:\n");
        for (int r = 3; r <= 5; ++r) {
            int64_t nt = 1 + ((int64_t)1 << r) - 1 + gauss_binom(r, 2);
            long double lm = 0;
            for (int i = 0; i < NQ; ++i) lm += log10l((long double)(nt + 13 - i));
            for (int i = 2; i <= NQ; ++i) lm -= log10l((long double)i);
            logline("    rank %d: %lld types, %s multisets\n", r, (long long)nt, log10_str(lm).c_str());
        }
        logline("\nConclusion: an exhaustive Stage-1 search at rank 8 is not possible.\n");
        return 0;
    }
    if (O.bbr > bb_max_rank()) {
        long double li, lim, lms; int64_t nt8; space_sizes(li, lim, lms, nt8);
        logline("\nThe column-multiset branch and bound is implemented for rank <= %d.\n",
                bb_max_rank());
        logline("At rank %d there are %lld column types and about %s multisets;\n",
                O.bbr, (long long)bb_type_count(O.bbr),
                O.bbr == R1 ? log10_str(lms).c_str() : "an astronomical number of");
        logline("no amount of pruning brings that into reach.  Use --mode hybrid.\n");
        return 3;
    }
    if (O.bbr > 0) {
        BB_ABORTED = false;
        if (O.bbtime > 0) BB_DEADLINE = elapsed() + O.bbtime;
        BBTypes T; build_types(O.bbr, T);
        logline("\nBRANCH AND BOUND  rank %d, generator weight %d, %d column types\n",
                O.bbr, O.bbw, T.ntypes);
        BBResult R; bb_search(O.bbr, O.bbw, O.symDepth, O.threads, 8, R);
        logline("  nodes %llu  leaves %llu  pruned %llu  %.2f s\n",
                (unsigned long long)R.nodes, (unsigned long long)R.leaves,
                (unsigned long long)R.pruned, R.seconds);
        logline("  tree completed           : %s\n", R.completed ? "YES" : "NO (cut short)");
        logline("  best sum P               : %lld\n", (long long)R.bestSum);
        logline("  |U| = %lld   M = %lld / %lld\n", (long long)R.U, (long long)R.M,
                (long long)TOT_ERR);
        logline("  raw optimal leaves       : %lld\n", (long long)R.rawOptimal);
        logline("  status                   : %s\n",
                R.completed ? "PROVEN GLOBAL OPTIMUM for this rank and weight" : "BEST KNOWN");
        for (size_t k = 0; k < R.gens.size(); ++k) {
            logline("  solution %zu:\n", k + 1);
            for (int i = 0; i < O.bbr; ++i)
                logline("    g%d = %s  (weight %d)\n", i + 1,
                        pstr(R.gens[k][i]).c_str(), pwt(R.gens[k][i]));
        }
        return R.completed ? 0 : 2;
    }
    if (O.selftest) { bool ok = selftest(O.selftest, O.threads, O.bbtime > 0 ? O.bbtime : 1800.0);
        return ok ? 0 : 1; }
    if (O.calibrate == 1) { calibrate(O); return 0; }

    if (O.mode == "exact" && !O.forceExact) {
        long double li, lim, lms; int64_t nt8; space_sizes(li, lim, lms, nt8);
        logline("\n--mode exact requests an exhaustive Stage-1 branch and bound at rank 8.\n");
        logline("That tree has about %s states before pruning; the engine measures\n",
                log10_str(lms).c_str());
        logline("a few million nodes per second, so it cannot finish in any human timeframe.\n");
        logline("Re-run with --force-exact if you want it started anyway, or use --mode hybrid.\n");
        return 3;
    }

    BoundResult B; certified_bound(B, 40, O.seed);
    logline("\nCERTIFIED LOWER BOUND on |U_8| over every 8-generator weight-10 code: %lld\n",
            (long long)B.intBound);
    logline("  (weak-duality certificate, value %.4f -- valid for any multipliers, so this\n",
            B.value);
    logline("   number cannot be wrong, only weak)\n");

    int64_t restartsDone = 0;
    if (O.resume && load_checkpoint(restartsDone))
        logline("\nRESUMED from checkpoint: %lld restarts, best |U_8| = %lld, %zu classes\n",
                (long long)restartsDone, (long long)(BEST_SUM / S1SZ), FAMILY.size());

    if (O.seedbb > 0 && O.seedbb <= 6) {
        BB_ABORTED = false;
        if (O.bbtime > 0) BB_DEADLINE = elapsed() + O.bbtime;
        logline("\n[phase 1] seeding from rank-%d weight-%d configurations\n", O.seedbb, W1);
        BBResult R; bb_search(O.seedbb, W1, O.symDepth, O.threads, 4000, R);
        logline("  nodes %llu, %.2f s, completed %s, |U| = %lld, %zu configurations kept\n",
                (unsigned long long)R.nodes, R.seconds, R.completed ? "yes" : "no",
                (long long)R.U, R.gens.size());
        for (auto& g : R.gens) SEEDPOOL.push_back(g);
        BB_DEADLINE = 0; BB_ABORTED = false;
    }

    if (O.restarts <= 0 && O.tlimit <= 0) O.restarts = 2000;
    phase_multistart(O);
    if (O.twoopt) phase_twoopt(O);
    if (O.threeopt) phase_threeopt(O);
    phase_stage2(O);
    write_reports(O, B, O.restarts <= 0);
    logline("\nfinished %s   total elapsed %s\n", now_stamp().c_str(), hms(elapsed()).c_str());
    if (LOGF) fclose(LOGF);
    return 0;
}
