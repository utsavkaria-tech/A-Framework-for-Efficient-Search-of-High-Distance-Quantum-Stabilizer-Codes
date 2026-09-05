// ===========================================================================================
//  stab14_refine.cpp -- TARGETED refinement around a known 3+8 stabilizer solution.
//
//  Stage-1 stabilizer S3 (frozen, part of the proven-optimal family, orbit 0006) and an
//  8-generator weight-10 extension with M_8 = 11181 / 11186 (5 residual errors) are given.
//  This program searches the neighbourhood of that exact incumbent for a smaller residual.
//
//  It is a SEPARATE program: stab14_2stage.cpp and every existing result file are left
//  untouched, and all output goes under refine/ .
//
//  =========================================================================================
//  THE SEARCH DESIGN, AND WHY IT IS NOT "THE OLD HEURISTIC WITH A NEW SEED"
//  =========================================================================================
//  Everything lives in V = C(S3)/S3, dim 25-3 = 22.  The eight generators span an
//  8-dimensional totally isotropic D <= V, and (Poisson summation, re-verified at run time)
//
//        covered(D) = (1/128) * SUM_{v in D\{0}} cov1(v),    cov1(v) = #{E in U_3 : <v,E>=1}.
//
//  So the objective is a SUM OVER THE 255 NON-ZERO ELEMENTS OF THE SUBSPACE.  That is what
//  makes an exact k-generator-replacement search affordable:
//
//    Let K = the (8-k)-dimensional subspace spanned by the generators we KEEP, and
//        f_K(w) = SUM_{u in K} cov1(w ^ u)        (the coset sum at w).
//    Then for any new generators a_1..a_k,
//        S(D') = S(K) + SUM over the 2^k - 1 non-zero combinations c of  f_K(c).
//
//    f_K is a XOR-convolution of cov1 with the indicator of K, so it can be tabulated for
//    ALL 2^22 classes at once by (8-k) in-place butterfly passes -- 4.2M additions each,
//    about 20 ms -- instead of 2^(8-k) lookups per candidate.  This is the whole trick.
//
//    k = 1:  S = S(K) + f(a)                                  -> exhaustive, ~16K candidates
//    k = 2:  S = S(K) + f(a) + f(b) + f(a^b)                   -> exhaustive over ~32K^2/2
//            pairs, with the bound f(a^b) <= max f giving a hard break once the sorted-by-f
//            candidate list drops below the threshold.
//    k = 3:  fix a, put g(w) = f(w) + f(w^a); then
//            S = S(K) + f(a) + g(b) + g(c) + g(b^c)
//            -- the SAME three-term pair form.  One tabulation pass builds g, and the k=2
//            routine is reused verbatim.
//
//  Consequences that matter scientifically:
//    * the 1-opt and 2-opt neighbourhoods are searched EXHAUSTIVELY, so if nothing is found
//      the incumbent is CERTIFIED optimal against every replacement of one or two of its
//      eight generators (over all weight-10 classes, not a sample);
//    * the 3-opt neighbourhood is searched exhaustively in its last two coordinates for each
//      first choice, so it is far stronger than a randomized rebuild;
//    * residual-targeted restriction: to kill a residual r we need some new generator v with
//      <v,r> = 1, which simply halves the candidate list -- a cheap, exact way to aim the
//      k = 3 search at specific residual errors instead of hoping.
//
//  The union objective is never proxied: the score above IS the exact union coverage, and
//  every improvement is re-verified from scratch against all 91770 errors.
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
#include <algorithm>
#include <numeric>
#include <random>
#include <atomic>
#include <mutex>
#include <thread>
#include <chrono>
#include <bit>
#include <unordered_set>
#include <filesystem>

constexpr int NQ    = 14;
constexpr int R1    = 3;                      // frozen Stage-1 generators
constexpr int R2    = 8;                      // Stage-2 generators
constexpr int MAXW  = 4;
constexpr int W1    = 8, W2 = 10;
constexpr uint32_t QM   = (1u << NQ) - 1u;
constexpr uint32_t NONE = 0xFFFFFFFFu;

constexpr int VDIM  = 2 * NQ - 2 * R1;        // 22
constexpr int VHALF = VDIM / 2;               // 11
constexpr int CDIM  = 2 * NQ - R1;            // 25
constexpr uint32_t VLO = (1u << VHALF) - 1u;
constexpr int VSZ   = 1 << VDIM;              // 4194304

static inline uint32_t pk(uint32_t x, uint32_t z) { return x | (z << 16); }
static inline uint32_t px(uint32_t p) { return p & QM; }
static inline uint32_t pz(uint32_t p) { return (p >> 16) & QM; }
static inline int pwt(uint32_t p) { return std::popcount((p | (p >> 16)) & QM); }
static inline int par(unsigned a) { return std::popcount(a) & 1; }
static inline int symp(uint32_t a, uint32_t b) { return par((px(a) & pz(b)) ^ (pz(a) & px(b))); }
static inline int csymp(uint32_t a, uint32_t b) {
    return std::popcount(((a & VLO) & (b >> VHALF)) ^ ((a >> VHALF) & (b & VLO))) & 1;
}
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
            default: fprintf(stderr, "bad Pauli char '%c' at %d\n", s[j], j); exit(1);
        }
    }
    return pk(x, z);
}

// ---------------------------------------------------------------- the frozen incumbent
static const char* G3_STR[R1] = {
    "IIIIXXXXIIZZXX",
    "IIIIXXIIXXXXZZ",
    "IIIIIIZZZZZZZZ"
};
static const char* H8_STR[R2] = {
    "IXIXXXIZYZIXXX",
    "IZXIXXYIXYIYXX",
    "IIIZYXZXYXZIZY",
    "IIIYXZYIYZZYZX",
    "IIZIXZXYZXZYIZ",
    "XIIIXYXZXYIXYX",
    "IIYIZXZXYXXIYX",
    "ZIIIZXIZZYZXYX"
};
static uint32_t G3[R1], H0[R2];

static std::chrono::steady_clock::time_point T0;
static double elapsed() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count(); }

static std::string RB;
static void ap(const char* fmt, ...) {
    char buf[8192]; va_list a; va_start(a, fmt); vsnprintf(buf, sizeof(buf), fmt, a); va_end(a); RB += buf;
}
static void dump(const std::string& path, bool echo = false) {
    if (echo) { fputs(RB.c_str(), stdout); fflush(stdout); }
    FILE* f = fopen(path.c_str(), "w"); if (f) { fputs(RB.c_str(), f); fclose(f); }
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

static std::vector<uint32_t> U3;          // Stage-1 undetected errors (the targets)
static int64_t TOT_ERR = 0, TOTW4 = 0;
static int64_t U3_BY_W[MAXW + 1] = {0};

static void build_targets() {
    U3.clear(); TOT_ERR = 0; TOTW4 = 0;
    for (int w = 0; w <= MAXW; ++w) U3_BY_W[w] = 0;
    for_each_error([&](uint32_t E, int w) {
        ++TOT_ERR; if (w == MAXW) ++TOTW4;
        int s = 0; for (int i = 0; i < R1; ++i) s |= symp(G3[i], E) << i;
        if (!s) { U3.push_back(E); U3_BY_W[w]++; }
    });
}

// ---------------------------------------------------------------- quotient V = C(S3)/S3
static uint32_t UB[VHALF], WB[VHALF], S3EL[1 << R1];
static inline uint32_t cls_of(uint32_t v) {
    uint32_t c = 0;
    for (int i = 0; i < VHALF; ++i) {
        c |= uint32_t(symp(WB[i], v)) << i;
        c |= uint32_t(symp(UB[i], v)) << (i + VHALF);
    }
    return c;
}
static void build_basis() {
    std::vector<uint32_t> cur;
    {
        std::vector<uint32_t> b;
        for (int j = 0; j < NQ; ++j) { b.push_back(pk(1u << j, 0)); b.push_back(pk(0, 1u << j)); }
        for (int g = 0; g < R1; ++g) {
            int piv = -1;
            for (size_t i = 0; i < b.size(); ++i) if (symp(G3[g], b[i])) { piv = (int)i; break; }
            if (piv < 0) { fprintf(stderr, "Stage-1 generators dependent\n"); exit(1); }
            uint32_t pv = b[piv]; b.erase(b.begin() + piv);
            for (auto& q : b) if (symp(G3[g], q)) q ^= pv;
        }
        cur = b;
    }
    if ((int)cur.size() != CDIM) { fprintf(stderr, "dim C(S3) wrong\n"); exit(1); }
    int np = 0;
    while (np < VHALF) {
        int pi = -1, qi = -1;
        for (size_t i = 0; i < cur.size() && pi < 0; ++i)
            for (size_t j = i + 1; j < cur.size(); ++j)
                if (symp(cur[i], cur[j])) { pi = (int)i; qi = (int)j; break; }
        if (pi < 0) break;
        uint32_t u = cur[pi], w = cur[qi];
        cur.erase(cur.begin() + qi); cur.erase(cur.begin() + pi);
        for (auto& r : cur) { if (symp(r, w)) r ^= u; if (symp(r, u)) r ^= w; }
        UB[np] = u; WB[np] = w; ++np;
    }
    if (np != VHALF || (int)cur.size() != R1) { fprintf(stderr, "symplectic reduction failed\n"); exit(1); }
    for (int m = 0; m < (1 << R1); ++m) { uint32_t t = 0;
        for (int g = 0; g < R1; ++g) if ((m >> g) & 1) t ^= G3[g];
        S3EL[m] = t; }
    for (uint32_t r : cur) { bool in = false;
        for (int m = 0; m < (1 << R1); ++m) if (S3EL[m] == r) in = true;
        if (!in) { fprintf(stderr, "radical is not S3\n"); exit(1); } }
}

static std::vector<uint16_t> COV1;
static std::vector<uint32_t> REPR10;
static std::vector<uint32_t> TCLS;
static std::vector<uint32_t> CAND;            // weight-10-realisable classes
static int COVMAX = 0;
static int64_t NCLS10 = 0, NW10 = 0;

static void build_cov1() {
    TCLS.resize(U3.size());
    std::vector<int32_t> f(VSZ, 0);
    for (size_t j = 0; j < U3.size(); ++j) {
        uint32_t c = cls_of(U3[j]); TCLS[j] = c;
        f[((c >> VHALF) & VLO) | ((c & VLO) << VHALF)] += 1;
    }
    for (int b = 0; b < VDIM; ++b) {
        int step = 1 << b;
        for (int i = 0; i < VSZ; i += step << 1)
            for (int k = i; k < i + step; ++k) {
                int32_t a = f[k], d = f[k + step]; f[k] = a + d; f[k + step] = a - d; }
    }
    const int32_t N = (int32_t)U3.size();
    COV1.resize(VSZ);
    for (int v = 0; v < VSZ; ++v) { int32_t c = (N - f[v]) / 2; COV1[v] = uint16_t(c);
        if (v && c > COVMAX) COVMAX = c; }
}
static void build_repr10() {
    REPR10.assign(VSZ, NONE);
    uint32_t basis[CDIM], bcls[CDIM];
    for (int g = 0; g < R1; ++g) { basis[g] = G3[g]; bcls[g] = 0; }
    for (int i = 0; i < VHALF; ++i) { basis[R1 + i] = UB[i]; bcls[R1 + i] = 1u << i; }
    for (int i = 0; i < VHALF; ++i) { basis[R1 + VHALF + i] = WB[i]; bcls[R1 + VHALF + i] = 1u << (i + VHALF); }
    uint32_t cur = 0, ccl = 0;
    for (uint64_t i = 1; i < (1ull << CDIM); ++i) {
        int b = std::countr_zero(i);
        cur ^= basis[b]; ccl ^= bcls[b];
        if (pwt(cur) == W2) { ++NW10;
            if (ccl && (REPR10[ccl] == NONE || cur < REPR10[ccl])) {
                if (REPR10[ccl] == NONE) ++NCLS10; REPR10[ccl] = cur; } }
    }
    CAND.clear(); CAND.reserve(NCLS10);
    for (uint32_t v = 1; v < (uint32_t)VSZ; ++v) if (REPR10[v] != NONE) CAND.push_back(v);
}

// ===========================================================================================
//  SOLUTION OBJECT, EXACT SCORE, RESIDUAL ANALYSIS, LEXICOGRAPHIC ORDER
// ===========================================================================================
struct Sol {
    uint32_t cls[R2]{};                 // the eight classes (they span D)
    uint32_t pau[R2]{};                 // a concrete weight-10 Pauli per class
    int64_t  span_sum = -1;             // SUM over the 255 non-zero elements of D
    int64_t  covered = -1, residual = -1;
    int64_t  res_w[MAXW + 1] = {0};
    bool     valid = false;
};

// preferred Pauli representative per class: the incumbent's own strings are kept so that an
// unchanged generator is reported exactly as the user supplied it
static std::vector<std::pair<uint32_t, uint32_t>> PREF;
static uint32_t pauli_for(uint32_t c) {
    for (auto& p : PREF) if (p.first == c) return p.second;
    return REPR10[c];
}

static inline int64_t span_sum_of(const uint32_t* cls, int d) {
    int64_t s = 0;
    for (int m = 1; m < (1 << d); ++m) {
        uint32_t v = 0; int mm = m;
        while (mm) { int b = std::countr_zero(unsigned(mm)); mm &= mm - 1; v ^= cls[b]; }
        s += COV1[v];
    }
    return s;
}
static bool independent_isotropic(const uint32_t* cls) {
    for (int a = 0; a < R2; ++a) {
        if (!cls[a]) return false;
        for (int b = a + 1; b < R2; ++b) if (csymp(cls[a], cls[b])) return false;
    }
    uint32_t piv[VDIM] = {0};
    for (int a = 0; a < R2; ++a) {
        uint32_t v = cls[a];
        while (v) { int b = 31 - std::countl_zero(v);
            if (piv[b]) v ^= piv[b]; else { piv[b] = v; break; } }
        if (!v) return false;
    }
    return true;
}
// residual list: targets whose class is orthogonal to every generator
static void residual_of(const uint32_t* cls, std::vector<uint32_t>& out) {
    out.clear();
    for (size_t j = 0; j < U3.size(); ++j) {
        uint32_t t = TCLS[j];
        int hit = 0;
        for (int a = 0; a < R2; ++a) if (csymp(cls[a], t)) { hit = 1; break; }
        if (!hit) out.push_back(U3[j]);
    }
}
static Sol make_sol(const uint32_t* cls) {
    Sol S; std::copy(cls, cls + R2, S.cls); std::sort(S.cls, S.cls + R2);
    S.span_sum = span_sum_of(S.cls, R2);
    S.covered  = S.span_sum / 128;
    S.residual = (int64_t)U3.size() - S.covered;
    std::vector<uint32_t> res; residual_of(S.cls, res);
    for (int w = 0; w <= MAXW; ++w) S.res_w[w] = 0;
    for (uint32_t E : res) S.res_w[pwt(E)]++;
    for (int a = 0; a < R2; ++a) S.pau[a] = pauli_for(S.cls[a]);
    S.valid = true;
    return S;
}
// Lexicographic objective: fewest residual errors, then fewest weight-4, weight-3, weight-2,
// weight-1 residuals.  Deliberately NOT a total order: two solutions with the same residual
// count AND the same weight profile are considered equal, so a mere relabelling of the same
// quality never counts as an improvement, never triggers a NEW BEST report, and cannot make
// the local search cycle.  Note a tie in the residual COUNT can still be a strict improvement
// in the weight profile, which is why equal-score candidates are still enumerated.
static bool better(const Sol& a, const Sol& b) {
    if (!b.valid) return a.valid;
    if (a.residual != b.residual) return a.residual < b.residual;
    for (int w = MAXW; w >= 1; --w) if (a.res_w[w] != b.res_w[w]) return a.res_w[w] < b.res_w[w];
    return false;
}

struct ResidualInfo {
    int64_t n = 0, by_w[MAXW + 1] = {0};
    bool subgroup = false, all_x = false;
    int  dim = 0;
    std::vector<uint32_t> gens, list, qubits;
};
static ResidualInfo analyse(const std::vector<uint32_t>& R) {
    ResidualInfo r; r.n = (int64_t)R.size(); r.list = R;
    for (uint32_t E : R) r.by_w[pwt(E)]++;
    uint32_t piv[2 * NQ] = {0};
    for (uint32_t E : R) {
        uint32_t v = px(E) | (pz(E) << NQ);
        while (v) { int b = 31 - std::countl_zero(v);
            if (piv[b]) v ^= piv[b];
            else { piv[b] = v; r.gens.push_back(pk(v & QM, (v >> NQ) & QM)); break; } }
    }
    r.dim = (int)r.gens.size();
    std::unordered_set<uint32_t> S; S.insert(0u); for (uint32_t E : R) S.insert(E);
    r.subgroup = true;
    for (uint32_t a : S) { for (uint32_t b : S) if (!S.count(a ^ b)) { r.subgroup = false; break; }
        if (!r.subgroup) break; }
    if (r.subgroup && (int64_t)S.size() != (1LL << r.dim)) r.subgroup = false;
    r.all_x = true; uint32_t touched = 0;
    for (uint32_t E : R) { if (pz(E)) r.all_x = false; touched |= px(E) | pz(E); }
    for (int j = 0; j < NQ; ++j) if ((touched >> j) & 1) r.qubits.push_back((uint32_t)j);
    return r;
}

// ---------------------------------------------------------------- independent verification
struct Verify {
    bool wt_ok = true, cg_ok = true, ch_ok = true, rank_ok = true, u3_ok = true, score_ok = true;
    int weights[R2]{}, rank = 0, comm[R1 + R2][R1 + R2]{};
    int64_t u3size = 0, covered = 0, per_gen[R2]{}, cov_by_w[MAXW + 1]{}, tgt_by_w[MAXW + 1]{};
    int64_t w4_total = 0, w4_det = 0, w4_und = 0;
    int64_t all_total = 0, orig_det = 0, det11 = 0;
    std::vector<uint32_t> leftover;
};
// Rebuilds the entire 91770-error universe from scratch and recomputes everything directly.
static Verify verify(const Sol& S) {
    Verify V;
    for (int a = 0; a < R2; ++a) { V.weights[a] = pwt(S.pau[a]); if (V.weights[a] != W2) V.wt_ok = false; }
    uint32_t all[R1 + R2];
    for (int i = 0; i < R1; ++i) all[i] = G3[i];
    for (int a = 0; a < R2; ++a) all[R1 + a] = S.pau[a];
    for (int i = 0; i < R1 + R2; ++i) for (int j = 0; j < R1 + R2; ++j) {
        V.comm[i][j] = symp(all[i], all[j]);
        if (i != j && V.comm[i][j]) { if (i < R1 || j < R1) V.cg_ok = false; else V.ch_ok = false; } }
    { uint32_t rows[R1 + R2];
      for (int i = 0; i < R1 + R2; ++i) rows[i] = px(all[i]) | (pz(all[i]) << NQ);
      int r = 0;
      for (int b = 0; b < 2 * NQ && r < R1 + R2; ++b) {
          int p = -1; for (int i = r; i < R1 + R2; ++i) if ((rows[i] >> b) & 1) { p = i; break; }
          if (p < 0) continue;
          std::swap(rows[r], rows[p]);
          for (int i = 0; i < R1 + R2; ++i) if (i != r && ((rows[i] >> b) & 1)) rows[i] ^= rows[r];
          ++r; }
      V.rank = r; V.rank_ok = (r == R1 + R2); }
    for_each_error([&](uint32_t E, int w) {
        int s0 = 0; for (int i = 0; i < R1; ++i) s0 |= symp(G3[i], E) << i;
        int s1 = 0; for (int a = 0; a < R2; ++a) s1 |= symp(S.pau[a], E) << a;
        V.all_total++;
        if (s0) V.orig_det++;
        if (s0 || s1) V.det11++;
        if (w == MAXW) { V.w4_total++; if (s0 || s1) V.w4_det++; else V.w4_und++; }
        if (!s0) {
            V.u3size++; V.tgt_by_w[w]++;
            if (s1) { V.covered++; V.cov_by_w[w]++;
                for (int a = 0; a < R2; ++a) if ((s1 >> a) & 1) V.per_gen[a]++; }
            else V.leftover.push_back(E);
        }
    });
    V.u3_ok = (V.u3size == (int64_t)U3.size());
    V.score_ok = (V.covered == S.covered);
    return V;
}

// ===========================================================================================
//  THE k-OPT ENGINE
//  f_K(w) = SUM_{u in K} cov1(w ^ u) is a XOR-convolution of cov1 with the indicator of the
//  subspace K.  Because XOR by a fixed b is an involution, one basis vector is absorbed by a
//  single in-place pass:  F(w) and F(w^b) both become F(w) + F(w^b).  dim K passes therefore
//  tabulate f_K over ALL 2^22 classes in ~20 ms, after which any k-generator replacement is
//  scored with 2^k - 1 array lookups instead of 2^(8-k) per candidate.
// ===========================================================================================
static std::vector<uint32_t> FTAB, GTAB;

static void build_ftab(std::vector<uint32_t>& T, const uint32_t* basis, int k) {
    T.resize(VSZ);
    std::copy(COV1.begin(), COV1.end(), T.begin());        // widen uint16 -> uint32
    for (int t = 0; t < k; ++t) {
        const uint32_t b = basis[t];
        const int p = std::countr_zero(b);
        for (uint32_t w = 0; w < (uint32_t)VSZ; ++w) {
            if ((w >> p) & 1) continue;
            const uint32_t w2 = w ^ b;
            const uint32_t s = T[w] + T[w2];
            T[w] = s; T[w2] = s;                            // f_K is constant on K-cosets
        }
    }
}
static uint32_t table_max(const std::vector<uint32_t>& T) {
    uint32_t m = 0; for (uint32_t v : T) if (v > m) m = v; return m;
}

struct Hit { uint32_t a, b; int64_t score; };

// Exhaustive over unordered pairs from `cls` (sorted so that `val` is non-increasing).
// Records every pair whose exact score exceeds `need`.  The bound score <= base + val[i] +
// val[j] + tabmax turns the sorted order into a hard break, so the swept region shrinks fast
// as the incumbent improves.
static void search_pairs(const std::vector<uint32_t>& cls, const std::vector<uint32_t>& val,
                         const std::vector<uint32_t>& tab, int64_t base, int64_t need,
                         uint32_t tabmax, const std::vector<uint32_t>& kspan,
                         int nthreads, std::vector<Hit>& out, size_t cap) {
    const size_t n = cls.size();
    if (n < 2) return;
    std::atomic<size_t> next{0};
    std::mutex mtx;
    std::vector<std::thread> th;
    for (int t = 0; t < nthreads; ++t) th.emplace_back([&]() {
        std::vector<Hit> local;
        for (;;) {
            size_t i = next.fetch_add(1);
            if (i + 1 >= n) break;
            const int64_t vi = val[i];
            if (base + 2 * vi + (int64_t)tabmax <= need) continue;      // and all larger i too
            const uint32_t ci = cls[i];
            for (size_t j = i + 1; j < n; ++j) {
                const int64_t vj = val[j];
                if (base + vi + vj + (int64_t)tabmax <= need) break;
                const uint32_t cj = cls[j];
                if (csymp(ci, cj)) continue;
                const uint32_t x = ci ^ cj;
                const int64_t s = base + vi + vj + (int64_t)tab[x];
                if (s > need) {
                    if (std::binary_search(kspan.begin(), kspan.end(), x)) continue;  // dependent
                    local.push_back({ci, cj, s});
                    if (local.size() > cap) { std::sort(local.begin(), local.end(),
                        [](const Hit& p, const Hit& q) { return p.score > q.score; });
                        local.resize(cap / 2); }
                }
            }
        }
        if (!local.empty()) { std::lock_guard<std::mutex> lk(mtx);
            out.insert(out.end(), local.begin(), local.end()); }
    });
    for (auto& x : th) x.join();
}

// ---------------------------------------------------------------- incumbent bookkeeping
static Sol BEST, CUR;
static std::mutex BESTMTX;
static int IMPROVEMENTS = 0;
static std::string OUTROOT = "refine";
static std::vector<uint8_t> OMASK;             // per candidate: commutes-with-generator bits

static void rebuild_omask(const Sol& S, int nthreads) {
    OMASK.assign(CAND.size(), 0);
    std::atomic<size_t> next{0};
    std::vector<std::thread> th;
    for (int t = 0; t < nthreads; ++t) th.emplace_back([&]() {
        for (;;) {
            size_t lo = next.fetch_add(65536);
            if (lo >= CAND.size()) break;
            size_t hi = std::min(lo + 65536, CAND.size());
            for (size_t i = lo; i < hi; ++i) {
                uint8_t m = 0;
                for (int a = 0; a < R2; ++a) if (!csymp(CAND[i], S.cls[a])) m |= uint8_t(1u << a);
                OMASK[i] = m;
            }
        }
    });
    for (auto& x : th) x.join();
}

// span of the kept generators, sorted (for the independence test)
static void kept_span(const Sol& S, uint8_t dropmask, std::vector<uint32_t>& keep,
                      std::vector<uint32_t>& span) {
    keep.clear();
    for (int a = 0; a < R2; ++a) if (!((dropmask >> a) & 1)) keep.push_back(S.cls[a]);
    span.assign(1, 0u);
    for (uint32_t k : keep) { size_t n = span.size();
        for (size_t i = 0; i < n; ++i) span.push_back(span[i] ^ k); }
    std::sort(span.begin(), span.end());
}
static int64_t span_sum_list(const std::vector<uint32_t>& span) {
    int64_t s = 0; for (uint32_t v : span) if (v) s += COV1[v]; return s;
}

static void report_new_best(const Sol& prev, const Sol& cur, const Verify& V, double t);

// try to install a candidate solution; returns true if it became the new incumbent
// Accept a move into the working point CUR; promote to the global incumbent BEST only after
// the candidate has passed a from-scratch re-verification over all 91770 errors.
static bool try_move(const uint32_t* cls) {
    if (!independent_isotropic(cls)) return false;
    Sol S = make_sol(cls);
    if (!better(S, CUR)) return false;
    bool isbest;
    { std::lock_guard<std::mutex> lk(BESTMTX); isbest = better(S, BEST); }
    if (isbest) {
        Verify V = verify(S);
        if (!V.wt_ok || !V.cg_ok || !V.ch_ok || !V.rank_ok || !V.u3_ok || !V.score_ok) {
            printf("*** candidate failed independent verification, rejected ***\n");
            return false;
        }
        std::lock_guard<std::mutex> lk(BESTMTX);
        if (better(S, BEST)) { Sol prev = BEST; BEST = S; CUR = S;
            report_new_best(prev, S, V, elapsed()); return true; }
    }
    CUR = S;
    return true;
}

// ---------------------------------------------------------------- 1-opt (exhaustive)
static bool opt1(int nthreads) {
    bool improved = false;
    for (;;) {
        bool any = false;
        for (int drop = 0; drop < R2; ++drop) {
            std::vector<uint32_t> keep, span;
            kept_span(CUR, uint8_t(1u << drop), keep, span);
            build_ftab(FTAB, keep.data(), (int)keep.size());
            const int64_t base = span_sum_list(span);
            const uint8_t need_bits = uint8_t(((1u << R2) - 1u) & ~(1u << drop));
            uint32_t bestv = 0; int64_t bests = CUR.span_sum;
            std::vector<uint32_t> tie;
            for (size_t i = 0; i < CAND.size(); ++i) {
                if ((OMASK[i] & need_bits) != need_bits) continue;
                const uint32_t c = CAND[i];
                if (std::binary_search(span.begin(), span.end(), c)) continue;
                const int64_t s = base + (int64_t)FTAB[c];
                if (s > bests) { bests = s; bestv = c; tie.clear(); tie.push_back(c); }
                else if (s == bests && tie.size() < 4096) tie.push_back(c);
            }
            if (bestv || !tie.empty()) {
                uint32_t trial[R2];
                for (uint32_t c : tie) {
                    std::copy(CUR.cls, CUR.cls + R2, trial); trial[drop] = c;
                    if (try_move(trial)) { any = true; improved = true; break; }
                }
            }
            if (any) break;
        }
        if (!any) break;
        rebuild_omask(CUR, nthreads);
    }
    return improved;
}

// ---------------------------------------------------------------- 2-opt (exhaustive)
static bool opt2(int nthreads, bool verbose) {
    bool improved = false;
    for (;;) {
        bool any = false;
        for (int d1 = 0; d1 < R2 && !any; ++d1) for (int d2 = d1 + 1; d2 < R2 && !any; ++d2) {
            const uint8_t dm = uint8_t((1u << d1) | (1u << d2));
            std::vector<uint32_t> keep, span;
            kept_span(CUR, dm, keep, span);
            build_ftab(FTAB, keep.data(), (int)keep.size());
            const uint32_t fmax = table_max(FTAB);
            const int64_t base = span_sum_list(span);
            const uint8_t need_bits = uint8_t(((1u << R2) - 1u) & ~dm);
            std::vector<uint32_t> cls, val;
            cls.reserve(CAND.size() >> 5);
            for (size_t i = 0; i < CAND.size(); ++i) {
                if ((OMASK[i] & need_bits) != need_bits) continue;
                const uint32_t c = CAND[i];
                if (std::binary_search(span.begin(), span.end(), c)) continue;
                cls.push_back(c);
            }
            std::sort(cls.begin(), cls.end(), [](uint32_t a, uint32_t b) {
                if (FTAB[a] != FTAB[b]) return FTAB[a] > FTAB[b]; return a < b; });
            val.resize(cls.size());
            for (size_t i = 0; i < cls.size(); ++i) val[i] = FTAB[cls[i]];
            std::vector<Hit> hits;
            search_pairs(cls, val, FTAB, base, CUR.span_sum - 1, fmax, span, nthreads, hits, 20000);
            if (verbose)
                printf("    2-opt drop(%d,%d): %zu candidates, %zu hits\n",
                       d1, d2, cls.size(), hits.size());
            std::sort(hits.begin(), hits.end(), [](const Hit& p, const Hit& q) { return p.score > q.score; });
            uint32_t trial[R2];
            for (const Hit& h : hits) {
                std::copy(CUR.cls, CUR.cls + R2, trial);
                trial[d1] = h.a; trial[d2] = h.b;
                if (try_move(trial)) { any = true; improved = true; break; }
            }
        }
        if (!any) break;
        rebuild_omask(CUR, nthreads);
    }
    return improved;
}

// ---------------------------------------------------------------- 3-opt
//  For a fixed first replacement a, g(w) = f(w) + f(w^a) turns the remaining problem into
//  exactly the same three-term pair search:  S = base + f(a) + g(b) + g(c) + g(b^c).
//  `target` (a residual class, or 0) restricts the first replacement to generators that
//  actually anticommute with that residual error -- an exact way to aim the search.
static bool opt3(int nthreads, int topM, uint32_t target, double deadline, bool verbose) {
    bool improved = false;
    for (int d1 = 0; d1 < R2; ++d1) for (int d2 = d1 + 1; d2 < R2; ++d2) for (int d3 = d2 + 1; d3 < R2; ++d3) {
        if (elapsed() > deadline) return improved;
        const uint8_t dm = uint8_t((1u << d1) | (1u << d2) | (1u << d3));
        std::vector<uint32_t> keep, span;
        kept_span(CUR, dm, keep, span);
        build_ftab(FTAB, keep.data(), (int)keep.size());
        const int64_t base = span_sum_list(span);
        const uint8_t need_bits = uint8_t(((1u << R2) - 1u) & ~dm);
        std::vector<uint32_t> pool;
        for (size_t i = 0; i < CAND.size(); ++i) {
            if ((OMASK[i] & need_bits) != need_bits) continue;
            const uint32_t c = CAND[i];
            if (std::binary_search(span.begin(), span.end(), c)) continue;
            if (target && !csymp(c, target)) continue;         // must hit the targeted residual
            pool.push_back(c);
        }
        std::sort(pool.begin(), pool.end(), [](uint32_t a, uint32_t b) {
            if (FTAB[a] != FTAB[b]) return FTAB[a] > FTAB[b]; return a < b; });
        const int M = std::min<int>(topM, (int)pool.size());
        for (int ai = 0; ai < M; ++ai) {
            if (elapsed() > deadline) return improved;
            const uint32_t a = pool[ai];
            GTAB.resize(VSZ);
            for (uint32_t w = 0; w < (uint32_t)VSZ; ++w) GTAB[w] = FTAB[w] + FTAB[w ^ a];
            const uint32_t gmax = table_max(GTAB);
            std::vector<uint32_t> span2(span);
            { size_t n = span2.size(); for (size_t i = 0; i < n; ++i) span2.push_back(span2[i] ^ a);
              std::sort(span2.begin(), span2.end()); }
            std::vector<uint32_t> cls, val;
            for (uint32_t c : pool) {
                if (c == a || csymp(c, a)) continue;
                if (std::binary_search(span2.begin(), span2.end(), c)) continue;
                cls.push_back(c);
            }
            std::sort(cls.begin(), cls.end(), [](uint32_t p, uint32_t q) {
                if (GTAB[p] != GTAB[q]) return GTAB[p] > GTAB[q]; return p < q; });
            val.resize(cls.size());
            for (size_t i = 0; i < cls.size(); ++i) val[i] = GTAB[cls[i]];
            std::vector<Hit> hits;
            search_pairs(cls, val, GTAB, base + (int64_t)FTAB[a], CUR.span_sum - 1, gmax,
                         span2, nthreads, hits, 20000);
            if (hits.empty()) continue;
            std::sort(hits.begin(), hits.end(), [](const Hit& p, const Hit& q) { return p.score > q.score; });
            uint32_t trial[R2];
            for (const Hit& h : hits) {
                std::copy(CUR.cls, CUR.cls + R2, trial);
                trial[d1] = a; trial[d2] = h.a; trial[d3] = h.b;
                if (try_move(trial)) { improved = true; break; }
            }
            if (improved) { rebuild_omask(CUR, nthreads);
                if (verbose) printf("    3-opt improved at drop(%d,%d,%d)\n", d1, d2, d3);
                return improved; }
        }
    }
    return improved;
}

// ===========================================================================================
//  OUTPUT
// ===========================================================================================
static void residual_block(const ResidualInfo& r) {
    ap("Remaining undetected errors: %lld\n\n", (long long)r.n);
    for (int w = 1; w <= MAXW; ++w) {
        ap("Weight %d:\n", w); int n = 0;
        for (uint32_t E : r.list) if (pwt(E) == w) { ap("  %s = %s\n", pstr(E).c_str(), supp(E).c_str()); ++n; }
        if (!n) ap("  (none)\n");
        ap("\n");
    }
    ap("Residual weight distribution: %lld / %lld / %lld / %lld  (weights 1/2/3/4)\n",
       (long long)r.by_w[1], (long long)r.by_w[2], (long long)r.by_w[3], (long long)r.by_w[4]);
    ap("Residual set forms a subgroup (with the identity): %s\n", r.subgroup ? "YES" : "NO");
    ap("Residual span dimension: %d\n", r.dim);
    ap("Generators:\n");
    for (uint32_t g : r.gens) ap("  %s = %s\n", pstr(g).c_str(), supp(g).c_str());
    ap("All residual errors are X-only: %s\n", r.all_x ? "YES" : "NO");
    ap("Physical qubits involved:");
    for (uint32_t q : r.qubits) ap(" %u", q);
    ap("\n");
}

static void write_sol_dir(const std::string& dir, const Sol& S, const Verify& V) {
    std::error_code ec; std::filesystem::create_directories(dir, ec);
    ResidualInfo r = analyse(V.leftover);
    uint32_t all[R1 + R2];
    for (int i = 0; i < R1; ++i) all[i] = G3[i];
    for (int a = 0; a < R2; ++a) all[R1 + a] = S.pau[a];

    RB.clear();
    ap("Stage-1 stabilizer (frozen, orbit 0006 of the proven-optimal family)\n\n");
    for (int i = 0; i < R1; ++i) ap("g%d = %s   (weight %d)\n", i + 1, pstr(G3[i]).c_str(), pwt(G3[i]));
    ap("\nM_3 = %lld / %lld\n|U_3| = %zu   (weights 1/2/3/4 = %lld/%lld/%lld/%lld)\n",
       (long long)V.orig_det, (long long)V.all_total, U3.size(),
       (long long)U3_BY_W[1], (long long)U3_BY_W[2], (long long)U3_BY_W[3], (long long)U3_BY_W[4]);
    dump(dir + "/stage1.txt");

    RB.clear();
    ap("Stage-2: eight weight-10 generators\n\n");
    for (int a = 0; a < R2; ++a) ap("h%d = %s   (weight %d)\n", a + 1, pstr(S.pau[a]).c_str(), V.weights[a]);
    ap("\nM_8 = %lld / %zu\nresidual = %lld\n", (long long)S.covered, U3.size(), (long long)S.residual);
    ap("M_final = %lld / %lld  (%.8f %%)\n", (long long)V.det11, (long long)V.all_total,
       100.0 * double(V.det11) / double(V.all_total));
    ap("\nTarget-error detection by each new generator:\n");
    for (int a = 0; a < R2; ++a)
        ap("h%d: %lld   (cov1 = %d %s)\n", a + 1, (long long)V.per_gen[a], (int)COV1[S.cls[a]],
           V.per_gen[a] == (int64_t)COV1[S.cls[a]] ? "OK" : "*** MISMATCH ***");
    ap("\nDetected Stage-1 targets by weight:\n");
    for (int w = 1; w <= MAXW; ++w)
        ap("  weight %d : %lld / %lld\n", w, (long long)V.cov_by_w[w], (long long)V.tgt_by_w[w]);
    ap("\nWeight-4 errors (total recomputed, not hard-coded): %lld\n", (long long)V.w4_total);
    ap("  detected   : %lld\n  undetected : %lld\n", (long long)V.w4_det, (long long)V.w4_und);
    ap("  ALL WEIGHT-4 ERRORS DETECTED: %s\n", V.w4_und == 0 ? "YES" : "NO");
    ap("\nCommutation (rows/cols g1..g3,h1..h8; 0 = commute):\n");
    for (int i = 0; i < R1 + R2; ++i) { ap("  ");
        for (int j = 0; j < R1 + R2; ++j) ap("%d ", V.comm[i][j]); ap("\n"); }
    ap("\nINDEPENDENT VERIFICATION (all 91770 errors rebuilt from scratch)\n");
    ap("  U_3 size                 : %lld %s\n", (long long)V.u3size, V.u3_ok ? "OK" : "MISMATCH");
    ap("  all weights == 10        : %s\n", V.wt_ok ? "OK" : "FAIL");
    ap("  [h_a , g_i] = 0          : %s\n", V.cg_ok ? "OK" : "FAIL");
    ap("  [h_a , h_b] = 0          : %s\n", V.ch_ok ? "OK" : "FAIL");
    ap("  rank = 11                : %s (rank %d)\n", V.rank_ok ? "OK" : "FAIL", V.rank);
    ap("  union coverage recount   : %lld %s\n", (long long)V.covered,
       V.score_ok ? "(matches the algebraic score)" : "*** MISMATCH ***");
    dump(dir + "/stage2_best.txt");

    RB.clear(); residual_block(r); dump(dir + "/remaining_errors.txt");

    RB.clear();
    ap("Complete 11-generator stabilizer  [[14,3]]\n\n");
    for (int i = 0; i < R1; ++i) ap("g%d = %s\n", i + 1, pstr(G3[i]).c_str());
    for (int a = 0; a < R2; ++a) ap("h%d = %s\n", a + 1, pstr(S.pau[a]).c_str());
    ap("\nParity-check matrix  H = [ X | Z ]   (11 x 28):\n");
    for (int i = 0; i < R1 + R2; ++i)
        ap("  %s | %s\n", bstr(px(all[i])).c_str(), bstr(pz(all[i])).c_str());
    ap("\nrank = %d,  weights =", V.rank);
    for (int i = 0; i < R1 + R2; ++i) ap(" %d", pwt(all[i]));
    ap("\n");
    dump(dir + "/final_matrix.txt");

    RB.clear();
    ap("Fixed Stage-1 orbit: 0006\n\n");
    ap("Stage-1:\n  M3 = %lld / %lld\n  U3 = %zu\n\n",
       (long long)V.orig_det, (long long)V.all_total, U3.size());
    ap("Stage-2 residual: %lld\n", (long long)S.residual);
    ap("M8:      %lld / %zu\n", (long long)S.covered, U3.size());
    ap("M_final: %lld / %lld  (%.8f %%)\n\n", (long long)V.det11, (long long)V.all_total,
       100.0 * double(V.det11) / double(V.all_total));
    ap("Residual weight distribution:\n  w1 = %lld\n  w2 = %lld\n  w3 = %lld\n  w4 = %lld\n\n",
       (long long)S.res_w[1], (long long)S.res_w[2], (long long)S.res_w[3], (long long)S.res_w[4]);
    ap("All weight-4 errors detected: %s\n\n", V.w4_und == 0 ? "YES" : "NO");
    residual_block(r);
    ap("\nStage-1 generators:\n");
    for (int i = 0; i < R1; ++i) ap("g%d = %s\n", i + 1, pstr(G3[i]).c_str());
    ap("\nStage-2 generators:\n");
    for (int a = 0; a < R2; ++a) ap("h%d = %s\n", a + 1, pstr(S.pau[a]).c_str());
    dump(dir + "/summary.txt");
}

static Sol START;                       // the incumbent we started from
static void report_new_best(const Sol& prev, const Sol& cur, const Verify& V, double t) {
    ++IMPROVEMENTS;
    ResidualInfo r = analyse(V.leftover);
    RB.clear();
    ap("============================================================\n");
    ap("NEW BEST\n");
    ap("============================================================\n\n");
    ap("Previous residual:\n%lld\n\n", (long long)prev.residual);
    ap("New residual:\n%lld\n\n", (long long)cur.residual);
    ap("M8:\n%lld / %zu\n\n", (long long)cur.covered, U3.size());
    ap("M_final:\n%lld / %lld\n\n", (long long)V.det11, (long long)V.all_total);
    ap("Residual weight distribution:\nw1 = %lld  w2 = %lld  w3 = %lld  w4 = %lld\n\n",
       (long long)cur.res_w[1], (long long)cur.res_w[2], (long long)cur.res_w[3], (long long)cur.res_w[4]);
    ap("Weight-4 residual:\n%lld   (all weight-4 detected: %s)\n\n",
       (long long)V.w4_und, V.w4_und == 0 ? "YES" : "NO");
    ap("Generators:\n");
    for (int a = 0; a < R2; ++a) ap("h%d = %s\n", a + 1, pstr(cur.pau[a]).c_str());
    ap("\n");
    residual_block(r);
    ap("\nTime:\n%.2f s\n", t);
    ap("============================================================\n");
    fputs(RB.c_str(), stdout); fflush(stdout);
    char dir[64]; snprintf(dir, sizeof(dir), "%s/improved_%04d", OUTROOT.c_str(), IMPROVEMENTS);
    std::string d(dir);
    { std::string keep = RB; std::error_code ec; std::filesystem::create_directories(d, ec);
      RB = keep; dump(d + "/new_best_report.txt"); }
    write_sol_dir(d, cur, V);
    write_sol_dir(OUTROOT + "/best_known", cur, V);      // only after full verification
    if (cur.residual == 0) {
        RB.clear();
        ap("============================================================\n");
        ap("FULL COVERAGE FOUND\n");
        ap("============================================================\n\n");
        ap("All %zu Stage-1-undetected errors are detected.\n\n", U3.size());
        ap("Final detection:\n%lld / %lld\n\n", (long long)V.det11, (long long)V.all_total);
        ap("All Pauli errors of weight <= 4 are detected.\n\n");
        ap("============================================================\n");
        dump(OUTROOT + "/FULL_COVERAGE.txt", true);
    }
    if (V.w4_und == 0) {
        printf("============================================================\n"
               "ALL WEIGHT-4 ERRORS DETECTED   %lld / %lld\n"
               "============================================================\n",
               (long long)V.w4_det, (long long)V.w4_total);
    }
}

// ===========================================================================================
//  ITERATED LOCAL SEARCH DRIVER
//  exhaustive 1-opt -> exhaustive 2-opt -> targeted/general 3-opt -> kick, repeat.
//  The kick deliberately accepts a worse working point (CUR) so the next local search starts
//  in a different basin; BEST is only ever replaced by a fully verified improvement.
// ===========================================================================================
// Ruin-and-recreate kick.  A uniformly random replacement is far too destructive here (a
// random weight-10 class has cov1 ~ 5593 against ~5790 for a good one, and the local search
// never climbs back), so the dropped generators are rebuilt greedily-randomly: each one is
// drawn uniformly from the topR classes by coset gain against what is already kept.  That
// lands in a different basin of comparable quality instead of a much worse one.
static bool kick(std::mt19937_64& rng, int k, int topR, int nthreads) {
    int idx[R2]; std::iota(idx, idx + R2, 0); std::shuffle(idx, idx + R2, rng);
    uint32_t trial[R2]; std::copy(CUR.cls, CUR.cls + R2, trial);
    std::vector<uint32_t> kb;
    for (int a = 0; a < R2; ++a) {
        bool dropped = false;
        for (int t = 0; t < k; ++t) if (idx[t] == a) dropped = true;
        if (!dropped) kb.push_back(CUR.cls[a]);
    }
    for (int t = 0; t < k; ++t) {
        build_ftab(FTAB, kb.data(), (int)kb.size());
        std::vector<uint32_t> span(1, 0u);
        for (uint32_t v : kb) { size_t n = span.size();
            for (size_t i = 0; i < n; ++i) span.push_back(span[i] ^ v); }
        std::sort(span.begin(), span.end());
        std::vector<uint32_t> pool;
        for (uint32_t c : CAND) {
            bool ok = true;
            for (uint32_t v : kb) if (csymp(c, v)) { ok = false; break; }
            if (!ok || std::binary_search(span.begin(), span.end(), c)) continue;
            pool.push_back(c);
        }
        if (pool.empty()) return false;
        const int R = std::min<int>(topR, (int)pool.size());
        std::partial_sort(pool.begin(), pool.begin() + R, pool.end(),
            [](uint32_t a, uint32_t b) { return FTAB[a] > FTAB[b]; });
        const uint32_t pick = pool[rng() % (size_t)R];
        trial[idx[t]] = pick; kb.push_back(pick);
    }
    if (!independent_isotropic(trial)) return false;
    CUR = make_sol(trial);
    rebuild_omask(CUR, nthreads);
    return true;
}

// ===========================================================================================
//  MAIN
// ===========================================================================================
static void usage() {
    printf("stab14_refine -- targeted refinement around the known 3+8 incumbent (orbit 0006)\n"
           "  --seed S            RNG seed for the kicks (default 1)\n"
           "  --time-limit S      wall-clock budget in seconds (default 3600)\n"
           "  --threads N         default: all cores\n"
           "  --deterministic     no kicks; exhaustive 1-opt / 2-opt / 3-opt only\n"
           "  --top-m N           first-replacement breadth of the 3-opt phase (default 48)\n"
           "  --no-3opt           stop after the exhaustive 1-opt and 2-opt certification\n"
           "  --out DIR           output root (default refine)\n"
           "  --verbose\n");
}

int main(int argc, char** argv) {
    int nthreads = (int)std::max(1u, std::thread::hardware_concurrency());
    uint64_t seed = 1;
    double tlimit = 3600.0;
    int topM = 48;
    bool deterministic = false, no3 = false, verbose = false, useKick = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto nxt = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : std::string(); };
        if      (a == "--seed")          seed = strtoull(nxt().c_str(), nullptr, 10);
        else if (a == "--time-limit")    tlimit = atof(nxt().c_str());
        else if (a == "--threads")       nthreads = std::max(1, atoi(nxt().c_str()));
        else if (a == "--deterministic") deterministic = true;
        else if (a == "--top-m")         topM = std::max(1, atoi(nxt().c_str()));
        else if (a == "--no-3opt")       no3 = true;
        else if (a == "--kick")          useKick = true;
        else if (a == "--out")           OUTROOT = nxt();
        else if (a == "--verbose")       verbose = true;
        else { usage(); return (a == "--help" || a == "-h") ? 0 : 1; }
    }
    T0 = std::chrono::steady_clock::now();
    std::error_code ec; std::filesystem::create_directories(OUTROOT, ec);

    for (int i = 0; i < R1; ++i) G3[i] = parse_pauli(G3_STR[i]);
    for (int a = 0; a < R2; ++a) H0[a] = parse_pauli(H8_STR[a]);

    printf("=== targeted refinement around the known 3+8 solution (Stage-1 orbit 0006) ===\n");
    for (int i = 0; i < R1; ++i) {
        printf("g%d = %s  weight %d\n", i + 1, pstr(G3[i]).c_str(), pwt(G3[i]));
        if (pwt(G3[i]) != W1) { printf("FATAL: Stage-1 weight != %d\n", W1); return 2; }
    }
    for (int i = 0; i < R1; ++i) for (int j = i + 1; j < R1; ++j)
        if (symp(G3[i], G3[j])) { printf("FATAL: Stage-1 generators do not commute\n"); return 2; }

    build_targets();
    printf("\n|U_3| = %zu   (weights 1/2/3/4 = %lld/%lld/%lld/%lld),  total errors %lld,"
           "  weight-4 total %lld\n", U3.size(), (long long)U3_BY_W[1], (long long)U3_BY_W[2],
           (long long)U3_BY_W[3], (long long)U3_BY_W[4], (long long)TOT_ERR, (long long)TOTW4);
    if (U3.size() != 11186) { printf("FATAL: |U_3| != 11186\n"); return 3; }

    build_basis(); build_cov1(); build_repr10();
    printf("V = C(S3)/S3 dim %d (%d classes);  %lld weight-10 Paulis reach %lld classes;"
           "  max cov1 = %d\n", VDIM, VSZ, (long long)NW10, (long long)NCLS10, COVMAX);

    uint32_t cls0[R2];
    for (int a = 0; a < R2; ++a) {
        cls0[a] = cls_of(H0[a]);
        PREF.push_back({cls0[a], H0[a]});
        if (pwt(H0[a]) != W2) { printf("FATAL: h%d does not have weight %d\n", a + 1, W2); return 4; }
        for (int i = 0; i < R1; ++i) if (symp(H0[a], G3[i])) {
            printf("FATAL: h%d does not commute with g%d\n", a + 1, i + 1); return 4; }
    }
    if (!independent_isotropic(cls0)) { printf("FATAL: incumbent is not isotropic/independent\n"); return 4; }

    Sol S0 = make_sol(cls0);
    Verify V0 = verify(S0);
    printf("\nincumbent: M_8 = %lld / %zu, residual %lld (w %lld/%lld/%lld/%lld),"
           " M_final = %lld / %lld\n", (long long)S0.covered, U3.size(), (long long)S0.residual,
           (long long)S0.res_w[1], (long long)S0.res_w[2], (long long)S0.res_w[3],
           (long long)S0.res_w[4], (long long)V0.det11, (long long)V0.all_total);
    printf("verification: weights %s, [h,g] %s, [h,h] %s, rank %d %s, brute-force coverage %lld %s\n",
           V0.wt_ok ? "OK" : "FAIL", V0.cg_ok ? "OK" : "FAIL", V0.ch_ok ? "OK" : "FAIL",
           V0.rank, V0.rank_ok ? "OK" : "FAIL", (long long)V0.covered,
           V0.score_ok ? "OK (matches the algebraic identity)" : "*** MISMATCH ***");
    if (!V0.wt_ok || !V0.cg_ok || !V0.ch_ok || !V0.rank_ok || !V0.score_ok) {
        printf("FATAL: the supplied incumbent does not verify\n"); return 5; }
    printf("weight-4 errors: %lld detected, %lld undetected (total %lld)\n",
           (long long)V0.w4_det, (long long)V0.w4_und, (long long)V0.w4_total);

    START = S0; BEST = S0; CUR = S0;
    write_sol_dir(OUTROOT + "/best_known", S0, V0);         // preserved before anything else
    rebuild_omask(CUR, nthreads);

    // ---------------------------------------------------------------- exhaustive phases
    printf("\n[phase 1] exhaustive 1-opt (every replacement of one generator by any of the"
           " %zu weight-10 classes) ...\n", CAND.size());
    fflush(stdout);
    bool imp1 = opt1(nthreads);
    printf("[phase 1] %s   residual now %lld   (%.1f s)\n",
           imp1 ? "IMPROVED" : "no improvement -- incumbent is 1-opt optimal",
           (long long)BEST.residual, elapsed());
    fflush(stdout);

    printf("\n[phase 2] exhaustive 2-opt over all %d pairs of dropped generators ...\n", 28);
    fflush(stdout);
    bool imp2 = opt2(nthreads, verbose);
    printf("[phase 2] %s   residual now %lld   (%.1f s)\n",
           imp2 ? "IMPROVED" : "no improvement -- incumbent is 1-opt AND 2-opt optimal",
           (long long)BEST.residual, elapsed());
    fflush(stdout);

    bool anyimp = imp1 || imp2;
    if (!no3) {
        printf("\n[phase 3] 3-opt: 56 drop-triples, first replacement over the top-%d classes"
               " by coset gain, last two exhaustive.  Residual-targeted rounds alternate with"
               " untargeted ones.\n", topM);
        fflush(stdout);
        // ANCHORED ITERATIVE BROADENING.  Kicking away from a 1-opt/2-opt-certified incumbent
        // only lands in worse basins, so instead every sweep restarts from BEST and the only
        // thing that grows is the breadth M of the first replacement.  One pass cycles the
        // residual target (untargeted, then one sweep aimed at each residual error in turn);
        // when the whole cycle fails, M doubles.  M = |pool| would be exhaustive 3-opt, so
        // this is a monotone ladder towards it rather than a random walk.
        std::mt19937_64 rng(seed);
        int M = topM, cyc = 0;
        while (elapsed() < tlimit) {
            CUR = BEST; rebuild_omask(CUR, nthreads);
            std::vector<uint32_t> res; residual_of(BEST.cls, res);
            std::vector<uint32_t> targets; targets.push_back(0);
            for (uint32_t E : res) targets.push_back(cls_of(E));
            const uint32_t tgt = targets[(size_t)cyc % targets.size()];
            if (verbose)
                printf("  sweep: breadth M=%d, target %s, residual %lld, %.0f s left\n",
                       M, tgt ? supp(res[(cyc % targets.size()) - 1]).c_str() : "none",
                       (long long)BEST.residual, tlimit - elapsed());
            fflush(stdout);
            bool imp = opt3(nthreads, M, tgt, tlimit, verbose);
            if (imp) { anyimp = true; opt1(nthreads); opt2(nthreads, false);
                       M = topM; cyc = 0; continue; }
            ++cyc;
            if ((size_t)cyc >= targets.size()) {                 // whole cycle exhausted
                cyc = 0;
                if (useKick) { if (!kick(rng, 2 + int(rng() % 2), 64, nthreads)) break;
                               opt1(nthreads); opt2(nthreads, false);
                               if (better(CUR, BEST)) { anyimp = true; } }
                else M = std::min<int>(M * 2, 1 << 20);
                if (deterministic && !useKick && M >= (1 << 20)) break;
            }
        }
    }

    // ---------------------------------------------------------------- final report
    Verify VB = verify(BEST);
    ResidualInfo rb = analyse(VB.leftover);
    RB.clear();
    ap("============================================================\n");
    ap("TARGETED REFINEMENT RESULT\n");
    ap("============================================================\n\n");
    ap("Fixed Stage-1 orbit:\n0006\n\n");
    ap("Stage-1:\nM3 = %lld / %lld\nU3 = %zu\n\n",
       (long long)VB.orig_det, (long long)VB.all_total, U3.size());
    ap("Starting Stage-2 residual:\n%lld\n\n", (long long)START.residual);
    ap("Best Stage-2 residual:\n%lld\n\n", (long long)BEST.residual);
    ap("Improvement:\n%lld\n\n", (long long)(START.residual - BEST.residual));
    ap("M8:\n%lld / %zu\n\n", (long long)BEST.covered, U3.size());
    ap("M_final:\n%lld / %lld\n\n", (long long)VB.det11, (long long)VB.all_total);
    ap("Residual weight distribution:\nw1 = %lld\nw2 = %lld\nw3 = %lld\nw4 = %lld\n\n",
       (long long)BEST.res_w[1], (long long)BEST.res_w[2], (long long)BEST.res_w[3],
       (long long)BEST.res_w[4]);
    ap("All weight-4 errors detected:\n%s   (%lld / %lld)\n\n",
       VB.w4_und == 0 ? "YES" : "NO", (long long)VB.w4_det, (long long)VB.w4_total);
    residual_block(rb);
    ap("\nStage-1 generators:\n");
    for (int i = 0; i < R1; ++i) ap("g%d = %s\n", i + 1, pstr(G3[i]).c_str());
    ap("\nStage-2 generators:\n");
    for (int a = 0; a < R2; ++a) ap("h%d = %s\n", a + 1, pstr(BEST.pau[a]).c_str());
    ap("\nSearch time:\n%.2f s\n\nSeeds:\n%llu\n\n", elapsed(), (unsigned long long)seed);
    ap("Neighbourhoods searched EXHAUSTIVELY (over all %zu weight-10 classes):\n"
       "  1-opt : every replacement of one of the eight generators\n"
       "  2-opt : every simultaneous replacement of two of them\n", CAND.size());
    ap("  3-opt : all 56 drop-triples, first replacement over the top-%d classes by coset\n"
       "          gain, the remaining two exhaustive%s\n", topM, no3 ? " (SKIPPED)" : "");
    ap("\nStatus:\n%s\n", BEST.residual == 0 ? "FULL COVERAGE"
                        : (anyimp ? "IMPROVED" : "BEST KNOWN"));
    if (!anyimp) {
        ap("\nNo improvement over the incumbent was found.\n"
           "The %lld-residual solution remains BEST KNOWN.\n"
           "No impossibility claim is made.\n", (long long)START.residual);
        ap("What IS established: the incumbent is optimal against every 1-generator and every\n"
           "2-generator replacement, exhaustively over all %zu weight-10 classes.\n", CAND.size());
    }
    ap("============================================================\n");
    dump(OUTROOT + "/refinement_report.txt", true);
    write_sol_dir(OUTROOT + "/best_known", BEST, VB);
    return 0;
}
