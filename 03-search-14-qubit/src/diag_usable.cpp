// ===========================================================================================
//  stab14_kopt.cpp -- exhaustive 3-opt and targeted 4-opt around the best-known 3+8 solutions.
//
//  NEW DIRECTORY, NEW SOURCE.  Nothing under stabilizer_14q/ outside targeted_3opt4opt/ is
//  written; the seven orbit incumbents below were READ from stage1_orbit_*/summary.txt.
//
//  =========================================================================================
//  THE REFORMULATION THAT MAKES EXHAUSTIVE 3-OPT POSSIBLE
//  =========================================================================================
//  With S3 frozen, the eight Stage-2 generators span an 8-dimensional totally isotropic
//  D <= V = C(S3)/S3 (dim 22), and
//        M_8(D) = (1/128) SUM_{v in D\{0}} cov1(v),   cov1(v) = #{E in U_3 : <v,E> = 1}.
//
//  A k-opt move keeps K = span of the 8-k retained generators (isotropic, dim 8-k) and asks
//  for an 8-dimensional isotropic D' with K <= D' <= K^perp.  Those D' are in bijection with
//  the k-dimensional TOTALLY ISOTROPIC SUBSPACES OF THE QUOTIENT
//
//        Q = K^perp / K ,    dim Q = 22 - 2(8-k) = 6 + 2k ,
//
//  which carries a non-degenerate symplectic form.  So
//        k = 3  ->  dim Q = 12,  |Q| =  4096   (score table 16 KB: L1-resident)
//        k = 4  ->  dim Q = 14,  |Q| = 16384   (64 KB)
//  instead of the 2^22 = 4.2M-element, 16 MB tables the previous refiner probed.  f_K is
//  constant on K-cosets, so it descends to phi: Q -> Z, phi(q) = SUM_{u in K} cov1(lift(q)^u),
//  and for a k-dimensional D'/K spanned by q_1..q_k
//
//        S(D') = S(K) + SUM over the 2^k - 1 non-zero q in D'/K of phi(q).
//
//  Consequences:
//    * the ENTIRE 3-opt neighbourhood of one drop-triple is a search for 3-dimensional
//      isotropic subspaces of a 12-dimensional symplectic space -- about 5.1e7 of them --
//      so all 56 drop-triples can be swept EXHAUSTIVELY.  No breadth limit is needed, and
//      the requested M = 256 is strictly subsumed.
//    * fixing q1 and putting psi(w) = phi(w) + phi(w^q1) turns the rest into the three-term
//      pair form  S = S(K) + phi(q1) + psi(q2) + psi(q3) + psi(q2^q3), so one table rebuild
//      per q1 (4096 adds) serves the whole inner double loop.
//    * q1 is forced to be the MINIMUM of the 7 non-zero span elements, which removes most of
//      the 168 bases per subspace without conflicting with the psi-sorted pruning order.
//
//  Weight-10 realisability: a generator is a class v in A10 (the classes reachable by a
//  weight-10 Pauli of C(S3)).  Replacing generators by q in Q is legitimate iff the K-coset
//  of q contains an A10 class; the program records one such representative per q.  Because
//  ~48.6% of all classes are weight-10 realisable and a coset holds 2^(8-k) >= 16 of them,
//  essentially every q is usable -- but the test is done exactly, never assumed.
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
constexpr int R1    = 3;
constexpr int R2    = 8;
constexpr int MAXW  = 4;
constexpr int W1    = 8, W2 = 10;
constexpr uint32_t QM   = (1u << NQ) - 1u;
constexpr uint32_t NONE = 0xFFFFFFFFu;

constexpr int VDIM  = 2 * NQ - 2 * R1;        // 22
constexpr int VHALF = VDIM / 2;               // 11
constexpr int CDIM  = 2 * NQ - R1;            // 25
constexpr uint32_t VLO = (1u << VHALF) - 1u;
constexpr int VSZ   = 1 << VDIM;

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
            default: fprintf(stderr, "bad Pauli char '%c'\n", s[j]); exit(1);
        }
    }
    return pk(x, z);
}

// ---------------------------------------------------------------- the seven orbit incumbents
struct OrbitDef { const char* name; const char* g[R1]; const char* h[R2]; int start_m8; };
static const OrbitDef ORBITS[] = {
{"0006",{"IIIIXXXXIIZZXX","IIIIXXIIXXXXZZ","IIIIIIZZZZZZZZ"},
        {"IXIXXXIZYZIXXX","IZXIXXYIXYIYXX","IIIZYXZXYXZIZY","IIIYXZYIYZZYZX",
         "IIZIXZXYZXZYIZ","XIIIXYXZXYIXYX","IIYIZXZXYXXIYX","ZIIIZXIZZYZXYX"},11181},
{"0002",{"IIIIXIXIYYXXZX","IIIIIXIXXXXXXZ","IIIIIIZZZZZZZZ"},
        {"ZXIIIIYZYXXZXX","IIIXZXXZZYYZIX","IIIXIYYYZXXYZX","IIXIZXXYZXZXZI",
         "IIXIZXZXXYZIZX","XIIIIYXZXYZYXX","YIIIXZXZXZZXXI","XIIIZXYXYZXYIX"},11179},
{"0012",{"IIIXXIXIIYXXXZ","IIIIIXIXXXXXXX","IIIIIIXXZZZZZZ"},
        {"XIIIZYXIYZXZXX","IXIIXIYXYYYZXX","IIIIXYZZXZYYXZ","XIIIXXXZZXZYIX",
         "IZIXIYXXYZYXXI","IYIXYXXXIXXZZI","IZIIXXYXZYYXIZ","IIXXIZXXYYXZIX"},11179},
{"0023",{"IIIXIXXIXIXXXX","IIIIXIIXZXXXXZ","IIIIIXXXIZZZZZ"},
        {"XIIYIZIXXXZZXX","IXIIZIXXXYYZXX","IIXZXIXXXZIZXX","IZIXXXIZZZIYZZ",
         "XIIXIYXXZYYIZX","IYIIXXIZXZZYXX","IIIYYXIYYXZYYZ","IXIXXYXXIYXZXI"},11179},
{"0031",{"IIXXIIIXIXYYXX","IIIIXXIIXXXXXX","IIIIIIXXXXZZZZ"},
        {"IIXYIXZYXIYYZZ","IIXYXXIYZYXIYY","IIYXXIZXXYYIYZ","XXIZIXIXXXZYXI",
         "ZZYIIXIXXXYZIX","IIZXIXXXYYYXZI","IIZZIXXXYXXIXY","IIYYXZXXXIIXZX"},11179},
{"0036",{"IIXXIXXXIXIIXX","IIIIXXIIXXXXXX","IIIIIIXXXXZZZZ"},
        {"IIXZYZXIXIYYXX","IIIZYYXXYIXZXX","IIZYXXIXXYYZIZ","XXIYIXIXXXIZYX",
         "ZZYIXIXIXXXYZI","IIYZXXXZIXZXZI","IIZXXXXIZZZYIX","IIYXYXIXZIYXZX"},11179},
{"0040",{"IIXIIXXIXXXXIY","IIIXIXIXXXZIXX","IIIIXIXXXXIZZZ"},
        {"IIZYXYXIXIXYXZ","XIYXIXXIXIYYZX","IXIXXIZXIYXZXX","IXXXXXYXIXZXII",
         "XIIXXIXYIXXYXX","IIYZXIZXXIXXYY","IIXXXYYZIXXZIZ","ZIIXXZXXXIXYXI"},11179}
};
static constexpr int NORB = (int)(sizeof(ORBITS) / sizeof(ORBITS[0]));

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

// ---------------------------------------------------------------- errors and targets
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
static std::vector<uint32_t> U3;
static int64_t TOT_ERR = 0, TOTW4 = 0, U3_BY_W[MAXW + 1];
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
static uint32_t UBS[VHALF], WBS[VHALF], S3EL[1 << R1];
static inline uint32_t cls_of(uint32_t v) {
    uint32_t c = 0;
    for (int i = 0; i < VHALF; ++i) {
        c |= uint32_t(symp(WBS[i], v)) << i;
        c |= uint32_t(symp(UBS[i], v)) << (i + VHALF);
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
        UBS[np] = u; WBS[np] = w; ++np;
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
static std::vector<uint32_t> REPR10, TCLS;
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
    COV1.resize(VSZ); COVMAX = 0;
    for (int v = 0; v < VSZ; ++v) { int32_t c = (N - f[v]) / 2; COV1[v] = uint16_t(c);
        if (v && c > COVMAX) COVMAX = c; }
}
static void build_repr10() {
    REPR10.assign(VSZ, NONE); NW10 = 0; NCLS10 = 0;
    uint32_t basis[CDIM], bcls[CDIM];
    for (int g = 0; g < R1; ++g) { basis[g] = G3[g]; bcls[g] = 0; }
    for (int i = 0; i < VHALF; ++i) { basis[R1 + i] = UBS[i]; bcls[R1 + i] = 1u << i; }
    for (int i = 0; i < VHALF; ++i) { basis[R1 + VHALF + i] = WBS[i]; bcls[R1 + VHALF + i] = 1u << (i + VHALF); }
    uint32_t cur = 0, ccl = 0;
    for (uint64_t i = 1; i < (1ull << CDIM); ++i) {
        int b = std::countr_zero(i);
        cur ^= basis[b]; ccl ^= bcls[b];
        if (pwt(cur) == W2) { ++NW10;
            if (ccl && (REPR10[ccl] == NONE || cur < REPR10[ccl])) {
                if (REPR10[ccl] == NONE) ++NCLS10; REPR10[ccl] = cur; } }
    }
}

// ===========================================================================================
//  SOLUTION, EXACT SCORE, RESIDUAL ANALYSIS, INDEPENDENT VERIFICATION
// ===========================================================================================
struct Sol {
    uint32_t cls[R2]{}, pau[R2]{};
    int64_t span_sum = -1, covered = -1, residual = -1, res_w[MAXW + 1] = {0};
    bool valid = false;
};
static std::vector<std::pair<uint32_t, uint32_t>> PREF;   // class -> preferred Pauli
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
static void residual_of(const uint32_t* cls, std::vector<uint32_t>& out) {
    out.clear();
    for (size_t j = 0; j < U3.size(); ++j) {
        uint32_t t = TCLS[j]; int hit = 0;
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
    S.valid = true; return S;
}
// primary: fewest residual; then fewest weight-4, weight-3, weight-2, weight-1 residuals.
// Deliberately not a total order: equal count AND equal profile compare equal, so a mere
// relabelling never counts as an improvement and the search cannot cycle.
static bool better(const Sol& a, const Sol& b) {
    if (!b.valid) return a.valid;
    if (a.residual != b.residual) return a.residual < b.residual;
    for (int w = MAXW; w >= 1; --w) if (a.res_w[w] != b.res_w[w]) return a.res_w[w] < b.res_w[w];
    return false;
}
struct ResidualInfo {
    int64_t n = 0, by_w[MAXW + 1] = {0};
    bool subgroup = false, all_x = false;
    int dim = 0;
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
struct Verify {
    bool wt_ok = true, cg_ok = true, ch_ok = true, rank_ok = true, u3_ok = true, score_ok = true;
    int weights[R2]{}, rank = 0, comm[R1 + R2][R1 + R2]{};
    int64_t u3size = 0, covered = 0, per_gen[R2]{}, cov_by_w[MAXW + 1]{}, tgt_by_w[MAXW + 1]{};
    int64_t w4_total = 0, w4_det = 0, w4_und = 0, all_total = 0, orig_det = 0, det11 = 0;
    std::vector<uint32_t> leftover;
};
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
//  THE QUOTIENT Q = K^perp / K FOR A GIVEN SET OF KEPT GENERATORS
// ===========================================================================================
struct QSpace {
    int k = 0, dim = 0, size = 0;           // dim = 6 + 2k
    std::vector<uint32_t> lift;             // Q -> V (a coset representative in K^perp)
    std::vector<uint32_t> phi;              // coset sum of cov1
    std::vector<uint32_t> rep10;            // an A10 class in the coset, or NONE
    std::vector<uint16_t> mq;               // mq[q] such that <a,b> = parity(a & mq[b])
    std::vector<uint32_t> kspan;            // the 2^(8-k) elements of K, sorted
    int64_t sK = 0;                         // SUM over K\{0} of cov1
    uint32_t phimax = 0;
};
// Build the quotient for the kept classes.  Everything is exact: the basis of K^perp is found
// by elimination against the symplectic form, and the Q-form is read off the lifts.
static bool build_qspace(const uint32_t* keep, int nkeep, QSpace& Q) {
    Q.k = R2 - nkeep; Q.dim = 6 + 2 * Q.k; Q.size = 1 << Q.dim;
    // K span
    Q.kspan.assign(1, 0u);
    for (int i = 0; i < nkeep; ++i) { size_t n = Q.kspan.size();
        for (size_t j = 0; j < n; ++j) Q.kspan.push_back(Q.kspan[j] ^ keep[i]); }
    Q.sK = 0; for (uint32_t v : Q.kspan) if (v) Q.sK += COV1[v];
    std::sort(Q.kspan.begin(), Q.kspan.end());
    // basis of K^perp inside V
    std::vector<uint32_t> b;
    for (int i = 0; i < VDIM; ++i) b.push_back(1u << i);
    for (int i = 0; i < nkeep; ++i) {
        int piv = -1;
        for (size_t j = 0; j < b.size(); ++j) if (csymp(keep[i], b[j])) { piv = (int)j; break; }
        if (piv < 0) return false;
        uint32_t pv = b[piv]; b.erase(b.begin() + piv);
        for (auto& q : b) if (csymp(keep[i], q)) q ^= pv;
    }
    if ((int)b.size() != VDIM - nkeep) return false;
    // quotient by K: reduce the K^perp basis modulo K, keep dim Q independent survivors
    std::vector<uint32_t> qb;
    {
        uint32_t piv[VDIM] = {0};
        for (uint32_t v : Q.kspan) if (v) { uint32_t t = v;
            while (t) { int p = 31 - std::countl_zero(t);
                if (piv[p]) t ^= piv[p]; else { piv[p] = t; break; } } }
        for (uint32_t v : b) { uint32_t t = v;
            while (t) { int p = 31 - std::countl_zero(t);
                if (piv[p]) t ^= piv[p]; else { piv[p] = t; qb.push_back(v); break; } }
            if ((int)qb.size() == Q.dim) break; }
    }
    if ((int)qb.size() != Q.dim) return false;
    // lifts, phi, usable representatives
    Q.lift.assign(Q.size, 0u); Q.phi.assign(Q.size, 0u); Q.rep10.assign(Q.size, NONE);
    for (int q = 0; q < Q.size; ++q) {
        uint32_t v = 0, m = (uint32_t)q;
        while (m) { int i = std::countr_zero(m); m &= m - 1; v ^= qb[i]; }
        Q.lift[q] = v;
        uint32_t s = 0, best = NONE;
        for (uint32_t u : Q.kspan) { uint32_t c = v ^ u; s += COV1[c];
            if (REPR10[c] != NONE && c < best) best = c; }
        Q.phi[q] = s; Q.rep10[q] = best;
        if (q && s > Q.phimax) Q.phimax = s;
    }
    // symplectic form on Q, as mq[q] with <a,b> = parity(a & mq[b])
    Q.mq.assign(Q.size, 0);
    for (int q = 0; q < Q.size; ++q) {
        uint16_t m = 0;
        for (int i = 0; i < Q.dim; ++i) if (csymp(qb[i], Q.lift[q])) m |= uint16_t(1u << i);
        Q.mq[q] = m;
    }
    return true;
}
static inline int qsymp(const QSpace& Q, uint32_t a, uint32_t b) {
    return std::popcount((unsigned)(a & Q.mq[b])) & 1;
}

// ===========================================================================================
//  INCUMBENT, REPORTING, CHECKPOINTING
// ===========================================================================================
static Sol BEST;
static std::mutex BESTMTX;
static int IMPROVEMENTS = 0;
static std::string OUTDIR = "targeted_3opt4opt";
static std::string ORBNAME = "0006";
static std::string PHASENAME = "3opt";
static Sol START;

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
    ap("Residual F_2 span dimension: %d\n", r.dim);
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
    ap("Stage-1 stabilizer (frozen), orbit %s\n\n", ORBNAME.c_str());
    for (int i = 0; i < R1; ++i) ap("g%d = %s   (weight %d)\n", i + 1, pstr(G3[i]).c_str(), pwt(G3[i]));
    ap("\nM_3 = %lld / %lld\n|U_3| = %zu\n", (long long)V.orig_det, (long long)V.all_total, U3.size());
    dump(dir + "/stage1.txt");
    RB.clear();
    ap("Orbit %s, phase %s\n\nStage-2 generators (all weight 10):\n", ORBNAME.c_str(), PHASENAME.c_str());
    for (int a = 0; a < R2; ++a) ap("h%d = %s   (weight %d)\n", a + 1, pstr(S.pau[a]).c_str(), V.weights[a]);
    ap("\nM_8 = %lld / %zu\nresidual = %lld\n", (long long)S.covered, U3.size(), (long long)S.residual);
    ap("M_final = %lld / %lld  (%.8f %%)\n", (long long)V.det11, (long long)V.all_total,
       100.0 * double(V.det11) / double(V.all_total));
    ap("\nWeight-4 total: %lld\nWeight-4 detected: %lld\nWeight-4 undetected: %lld\n",
       (long long)V.w4_total, (long long)V.w4_det, (long long)V.w4_und);
    if (V.w4_und == 0) ap("ALL WEIGHT-4 ERRORS DETECTED\n");
    ap("\nTarget detection per generator:\n");
    for (int a = 0; a < R2; ++a)
        ap("h%d: %lld   (cov1 = %d %s)\n", a + 1, (long long)V.per_gen[a], (int)COV1[S.cls[a]],
           V.per_gen[a] == (int64_t)COV1[S.cls[a]] ? "OK" : "*** MISMATCH ***");
    ap("\nCommutation (g1..g3,h1..h8; 0 = commute):\n");
    for (int i = 0; i < R1 + R2; ++i) { ap("  ");
        for (int j = 0; j < R1 + R2; ++j) ap("%d ", V.comm[i][j]); ap("\n"); }
    ap("\nINDEPENDENT VERIFICATION (all 91770 errors regenerated from scratch)\n");
    ap("  |U_3| rebuilt            : %lld %s\n", (long long)V.u3size, V.u3_ok ? "OK" : "MISMATCH");
    ap("  all weights == 10        : %s\n", V.wt_ok ? "OK" : "FAIL");
    ap("  [h_a , g_i] = 0          : %s\n", V.cg_ok ? "OK" : "FAIL");
    ap("  [h_a , h_b] = 0          : %s\n", V.ch_ok ? "OK" : "FAIL");
    ap("  rank = 11                : %s (rank %d)\n", V.rank_ok ? "OK" : "FAIL", V.rank);
    ap("  brute-force union recount: %lld %s\n", (long long)V.covered,
       V.score_ok ? "(matches the algebraic score)" : "*** MISMATCH ***");
    dump(dir + "/stage2_best.txt");
    RB.clear(); residual_block(r); dump(dir + "/remaining_errors.txt");
    RB.clear();
    ap("Complete 11-generator stabilizer, orbit %s\n\n", ORBNAME.c_str());
    for (int i = 0; i < R1; ++i) ap("g%d = %s\n", i + 1, pstr(G3[i]).c_str());
    for (int a = 0; a < R2; ++a) ap("h%d = %s\n", a + 1, pstr(S.pau[a]).c_str());
    ap("\nH = [ X | Z ]  (11 x 28):\n");
    for (int i = 0; i < R1 + R2; ++i) ap("  %s | %s\n", bstr(px(all[i])).c_str(), bstr(pz(all[i])).c_str());
    ap("\nrank = %d, weights =", V.rank);
    for (int i = 0; i < R1 + R2; ++i) ap(" %d", pwt(all[i]));
    ap("\n");
    dump(dir + "/final_matrix.txt");
    RB.clear();
    ap("Orbit: %s\nPhase: %s\n\nM8: %lld / %zu\nResidual: %lld\n",
       ORBNAME.c_str(), PHASENAME.c_str(), (long long)S.covered, U3.size(), (long long)S.residual);
    ap("M_final: %lld / %lld\n\nResidual weight distribution: %lld / %lld / %lld / %lld\n\n",
       (long long)V.det11, (long long)V.all_total, (long long)S.res_w[1], (long long)S.res_w[2],
       (long long)S.res_w[3], (long long)S.res_w[4]);
    ap("Weight-4 total: %lld\nWeight-4 detected: %lld\nWeight-4 undetected: %lld\n\n",
       (long long)V.w4_total, (long long)V.w4_det, (long long)V.w4_und);
    residual_block(r);
    dump(dir + "/summary.txt");
}
static void report_new_best(const Sol& prev, const Sol& cur, const Verify& V) {
    ++IMPROVEMENTS;
    ResidualInfo r = analyse(V.leftover);
    RB.clear();
    ap("============================================================\nNEW BEST\n");
    ap("============================================================\n\n");
    ap("Orbit %s, phase %s\n\n", ORBNAME.c_str(), PHASENAME.c_str());
    ap("Previous residual: %lld\nNew residual:      %lld\n\n",
       (long long)prev.residual, (long long)cur.residual);
    ap("M8:      %lld / %zu\nM_final: %lld / %lld\n\n",
       (long long)cur.covered, U3.size(), (long long)V.det11, (long long)V.all_total);
    ap("Residual weight distribution: w1=%lld w2=%lld w3=%lld w4=%lld\n",
       (long long)cur.res_w[1], (long long)cur.res_w[2], (long long)cur.res_w[3], (long long)cur.res_w[4]);
    ap("Weight-4 total: %lld\nWeight-4 detected: %lld\nWeight-4 undetected: %lld\n",
       (long long)V.w4_total, (long long)V.w4_det, (long long)V.w4_und);
    if (V.w4_und == 0) ap("\nALL WEIGHT-4 ERRORS DETECTED\n");
    ap("\nGenerators:\n");
    for (int a = 0; a < R2; ++a) ap("h%d = %s\n", a + 1, pstr(cur.pau[a]).c_str());
    ap("\n"); residual_block(r);
    ap("\nTime: %.2f s\n============================================================\n", elapsed());
    fputs(RB.c_str(), stdout); fflush(stdout);
    char d[256]; snprintf(d, sizeof(d), "%s/orbit%s_%s/improved_%04d",
                          OUTDIR.c_str(), ORBNAME.c_str(), PHASENAME.c_str(), IMPROVEMENTS);
    std::string dir(d);
    { std::string keep = RB; std::error_code ec; std::filesystem::create_directories(dir, ec);
      RB = keep; dump(dir + "/new_best_report.txt"); }
    write_sol_dir(dir, cur, V);
    char b[256]; snprintf(b, sizeof(b), "%s/orbit%s_%s/best_known", OUTDIR.c_str(), ORBNAME.c_str(), PHASENAME.c_str());
    write_sol_dir(std::string(b), cur, V);
    if (cur.residual == 0) {
        RB.clear();
        ap("============================================================\nFULL COVERAGE FOUND\n");
        ap("============================================================\n\n");
        ap("All %zu Stage-1-undetected errors are detected.\n\n", U3.size());
        ap("Final detection:\n%lld / %lld\n\n", (long long)V.det11, (long long)V.all_total);
        ap("All Pauli errors of weight <= 4 are detected.\n\n");
        ap("Orbit %s.\n============================================================\n", ORBNAME.c_str());
        dump(OUTDIR + "/FULL_COVERAGE.txt", true);
    }
}
// Accept only after a from-scratch re-verification over all 91770 errors.
static bool offer(const uint32_t* cls) {
    if (!independent_isotropic(cls)) return false;
    Sol S = make_sol(cls);
    { std::lock_guard<std::mutex> lk(BESTMTX); if (!better(S, BEST)) return false; }
    Verify V = verify(S);
    if (!V.wt_ok || !V.cg_ok || !V.ch_ok || !V.rank_ok || !V.u3_ok || !V.score_ok) {
        printf("*** candidate failed independent verification, rejected ***\n"); return false; }
    std::lock_guard<std::mutex> lk(BESTMTX);
    if (!better(S, BEST)) return false;
    Sol prev = BEST; BEST = S;
    report_new_best(prev, S, V);
    return true;
}
static void save_checkpoint(int orbitIdx, int phase, int dropIdx, double secs) {
    RB.clear();
    ap("phase %d\norbit %d\ndrop %d\nelapsed %.2f\nresidual %lld\ncovered %lld\n",
       phase, orbitIdx, dropIdx, secs, (long long)BEST.residual, (long long)BEST.covered);
    ap("classes");
    for (int a = 0; a < R2; ++a) ap(" %u", BEST.cls[a]);
    ap("\ngenerators\n");
    for (int a = 0; a < R2; ++a) ap("%s\n", pstr(BEST.pau[a]).c_str());
    ap("residual_list\n");
    { std::vector<uint32_t> res; residual_of(BEST.cls, res);
      for (uint32_t E : res) ap("%s\n", pstr(E).c_str()); }
    dump(OUTDIR + "/checkpoint.txt");
}

// ===========================================================================================
//  EXHAUSTIVE 3-OPT  (complete sweep of the whole 3-opt neighbourhood, no breadth limit)
//  For every one of the C(8,3) = 56 drop-triples, every 3-dimensional isotropic subspace of
//  Q = K^perp/K (dim 12) is enumerated.  q1 is forced to be the MINIMUM of the seven non-zero
//  span elements, via P(x) = (x > q1) && ((q1^x) > q1), which removes most of the 168 bases
//  per subspace at the cost of three integer comparisons.  The remaining pair {q2,q3} is
//  swept with the psi-sorted bound  S <= base + psi(q2) + psi(q3) + max psi.
// ===========================================================================================
static std::atomic<uint64_t> LEAVES{0};
static bool opt3_drop(const QSpace& Q, const uint32_t* keep, const int* dropPos, int nthreads,
                      int64_t need, std::vector<std::array<uint32_t,3>>& hits, size_t cap) {
    (void)keep; (void)dropPos;
    std::atomic<int> next{0};
    std::mutex mtx;
    std::vector<std::thread> th;
    const int n = Q.size;
    for (int t = 0; t < nthreads; ++t) th.emplace_back([&]() {
        std::vector<uint32_t> psi((size_t)n), L, val;
        std::vector<std::array<uint32_t,3>> local;
        uint64_t leaves = 0;
        for (;;) {
            int q1 = next.fetch_add(1);
            if (q1 >= n) break;
            if (q1 == 0 || Q.rep10[q1] == NONE) continue;
            const int64_t base = Q.sK + (int64_t)Q.phi[q1];
            uint32_t pmax = 0;
            for (int w = 0; w < n; ++w) { uint32_t s = Q.phi[w] + Q.phi[w ^ (uint32_t)q1];
                psi[w] = s; if (s > pmax) pmax = s; }
            L.clear();
            for (int q = q1 + 1; q < n; ++q) {
                if (Q.rep10[q] == NONE) continue;
                if ((uint32_t)(q1 ^ q) <= (uint32_t)q1) continue;      // P(q)
                if (qsymp(Q, (uint32_t)q, (uint32_t)q1)) continue;
                L.push_back((uint32_t)q);
            }
            std::sort(L.begin(), L.end(), [&](uint32_t a, uint32_t b) {
                if (psi[a] != psi[b]) return psi[a] > psi[b]; return a < b; });
            const size_t m = L.size();
            val.resize(m);
            for (size_t i = 0; i < m; ++i) val[i] = psi[L[i]];
            for (size_t i = 0; i + 1 < m; ++i) {
                const int64_t vi = val[i];
                if (base + 2 * vi + (int64_t)pmax <= need) break;
                const uint32_t a = L[i];
                for (size_t j = i + 1; j < m; ++j) {
                    const int64_t vj = val[j];
                    if (base + vi + vj + (int64_t)pmax <= need) break;
                    const uint32_t b = L[j];
                    if (qsymp(Q, a, b)) continue;
                    const uint32_t x = a ^ b;
                    if ((x <= (uint32_t)q1) || ((uint32_t)(q1 ^ x) <= (uint32_t)q1)) continue;  // P(a^b)
                    ++leaves;
                    const int64_t s = base + vi + vj + (int64_t)psi[x];
                    if (s > need) {
                        local.push_back({(uint32_t)q1, a, b});
                        if (local.size() > cap) local.resize(cap);
                    }
                }
            }
        }
        LEAVES += leaves;
        if (!local.empty()) { std::lock_guard<std::mutex> lk(mtx);
            hits.insert(hits.end(), local.begin(), local.end()); }
    });
    for (auto& x : th) x.join();
    return !hits.empty();
}

// ---------------------------------------------------------------- 3-opt over all 56 drops
static bool run_3opt(int nthreads, bool verbose, double deadline, int orbitIdx, int phaseId) {
    bool improved = false;
    for (;;) {
        bool any = false;
        int dropIdx = 0;
        for (int d1 = 0; d1 < R2 && !any; ++d1)
        for (int d2 = d1 + 1; d2 < R2 && !any; ++d2)
        for (int d3 = d2 + 1; d3 < R2 && !any; ++d3, ++dropIdx) {
            if (elapsed() > deadline) return improved;
            uint32_t keep[R2 - 3]; int nk = 0, pos[3] = {d1, d2, d3};
            for (int a = 0; a < R2; ++a) if (a != d1 && a != d2 && a != d3) keep[nk++] = BEST.cls[a];
            QSpace Q;
            if (!build_qspace(keep, nk, Q)) { printf("  quotient construction failed\n"); continue; }
            std::vector<std::array<uint32_t,3>> hits;
            double t0 = elapsed();
            opt3_drop(Q, keep, pos, nthreads, BEST.span_sum, hits, 4096);
            if (verbose)
                printf("    drop(%d,%d,%d): |Q|=%d, sK=%lld, phimax=%u, %zu hits, %.1f s\n",
                       d1, d2, d3, Q.size, (long long)Q.sK, Q.phimax, hits.size(), elapsed() - t0);
            fflush(stdout);
            for (auto& h : hits) {
                uint32_t trial[R2]; std::copy(BEST.cls, BEST.cls + R2, trial);
                trial[d1] = Q.rep10[h[0]]; trial[d2] = Q.rep10[h[1]]; trial[d3] = Q.rep10[h[2]];
                if (trial[d1] == NONE || trial[d2] == NONE || trial[d3] == NONE) continue;
                if (offer(trial)) { any = true; improved = true; break; }
            }
            save_checkpoint(orbitIdx, phaseId, dropIdx, elapsed());
        }
        if (!any) break;
    }
    return improved;
}

// ===========================================================================================
//  TARGETED 4-OPT.  dim Q = 14 (16384 classes); the number of 4-dimensional isotropic
//  subspaces per drop-quadruple is about 5.6e10, so an exhaustive sweep is out of reach and
//  the first TWO basis elements are breadth-limited (top-M by coset gain), while the last two
//  are swept exhaustively with the same psi-sorted bound.  psi2(w) = psi1(w) + psi1(w^q2)
//  makes the bottom level the identical three-term pair problem.
// ===========================================================================================
static bool run_4opt(int nthreads, int M1, int M2, bool verbose, double deadline,
                     int orbitIdx, int phaseId) {
    bool improved = false;
    int dropIdx = 0;
    for (int d1 = 0; d1 < R2; ++d1)
    for (int d2 = d1 + 1; d2 < R2; ++d2)
    for (int d3 = d2 + 1; d3 < R2; ++d3)
    for (int d4 = d3 + 1; d4 < R2; ++d4, ++dropIdx) {
        if (elapsed() > deadline) return improved;
        uint32_t keep[R2 - 4]; int nk = 0;
        for (int a = 0; a < R2; ++a) if (a != d1 && a != d2 && a != d3 && a != d4) keep[nk++] = BEST.cls[a];
        QSpace Q;
        if (!build_qspace(keep, nk, Q)) continue;
        const int n = Q.size;
        std::vector<uint32_t> ord;
        for (int q = 1; q < n; ++q) if (Q.rep10[q] != NONE) ord.push_back((uint32_t)q);
        std::sort(ord.begin(), ord.end(), [&](uint32_t a, uint32_t b) {
            if (Q.phi[a] != Q.phi[b]) return Q.phi[a] > Q.phi[b]; return a < b; });
        const int m1 = (M1 <= 0) ? (int)ord.size() : std::min<int>(M1, (int)ord.size());
        double t0 = elapsed();
        std::atomic<int> nexti{0};
        std::mutex mtx;
        std::vector<std::array<uint32_t,4>> hits;
        std::vector<std::thread> th;
        for (int t = 0; t < nthreads; ++t) th.emplace_back([&]() {
            std::vector<uint32_t> psi1((size_t)n), psi2((size_t)n), L1, L2, val;
            std::vector<std::array<uint32_t,4>> local;
            for (;;) {
                int ii = nexti.fetch_add(1);
                if (ii >= m1) break;
                if (elapsed() > deadline) break;
                const uint32_t q1 = ord[ii];
                for (int w = 0; w < n; ++w) psi1[w] = Q.phi[w] + Q.phi[w ^ q1];
                L1.clear();
                for (uint32_t q : ord) { if (q == q1) continue;
                    if (qsymp(Q, q, q1)) continue; L1.push_back(q); }
                std::sort(L1.begin(), L1.end(), [&](uint32_t a, uint32_t b) {
                    if (psi1[a] != psi1[b]) return psi1[a] > psi1[b]; return a < b; });
                const int m2 = (M2 <= 0) ? (int)L1.size() : std::min<int>(M2, (int)L1.size());
                for (int jj = 0; jj < m2; ++jj) {
                    const uint32_t q2 = L1[jj];
                    uint32_t p2max = 0;
                    for (int w = 0; w < n; ++w) { uint32_t s = psi1[w] + psi1[w ^ q2];
                        psi2[w] = s; if (s > p2max) p2max = s; }
                    const int64_t base = Q.sK + (int64_t)Q.phi[q1] + (int64_t)psi1[q2];
                    const uint32_t sp[4] = {0u, q1, q2, q1 ^ q2};
                    // A qualifying pair needs psi2(a)+psi2(b)+psi2(a^b) > need-base, and both
                    // psi2(b) and psi2(a^b) are at most p2max, so BOTH members must satisfy
                    // psi2 > need - base - 2*p2max.  Filtering on that first shrinks L2 from
                    // several thousand to a handful, which is what makes the per-(q1,q2) sort
                    // affordable -- sorting the unfiltered list was the real bottleneck.
                    const int64_t tau = BEST.span_sum - base - 2 * (int64_t)p2max;
                    L2.clear();
                    for (uint32_t q : L1) {
                        if ((int64_t)psi2[q] <= tau) continue;
                        if (qsymp(Q, q, q2)) continue;
                        if (q == sp[1] || q == sp[2] || q == sp[3]) continue;
                        L2.push_back(q);
                    }
                    if (L2.size() < 2) continue;
                    std::sort(L2.begin(), L2.end(), [&](uint32_t a, uint32_t b) {
                        if (psi2[a] != psi2[b]) return psi2[a] > psi2[b]; return a < b; });
                    const size_t mm = L2.size();
                    val.resize(mm);
                    for (size_t i = 0; i < mm; ++i) val[i] = psi2[L2[i]];
                    for (size_t i = 0; i + 1 < mm; ++i) {
                        const int64_t vi = val[i];
                        if (base + 2 * vi + (int64_t)p2max <= BEST.span_sum) break;
                        const uint32_t a = L2[i];
                        for (size_t j = i + 1; j < mm; ++j) {
                            const int64_t vj = val[j];
                            if (base + vi + vj + (int64_t)p2max <= BEST.span_sum) break;
                            const uint32_t b = L2[j];
                            if (qsymp(Q, a, b)) continue;
                            const uint32_t x = a ^ b;
                            if (x == sp[1] || x == sp[2] || x == sp[3]) continue;   // independence
                            const int64_t s = base + vi + vj + (int64_t)psi2[x];
                            if (s > BEST.span_sum) { local.push_back({q1, q2, a, b});
                                if (local.size() > 2048) local.resize(2048); }
                        }
                    }
                }
            }
            if (!local.empty()) { std::lock_guard<std::mutex> lk(mtx);
                hits.insert(hits.end(), local.begin(), local.end()); }
        });
        for (auto& x : th) x.join();
        if (verbose)
            printf("    drop(%d,%d,%d,%d): |Q|=%d, M1=%d, %zu hits, %.1f s\n",
                   d1, d2, d3, d4, n, m1, hits.size(), elapsed() - t0);
        fflush(stdout);
        for (auto& h : hits) {
            uint32_t trial[R2]; std::copy(BEST.cls, BEST.cls + R2, trial);
            trial[d1] = Q.rep10[h[0]]; trial[d2] = Q.rep10[h[1]];
            trial[d3] = Q.rep10[h[2]]; trial[d4] = Q.rep10[h[3]];
            bool bad = false;
            for (int a = 0; a < R2; ++a) if (trial[a] == NONE) bad = true;
            if (bad) continue;
            if (offer(trial)) { improved = true; break; }
        }
        save_checkpoint(orbitIdx, phaseId, dropIdx, elapsed());
        if (improved) return improved;
    }
    return improved;
}

// ---------------------------------------------------------------- per-orbit setup
static bool load_orbit(int idx, bool quiet) {
    const OrbitDef& O = ORBITS[idx];
    ORBNAME = O.name;
    for (int i = 0; i < R1; ++i) G3[i] = parse_pauli(O.g[i]);
    for (int a = 0; a < R2; ++a) H0[a] = parse_pauli(O.h[a]);
    for (int i = 0; i < R1; ++i) if (pwt(G3[i]) != W1) { printf("orbit %s: bad g weight\n", O.name); return false; }
    for (int i = 0; i < R1; ++i) for (int j = i + 1; j < R1; ++j)
        if (symp(G3[i], G3[j])) { printf("orbit %s: g do not commute\n", O.name); return false; }
    build_targets();
    if (U3.size() != 11186) { printf("orbit %s: |U_3| = %zu != 11186\n", O.name, U3.size()); return false; }
    build_basis(); build_cov1(); build_repr10();
    PREF.clear();
    uint32_t cls0[R2];
    for (int a = 0; a < R2; ++a) {
        cls0[a] = cls_of(H0[a]);
        PREF.push_back({cls0[a], H0[a]});
        if (pwt(H0[a]) != W2) { printf("orbit %s: h%d weight != 10\n", O.name, a + 1); return false; }
        for (int i = 0; i < R1; ++i) if (symp(H0[a], G3[i])) {
            printf("orbit %s: h%d does not commute with g%d\n", O.name, a + 1, i + 1); return false; }
    }
    if (!independent_isotropic(cls0)) { printf("orbit %s: incumbent not isotropic/independent\n", O.name); return false; }
    Sol S0 = make_sol(cls0);
    Verify V0 = verify(S0);
    if (!V0.wt_ok || !V0.cg_ok || !V0.ch_ok || !V0.rank_ok || !V0.score_ok) {
        printf("orbit %s: incumbent FAILED verification\n", O.name); return false; }
    if (S0.covered != O.start_m8)
        printf("  note: orbit %s incumbent scores %lld, the table said %d\n",
               O.name, (long long)S0.covered, O.start_m8);
    BEST = S0; START = S0; IMPROVEMENTS = 0;
    if (!quiet) {
        printf("orbit %s: M_8 = %lld / %zu, residual %lld (w %lld/%lld/%lld/%lld), M_final = %lld / %lld\n",
               O.name, (long long)S0.covered, U3.size(), (long long)S0.residual,
               (long long)S0.res_w[1], (long long)S0.res_w[2], (long long)S0.res_w[3],
               (long long)S0.res_w[4], (long long)V0.det11, (long long)V0.all_total);
        printf("         verification: weights OK, [h,g] OK, [h,h] OK, rank %d, brute-force %lld OK;"
               "  weight-4 undetected %lld / %lld\n", V0.rank, (long long)V0.covered,
               (long long)V0.w4_und, (long long)V0.w4_total);
    }
    return true;
}

// ===========================================================================================
//  MAIN
// ===========================================================================================
struct Row { std::string orbit; int64_t start_m8, best_m8, residual, w4res; bool improved; };
static std::vector<Row> TABLE;

static void usage() {
    printf("stab14_kopt -- exhaustive 3-opt and targeted 4-opt around the best-known 3+8 solutions\n"
           "  --phase 3opt|4opt|calibrate     --pipeline    (phase 1 -> 2 -> 3)\n"
           "  --orbit 0006 | --all-top-orbits\n"
           "  --M1 N  --M2 N     4-opt breadth of the first two replacements (default 128/128)\n"
           "  --threads N  --seed N  --time-limit N  --deterministic  --selftest N\n"
           "  --out DIR  --verbose  --resume  --checkpoint N\n");
}
static int orbit_index(const std::string& nm) {
    for (int i = 0; i < NORB; ++i) if (nm == ORBITS[i].name) return i;
    return -1;
}
// One drop-triple / drop-quadruple timed, then extrapolated to the whole phase.
static void calibrate(int nthreads, int M1, int M2) {
    printf("\n=== CALIBRATION (orbit %s) ===\n", ORBNAME.c_str());
    { uint32_t keep[R2 - 3]; int nk = 0;
      for (int a = 3; a < R2; ++a) keep[nk++] = BEST.cls[a];
      QSpace Q; double t0 = elapsed();
      if (!build_qspace(keep, nk, Q)) { printf("quotient build failed\n"); return; }
      double tb = elapsed() - t0;
      std::vector<std::array<uint32_t,3>> hits;
      t0 = elapsed();
      LEAVES = 0;
      opt3_drop(Q, keep, nullptr, nthreads, BEST.span_sum, hits, 4096);
      double ts = elapsed() - t0;
      printf("3-opt: |Q| = %d (dim %d),  quotient build %.2f s,  one drop-triple %.2f s,"
             "  %llu leaves\n", Q.size, Q.dim, tb, ts, (unsigned long long)LEAVES.load());
      printf("  -> 56 drop-triples ~ %.1f s per orbit;  7 orbits ~ %.1f min"
             " (plus ~%.0f s setup per orbit)\n", 56 * (ts + tb), 7 * 56 * (ts + tb) / 60.0, 20.0);
    }
    { uint32_t keep[R2 - 4]; int nk = 0;
      for (int a = 4; a < R2; ++a) keep[nk++] = BEST.cls[a];
      QSpace Q; double t0 = elapsed();
      if (!build_qspace(keep, nk, Q)) { printf("quotient build failed\n"); return; }
      double tb = elapsed() - t0;
      printf("4-opt: |Q| = %d (dim %d), quotient build %.2f s\n", Q.size, Q.dim, tb);
      // time a single (q1,q2) inner sweep and scale
      std::vector<uint32_t> ord;
      for (int q = 1; q < Q.size; ++q) if (Q.rep10[q] != NONE) ord.push_back((uint32_t)q);
      std::sort(ord.begin(), ord.end(), [&](uint32_t a, uint32_t b) { return Q.phi[a] > Q.phi[b]; });
      std::vector<uint32_t> psi1(Q.size), psi2(Q.size);
      const uint32_t q1 = ord[0];
      for (int w = 0; w < Q.size; ++w) psi1[w] = Q.phi[w] + Q.phi[w ^ q1];
      std::vector<uint32_t> L1;
      for (uint32_t q : ord) if (q != q1 && !qsymp(Q, q, q1)) L1.push_back(q);
      std::sort(L1.begin(), L1.end(), [&](uint32_t a, uint32_t b) { return psi1[a] > psi1[b]; });
      t0 = elapsed();
      uint64_t cnt = 0;
      const uint32_t q2 = L1[0];
      uint32_t p2max = 0;
      for (int w = 0; w < Q.size; ++w) { uint32_t s = psi1[w] + psi1[w ^ q2];
          psi2[w] = s; if (s > p2max) p2max = s; }
      const int64_t base = Q.sK + (int64_t)Q.phi[q1] + (int64_t)psi1[q2];
      std::vector<uint32_t> L2;
      for (uint32_t q : L1) if (!qsymp(Q, q, q2) && q != q2 && q != (q1 ^ q2)) L2.push_back(q);
      std::sort(L2.begin(), L2.end(), [&](uint32_t a, uint32_t b) { return psi2[a] > psi2[b]; });
      for (size_t i = 0; i + 1 < L2.size(); ++i) {
          if (base + 2 * (int64_t)psi2[L2[i]] + (int64_t)p2max <= BEST.span_sum) break;
          for (size_t j = i + 1; j < L2.size(); ++j) {
              if (base + (int64_t)psi2[L2[i]] + (int64_t)psi2[L2[j]] + (int64_t)p2max <= BEST.span_sum) break;
              ++cnt; }
      }
      double ts = elapsed() - t0;
      printf("  one (q1,q2) inner sweep: |L2| = %zu, %llu pairs surviving the bound, %.4f s"
             " (single-threaded)\n", L2.size(), (unsigned long long)cnt, ts);
      double per = ts / std::max(1.0, (double)nthreads);
      printf("  -> 70 drop-quadruples x M1=%d x M2=%d ~ %.1f min at %d threads\n",
             M1, M2, 70.0 * M1 * M2 * per / 60.0, nthreads);
    }
    printf("=== END CALIBRATION ===\n\n");
    fflush(stdout);
}

// POSITIVE CONTROL.  Deliberately replace `n` generators of the loaded incumbent by random
// valid weight-10 classes, making the solution worse, then let the search run.  If the k-opt
// engine is sound it must climb back to (at least) the original score; "no improvement found"
// on the real incumbent is only meaningful once this test passes.
static bool degrade(int n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    uint32_t trial[R2]; std::copy(BEST.cls, BEST.cls + R2, trial);
    int idx[R2]; std::iota(idx, idx + R2, 0); std::shuffle(idx, idx + R2, rng);
    for (int t = 0; t < n; ++t) {
        const int pos = idx[t];
        bool ok = false;
        for (int tries = 0; tries < 4000000 && !ok; ++tries) {
            uint32_t c = (uint32_t)(rng() % (uint64_t)VSZ);
            if (!c || REPR10[c] == NONE) continue;
            bool good = true;
            for (int a = 0; a < R2 && good; ++a) if (a != pos && csymp(c, trial[a])) good = false;
            if (!good) continue;
            uint32_t save = trial[pos]; trial[pos] = c;
            if (independent_isotropic(trial)) ok = true; else trial[pos] = save;
        }
        if (!ok) return false;
    }
    BEST = make_sol(trial);
    return true;
}

int kopt_main(int argc, char** argv) {
    int nthreads = (int)std::max(1u, std::thread::hardware_concurrency());
    int M1 = 128, M2 = 128, selftestN = 0, ckpt = 1, degradeN = 0;
    uint64_t seed = 1;
    double tlimit = 1e18;
    std::string phase = "3opt", orbitArg = "0006";
    bool allOrbits = false, pipeline = false, verbose = false, deterministic = false, resume = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto nxt = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : std::string(); };
        if      (a == "--phase")            phase = nxt();
        else if (a == "--orbit")            orbitArg = nxt();
        else if (a == "--all-top-orbits")   allOrbits = true;
        else if (a == "--pipeline")         pipeline = true;
        else if (a == "--M1")               M1 = atoi(nxt().c_str());
        else if (a == "--M2")               M2 = atoi(nxt().c_str());
        else if (a == "--M")              { M1 = M2 = atoi(nxt().c_str()); }
        else if (a == "--threads")          nthreads = std::max(1, atoi(nxt().c_str()));
        else if (a == "--seed")             seed = strtoull(nxt().c_str(), nullptr, 10);
        else if (a == "--time-limit")       tlimit = atof(nxt().c_str());
        else if (a == "--deterministic")    deterministic = true;
        else if (a == "--selftest")         selftestN = atoi(nxt().c_str());
        else if (a == "--out")              OUTDIR = nxt();
        else if (a == "--verbose")          verbose = true;
        else if (a == "--resume")           resume = true;
        else if (a == "--checkpoint")       ckpt = std::max(1, atoi(nxt().c_str()));
        else if (a == "--degrade")          degradeN = atoi(nxt().c_str());
        else { usage(); return (a == "--help" || a == "-h") ? 0 : 1; }
    }
    (void)seed; (void)deterministic; (void)resume; (void)ckpt; (void)selftestN;
    T0 = std::chrono::steady_clock::now();
    std::error_code ec; std::filesystem::create_directories(OUTDIR, ec);

    printf("=== targeted 3-opt / 4-opt around the best-known 3+8 solutions ===\n");
    printf("threads %d,  output root %s\n", nthreads, OUTDIR.c_str());

    std::vector<int> orbitList;
    if (allOrbits || pipeline) { for (int i = 0; i < NORB; ++i) orbitList.push_back(i); }
    else { int k = orbit_index(orbitArg);
           if (k < 0) { printf("unknown orbit %s\n", orbitArg.c_str()); return 1; }
           orbitList.push_back(k); }

    if (phase == "calibrate") {
        if (!load_orbit(orbitList[0], false)) return 2;
        calibrate(nthreads, M1, M2);
        return 0;
    }

    // ---------------------------------------------------------------- PHASE 1 / 2 : 3-opt
    int64_t globalBest = -1; std::string globalOrbit;
    Sol globalSol; Verify globalV;
    for (size_t oi = 0; oi < orbitList.size(); ++oi) {
        if (elapsed() > tlimit) { printf("\ntime limit reached\n"); break; }
        PHASENAME = "3opt";
        printf("\n--- orbit %s : exhaustive 3-opt (complete sweep of all 56 drop-triples) ---\n",
               ORBITS[orbitList[oi]].name);
        fflush(stdout);
        if (!load_orbit(orbitList[oi], false)) continue;
        if (degradeN > 0) {
            if (!degrade(degradeN, seed)) { printf("  degradation failed\n"); continue; }
            printf("  [positive control] degraded %d generators -> M_8 = %lld, residual %lld;"
                   " the search must climb back\n", degradeN, (long long)BEST.covered,
                   (long long)BEST.residual);
            fflush(stdout);
        }
        const int64_t start = BEST.covered;
        double t0 = elapsed();
        bool imp = run_3opt(nthreads, verbose, tlimit, orbitList[oi], 1);
        Verify V = verify(BEST);
        printf("  orbit %s: %s   M_8 %lld -> %lld,  residual %lld (w4 %lld),  %.1f s\n",
               ORBNAME.c_str(), imp ? "IMPROVED" : "no improvement -- 3-opt neighbourhood EXHAUSTED",
               (long long)start, (long long)BEST.covered, (long long)BEST.residual,
               (long long)V.w4_und, elapsed() - t0);
        fflush(stdout);
        char dir[256]; snprintf(dir, sizeof(dir), "%s/orbit%s_3opt/best_known", OUTDIR.c_str(), ORBNAME.c_str());
        write_sol_dir(std::string(dir), BEST, V);
        TABLE.push_back({ORBNAME, start, BEST.covered, BEST.residual, V.w4_und, imp});
        if (BEST.covered > globalBest) { globalBest = BEST.covered; globalOrbit = ORBNAME;
                                         globalSol = BEST; globalV = V; }
        if (BEST.residual == 0) { printf("\nFULL COVERAGE -- stopping\n"); break; }
        if (!pipeline && !allOrbits) break;
    }

    // ---------------------------------------------------------------- comparison table
    RB.clear();
    ap("Orbit    Starting M8    Best M8    Improvement    Residual    W4 residual\n");
    ap("---------------------------------------------------------------------------\n");
    for (const Row& r : TABLE)
        ap("%-8s %-14lld %-10lld %-14lld %-11lld %lld\n", r.orbit.c_str(),
           (long long)r.start_m8, (long long)r.best_m8, (long long)(r.best_m8 - r.start_m8),
           (long long)r.residual, (long long)r.w4res);
    ap("\nGlobal best after 3-opt: orbit %s, M_8 = %lld / %zu\n",
       globalOrbit.c_str(), (long long)globalBest, U3.size());
    dump(OUTDIR + "/phase12_table.txt", true);

    // ---------------------------------------------------------------- PHASE 3 : 4-opt
    bool did4 = false;
    if ((pipeline && globalBest <= 11181) || phase == "4opt") {
        did4 = true;
        int gi = orbit_index(globalOrbit.empty() ? orbitArg : globalOrbit);
        if (gi < 0) gi = 0;
        PHASENAME = "4opt";
        printf("\n--- orbit %s : targeted 4-opt (M1 = %d, M2 = %d; last two replacements exhaustive) ---\n",
               ORBITS[gi].name, M1, M2);
        fflush(stdout);
        if (load_orbit(gi, false)) {
            double t0 = elapsed();
            bool imp = run_4opt(nthreads, M1, M2, verbose, tlimit, gi, 3);
            Verify V = verify(BEST);
            printf("  4-opt: %s   M_8 = %lld,  residual %lld,  %.1f s\n",
                   imp ? "IMPROVED" : "no improvement at this breadth",
                   (long long)BEST.covered, (long long)BEST.residual, elapsed() - t0);
            char dir[256]; snprintf(dir, sizeof(dir), "%s/final_4opt/best_known", OUTDIR.c_str());
            write_sol_dir(std::string(dir), BEST, V);
            if (BEST.covered > globalBest) { globalBest = BEST.covered; globalOrbit = ORBITS[gi].name;
                                             globalSol = BEST; globalV = V; }
        }
    }

    // ---------------------------------------------------------------- final report
    RB.clear();
    ap("============================================================\n");
    ap("TARGETED 3-OPT / 4-OPT RESULT\n");
    ap("============================================================\n\n");
    ap("PROVEN STATEMENTS\n");
    ap("  Stage-1 optimum:  M3 = 80584 / 91770   PROVEN GLOBAL OPTIMUM\n");
    ap("                    |U_3| = 11186 for every orbit of the optimal family\n");
    ap("  Orbit 0006 incumbent:  1-opt optimality PROVEN, 2-opt optimality PROVEN\n");
    ap("                         (established by the earlier exhaustive refiner)\n");
    for (const Row& r : TABLE)
        if (!r.improved)
            ap("  Orbit %-5s 3-opt neighbourhood EXHAUSTED with no improvement -> the incumbent\n"
               "             is 3-opt optimal (every replacement of any three of its eight\n"
               "             generators, over every 3-dim isotropic subspace of K^perp/K)\n", r.orbit.c_str());
    ap("\nSEARCH RESULTS\n");
    ap("  3-opt: BEST FOUND = %lld / %zu  (orbit %s)\n", (long long)globalBest, U3.size(), globalOrbit.c_str());
    if (did4) ap("  4-opt: BEST FOUND = %lld / %zu  (breadth-limited M1=%d M2=%d, NOT exhaustive)\n",
                 (long long)globalBest, U3.size(), M1, M2);
    else      ap("  4-opt: not run\n");
    ap("\n");
    for (const Row& r : TABLE)
        ap("  orbit %-5s  M_8 %lld -> %lld   residual %lld   weight-4 residual %lld\n",
           r.orbit.c_str(), (long long)r.start_m8, (long long)r.best_m8,
           (long long)r.residual, (long long)r.w4res);
    if (globalSol.valid) {
        ap("\nGLOBAL BEST (orbit %s)\n", globalOrbit.c_str());
        ap("  M_8 = %lld / %zu,  M_final = %lld / %lld\n", (long long)globalSol.covered, U3.size(),
           (long long)globalV.det11, (long long)globalV.all_total);
        ap("  Weight-4 total: %lld\n  Weight-4 detected: %lld\n  Weight-4 undetected: %lld\n",
           (long long)globalV.w4_total, (long long)globalV.w4_det, (long long)globalV.w4_und);
        if (globalV.w4_und == 0) ap("  ALL WEIGHT-4 ERRORS DETECTED\n");
        ap("\n");
        ResidualInfo r = analyse(globalV.leftover);
        residual_block(r);
    }
    ap("\nNO GLOBAL OPTIMALITY IS CLAIMED.  The 3-opt sweeps are exhaustive over their own\n"
       "neighbourhood only; the 4-opt sweep is breadth-limited and therefore heuristic.\n");
    ap("Elapsed: %.1f s\n", elapsed());
    ap("============================================================\n");
    dump(OUTDIR + "/final_report.txt", true);
    return 0;
}

// Diagnostic: is the "weight-10 usable" flag on quotient elements structural or a bug?
int main() {
    T0 = std::chrono::steady_clock::now();
    if (!load_orbit(0, true)) return 1;
    printf("orbit 0006 loaded: |A10| = %lld of %d classes (%.2f%%)\n",
           (long long)NCLS10, VSZ - 1, 100.0 * double(NCLS10) / double(VSZ - 1));
    for (int k = 3; k <= 4; ++k) {
        int nk = R2 - k;
        std::vector<uint32_t> keep;
        for (int a = k; a < R2; ++a) keep.push_back(BEST.cls[a]);
        QSpace Q;
        if (!build_qspace(keep.data(), nk, Q)) { printf("build failed\n"); return 1; }
        // histogram: how many A10 classes does each K-coset contain?
        std::vector<int> hist(Q.kspan.size() + 1, 0);
        int usable = 0;
        for (int q = 1; q < Q.size; ++q) {
            int c = 0;
            for (uint32_t u : Q.kspan) if (REPR10[Q.lift[q] ^ u] != NONE) ++c;
            hist[c]++;
            if (Q.rep10[q] != NONE) ++usable;
        }
        printf("\nk = %d: |Q| = %d, |K| = %zu, usable q = %d (%.1f%%)\n",
               k, Q.size, Q.kspan.size(), usable, 100.0 * usable / (Q.size - 1));
        printf("  histogram of A10 members per K-coset:\n");
        for (size_t c = 0; c < hist.size(); ++c) if (hist[c]) printf("    %2zu members : %d cosets\n", c, hist[c]);
    }
    return 0;
}
