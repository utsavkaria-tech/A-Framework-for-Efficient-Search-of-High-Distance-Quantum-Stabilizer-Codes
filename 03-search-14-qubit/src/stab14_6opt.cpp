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
#include <ctime>

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
    std::vector<uint32_t> mq;               // mq[q] such that <a,b> = parity(a & mq[b])  (18 bits for k=6)
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
        uint32_t m = 0;
        for (int i = 0; i < Q.dim; ++i) if (csymp(qb[i], Q.lift[q])) m |= uint32_t(1u << i);
        Q.mq[q] = m;
    }
    return true;
}
static inline int qsymp(const QSpace& Q, uint32_t a, uint32_t b) {
    return std::popcount((unsigned)(a & Q.mq[b])) & 1;
}

// ===========================================================================================
//  6-OPT ENGINE, LOGGING, CHECKPOINTING
//
//  A 6-opt move retains exactly two of the eight weight-10 generators and replaces the other
//  six.  With K = span of the two retained classes (dim 2), the admissible results are the
//  8-dimensional isotropic D' with K <= D' <= K^perp, i.e. the 6-dimensional totally
//  isotropic subspaces of  Q = K^perp/K,  dim Q = 22 - 2*2 = 18,  |Q| = 262144.
//
//  Score telescoping (exact, no proxy).  With phi(q) = SUM_{u in K} cov1(lift(q)^u) and
//        psi1(w) = phi(w)  + phi(w^q1)
//        psi2(w) = psi1(w) + psi1(w^q2)
//        psi3(w) = psi2(w) + psi2(w^q3)
//        psi4(w) = psi3(w) + psi3(w^q4)
//  the span sum of D' = K + span(q1..q6) is exactly
//        S = S(K) + phi(q1) + psi1(q2) + psi2(q3) + psi3(q4)
//                 + psi4(q5) + psi4(q6) + psi4(q5^q6).
//  So the bottom two levels are again the three-term pair problem, swept exhaustively with
//  the tau filter (both members of a qualifying pair must have psi4 > need - base - 2*psi4max).
//
//  SIZE OF THE EXHAUSTIVE NEIGHBOURHOOD (why --M 0 is not runnable, see README):
//        6-dim isotropic subspaces of Sp(18,2) = 4.887e17 per drop-pair
//        x C(8,2) = 28 retained pairs                = 1.368e19 subspaces
//        a nested DFS walks 1.024e17 four-dimensional prefixes.
// ===========================================================================================
static const int RETAIN = 2, NREPL = 6;

static FILE* LOGF = nullptr;
static std::mutex LOGMTX;
static void logline(const char* fmt, ...) {
    char buf[8192]; va_list a; va_start(a, fmt); vsnprintf(buf, sizeof(buf), fmt, a); va_end(a);
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

// ---------------------------------------------------------------- global search state
static Sol BEST, START, ANCHOR;
static uint32_t KEEPMASK = 0;
static std::mutex BESTMTX;
static int IMPROVEMENTS = 0;
static std::string OUTDIR = ".";
static std::atomic<uint64_t> PREFIXES{0}, LEAVES{0}, PRUNED{0}, UNITS_DONE{0};
static uint64_t UNITS_TOTAL = 0;
static std::atomic<int> CUR_DROP{0};
static std::atomic<int> CUR_I1{0};
static double LAST_REPORT = 0;
static int M1 = 64, M2 = 64, M3 = 64, M4 = 64;
static int RESUME_DROP = 0, RESUME_I1 = 0;

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
static std::string residual_oneline(const Sol& S) {
    std::vector<uint32_t> res; residual_of(S.cls, res);
    std::string s;
    for (size_t i = 0; i < res.size(); ++i) { if (i) s += " | "; s += supp(res[i]); }
    return s.empty() ? std::string("(none)") : s;
}

static void write_sol_dir(const std::string& dir, const Sol& S, const Verify& V) {
    std::error_code ec; std::filesystem::create_directories(dir, ec);
    ResidualInfo r = analyse(V.leftover);
    uint32_t all[R1 + R2];
    for (int i = 0; i < R1; ++i) all[i] = G3[i];
    for (int a = 0; a < R2; ++a) all[R1 + a] = S.pau[a];
    RB.clear();
    ap("Eight additional weight-10 generators (orbit 0006, 6-opt search)\n\n");
    for (int a = 0; a < R2; ++a) ap("h%d = %s   (weight %d)\n", a + 1, pstr(S.pau[a]).c_str(), V.weights[a]);
    ap("\nM_8 = %lld / %zu\nresidual = %lld\n", (long long)S.covered, U3.size(), (long long)S.residual);
    ap("M_final = %lld / %lld  (%.8f %%)\n", (long long)V.det11, (long long)V.all_total,
       100.0 * double(V.det11) / double(V.all_total));
    dump(dir + "/generators.txt");
    RB.clear();
    ap("Complete 11-generator stabilizer  [[14,3]]\n\n");
    for (int i = 0; i < R1; ++i) ap("g%d = %s   (weight %d)\n", i + 1, pstr(G3[i]).c_str(), pwt(G3[i]));
    for (int a = 0; a < R2; ++a) ap("h%d = %s   (weight %d)\n", a + 1, pstr(S.pau[a]).c_str(), V.weights[a]);
    ap("\nH = [ X | Z ]   (11 x 28):\n");
    for (int i = 0; i < R1 + R2; ++i) ap("  %s | %s\n", bstr(px(all[i])).c_str(), bstr(pz(all[i])).c_str());
    ap("\nrank = %d\n", V.rank);
    dump(dir + "/final_matrix.txt");
    RB.clear(); residual_block(r); dump(dir + "/remaining_errors.txt");
    RB.clear();
    ap("VERIFICATION (all %lld weight-<=4 errors regenerated from scratch)\n\n", (long long)V.all_total);
    ap("  1. every new generator has weight exactly 10 : %s\n", V.wt_ok ? "OK" : "FAIL");
    ap("  2. [h_a , g_i] = 0 for all a,i              : %s\n", V.cg_ok ? "OK" : "FAIL");
    ap("  3. [h_a , h_b] = 0 for all a,b              : %s\n", V.ch_ok ? "OK" : "FAIL");
    ap("  4. rank of the 11-generator stabilizer      : %d %s\n", V.rank, V.rank_ok ? "OK" : "FAIL");
    ap("  5. Stage-1-undetected set rebuilt           : %lld %s\n", (long long)V.u3size,
       V.u3_ok ? "OK" : "MISMATCH");
    ap("  6. all %zu targets tested explicitly     : done\n", U3.size());
    ap("  7. union coverage recomputed independently  : %lld %s\n", (long long)V.covered,
       V.score_ok ? "OK (matches the algebraic score)" : "*** MISMATCH ***");
    ap("  8. residual errors listed                   : %lld (see remaining_errors.txt)\n",
       (long long)V.leftover.size());
    ap("\nWeight-4 total: %lld\nWeight-4 detected: %lld\nWeight-4 undetected: %lld\n",
       (long long)V.w4_total, (long long)V.w4_det, (long long)V.w4_und);
    if (V.w4_und == 0) ap("ALL WEIGHT-4 ERRORS DETECTED\n");
    if (S.residual == 0) {
        ap("\n  9. FULL COVERAGE: every one of the %lld weight-<=4 Pauli errors has a non-zero\n"
           "     11-bit syndrome, so N(L) contains no non-zero Pauli of weight <= 4 and the\n"
           "     code has PURE DISTANCE >= 5:  a [[14,3,5]] construction.\n", (long long)V.all_total);
    }
    dump(dir + "/verification.txt");
    RB.clear();
    ap("search: exhaustive-in-the-last-two 6-opt, orbit 0006\n");
    ap("breadth M1=%d M2=%d M3=%d M4=%d   (0 = unlimited)\n", M1, M2, M3, M4);
    ap("prefixes examined: %llu\nleaves scored: %llu\npruned: %llu\n",
       (unsigned long long)PREFIXES.load(), (unsigned long long)LEAVES.load(),
       (unsigned long long)PRUNED.load());
    ap("elapsed at discovery: %.1f s\ntimestamp: %s\n", elapsed(), now_stamp().c_str());
    ap("starting incumbent: M_8 = %lld ; this solution: M_8 = %lld ; improvement %lld\n",
       (long long)START.covered, (long long)S.covered, (long long)(S.covered - START.covered));
    dump(dir + "/metadata.txt");
}

// checkpoint: enough to know exactly which outer blocks are finished
static void save_checkpoint(int dropIdx, int i1done) {
    RB.clear();
    ap("# stab14_6opt checkpoint -- resume with --resume\n");
    ap("version 1\ntimestamp %s\nelapsed %.2f\n", now_stamp().c_str(), elapsed());
    ap("breadth %d %d %d %d\n", M1, M2, M3, M4);
    ap("drop %d\ni1_completed %d\n", dropIdx, i1done);
    ap("units_done %llu\nunits_total %llu\n",
       (unsigned long long)UNITS_DONE.load(), (unsigned long long)UNITS_TOTAL);
    ap("prefixes %llu\nleaves %llu\npruned %llu\n", (unsigned long long)PREFIXES.load(),
       (unsigned long long)LEAVES.load(), (unsigned long long)PRUNED.load());
    ap("best_covered %lld\nbest_residual %lld\n", (long long)BEST.covered, (long long)BEST.residual);
    ap("best_classes");
    for (int a = 0; a < R2; ++a) ap(" %u", BEST.cls[a]);
    ap("\nbest_generators\n");
    for (int a = 0; a < R2; ++a) ap("%s\n", pstr(BEST.pau[a]).c_str());
    ap("best_residual_errors\n");
    { std::vector<uint32_t> res; residual_of(BEST.cls, res);
      for (uint32_t E : res) ap("%s  %s\n", pstr(E).c_str(), supp(E).c_str()); }
    dump(OUTDIR + "/checkpoint.txt");
}
static bool load_checkpoint() {
    FILE* f = fopen((OUTDIR + "/checkpoint.txt").c_str(), "r");
    if (!f) return false;
    char line[512]; int drop = 0, i1 = 0; bool got = false;
    uint32_t cls[R2]; int ncls = 0;
    while (fgets(line, sizeof(line), f)) {
        if (!strncmp(line, "drop ", 5)) { drop = atoi(line + 5); got = true; }
        else if (!strncmp(line, "i1_completed ", 13)) i1 = atoi(line + 13);
        // restore the counters so progress/ETA stay continuous across a restart
        else if (!strncmp(line, "units_done ", 11)) UNITS_DONE.store(strtoull(line + 11, nullptr, 10));
        else if (!strncmp(line, "prefixes ", 9))    PREFIXES.store(strtoull(line + 9, nullptr, 10));
        else if (!strncmp(line, "leaves ", 7))      LEAVES.store(strtoull(line + 7, nullptr, 10));
        else if (!strncmp(line, "pruned ", 7))      PRUNED.store(strtoull(line + 7, nullptr, 10));
        else if (!strncmp(line, "best_classes", 12)) {
            char* p = line + 12; ncls = 0;
            while (ncls < R2) { while (*p == ' ') ++p; if (!*p || *p == '\n') break;
                cls[ncls++] = (uint32_t)strtoul(p, &p, 10); }
        }
    }
    fclose(f);
    if (!got) return false;
    RESUME_DROP = drop; RESUME_I1 = i1;
    if (ncls == R2 && independent_isotropic(cls)) {
        Sol S = make_sol(cls);
        if (better(S, BEST)) BEST = S;
    }
    return true;
}

// ---------------------------------------------------------------- accept an improvement
static bool offer(const uint32_t* cls) {
    if (!independent_isotropic(cls)) return false;
    Sol S = make_sol(cls);
    { std::lock_guard<std::mutex> lk(BESTMTX); if (!better(S, BEST)) return false; }
    Verify V = verify(S);                                  // full 91770-error recheck
    if (!V.wt_ok || !V.cg_ok || !V.ch_ok || !V.rank_ok || !V.u3_ok || !V.score_ok) {
        logline("*** candidate FAILED independent verification -- rejected ***\n");
        return false;
    }
    std::lock_guard<std::mutex> lk(BESTMTX);
    if (!better(S, BEST)) return false;
    Sol prev = BEST; BEST = S; ++IMPROVEMENTS;
    ResidualInfo r = analyse(V.leftover);
    logline("\n============================================================\n"
            "NEW BEST  (%s)\n============================================================\n"
            "previous residual %lld -> new residual %lld\n"
            "M_8      = %lld / %zu   (improvement %lld over the starting %lld)\n"
            "M_final  = %lld / %lld\n"
            "residual weights w1=%lld w2=%lld w3=%lld w4=%lld\n"
            "weight-4: %lld detected, %lld undetected of %lld\n"
            "residual: %s\n"
            "============================================================\n\n",
            now_stamp().c_str(), (long long)prev.residual, (long long)S.residual,
            (long long)S.covered, U3.size(), (long long)(S.covered - START.covered),
            (long long)START.covered, (long long)V.det11, (long long)V.all_total,
            (long long)S.res_w[1], (long long)S.res_w[2], (long long)S.res_w[3], (long long)S.res_w[4],
            (long long)V.w4_det, (long long)V.w4_und, (long long)V.w4_total,
            residual_oneline(S).c_str());
    char d[256];
    snprintf(d, sizeof(d), "%s/best_known/improvement_%03d_M8_%lld", OUTDIR.c_str(),
             IMPROVEMENTS, (long long)S.covered);
    write_sol_dir(std::string(d), S, V);                   // kept, never overwritten
    write_sol_dir(OUTDIR + "/best_known/current", S, V);
    if (V.w4_und == 0)
        logline("ALL WEIGHT-4 ERRORS DETECTED  (%lld / %lld)\n",
                (long long)V.w4_det, (long long)V.w4_total);
    if (S.residual == 0) {
        RB.clear();
        ap("============================================================\n");
        ap("FULL COVERAGE -- [[14,3,5]] CONSTRUCTION FOUND\n");
        ap("============================================================\n\n");
        ap("M_8 = %lld / %zu : every Stage-1-undetected error is detected.\n",
           (long long)S.covered, U3.size());
        ap("M_final = %lld / %lld : every weight-<=4 Pauli error has a non-zero 11-bit\n"
           "syndrome, so N(L) contains no non-zero Pauli of weight <= 4 and the code has\n"
           "PURE DISTANCE >= 5.\n\n", (long long)V.det11, (long long)V.all_total);
        for (int i = 0; i < R1; ++i) ap("g%d = %s\n", i + 1, pstr(G3[i]).c_str());
        for (int a = 0; a < R2; ++a) ap("h%d = %s\n", a + 1, pstr(S.pau[a]).c_str());
        ap("\n============================================================\n");
        dump(OUTDIR + "/FULL_COVERAGE_[[14,3,5]].txt", true);
        if (LOGF) { fputs(RB.c_str(), LOGF); fflush(LOGF); }
    }
    return true;
}

// ---------------------------------------------------------------- 30-minute status report
static void status_report(bool force, const char* where) {
    double e = elapsed();
    { std::lock_guard<std::mutex> lk(LOGMTX);
      if (!force && e - LAST_REPORT < 1800.0) return;
      LAST_REPORT = e; }
    const uint64_t done = UNITS_DONE.load();
    const double frac = UNITS_TOTAL ? double(done) / double(UNITS_TOTAL) : 0.0;
    const double rate = done > 0 ? e / double(done) : 0.0;           // seconds per unit
    const double total_est = rate * double(UNITS_TOTAL);
    const double remain = total_est - e;
    const uint64_t pf = PREFIXES.load();
    logline(
      "\n---------------- STATUS %s ----------------\n"
      "elapsed                : %s (%.1f s)\n"
      "position               : drop-pair %d / 28, outer q1 index %d, %s\n"
      "progress               : %llu / %llu outer blocks = %.4f %% of the CONFIGURED search\n"
      "prefixes examined      : %llu\n"
      "leaves scored          : %llu\n"
      "pruned (bound/filter)  : %llu\n"
      "current best M_8       : %lld / %zu\n"
      "improvement over start : %lld  (start was %lld)\n"
      "residual errors        : %lld  [%s]\n"
      "throughput             : %.2f outer blocks/h, %.3g prefixes/s\n"
      "estimated total run    : %s\n"
      "estimated remaining    : %s\n"
      "checkpoint             : %s/checkpoint.txt\n"
      "-------------------------------------------------------------\n\n",
      now_stamp().c_str(), hms(e).c_str(), e,
      CUR_DROP.load(), CUR_I1.load(), where,
      (unsigned long long)done, (unsigned long long)UNITS_TOTAL, 100.0 * frac,
      (unsigned long long)pf, (unsigned long long)LEAVES.load(), (unsigned long long)PRUNED.load(),
      (long long)BEST.covered, U3.size(),
      (long long)(BEST.covered - START.covered), (long long)START.covered,
      (long long)BEST.residual, residual_oneline(BEST).c_str(),
      rate > 0 ? 3600.0 / rate : 0.0, e > 0 ? double(pf) / e : 0.0,
      rate > 0 ? hms(total_est).c_str() : "(measuring)",
      rate > 0 ? hms(remain).c_str() : "(measuring)",
      OUTDIR.c_str());
}

// ---------------------------------------------------------------- the 6-opt sweep
//  Breadth M1..M4 limit the first four replacements (0 = unlimited); the last two are ALWAYS
//  swept exhaustively with the tau filter.  Redundant bases of the same subspace are not
//  filtered out: doing so would conflict with the psi-sorted pruning order, and visiting a
//  subspace more than once is wasteful but never incorrect.
struct Worker {
    std::vector<uint32_t> psi1, psi2, psi3, psi4, L1, L2, L3, L4, val;
    void alloc(int n) { psi1.resize(n); psi2.resize(n); psi3.resize(n); psi4.resize(n); }
};
static void sweep_q1(const QSpace& Q, uint32_t q1, Worker& W, int64_t need) {
    const int n = Q.size;
    const int64_t base1 = Q.sK + (int64_t)Q.phi[q1];
    for (int w = 0; w < n; ++w) W.psi1[w] = Q.phi[w] + Q.phi[w ^ q1];
    W.L1.clear();
    for (int q = 1; q < n; ++q) {
        if (Q.rep10[q] == NONE || (uint32_t)q == q1) continue;
        if (qsymp(Q, (uint32_t)q, q1)) continue;
        W.L1.push_back((uint32_t)q);
    }
    std::sort(W.L1.begin(), W.L1.end(), [&](uint32_t a, uint32_t b) {
        if (W.psi1[a] != W.psi1[b]) return W.psi1[a] > W.psi1[b]; return a < b; });
    const int m2 = (M2 <= 0) ? (int)W.L1.size() : std::min<int>(M2, (int)W.L1.size());
    uint32_t sp2[2] = {0u, q1};
    for (int i2 = 0; i2 < m2; ++i2) {
        const uint32_t q2 = W.L1[i2];
        const int64_t base2 = base1 + (int64_t)W.psi1[q2];
        for (int w = 0; w < n; ++w) W.psi2[w] = W.psi1[w] + W.psi1[w ^ q2];
        uint32_t sp4[4] = {0u, q1, q2, q1 ^ q2};
        W.L2.clear();
        for (uint32_t q : W.L1) {
            if (qsymp(Q, q, q2)) continue;
            if (q == sp4[1] || q == sp4[2] || q == sp4[3]) continue;
            W.L2.push_back(q);
        }
        std::sort(W.L2.begin(), W.L2.end(), [&](uint32_t a, uint32_t b) {
            if (W.psi2[a] != W.psi2[b]) return W.psi2[a] > W.psi2[b]; return a < b; });
        const int m3 = (M3 <= 0) ? (int)W.L2.size() : std::min<int>(M3, (int)W.L2.size());
        for (int i3 = 0; i3 < m3; ++i3) {
            const uint32_t q3 = W.L2[i3];
            const int64_t base3 = base2 + (int64_t)W.psi2[q3];
            for (int w = 0; w < n; ++w) W.psi3[w] = W.psi2[w] + W.psi2[w ^ q3];
            uint32_t sp8[8];
            for (int m = 0; m < 8; ++m) { uint32_t v = 0;
                if (m & 1) v ^= q1; if (m & 2) v ^= q2; if (m & 4) v ^= q3; sp8[m] = v; }
            W.L3.clear();
            for (uint32_t q : W.L2) {
                if (qsymp(Q, q, q3)) continue;
                bool dep = false; for (int m = 1; m < 8; ++m) if (q == sp8[m]) dep = true;
                if (dep) continue;
                W.L3.push_back(q);
            }
            std::sort(W.L3.begin(), W.L3.end(), [&](uint32_t a, uint32_t b) {
                if (W.psi3[a] != W.psi3[b]) return W.psi3[a] > W.psi3[b]; return a < b; });
            const int m4 = (M4 <= 0) ? (int)W.L3.size() : std::min<int>(M4, (int)W.L3.size());
            for (int i4 = 0; i4 < m4; ++i4) {
                const uint32_t q4 = W.L3[i4];
                const int64_t base4 = base3 + (int64_t)W.psi3[q4];
                uint32_t p4max = 0;
                for (int w = 0; w < n; ++w) { uint32_t s = W.psi3[w] + W.psi3[w ^ q4];
                    W.psi4[w] = s; if (s > p4max) p4max = s; }
                uint32_t sp16[16];
                for (int m = 0; m < 16; ++m) { uint32_t v = 0;
                    if (m & 1) v ^= q1; if (m & 2) v ^= q2; if (m & 4) v ^= q3; if (m & 8) v ^= q4;
                    sp16[m] = v; }
                ++PREFIXES;
                // tau filter: both members of a qualifying pair need psi4 > need-base-2*p4max
                const int64_t tau = need - base4 - 2 * (int64_t)p4max;
                W.L4.clear();
                for (uint32_t q : W.L3) {
                    if ((int64_t)W.psi4[q] <= tau) continue;
                    if (qsymp(Q, q, q4)) continue;
                    bool dep = false; for (int m = 1; m < 16; ++m) if (q == sp16[m]) dep = true;
                    if (dep) continue;
                    W.L4.push_back(q);
                }
                if (W.L4.size() < 2) { ++PRUNED; continue; }
                std::sort(W.L4.begin(), W.L4.end(), [&](uint32_t a, uint32_t b) {
                    if (W.psi4[a] != W.psi4[b]) return W.psi4[a] > W.psi4[b]; return a < b; });
                const size_t mm = W.L4.size();
                W.val.resize(mm);
                for (size_t i = 0; i < mm; ++i) W.val[i] = W.psi4[W.L4[i]];
                for (size_t i = 0; i + 1 < mm; ++i) {
                    const int64_t vi = W.val[i];
                    if (base4 + 2 * vi + (int64_t)p4max <= need) break;
                    const uint32_t a = W.L4[i];
                    for (size_t j = i + 1; j < mm; ++j) {
                        const int64_t vj = W.val[j];
                        if (base4 + vi + vj + (int64_t)p4max <= need) break;
                        const uint32_t b = W.L4[j];
                        if (qsymp(Q, a, b)) continue;
                        const uint32_t x = a ^ b;
                        bool dep = false; for (int m = 1; m < 16; ++m) if (x == sp16[m]) dep = true;
                        if (dep) continue;
                        ++LEAVES;
                        const int64_t s = base4 + vi + vj + (int64_t)W.psi4[x];
                        if (s > need) {
                            uint32_t six[6] = {q1, q2, q3, q4, a, b};
                            uint32_t trial[R2];
                            int t = 0;
                            for (int r = 0; r < R2; ++r) if (KEEPMASK & (1u << r)) trial[r] = ANCHOR.cls[r];
                            for (int r = 0; r < R2; ++r) if (!(KEEPMASK & (1u << r))) {
                                trial[r] = Q.rep10[six[t]]; ++t; }
                            bool bad = false;
                            for (int r = 0; r < R2; ++r) if (trial[r] == NONE) bad = true;
                            if (!bad) offer(trial);
                        }
                    }
                }
            }
        }
    }
}


// ---------------------------------------------------------------- orbit load + verification
static bool load_orbit(int idx, bool quiet) {
    const OrbitDef& O = ORBITS[idx];
    for (int i = 0; i < R1; ++i) G3[i] = parse_pauli(O.g[i]);
    for (int a = 0; a < R2; ++a) H0[a] = parse_pauli(O.h[a]);
    for (int i = 0; i < R1; ++i) if (pwt(G3[i]) != W1) return false;
    for (int i = 0; i < R1; ++i) for (int j = i + 1; j < R1; ++j) if (symp(G3[i], G3[j])) return false;
    build_targets();
    if (U3.size() != 11186) { if (!quiet) printf("|U_3| = %zu != 11186\n", U3.size()); return false; }
    build_basis(); build_cov1(); build_repr10();
    PREF.clear();
    uint32_t cls0[R2];
    for (int a = 0; a < R2; ++a) {
        cls0[a] = cls_of(H0[a]);
        PREF.push_back({cls0[a], H0[a]});
        if (pwt(H0[a]) != W2) return false;
        for (int i = 0; i < R1; ++i) if (symp(H0[a], G3[i])) return false;
    }
    if (!independent_isotropic(cls0)) return false;
    BEST = make_sol(cls0);
    return true;
}
// ---------------------------------------------------------------- driver
static void run_droppair(const QSpace& Q, int dropIdx, int nthreads, double tlimit,
                         int startI1, bool& aborted) {
    const int n = Q.size;
    std::vector<uint32_t> ord;
    for (int q = 1; q < n; ++q) if (Q.rep10[q] != NONE) ord.push_back((uint32_t)q);
    std::sort(ord.begin(), ord.end(), [&](uint32_t a, uint32_t b) {
        if (Q.phi[a] != Q.phi[b]) return Q.phi[a] > Q.phi[b]; return a < b; });
    const int m1 = (M1 <= 0) ? (int)ord.size() : std::min<int>(M1, (int)ord.size());
    logline("[drop-pair %d/28] retained h%d,h%d ; |Q| = %d (dim %d) ; usable q = %zu ;"
            " outer blocks = %d ; starting at %d\n",
            dropIdx + 1, std::countr_zero(KEEPMASK) + 1,
            (int)std::bit_width(KEEPMASK), n, Q.dim, ord.size(), m1, startI1);

    std::atomic<int> next{startI1};
    std::mutex donemtx;
    std::vector<char> done((size_t)m1, 0);
    int completed = startI1;
    std::vector<std::thread> th;
    std::atomic<bool> stop{false};
    for (int t = 0; t < nthreads; ++t) th.emplace_back([&]() {
        Worker W; W.alloc(n);
        for (;;) {
            if (stop.load()) break;
            int i1 = next.fetch_add(1);
            if (i1 >= m1) break;
            CUR_I1.store(i1);
            sweep_q1(Q, ord[i1], W, BEST.span_sum);
            ++UNITS_DONE;
            {
                std::lock_guard<std::mutex> lk(donemtx);
                done[(size_t)i1] = 1;
                while (completed < m1 && done[(size_t)completed]) ++completed;
                save_checkpoint(dropIdx, completed);
            }
            status_report(false, "inside drop-pair sweep");
            if (elapsed() > tlimit) { stop.store(true); break; }
        }
    });
    for (auto& x : th) x.join();
    aborted = stop.load();
}

static void usage() {
    printf("stab14_6opt -- 6-opt search around the orbit-0006 M_8 = 11181 incumbent\n"
           "  --M N            breadth of ALL four leading replacements (0 = unlimited)\n"
           "  --M1/--M2/--M3/--M4 N   set them individually\n"
           "  --threads N      default: all cores\n"
           "  --time-limit S   wall-clock cap in seconds\n"
           "  --resume         continue from checkpoint.txt\n"
           "  --calibrate      measure real throughput and print runtime estimates, then exit\n"
           "  --selftest       verify the incumbent and the algebraic identity, then exit\n"
           "  --out DIR        output root (default: the directory of the executable)\n"
           "  --drop-pair N    run only drop-pair N (0..27)\n");
}

int main(int argc, char** argv) {
    int nthreads = (int)std::max(1u, std::thread::hardware_concurrency());
    double tlimit = 1e18;
    bool resume = false, calib = false, selftest = false;
    int onlyDrop = -1;
    OUTDIR = ".";
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto nxt = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : std::string(); };
        if      (a == "--M")            { int v = atoi(nxt().c_str()); M1 = M2 = M3 = M4 = v; }
        else if (a == "--M1")           M1 = atoi(nxt().c_str());
        else if (a == "--M2")           M2 = atoi(nxt().c_str());
        else if (a == "--M3")           M3 = atoi(nxt().c_str());
        else if (a == "--M4")           M4 = atoi(nxt().c_str());
        else if (a == "--threads")      nthreads = std::max(1, atoi(nxt().c_str()));
        else if (a == "--time-limit")   tlimit = atof(nxt().c_str());
        else if (a == "--resume")       resume = true;
        else if (a == "--calibrate")    calib = true;
        else if (a == "--selftest")     selftest = true;
        else if (a == "--out")          OUTDIR = nxt();
        else if (a == "--drop-pair")    onlyDrop = atoi(nxt().c_str());
        else { usage(); return (a == "--help" || a == "-h") ? 0 : 1; }
    }
    T0 = std::chrono::steady_clock::now();
    std::error_code ec; std::filesystem::create_directories(OUTDIR, ec);
    std::filesystem::create_directories(OUTDIR + "/best_known", ec);
    LOGF = fopen((OUTDIR + "/6opt_search.log").c_str(), resume ? "a" : "w");

    logline("=== stab14_6opt : 6-opt search around the orbit-0006 incumbent ===\n");
    logline("started %s ; threads %d ; output %s\n", now_stamp().c_str(), nthreads, OUTDIR.c_str());

    if (!load_orbit(0, true)) { logline("FATAL: orbit load/verification failed\n"); return 2; }
    START = BEST; ANCHOR = BEST;
    Verify V0 = verify(BEST);
    logline("incumbent verified: M_8 = %lld / %zu, residual %lld (w %lld/%lld/%lld/%lld),"
            " M_final = %lld / %lld\n", (long long)BEST.covered, U3.size(), (long long)BEST.residual,
            (long long)BEST.res_w[1], (long long)BEST.res_w[2], (long long)BEST.res_w[3],
            (long long)BEST.res_w[4], (long long)V0.det11, (long long)V0.all_total);
    logline("  weights OK=%d  [h,g]=%d  [h,h]=%d  rank=%d  brute-force coverage %lld (match=%d)\n",
            V0.wt_ok, V0.cg_ok, V0.ch_ok, V0.rank, (long long)V0.covered, (int)V0.score_ok);
    logline("  weight-4: %lld detected, %lld undetected of %lld\n",
            (long long)V0.w4_det, (long long)V0.w4_und, (long long)V0.w4_total);
    logline("  residual: %s\n", residual_oneline(BEST).c_str());
    if (!V0.wt_ok || !V0.cg_ok || !V0.ch_ok || !V0.rank_ok || !V0.score_ok) {
        logline("FATAL: incumbent does not verify\n"); return 3; }
    if (selftest) { logline("\nSELFTEST PASSED: incumbent verifies and the algebraic score matches"
                            " the brute-force union.\n"); return 0; }

    // ---------------------------------------------------------------- calibration
    if (calib) {
        logline("\n=== CALIBRATION ===\n");
        KEEPMASK = (1u << 6) | (1u << 7);
        uint32_t keep[RETAIN]; int nk = 0;
        for (int r = 0; r < R2; ++r) if (KEEPMASK & (1u << r)) keep[nk++] = ANCHOR.cls[r];
        QSpace Q;
        double t0 = elapsed();
        if (!build_qspace(keep, nk, Q)) { logline("quotient build failed\n"); return 4; }
        double tb = elapsed() - t0;
        int usable = 0; for (int q = 1; q < Q.size; ++q) if (Q.rep10[q] != NONE) ++usable;
        logline("quotient: dim %d, |Q| = %d, build %.2f s, usable q = %d (%.1f%%), phimax = %u,"
                " S(K) = %lld\n", Q.dim, Q.size, tb, usable, 100.0 * usable / (Q.size - 1),
                Q.phimax, (long long)Q.sK);
        std::vector<uint32_t> ord;
        for (int q = 1; q < Q.size; ++q) if (Q.rep10[q] != NONE) ord.push_back((uint32_t)q);
        std::sort(ord.begin(), ord.end(), [&](uint32_t a, uint32_t b) { return Q.phi[a] > Q.phi[b]; });
        // time ONE outer block at the configured breadth
        Worker W; W.alloc(Q.size);
        t0 = elapsed();
        uint64_t p0 = PREFIXES.load();
        sweep_q1(Q, ord[0], W, BEST.span_sum);
        double t1 = elapsed() - t0;
        uint64_t pf = PREFIXES.load() - p0;
        logline("one outer block (q1 = best by phi) at M2=%d M3=%d M4=%d : %.3f s,"
                " %llu prefixes, %llu leaves\n", M2, M3, M4, t1, (unsigned long long)pf,
                (unsigned long long)LEAVES.load());
        const double perPrefix = pf ? t1 / double(pf) : 0.0;
        logline("  -> %.1f us per (q1,q2,q3,q4) prefix, single-threaded\n", perPrefix * 1e6);
        const int m1 = (M1 <= 0) ? (int)ord.size() : std::min<int>(M1, (int)ord.size());
        double cfg = 28.0 * double(m1) * t1 / std::max(1, nthreads);
        logline("\nCONFIGURED RUN  (M1=%d M2=%d M3=%d M4=%d, %d threads)\n", M1, M2, M3, M4, nthreads);
        logline("  28 drop-pairs x %d outer blocks x %.3f s / %d threads = %s\n",
                m1, t1, nthreads, hms(cfg).c_str());
        // exhaustive extrapolation
        const double PREF_EXH = 1.024e17;      // 4-dim isotropic prefixes over all 28 drop-pairs
        const double SUB_EXH  = 1.368e19;      // 6-dim isotropic subspaces over all 28 drop-pairs
        double exh = PREF_EXH * perPrefix / std::max(1, nthreads);
        logline("\nEXHAUSTIVE 6-OPT (M = 0 everywhere)\n");
        logline("  6-dim isotropic subspaces of Q over all 28 drop-pairs : %.3g\n", SUB_EXH);
        logline("  4-dim prefixes a nested DFS must walk                 : %.3g\n", PREF_EXH);
        logline("  at the measured %.1f us/prefix on %d threads          : %.3g s = %.3g years\n",
                perPrefix * 1e6, nthreads, exh, exh / 3.15576e7);
        logline("  => EXHAUSTIVE 6-OPT IS NOT RUNNABLE.  Use a finite breadth.\n");
        double frac = double(m1) * double(M2 <= 0 ? 1 : M2) * double(M3 <= 0 ? 1 : M3)
                    * double(M4 <= 0 ? 1 : M4) * 28.0 / PREF_EXH;
        logline("  the configured run covers about %.3g of the exhaustive prefix space\n", frac);
        logline("=== END CALIBRATION ===\n");
        return 0;
    }

    // ---------------------------------------------------------------- the search
    int startDrop = 0, startI1 = 0;
    if (resume && load_checkpoint()) {
        startDrop = RESUME_DROP; startI1 = RESUME_I1;
        logline("RESUMED from checkpoint: drop-pair %d, %d outer blocks already complete;"
                " best M_8 = %lld\n", startDrop, startI1, (long long)BEST.covered);
    }
    // work accounting
    {
        KEEPMASK = (1u << 6) | (1u << 7);
        uint32_t keep[RETAIN]; int nk = 0;
        for (int r = 0; r < R2; ++r) if (KEEPMASK & (1u << r)) keep[nk++] = ANCHOR.cls[r];
        QSpace Q; build_qspace(keep, nk, Q);
        int usable = 0; for (int q = 1; q < Q.size; ++q) if (Q.rep10[q] != NONE) ++usable;
        int m1 = (M1 <= 0) ? usable : std::min(M1, usable);
        UNITS_TOTAL = (uint64_t)28 * (uint64_t)m1;
    }
    logline("\nbreadth M1=%d M2=%d M3=%d M4=%d (0 = unlimited); the last TWO replacements are\n"
            "always swept exhaustively.  Total outer blocks = %llu.\n",
            M1, M2, M3, M4, (unsigned long long)UNITS_TOTAL);
    status_report(true, "start");

    bool aborted = false;
    int dropIdx = 0;
    for (int d1 = 0; d1 < R2 && !aborted; ++d1)
    for (int d2 = d1 + 1; d2 < R2 && !aborted; ++d2, ++dropIdx) {
        if (dropIdx < startDrop) continue;
        if (onlyDrop >= 0 && dropIdx != onlyDrop) continue;
        CUR_DROP.store(dropIdx);
        KEEPMASK = (1u << d1) | (1u << d2);
        uint32_t keep[RETAIN]; int nk = 0;
        for (int r = 0; r < R2; ++r) if (KEEPMASK & (1u << r)) keep[nk++] = ANCHOR.cls[r];
        QSpace Q;
        if (!build_qspace(keep, nk, Q)) { logline("drop-pair %d: quotient build FAILED\n", dropIdx); continue; }
        int s1 = (dropIdx == startDrop) ? startI1 : 0;
        run_droppair(Q, dropIdx, nthreads, tlimit, s1, aborted);
        save_checkpoint(dropIdx + 1, 0);
        status_report(true, "drop-pair complete");
    }

    // ---------------------------------------------------------------- final report
    Verify VB = verify(BEST);
    ResidualInfo rb = analyse(VB.leftover);
    RB.clear();
    ap("============================================================\n");
    ap("6-OPT SEARCH RESULT\n");
    ap("============================================================\n\n");
    ap("finished %s after %s\n\n", now_stamp().c_str(), hms(elapsed()).c_str());
    ap("breadth: M1=%d M2=%d M3=%d M4=%d  (0 = unlimited); last two always exhaustive\n", M1, M2, M3, M4);
    ap("outer blocks completed: %llu / %llu\n", (unsigned long long)UNITS_DONE.load(),
       (unsigned long long)UNITS_TOTAL);
    ap("prefixes examined: %llu   leaves scored: %llu   pruned: %llu\n\n",
       (unsigned long long)PREFIXES.load(), (unsigned long long)LEAVES.load(),
       (unsigned long long)PRUNED.load());
    ap("starting M_8 : %lld / %zu\n", (long long)START.covered, U3.size());
    ap("best    M_8 : %lld / %zu\n", (long long)BEST.covered, U3.size());
    ap("improvement : %lld\n\n", (long long)(BEST.covered - START.covered));
    ap("M_final = %lld / %lld  (%.8f %%)\n", (long long)VB.det11, (long long)VB.all_total,
       100.0 * double(VB.det11) / double(VB.all_total));
    ap("Weight-4 total: %lld\nWeight-4 detected: %lld\nWeight-4 undetected: %lld\n",
       (long long)VB.w4_total, (long long)VB.w4_det, (long long)VB.w4_und);
    if (VB.w4_und == 0) ap("ALL WEIGHT-4 ERRORS DETECTED\n");
    ap("\n"); residual_block(rb); ap("\n");
    if (BEST.covered > START.covered) {
        ap("IMPROVEMENT FOUND: M_8 = %lld.\n", (long long)BEST.covered);
        if (BEST.residual == 0)
            ap("M_8 = %zu : FULL COVERAGE -- a [[14,3,5]] construction (see the certificate file).\n",
               U3.size());
    } else if (aborted || UNITS_DONE.load() < UNITS_TOTAL) {
        ap("Search INCOMPLETE (time limit or --drop-pair restriction): %llu of %llu outer\n"
           "blocks done.  M_8 = %lld stands as best known.  NO optimality claim.\n",
           (unsigned long long)UNITS_DONE.load(), (unsigned long long)UNITS_TOTAL,
           (long long)BEST.covered);
    } else if (M1 <= 0 && M2 <= 0 && M3 <= 0 && M4 <= 0) {
        ap("Exhaustive 6-opt search found no improvement over the incumbent.  M_8 = %lld.\n",
           (long long)BEST.covered);
        ap("This proves 6-OPTIMALITY OF THIS INCUMBENT ONLY.  It is NOT a proof that 11181 is\n"
           "globally optimal over all eight-generator configurations.\n");
    } else {
        ap("The configured (breadth-limited) 6-opt search found no improvement.  M_8 = %lld.\n",
           (long long)BEST.covered);
        ap("This is NOT an exhaustive 6-opt result: with M1=%d M2=%d M3=%d M4=%d only part of\n"
           "the 6-opt neighbourhood was visited.  NO optimality claim of any kind.\n",
           M1, M2, M3, M4);
    }
    ap("============================================================\n");
    dump(OUTDIR + "/final_report.txt", true);
    if (LOGF) { fputs(RB.c_str(), LOGF); fflush(LOGF); fclose(LOGF); }
    return 0;
}
