// ===========================================================================================
//  stab14_2stage.cpp -- TWO SEPARATE optimisations on 14 qubits.
//
//    STAGE 1 : three generators of Pauli weight 8, commuting, rank 3, maximising the number
//              of the 91770 Pauli errors of weight 1..4 with non-zero 3-bit syndrome.
//    STAGE 2 : the Stage-1 optimum is FROZEN; eight further generators of Pauli weight 10,
//              commuting with it and with each other, rank 11 overall, maximising coverage
//              of ONLY the errors U_3 that Stage 1 missed.
//
//  The two objectives are never mixed: Stage 2 scores against U_3 alone, and the final
//  11-generator numbers are reported, not optimised.
//
//  -----------------------------------------------------------------------------------------
//  COMMON ALGEBRA (used by both stages)
//  -----------------------------------------------------------------------------------------
//  Pauli mod phase = (x|z) in F_2^28, symplectic form <a,b> = x_a.z_b + z_a.x_b.
//  For an r-dimensional isotropic subspace L (a stabilizer group of rank r) the number of
//  weight-w errors with zero syndrome is given by the quaternary Krawtchouk transform
//
//        #{E : wt(E)=w, s(E)=0}  =  (1/2^r) SUM_{h in L} Ew[wt(h)][w],
//        Ew[m][w] = SUM_k C(m,k) (-1)^k C(14-m, w-k) 3^{w-k}.
//
//  Proof: the Walsh transform over the syndrome space of the per-qubit error-counting
//  function is F_j(h) = 3 when h acts trivially on qubit j and -1 otherwise, so the
//  generating function of the zero-syndrome count is the elementary symmetric polynomial
//  e_w(F_1,...,F_14), which depends only on wt(h).  Writing P(m) = sum_{w=1..4} Ew[m][w],
//
//        #undetected(weight 1..4)  =  (1/2^r) [ 91770 + SUM_{h in L, h != 0} P(wt h) ].
//
//  So BOTH stages are ultimately functions of the weight enumerator of the stabilizer group.
//  P(m), m = 0..14:
//     91770 57914 34154 18314 8474 2970 394 -406 -326 -6 170 74 -166 -166 714
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
#include <set>
#include <map>
#include <filesystem>

// ------------------------------------------------------------------ global problem constants
constexpr int NQ    = 14;                  // physical qubits
constexpr int R1    = 3;                   // Stage-1 generators
constexpr int R2    = 8;                   // Stage-2 generators
constexpr int MAXW  = 4;                   // error weight considered
constexpr int W1    = 8;                   // Stage-1 generator weight
constexpr int W2    = 10;                  // Stage-2 generator weight
constexpr uint32_t QM   = (1u << NQ) - 1u;
constexpr uint32_t NONE = 0xFFFFFFFFu;
constexpr int64_t  INFP = (int64_t)1 << 40;

// Pauli packed as x | (z << 16)
static inline uint32_t pk(uint32_t x, uint32_t z) { return x | (z << 16); }
static inline uint32_t px(uint32_t p) { return p & QM; }
static inline uint32_t pz(uint32_t p) { return (p >> 16) & QM; }
static inline int pwt(uint32_t p) { return std::popcount((p | (p >> 16)) & QM); }
static inline int par(unsigned a) { return std::popcount(a) & 1; }
static inline int symp(uint32_t a, uint32_t b) {
    return par((px(a) & pz(b)) ^ (pz(a) & px(b)));
}
static const char PAULI_CH[2][2] = { {'I','Z'}, {'X','Y'} };      // [xbit][zbit]
static std::string pstr(uint32_t p) {
    std::string s; s.reserve(NQ);
    for (int j = 0; j < NQ; ++j) s += PAULI_CH[(px(p) >> j) & 1][(pz(p) >> j) & 1];
    return s;
}
static std::string bstr(uint32_t v) {
    std::string s; for (int j = 0; j < NQ; ++j) s += char('0' + ((v >> j) & 1)); return s;
}

// ------------------------------------------------------------------ Krawtchouk / P tables
static int64_t Cbin[40][40], Ew[NQ + 1][MAXW + 1], Ptab[NQ + 1], PbyC[NQ + 1];
static int64_t TOTW[MAXW + 1], TOTAL_ERRORS = 0;
static void init_tables() {
    for (int n = 0; n < 40; ++n) { Cbin[n][0] = 1;
        for (int k = 1; k <= n; ++k) Cbin[n][k] = Cbin[n-1][k-1] + (k <= n-1 ? Cbin[n-1][k] : 0); }
    auto C = [](int n, int k) { return (k < 0 || k > n || n < 0) ? (int64_t)0 : Cbin[n][k]; };
    int64_t p3[MAXW + 1]; p3[0] = 1; for (int w = 1; w <= MAXW; ++w) p3[w] = p3[w-1] * 3;
    for (int w = 0; w <= MAXW; ++w) { TOTW[w] = C(NQ, w) * p3[w]; if (w) TOTAL_ERRORS += TOTW[w]; }
    for (int m = 0; m <= NQ; ++m) {
        Ptab[m] = 0;
        for (int w = 0; w <= MAXW; ++w) {
            int64_t s = 0;
            for (int k = 0; k <= w; ++k) { int64_t t = C(m,k) * C(NQ-m, w-k) * p3[w-k]; s += (k&1) ? -t : t; }
            Ew[m][w] = s; if (w) Ptab[m] += s;
        }
    }
    for (int c = 0; c <= NQ; ++c) PbyC[c] = Ptab[NQ - c];
}

// ------------------------------------------------------------------ error enumeration
// Calls f(E, w) for every non-identity Pauli error of qubit weight 1..4 (91770 of them).
template <class F> static void for_each_error(F&& f) {
    static const uint32_t PX3[3] = {1,1,0}, PZ3[3] = {0,1,1};   // X, Y, Z
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

static std::chrono::steady_clock::time_point T0;
static double elapsed() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count(); }

// small formatted-report buffer
static std::string RB;
static void ap(const char* fmt, ...) {
    char buf[8192]; va_list a; va_start(a, fmt); vsnprintf(buf, sizeof(buf), fmt, a); va_end(a);
    RB += buf;
}
static std::string OUTDIR = ".";
static void dump(const char* fname, bool echo = false) {
    if (echo) { fputs(RB.c_str(), stdout); fflush(stdout); }
    std::string path = OUTDIR + "/" + fname;
    FILE* f = fopen(path.c_str(), "w"); if (f) { fputs(RB.c_str(), f); fclose(f); }
}

// ===========================================================================================
//  ==========================  S T A G E   1  ==============================================
//  Three weight-8 generators, exhaustive and PROVABLY optimal.
//
//  COLUMN-TYPE REDUCTION (identical in spirit to the 4-generator solver, now over F_2^3).
//  Qubit j contributes (x_j, z_j) in F_2^3 x F_2^3.  A single-qubit Clifford acts as
//  GL(2,2) on that pair; it preserves the row weights, preserves the commutation form
//  exactly ((ax+bz)(cx+dz)^T + sym = (ad+bc)(xz^T+zx^T)), and permutes the three
//  single-qubit syndromes {z, x, x+z} setwise -- and the objective sums over all three
//  Paulis of every support qubit, so it is invariant.  The complete invariant of a column
//  is therefore the SUBSPACE W_j = span{x_j,z_j} <= F_2^3, dim 0,1,2: only 1+7+7 = 15 types.
//
//  With U_j = W_j^perp and c(chi) = #{j : chi in U_j}, the weight of the group element
//  g_chi is m(chi) = 14 - c(chi), and
//        weight 8 for the three generators  <=>  c(e_i) = 6, i = 1,2,3
//        mutual commutation                 <=>  XOR_j plucker(W_j) = 0 in F_2^3
//        rank 3                             <=>  c(chi) <= 13 for all chi != 0
//        objective                          <=>  minimise SUM_{chi != 0} P(14 - c(chi)).
//  Multisets of 14 types out of 15 number C(28,14) = 40 116 600, and the four hard
//  constraints plus the bound below cut that to a few thousand nodes: Stage 1 is exhausted
//  in milliseconds, so its optimum is PROVEN.
// ===========================================================================================
constexpr int NCHI1 = 1 << R1;                 // 8 syndromes
constexpr int NTYPE1_MAX = 15;
static int TARGET_C1 = NQ - W1;                // required c(e_i) = 6

struct Type1 {
    uint8_t wmask, umask;      // subsets of F_2^3 (8 bits)
    uint8_t dimW, sig, pl, degU;
    uint8_t bx, bz;
    uint8_t chi[7];
};
static Type1 TY1[NTYPE1_MAX];
static int NTYPE1 = 0;
static const int P1a[3] = {0,0,1}, P1b[3] = {1,2,2};      // the 3 generator pairs
static const int EPT1[R1] = {1, 2, 4};                    // e_1,e_2,e_3 as elements of F_2^3

static uint8_t plucker3(int x, int z) {
    uint8_t p = 0;
    for (int b = 0; b < 3; ++b) { int i = P1a[b], k = P1b[b];
        int v = (((x >> i) & 1) & ((z >> k) & 1)) ^ (((z >> i) & 1) & ((x >> k) & 1));
        if (v) p |= uint8_t(1u << b); }
    return p;
}
static void build_types1() {
    std::vector<uint8_t> masks; masks.push_back(1u);
    for (int v = 1; v < 8; ++v) masks.push_back(uint8_t(1u | (1u << v)));
    for (int v = 1; v < 8; ++v) for (int w = v + 1; w < 8; ++w) {
        uint8_t m = uint8_t(1u | (1u << v) | (1u << w) | (1u << (v ^ w)));
        if (std::find(masks.begin(), masks.end(), m) == masks.end()) masks.push_back(m); }
    std::vector<Type1> tmp;
    for (uint8_t wm : masks) {
        Type1 t{}; t.wmask = wm;
        t.dimW = uint8_t(std::countr_zero(unsigned(std::popcount(wm))));      // |W| = 2^dim
        uint8_t um = 0;
        for (int chi = 0; chi < 8; ++chi) { bool ok = true;
            for (int w = 0; w < 8; ++w) if ((wm >> w) & 1) if (par(unsigned(chi & w))) { ok = false; break; }
            if (ok) um = uint8_t(um | (1u << chi)); }
        t.umask = um; t.degU = uint8_t(std::popcount(um) - 1);
        t.sig = 0; for (int i = 0; i < R1; ++i) if ((um >> EPT1[i]) & 1) t.sig |= uint8_t(1u << i);
        int n = 0; for (int chi = 1; chi < 8; ++chi) if ((um >> chi) & 1) t.chi[n++] = uint8_t(chi);
        int b[2] = {0,0}, nb = 0;
        for (int v = 1; v < 8 && nb < 2; ++v) if ((wm >> v) & 1) b[nb++] = v;
        t.bx = uint8_t(b[0]); t.bz = uint8_t(b[1]);
        t.pl = (t.dimW == 2) ? plucker3(b[0], b[1]) : uint8_t(0);
        tmp.push_back(t);
    }
    std::sort(tmp.begin(), tmp.end(), [](const Type1& a, const Type1& b) {
        if (a.dimW != b.dimW) return a.dimW < b.dimW;
        int pa = std::popcount(a.sig), pb = std::popcount(b.sig);
        if (pa != pb) return pa > pb;
        return a.wmask < b.wmask; });
    NTYPE1 = (int)tmp.size();
    for (int i = 0; i < NTYPE1; ++i) TY1[i] = tmp[i];
}

// S_3 acting on the three generators permutes coordinates of F_2^3, hence the column types.
static int S3act[6][NTYPE1_MAX]; static int NS3 = 0;
static void build_s3() {
    int p[3] = {0,1,2}; NS3 = 0;
    do {
        int map8[8];
        for (int v = 0; v < 8; ++v) { int r = 0;
            for (int i = 0; i < 3; ++i) if ((v >> i) & 1) r |= 1 << p[i];
            map8[v] = r; }
        for (int t = 0; t < NTYPE1; ++t) {
            uint8_t nm = 0;
            for (int v = 0; v < 8; ++v) if ((TY1[t].wmask >> v) & 1) nm = uint8_t(nm | (1u << map8[v]));
            int found = -1; for (int u = 0; u < NTYPE1; ++u) if (TY1[u].wmask == nm) { found = u; break; }
            S3act[NS3][t] = found;
        }
        ++NS3;
    } while (std::next_permutation(p, p + 3));
}

// ---- bound tables (identical construction to the 4-generator solver) ----------------------
constexpr int64_t RSCALE = 1 << 20;
static int64_t Pmin1[NQ + 2][NQ + 1], Rmax1[NQ + 2][NQ + 1];
static uint8_t sigAvail1[NTYPE1_MAX + 1], sigPopMax1[NTYPE1_MAX + 1], degMax1[NTYPE1_MAX + 1];
static void init_bounds1() {
    for (int c = 0; c <= NQ; ++c) for (int r = 0; r <= NQ; ++r) {
        int64_t best = (c <= NQ - 1) ? PbyC[c] : INFP, rat = 0;
        int dmax = std::min(r, (NQ - 1) - c);
        for (int d = 1; d <= dmax; ++d) {
            best = std::min(best, PbyC[c + d]);
            int64_t gain = PbyC[c] - PbyC[c + d];
            if (gain > 0) rat = std::max(rat, (gain * RSCALE + d - 1) / d);
        }
        Pmin1[c][r] = best; Rmax1[c][r] = rat;
    }
    sigAvail1[NTYPE1] = 0; sigPopMax1[NTYPE1] = 0; degMax1[NTYPE1] = 0;
    for (int t = NTYPE1 - 1; t >= 0; --t) {
        sigAvail1[t]  = uint8_t(sigAvail1[t + 1] | TY1[t].sig);
        sigPopMax1[t] = uint8_t(std::max<int>(sigPopMax1[t + 1], std::popcount(TY1[t].sig)));
        degMax1[t]    = uint8_t(std::max<int>(degMax1[t + 1], TY1[t].degU));
    }
}
static inline int64_t sumP1(const uint8_t* c) {
    int64_t s = 0; for (int chi = 1; chi < NCHI1; ++chi) s += PbyC[c[chi]]; return s;
}
// Admissible lower bound on the final sum of P over the 7 non-identity group elements.
static inline int64_t node_bound1(const uint8_t* c, int r, int t0) {
    int need = 0;
    for (int i = 0; i < R1; ++i) {
        int d = TARGET_C1 - c[EPT1[i]];
        if (d < 0 || d > r) return INFP;
        if (d > 0 && !((sigAvail1[t0] >> i) & 1)) return INFP;
        need += d;
    }
    if (r > 0 && need > r * (int)sigPopMax1[t0]) return INFP;
    const int64_t base = R1 * PbyC[TARGET_C1];
    int64_t s1 = 0, s2 = 0, rat = 0;
    for (int chi = 1; chi < NCHI1; ++chi) {
        if (chi == 1 || chi == 2 || chi == 4) continue;          // the three generator points
        int cc = c[chi];
        if (cc > NQ - 1) return INFP;                            // rank would drop below 3
        int64_t pm = Pmin1[cc][r]; if (pm >= INFP) return INFP;
        s1 += pm; s2 += PbyC[cc]; rat = std::max(rat, Rmax1[cc][r]);
    }
    int64_t Bmax = int64_t(degMax1[t0]) * r - need;
    if (Bmax < 0) return INFP;
    int64_t lb2 = s2 - (Bmax * rat + RSCALE - 1) / RSCALE;
    return base + std::max(s1, lb2);
}

static inline void compute_c1(const uint8_t* types, int n, uint8_t* c, uint8_t& plx) {
    std::memset(c, 0, NCHI1); plx = 0;
    for (int j = 0; j < n; ++j) { const Type1& t = TY1[types[j]];
        for (int k = 0; k < t.degU; ++k) c[t.chi[k]]++;
        plx ^= t.pl; }
}
static void canon1(const uint8_t* types, int n, uint8_t* out) {
    uint8_t buf[NQ], best[NQ]; bool first = true;
    for (int p = 0; p < NS3; ++p) {
        const int* A = S3act[p];
        for (int j = 0; j < n; ++j) buf[j] = uint8_t(A[types[j]]);
        std::sort(buf, buf + n);
        if (first || std::lexicographical_compare(buf, buf + n, best, best + n)) {
            std::copy(buf, buf + n, best); first = false; } }
    std::copy(best, best + n, out);
}
static inline bool lexmin1(const uint8_t* types, int n) {
    uint8_t buf[NQ];
    for (int p = 1; p < NS3; ++p) { const int* A = S3act[p];
        for (int j = 0; j < n; ++j) buf[j] = uint8_t(A[types[j]]);
        std::sort(buf, buf + n);
        if (std::lexicographical_compare(buf, buf + n, types, types + n)) return false; }
    return true;
}

// ---- Stage-1 incumbent + exhaustive DFS ---------------------------------------------------
static int64_t  S1_bestSumP = INFP;
static uint8_t  S1_bestTypes[NQ];
static bool     S1_have = false;
static uint64_t S1_nodes = 0, S1_pruned = 0, S1_leaves = 0;

static int S1_bestIdle = NQ + 1;
// Every column multiset attaining the optimum, in S_3-canonical form.  Because the state is
// already a MULTISET of SUBSPACE types, this set is the exact quotient of the optimal
// Stage-1 solutions by  < S_14 qubit permutations , (S_3 local Cliffords)^14 , S_3 generator
// permutations > -- the largest group that provably preserves generator weight, mutual
// commutation, rank and the whole weight-<=4 detection objective.
static std::set<std::array<uint8_t, NQ>> S1_ALL;
static uint64_t S1_RAWOPT = 0;             // optimal multisets before the S_3 quotient
static bool S1_USE_SYM = true;             // apply the S_3 lex-min prefix test
// Tie-break among Stage-1 optima (all have the identical, optimal objective): prefer
// representatives that leave no physical qubit completely untouched -- an idle qubit is
// dead weight for the code that Stage 2 will build on -- then take the S_3-lex-minimum so
// that the frozen Stage-1 solution is reproducible.
static void s1_offer(const uint8_t* types) {
    uint8_t c[NCHI1], plx; compute_c1(types, NQ, c, plx);
    if (plx) return;
    for (int i = 0; i < R1; ++i) if (c[EPT1[i]] != TARGET_C1) return;
    for (int chi = 1; chi < NCHI1; ++chi) if (c[chi] > NQ - 1) return;
    int64_t sp = sumP1(c);
    uint8_t cn[NQ]; canon1(types, NQ, cn);
    int idle = 0; for (int j = 0; j < NQ; ++j) if (TY1[cn[j]].dimW == 0) ++idle;
    // ---- collect the whole optimal family, not just one representative -------------------
    // The DFS only prunes subtrees whose admissible lower bound is strictly worse than the
    // current incumbent, so no optimal leaf can ever be cut: this single pass sees them all.
    if (!S1_have || sp < S1_bestSumP) { S1_ALL.clear(); S1_RAWOPT = 0; }
    if (!S1_have || sp <= S1_bestSumP) {
        std::array<uint8_t, NQ> key; std::copy(cn, cn + NQ, key.begin());
        S1_ALL.insert(key); ++S1_RAWOPT;
    }
    if (S1_have) {
        if (sp > S1_bestSumP) return;
        if (sp == S1_bestSumP) {
            if (idle > S1_bestIdle) return;
            if (idle == S1_bestIdle &&
                !std::lexicographical_compare(cn, cn + NQ, S1_bestTypes, S1_bestTypes + NQ)) return;
        }
    }
    std::copy(cn, cn + NQ, S1_bestTypes); S1_bestSumP = sp; S1_bestIdle = idle; S1_have = true;
}
static uint8_t s1_types[NQ]; static uint8_t s1_c[NCHI1]; static uint8_t s1_plx;
static void s1_dfs(int depth, int t0) {
    ++S1_nodes;
    const int r = NQ - depth;
    int64_t b = node_bound1(s1_c, r, t0);
    // keep ties so that the reported representative is the S_3-canonical minimum
    if (b >= INFP || (S1_have && b > S1_bestSumP)) { ++S1_pruned; return; }
    if (r == 0) { ++S1_leaves; if (s1_plx == 0) s1_offer(s1_types); return; }
    for (int t = t0; t < NTYPE1; ++t) {
        const Type1& ti = TY1[t];
        bool bad = false;
        for (int k = 0; k < ti.degU; ++k) if (++s1_c[ti.chi[k]] > NQ - 1) bad = true;
        s1_plx = uint8_t(s1_plx ^ ti.pl); s1_types[depth] = uint8_t(t);
        if (!bad && (!S1_USE_SYM || lexmin1(s1_types, depth + 1))) s1_dfs(depth + 1, t);
        s1_plx = uint8_t(s1_plx ^ ti.pl);
        for (int k = 0; k < ti.degU; ++k) s1_c[ti.chi[k]]--;
    }
}

// ---- Stage-1 reconstruction, independent verification and report --------------------------
static uint32_t G1[R1];                      // the frozen Stage-1 generators
static std::vector<uint32_t> U3;             // Stage-1 undetected errors (the Stage-2 targets)
static int64_t U3_BY_W[MAXW + 1] = {0};

static void s1_build_matrix(const uint8_t* types, uint32_t* G) {
    for (int i = 0; i < R1; ++i) G[i] = 0;
    for (int j = 0; j < NQ; ++j) {
        const Type1& t = TY1[types[j]];
        for (int i = 0; i < R1; ++i) {
            uint32_t x = ((t.bx >> i) & 1u) << j, z = ((t.bz >> i) & 1u) << j;
            G[i] |= pk(x, z);
        }
    }
}

struct S1Verify {
    int  weights[R1]{}; int rank = 0; int comm[R1][R1]{};
    bool ok_w = true, ok_comm = true, ok_rank = true;
    int64_t synd[NCHI1]{}, det = 0, und = 0;
    int64_t det_w[MAXW + 1]{}, und_w[MAXW + 1]{};
};
// completely independent: rebuilds all 91770 errors and evaluates every 3-bit syndrome
static S1Verify s1_verify(const uint32_t* G, std::vector<uint32_t>* undetected) {
    S1Verify V;
    for (int i = 0; i < R1; ++i) V.weights[i] = pwt(G[i]);
    for (int i = 0; i < R1; ++i) if (V.weights[i] != W1) V.ok_w = false;
    for (int i = 0; i < R1; ++i) for (int k = 0; k < R1; ++k) {
        V.comm[i][k] = symp(G[i], G[k]);
        if (i != k && V.comm[i][k]) V.ok_comm = false; }
    { uint32_t rows[R1]; for (int i = 0; i < R1; ++i) rows[i] = px(G[i]) | (pz(G[i]) << NQ);
      int rk = 0;
      for (int b = 0; b < 2 * NQ && rk < R1; ++b) {
          int p = -1; for (int i = rk; i < R1; ++i) if ((rows[i] >> b) & 1) { p = i; break; }
          if (p < 0) continue;
          std::swap(rows[rk], rows[p]);
          for (int i = 0; i < R1; ++i) if (i != rk && ((rows[i] >> b) & 1)) rows[i] ^= rows[rk];
          ++rk; }
      V.rank = rk; V.ok_rank = (rk == R1); }
    if (undetected) undetected->clear();
    for_each_error([&](uint32_t E, int w) {
        int s = 0; for (int i = 0; i < R1; ++i) s |= symp(G[i], E) << i;
        V.synd[s]++;
        if (s) { V.det++; V.det_w[w]++; }
        else { V.und++; V.und_w[w]++; if (undetected) undetected->push_back(E); }
    });
    return V;
}

static void s1_report(const uint32_t* G, const S1Verify& V, bool proven) {
    RB.clear();
    ap("============================================================\n");
    ap("STAGE 1 - THREE WEIGHT-8 GENERATORS\n");
    ap("============================================================\n\n");
    ap("Status:\nGlobal optimum: %s\n\n", proven ? "PROVEN" : "NOT PROVEN");
    ap("Total weight <=4 errors: %lld\n\n", (long long)TOTAL_ERRORS);
    ap("Detected:   %lld\n", (long long)V.det);
    ap("Undetected: %lld\n\n", (long long)V.und);
    ap("Detection percentage: %.8f %%\n\n", 100.0 * double(V.det) / double(TOTAL_ERRORS));
    for (int w = 1; w <= MAXW; ++w) {
        ap("Weight %d:\nDetected / Total : %lld / %lld\nUndetected       : %lld\n\n",
           w, (long long)V.det_w[w], (long long)TOTW[w], (long long)V.und_w[w]);
    }
    ap("Generators:\n");
    for (int i = 0; i < R1; ++i) ap("g%d = %s\n", i + 1, pstr(G[i]).c_str());
    ap("\nWeights:\n");
    for (int i = 0; i < R1; ++i) ap("g%d: %d\n", i + 1, V.weights[i]);
    ap("\nRank:\n%d  (target %d)\n", V.rank, R1);
    ap("\nCommutation (0 = commute):\n");
    for (int i = 0; i < R1; ++i) { ap("  ");
        for (int k = 0; k < R1; ++k) ap("%d ", V.comm[i][k]); ap("\n"); }
    ap("\nX part / Z part:\n");
    for (int i = 0; i < R1; ++i) ap("  %s | %s\n", bstr(px(G[i])).c_str(), bstr(pz(G[i])).c_str());
    ap("\nStabilizer group weight enumerator (weights of the 7 non-identity elements):\n ");
    { uint8_t c[NCHI1], plx; uint8_t tt[NQ]; std::copy(S1_bestTypes, S1_bestTypes + NQ, tt);
      compute_c1(tt, NQ, c, plx);
      for (int chi = 1; chi < NCHI1; ++chi) ap(" %d", NQ - c[chi]); }
    ap("\n\nSyndrome distribution:\n");
    for (int s = 0; s < NCHI1; ++s)
        ap("%d%d%d %lld\n", (s>>2)&1, (s>>1)&1, s&1, (long long)V.synd[s]);
    ap("\nINDEPENDENT VERIFICATION\n");
    ap("  generator weights == %d          : %s\n", W1, V.ok_w ? "OK" : "FAIL");
    ap("  pairwise commutation             : %s\n", V.ok_comm ? "OK" : "FAIL");
    ap("  rank(H3) = 3                     : %s\n", V.ok_rank ? "OK" : "FAIL");
    { int64_t alg = TOTAL_ERRORS - ((91778 + S1_bestSumP) / 8 - 1);
      ap("  algebraic score (Krawtchouk)     : %lld %s\n", (long long)alg,
         alg == V.det ? "(matches brute force)" : "*** MISMATCH ***"); }
    ap("\nSearch nodes: %llu\nPruned:       %llu\nLeaves:       %llu\nElapsed:      %.3f s\n",
       (unsigned long long)S1_nodes, (unsigned long long)S1_pruned,
       (unsigned long long)S1_leaves, elapsed());
    ap("============================================================\n");
    dump("stage1.txt");
}

static void s1_save_undetected() {
    RB.clear();
    ap("Stage-1 undetected error set U_3\n");
    ap("size = %zu   (weights 1/2/3/4 = %lld/%lld/%lld/%lld)\n\n", U3.size(),
       (long long)U3_BY_W[1], (long long)U3_BY_W[2], (long long)U3_BY_W[3], (long long)U3_BY_W[4]);
    for (int w = 1; w <= MAXW; ++w) {
        ap("Weight %d:\n", w);
        int n = 0;
        for (uint32_t E : U3) if (pwt(E) == w) ap("  E%-6d = %s\n", ++n, pstr(E).c_str());
        if (!n) ap("  (none)\n");
        ap("\n");
    }
    dump("stage1_undetected.txt");
}


// ===========================================================================================
//  ==========================  S T A G E   2  ==============================================
//  The Stage-1 optimum S3 = <g1,g2,g3> is now FROZEN.  We look for eight weight-10 Paulis,
//  commuting with S3 and with each other, independent, maximising coverage of U_3 ONLY.
//
//  (A) Every admissible h and every target E lies in C(S3) = S3^perp, dim 28-3 = 25, and
//      <h,E> is unchanged if h or E is shifted by an element of S3.  So the entire problem
//      lives in the quotient  V = C(S3)/S3,  dim 25-3 = 22, on which the symplectic form is
//      non-degenerate (the radical of the form on C(S3) is exactly S3 -- verified at run
//      time).  A symplectic basis u_1..u_11, w_1..w_11 turns a class into a 22-bit integer
//              class(v) = ( <w_i,v> )_i | ( <u_i,v> )_i << 11
//      and the form into   parity( (a & b') ^ (b & a') ):  two ANDs and a popcount.
//
//  (B) Mutual commutation + independence  <=>  the eight classes span an 8-dimensional
//      TOTALLY ISOTROPIC subspace D <= V, and rank 11 for the whole stabilizer is dim D = 8.
//      A target is undetected iff its class lies in D^perp, so the objective is a function
//      of the SUBSPACE, not of the basis.
//
//  (C) THE COVERAGE IDENTITY.  With cov1(v) = #{E in U_3 : <v,E> = 1},
//
//          #covered  =  (1/2^{d-1}) * SUM_{v in D \ {0}} cov1(v)      for dim D = d,
//                    =  (1/128) * SUM over the 255 non-zero elements  for d = 8.
//
//      Derivation: #uncovered = sum_{p in D^perp} f(p); D^perp is the dot-product dual of
//      sigma(D) (sigma = swap of the two 11-bit halves), so Poisson summation gives
//      #uncovered = 2^{-d} sum_{v in sigma(D)} fhat(v) with fhat(sigma(w)) = |U_3| -
//      2 cov1(w).  Summing the 2^d terms leaves |U_3| - 2^{1-d} sum_{w != 0} cov1(w).
//      NOTE the normalisation is 1/2^{d-1} = 1/128 here, NOT the 1/64 of the d = 7 problem.
//      The identity is checked against the explicit |U_3|-bit union in the self-test.
//
//  (D) cov1 for all 2^22 classes comes from ONE Walsh-Hadamard transform; the classes that
//      a weight-10 Pauli can reach come from a Gray-code sweep of all 2^25 elements of C(S3).
// ===========================================================================================
constexpr int VDIM2 = 2 * NQ - 2 * R1;        // dim C(S3)/S3 = 22
constexpr int VHALF = VDIM2 / 2;              // 11
constexpr int CDIM2 = 2 * NQ - R1;            // dim C(S3) = 25
constexpr uint32_t VLO = (1u << VHALF) - 1u;
constexpr int VSZ2 = 1 << VDIM2;              // 4194304 classes

static uint32_t UB2[VHALF], WB2[VHALF];       // symplectic basis of V
static uint32_t S3EL[1 << R1];                // the 8 elements of S3

static inline uint32_t cls_of(uint32_t v) {
    uint32_t c = 0;
    for (int i = 0; i < VHALF; ++i) {
        c |= uint32_t(symp(WB2[i], v)) << i;
        c |= uint32_t(symp(UB2[i], v)) << (i + VHALF);
    }
    return c;
}
static inline uint32_t lift_of(uint32_t c) {
    uint32_t v = 0;
    for (int i = 0; i < VHALF; ++i) {
        if ((c >> i) & 1)            v ^= UB2[i];
        if ((c >> (i + VHALF)) & 1)  v ^= WB2[i];
    }
    return v;
}
static inline int csymp(uint32_t a, uint32_t b) {
    return std::popcount(((a & VLO) & (b >> VHALF)) ^ ((a >> VHALF) & (b & VLO))) & 1;
}

static void build_symplectic_basis2() {
    std::vector<uint32_t> cur;
    {
        std::vector<uint32_t> basis;
        for (int j = 0; j < NQ; ++j) { basis.push_back(pk(1u << j, 0)); basis.push_back(pk(0, 1u << j)); }
        for (int g = 0; g < R1; ++g) {
            int piv = -1;
            for (size_t i = 0; i < basis.size(); ++i) if (symp(G1[g], basis[i])) { piv = (int)i; break; }
            if (piv < 0) { fprintf(stderr, "Stage-1 generators are dependent\n"); exit(1); }
            uint32_t pv = basis[piv]; basis.erase(basis.begin() + piv);
            for (auto& b : basis) if (symp(G1[g], b)) b ^= pv;
        }
        cur = basis;
    }
    if ((int)cur.size() != CDIM2) { fprintf(stderr, "dim C(S3) = %zu != %d\n", cur.size(), CDIM2); exit(1); }
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
        UB2[np] = u; WB2[np] = w; ++np;
    }
    if (np != VHALF || (int)cur.size() != R1) {
        fprintf(stderr, "symplectic reduction failed (%d pairs, radical %zu)\n", np, cur.size()); exit(1); }
    for (int m = 0; m < (1 << R1); ++m) { uint32_t t = 0;
        for (int g = 0; g < R1; ++g) if ((m >> g) & 1) t ^= G1[g];
        S3EL[m] = t; }
    for (uint32_t r : cur) {                              // the radical must be exactly S3
        bool in = false; for (int m = 0; m < (1 << R1); ++m) if (S3EL[m] == r) in = true;
        if (!in) { fprintf(stderr, "radical of the form on C(S3) is not S3\n"); exit(1); } }
}

// ---- cov1 over all 2^22 classes, and weight-10 class representatives ----------------------
static std::vector<uint16_t> COV1;
static std::vector<uint32_t> REPR10;
static std::vector<uint32_t> TCLS;
static int COVMAX = 0;
static int64_t N_W10 = 0, NCLS10 = 0;

static void build_cov1() {
    TCLS.resize(U3.size());
    std::vector<int32_t> f(VSZ2, 0);
    for (size_t j = 0; j < U3.size(); ++j) {
        uint32_t c = cls_of(U3[j]); TCLS[j] = c;
        f[((c >> VHALF) & VLO) | ((c & VLO) << VHALF)] += 1;      // sigma(c)
    }
    for (int b = 0; b < VDIM2; ++b) {
        int step = 1 << b;
        for (int i = 0; i < VSZ2; i += step << 1)
            for (int k = i; k < i + step; ++k) {
                int32_t a = f[k], d = f[k + step]; f[k] = a + d; f[k + step] = a - d; }
    }
    const int32_t N = (int32_t)U3.size();
    COV1.resize(VSZ2);
    for (int v = 0; v < VSZ2; ++v) { int32_t c = (N - f[v]) / 2; COV1[v] = uint16_t(c);
        if (v && c > COVMAX) COVMAX = c; }
}
static void build_repr10() {
    REPR10.assign(VSZ2, NONE);
    uint32_t basis[CDIM2], bcls[CDIM2];
    for (int g = 0; g < R1; ++g) { basis[g] = G1[g]; bcls[g] = 0; }
    for (int i = 0; i < VHALF; ++i) { basis[R1 + i] = UB2[i]; bcls[R1 + i] = 1u << i; }
    for (int i = 0; i < VHALF; ++i) { basis[R1 + VHALF + i] = WB2[i]; bcls[R1 + VHALF + i] = 1u << (i + VHALF); }
    uint32_t cur = 0, ccl = 0;
    const uint64_t TOT = 1ull << CDIM2;
    for (uint64_t i = 1; i < TOT; ++i) {
        int b = std::countr_zero(i);
        cur ^= basis[b]; ccl ^= bcls[b];
        if (pwt(cur) == W2) { ++N_W10;
            if (ccl && (REPR10[ccl] == NONE || cur < REPR10[ccl])) {
                if (REPR10[ccl] == NONE) ++NCLS10; REPR10[ccl] = cur; } }
    }
}

// ---- bitset engine over U_3 (used only for independent verification) ----------------------
static int TGTN = 0, WORDS = 0;
static std::vector<uint64_t> MASKBASE;
static void build_maskbase() {
    TGTN = (int)U3.size(); WORDS = (TGTN + 63) / 64;
    MASKBASE.assign(size_t(VDIM2) * WORDS, 0);
    for (int b = 0; b < VDIM2; ++b) {
        uint64_t* m = &MASKBASE[size_t(b) * WORDS];
        uint32_t v = 1u << b;
        for (int j = 0; j < TGTN; ++j) if (csymp(v, TCLS[j])) m[j >> 6] |= 1ull << (j & 63);
    }
}
static void class_mask(uint32_t v, uint64_t* out) {
    std::memset(out, 0, size_t(WORDS) * 8);
    while (v) { int b = std::countr_zero(v); v &= v - 1;
        const uint64_t* m = &MASKBASE[size_t(b) * WORDS];
        for (int k = 0; k < WORDS; ++k) out[k] ^= m[k]; }
}

// ---- Stage-2 solution bookkeeping ---------------------------------------------------------
struct Sol2 {
    uint32_t cls[R2]{}; uint32_t pau[R2]{};
    int64_t coverage = -1, sumIndiv = 0;
    int minIndiv = 0, ncol1 = 0; bool valid = false;
};
static inline int64_t span_sum(const uint32_t* cls, int d) {
    int64_t s = 0;
    for (int m = 1; m < (1 << d); ++m) {
        uint32_t v = 0; int mm = m;
        while (mm) { int b = std::countr_zero(unsigned(mm)); mm &= mm - 1; v ^= cls[b]; }
        s += COV1[v];
    }
    return s;
}
static inline int64_t coverage_of(const uint32_t* cls, int d) { return span_sum(cls, d) / (1LL << (d - 1)); }

static std::vector<uint32_t> CAND;
static std::mutex BESTMTX; static Sol2 BEST; static std::atomic<int64_t> BESTCOV{-1};

static bool valid_set(const uint32_t* cls) {
    for (int a = 0; a < R2; ++a) { if (!cls[a]) return false;
        for (int b = a + 1; b < R2; ++b) if (csymp(cls[a], cls[b])) return false; }
    uint32_t piv[VDIM2] = {0};
    for (int a = 0; a < R2; ++a) {
        uint32_t v = cls[a];
        while (v) { int b = 31 - std::countl_zero(v);
            if (piv[b]) v ^= piv[b]; else { piv[b] = v; break; } }
        if (!v) return false;
    }
    return true;
}
// Number of physical qubits whose column of the complete 11-generator matrix spans less than
// two dimensions.  The Krawtchouk identity at w = 1 gives exactly
//        (undetected weight-1 targets) = 3*n0 + n1 = this count,
// and because N(L) is a group, each such qubit also drags in the products of its weight-1
// error with the other leftovers.  --full-columns makes the search drive it to zero first;
// the reported coverage is always the true objective, the flag only steers the heuristic.
static bool REQ_FULLCOL = false;
static int deficient_columns(const uint32_t* pau) {
    uint32_t all[R1 + R2];
    for (int i = 0; i < R1; ++i) all[i] = G1[i];
    for (int a = 0; a < R2; ++a) all[R1 + a] = pau[a];
    int bad = 0;
    for (int j = 0; j < NQ; ++j) {
        uint32_t xc = 0, zc = 0;
        for (int i = 0; i < R1 + R2; ++i) {
            xc |= ((px(all[i]) >> j) & 1u) << i;
            zc |= ((pz(all[i]) >> j) & 1u) << i;
        }
        if (xc == 0 || zc == 0 || xc == zc) ++bad;
    }
    return bad;
}
static Sol2 finalize_sol(const uint32_t* cls) {
    Sol2 S; std::copy(cls, cls + R2, S.cls); std::sort(S.cls, S.cls + R2);
    S.coverage = coverage_of(S.cls, R2);
    S.sumIndiv = 0; S.minIndiv = 1 << 30;
    for (int a = 0; a < R2; ++a) { S.pau[a] = REPR10[S.cls[a]];
        S.sumIndiv += COV1[S.cls[a]]; S.minIndiv = std::min<int>(S.minIndiv, COV1[S.cls[a]]); }
    S.ncol1 = REQ_FULLCOL ? deficient_columns(S.pau) : 0;
    S.valid = true; return S;
}
static bool better(const Sol2& a, const Sol2& b) {
    if (!b.valid) return a.valid;
    if (REQ_FULLCOL && a.ncol1 != b.ncol1) return a.ncol1 < b.ncol1;
    if (a.coverage != b.coverage) return a.coverage > b.coverage;
    if (a.sumIndiv != b.sumIndiv) return a.sumIndiv > b.sumIndiv;
    if (a.minIndiv != b.minIndiv) return a.minIndiv > b.minIndiv;
    return std::lexicographical_compare(a.cls, a.cls + R2, b.cls, b.cls + R2);
}
static bool offer2(const uint32_t* cls) {
    if (!valid_set(cls)) return false;
    Sol2 S = finalize_sol(cls);
    if (S.coverage < BESTCOV.load(std::memory_order_relaxed)) return false;
    std::lock_guard<std::mutex> lk(BESTMTX);
    if (!better(S, BEST)) return false;
    BEST = S; BESTCOV.store(S.coverage); return true;
}

// ---- Stage-2 heuristics: randomized greedy, steepest-descent 1-opt, ruin-and-recreate -----
// The greedy gain of adding class v to a d-dimensional span D is exactly
//     G(v) = SUM_{u in D} cov1(v ^ u),
// i.e. the true increment of the span sum, so greedy optimises the real objective at every
// step rather than a surrogate.
struct Builder {
    std::vector<uint32_t> allow, tmp;
    uint32_t span[1 << R2], spanSorted[1 << R2];
    int nspan = 1;
    void reset() { span[0] = 0; spanSorted[0] = 0; nspan = 1; }
    bool in_span(uint32_t v) const { return std::binary_search(spanSorted, spanSorted + nspan, v); }
    void add_to_span(uint32_t v) {
        for (int i = 0; i < nspan; ++i) span[nspan + i] = span[i] ^ v;
        nspan <<= 1;
        std::copy(span, span + nspan, spanSorted); std::sort(spanSorted, spanSorted + nspan);
    }
    int64_t gain(uint32_t v) const {
        int64_t g = 0; for (int i = 0; i < nspan; ++i) g += COV1[span[i] ^ v]; return g;
    }
};
using GPair = std::pair<int64_t, uint32_t>;
static inline bool gcmp(const GPair& a, const GPair& b) { return a.first > b.first; }

static bool greedy_from(const uint32_t* keep, int nk, std::mt19937_64& rng, int rndTop,
                        uint32_t* out, Builder& B) {
    B.reset();
    for (int i = 0; i < nk; ++i) { out[i] = keep[i]; B.add_to_span(keep[i]); }
    if (nk == 0) B.allow.assign(CAND.begin(), CAND.end());
    else {
        B.allow.clear();
        for (uint32_t v : CAND) {
            bool ok = true;
            for (int k = 0; k < nk && ok; ++k) if (csymp(v, keep[k])) ok = false;
            if (ok && !B.in_span(v)) B.allow.push_back(v);
        }
    }
    for (int d = nk; d < R2; ++d) {
        if (B.allow.empty()) return false;
        std::vector<GPair> top; top.reserve(size_t(rndTop) + 4);
        for (uint32_t v : B.allow) {
            int64_t g = B.gain(v);
            if ((int)top.size() < rndTop) { top.emplace_back(g, v); std::push_heap(top.begin(), top.end(), gcmp); }
            else if (g > top.front().first) { std::pop_heap(top.begin(), top.end(), gcmp);
                top.back() = {g, v}; std::push_heap(top.begin(), top.end(), gcmp); }
        }
        uint32_t pick = top[rng() % top.size()].second;
        out[d] = pick; B.add_to_span(pick);
        B.tmp.clear();
        for (uint32_t v : B.allow) if (!csymp(v, pick) && !B.in_span(v)) B.tmp.push_back(v);
        B.allow.swap(B.tmp);
    }
    return true;
}

struct LocalWS { std::vector<uint8_t> om; };
static void om_column(std::vector<uint8_t>& om, int a, uint32_t c) {
    const uint8_t bit = uint8_t(1u << a), clr = uint8_t(~bit);
    for (size_t i = 0; i < CAND.size(); ++i)
        om[i] = uint8_t((om[i] & clr) | (csymp(CAND[i], c) ? 0 : bit));
}
// Replace one generator by the best compatible candidate, repeatedly.  om[] caches which of
// the eight current generators each candidate commutes with; after a swap only that one
// column is stale, so a sweep costs eight linear passes over the candidate array.
static bool local_improve(uint32_t* cls, LocalWS& WS) {
    if (WS.om.size() != CAND.size()) WS.om.assign(CAND.size(), 0);
    for (int a = 0; a < R2; ++a) om_column(WS.om, a, cls[a]);
    Sol2 cur = finalize_sol(cls);
    bool improved = false;
    for (;;) {
        bool any = false;
        for (int a = 0; a < R2; ++a) {
            uint32_t sp[1 << (R2 - 1)]; int ns = 1; sp[0] = 0;
            for (int b = 0; b < R2; ++b) if (b != a) {
                for (int i = 0; i < ns; ++i) sp[ns + i] = sp[i] ^ cls[b];
                ns <<= 1; }
            uint32_t spS[1 << (R2 - 1)]; std::copy(sp, sp + ns, spS); std::sort(spS, spS + ns);
            const uint8_t need = uint8_t(((1u << R2) - 1u) & ~(1u << a));
            uint32_t trial[R2]; std::copy(cls, cls + R2, trial);
            uint32_t bestV = cls[a]; Sol2 bestS = cur;
            for (size_t i = 0; i < CAND.size(); ++i) {
                if ((WS.om[i] & need) != need) continue;
                uint32_t v = CAND[i];
                if (std::binary_search(spS, spS + ns, v)) continue;
                trial[a] = v;
                Sol2 T = finalize_sol(trial);
                if (better(T, bestS)) { bestS = T; bestV = v; }
            }
            if (bestV != cls[a]) { cls[a] = bestV; cur = bestS; any = true; improved = true;
                om_column(WS.om, a, bestV); }
        }
        if (!any) break;
    }
    return improved;
}

// ---- beam search over subspaces, de-duplicated by a hash of the sorted span ---------------
struct BeamState { uint32_t cls[R2]; int32_t last; int64_t S; uint64_t key; };
static std::vector<uint32_t> POOL;
static uint64_t span_key(const uint32_t* cls, int d) {
    uint32_t sp[1 << R2]; int ns = 1; sp[0] = 0;
    for (int b = 0; b < d; ++b) { for (int i = 0; i < ns; ++i) sp[ns + i] = sp[i] ^ cls[b]; ns <<= 1; }
    std::sort(sp, sp + ns);
    uint64_t h = 1469598103934665603ull;
    for (int i = 1; i < ns; ++i) { h ^= sp[i]; h *= 1099511628211ull; }
    return h;
}
static void beam_search(int K, int poolN, int nthreads) {
    if (K <= 0) return;
    POOL.assign(CAND.begin(), CAND.begin() + std::min<size_t>(CAND.size(), size_t(poolN)));
    std::vector<BeamState> cur, nxt;
    for (int i = 0; i < (int)POOL.size() && (int)cur.size() < K; ++i) {
        BeamState s{}; s.cls[0] = POOL[i]; s.last = i; s.S = COV1[POOL[i]];
        s.key = span_key(s.cls, 1); cur.push_back(s);
    }
    for (int d = 1; d < R2; ++d) {
        std::vector<std::vector<BeamState>> perThread(nthreads);
        std::atomic<size_t> next{0};
        std::vector<std::thread> th;
        for (int t = 0; t < nthreads; ++t) th.emplace_back([&, t]() {
            auto& out = perThread[t];
            uint32_t sp[1 << R2], spS[1 << R2];
            for (;;) {
                size_t si = next.fetch_add(1);
                if (si >= cur.size()) break;
                const BeamState& s = cur[si];
                int ns = 1; sp[0] = 0;
                for (int b = 0; b < d; ++b) { for (int i = 0; i < ns; ++i) sp[ns + i] = sp[i] ^ s.cls[b]; ns <<= 1; }
                std::copy(sp, sp + ns, spS); std::sort(spS, spS + ns);
                constexpr int KPER = 24;
                std::vector<std::pair<int64_t,int>> best; best.reserve(KPER + 1);
                auto cmp = [](const std::pair<int64_t,int>& a, const std::pair<int64_t,int>& b) { return a.first > b.first; };
                for (int j = s.last + 1; j < (int)POOL.size(); ++j) {
                    uint32_t v = POOL[j];
                    bool ok = true;
                    for (int b = 0; b < d && ok; ++b) if (csymp(v, s.cls[b])) ok = false;
                    if (!ok || std::binary_search(spS, spS + ns, v)) continue;
                    int64_t g = 0; for (int i = 0; i < ns; ++i) g += COV1[sp[i] ^ v];
                    if ((int)best.size() < KPER) { best.emplace_back(g, j); std::push_heap(best.begin(), best.end(), cmp); }
                    else if (g > best.front().first) { std::pop_heap(best.begin(), best.end(), cmp);
                        best.back() = {g, j}; std::push_heap(best.begin(), best.end(), cmp); }
                }
                for (auto& b : best) { BeamState n = s; n.cls[d] = POOL[b.second]; n.last = b.second;
                    n.S = s.S + b.first; n.key = span_key(n.cls, d + 1); out.push_back(n); }
            }
        });
        for (auto& x : th) x.join();
        nxt.clear();
        for (auto& v : perThread) nxt.insert(nxt.end(), v.begin(), v.end());
        if (nxt.empty()) return;
        std::sort(nxt.begin(), nxt.end(), [](const BeamState& a, const BeamState& b) {
            if (a.S != b.S) return a.S > b.S;
            return std::lexicographical_compare(a.cls, a.cls + R2, b.cls, b.cls + R2); });
        std::vector<BeamState> ded; ded.reserve(std::min<size_t>(nxt.size(), size_t(K)));
        std::unordered_set<uint64_t> seen; seen.reserve(nxt.size() * 2);
        for (const auto& s : nxt) { if (!seen.insert(s.key).second) continue;
            ded.push_back(s); if ((int)ded.size() >= K) break; }
        cur.swap(ded);
    }
    for (const auto& s : cur) offer2(s.cls);
}

// ---- exact branch and bound ---------------------------------------------------------------
//  Admissible bound.  With D_d fixed (2^d elements, span sum S_d), the next generator adds
//  exactly G(v) and the remaining ones add cosets totalling 256 - 2^{d+1} further elements,
//  each of cov1 at most COVMAX.  Relaxing isotropy, independence and weight-10 realisability
//  can only enlarge the feasible set, so
//        coverage <= ( S_d + max_{v allowed} G(v) + (256 - 2^{d+1}) * COVMAX ) / 128
//  is valid.  It can prune only while COVMAX <= 128*incumbent/255, which is reported.
static std::atomic<bool> ABORT{false};
static std::atomic<uint64_t> NODES{0}, PRUNED{0};
static double TIME_LIMIT = 1e18;

struct Exact {
    std::vector<uint32_t> allow[R2 + 1];
    std::vector<int64_t> gains;
    uint32_t chosen[R2]{}, span[1 << R2]{};
    uint64_t nodes = 0, pruned = 0;
    void dfs(int d, int nspan, int64_t S) {
        if (ABORT.load(std::memory_order_relaxed)) return;
        ++nodes;
        if ((nodes & 0xFFFF) == 0 && elapsed() > TIME_LIMIT) { ABORT.store(true); return; }
        if (d == R2) { offer2(chosen); return; }
        auto& A = allow[d];
        if (A.empty()) { ++pruned; return; }
        gains.resize(A.size());
        int64_t gbest = -1;
        for (size_t i = 0; i < A.size(); ++i) {
            uint32_t v = CAND[A[i]];
            int64_t g = 0; for (int k = 0; k < nspan; ++k) g += COV1[span[k] ^ v];
            gains[i] = g; if (g > gbest) gbest = g;
        }
        const int64_t rest = int64_t(256 - (1 << (d + 1))) * COVMAX;
        if ((S + gbest + rest) / 128 <= BESTCOV.load(std::memory_order_relaxed)) { ++pruned; return; }
        for (size_t i = 0; i < A.size(); ++i) {
            if (ABORT.load(std::memory_order_relaxed)) return;
            int64_t S2 = S + gains[i];
            if ((S2 + rest) / 128 <= BESTCOV.load(std::memory_order_relaxed)) { ++pruned; continue; }
            uint32_t v = CAND[A[i]];
            chosen[d] = v;
            int ns2 = nspan << 1;
            for (int k = 0; k < nspan; ++k) span[nspan + k] = span[k] ^ v;
            uint32_t sorted[1 << R2]; std::copy(span, span + ns2, sorted); std::sort(sorted, sorted + ns2);
            auto& B = allow[d + 1]; B.clear();
            for (size_t j = i + 1; j < A.size(); ++j) {
                uint32_t u = CAND[A[j]];
                if (csymp(u, v) || std::binary_search(sorted, sorted + ns2, u)) continue;
                B.push_back(A[j]);
            }
            dfs(d + 1, ns2, S2);
        }
    }
};
static void exact_search(int nthreads) {
    std::atomic<size_t> next{0};
    std::vector<std::thread> th;
    for (int t = 0; t < nthreads; ++t) th.emplace_back([&]() {
        Exact E;
        for (;;) {
            size_t i = next.fetch_add(1);
            if (i >= CAND.size() || ABORT.load()) break;
            uint32_t v = CAND[i];
            E.chosen[0] = v; E.span[0] = 0; E.span[1] = v;
            auto& B = E.allow[1]; B.clear();
            for (size_t j = i + 1; j < CAND.size(); ++j) {
                uint32_t u = CAND[j];
                if (csymp(u, v) || u == v) continue;
                B.push_back((uint32_t)j);
            }
            E.dfs(1, 2, COV1[v]);
        }
        NODES += E.nodes; PRUNED += E.pruned;
    });
    for (auto& x : th) x.join();
}

// ---- self-tests ---------------------------------------------------------------------------
static uint64_t splitmix(uint64_t x) {
    x += 0x9E3779B97F4A7C15ull; x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBull; return x ^ (x >> 31);
}
static bool self_test(int n, uint64_t seed) {
    std::mt19937_64 rng(seed); bool ok = true;
    for (int t = 0; t < n && ok; ++t) {                       // (a) cov1 table vs direct count
        uint32_t v = uint32_t(rng() % VSZ2);
        int64_t c = 0; for (int j = 0; j < TGTN; ++j) c += csymp(v, TCLS[j]);
        if (c != (int64_t)COV1[v]) { printf("SELFTEST FAIL cov1[%u]: %lld vs %d\n", v, (long long)c, (int)COV1[v]); ok = false; }
    }
    for (int t = 0; t < n && ok; ++t) {                       // (b) class map is a section
        uint32_t v = 0;
        for (int k = 0; k < CDIM2; ++k) if (rng() & 1)
            v ^= (k < R1 ? G1[k] : (k < R1 + VHALF ? UB2[k - R1] : WB2[k - R1 - VHALF]));
        uint32_t r = v ^ lift_of(cls_of(v));
        bool in = false; for (int m = 0; m < (1 << R1); ++m) if (S3EL[m] == r) in = true;
        if (!in) { printf("SELFTEST FAIL: class map is not a section of C(S3) -> V\n"); ok = false; }
    }
    Builder B; std::vector<uint64_t> uni(WORDS), tm(WORDS);   // (c) identity vs explicit union
    for (int t = 0; t < std::max(1, n / 20) && ok; ++t) {
        uint32_t cls[R2];
        if (!greedy_from(nullptr, 0, rng, 64, cls, B)) continue;
        std::fill(uni.begin(), uni.end(), 0ull);
        for (int a = 0; a < R2; ++a) { class_mask(cls[a], tm.data());
            for (int k = 0; k < WORDS; ++k) uni[k] |= tm[k]; }
        int64_t pc = 0; for (int k = 0; k < WORDS; ++k) pc += std::popcount(uni[k]);
        int64_t alg = coverage_of(cls, R2);
        if (pc != alg) { printf("SELFTEST FAIL: union %lld vs identity %lld\n", (long long)pc, (long long)alg); ok = false; }
    }
    return ok;
}

// ---- heuristic driver ---------------------------------------------------------------------
static bool DETERMINISTIC = false;
static void heuristic(int nthreads, uint64_t seed, int restarts, double tlimit) {
    std::atomic<int> next{0};
    std::vector<std::thread> th;
    for (int t = 0; t < nthreads; ++t) th.emplace_back([&]() {
        Builder B; LocalWS WS; uint32_t cls[R2], keep[R2];
        Sol2 local;
        for (;;) {
            int r = next.fetch_add(1);
            if (restarts > 0 && r >= restarts) break;
            if (elapsed() > tlimit) break;
            std::mt19937_64 rng(splitmix(seed ^ (uint64_t(r) * 0x1000193ull)));
            bool ok;
            if (DETERMINISTIC || !local.valid || (r % 4) == 0) {
                ok = greedy_from(nullptr, 0, rng, (r == 0) ? 1 : 1 + int(r % 12), cls, B);
            } else {                                       // ruin and recreate
                int nk = R2 - (1 + int(rng() % 3));
                int idx[R2]; std::iota(idx, idx + R2, 0);
                std::shuffle(idx, idx + R2, rng);
                for (int i = 0; i < nk; ++i) keep[i] = local.cls[idx[i]];
                ok = greedy_from(keep, nk, rng, 1 + int(rng() % 8), cls, B);
            }
            if (!ok) continue;
            local_improve(cls, WS);
            Sol2 S = finalize_sol(cls);
            if (valid_set(S.cls) && (!local.valid || better(S, local))) local = S;
            offer2(cls);
            if (!DETERMINISTIC && (r % 32) == 31) {
                std::lock_guard<std::mutex> lk(BESTMTX);
                if (BEST.valid && (!local.valid || better(BEST, local))) local = BEST;
            }
        }
    });
    for (auto& x : th) x.join();
}

// ---- Stage-2 independent verification -----------------------------------------------------
struct S2Verify {
    bool wt_ok = true, comm_g_ok = true, comm_h_ok = true, rank_ok = true, u3_ok = true;
    int  weights[R2]{}, rank = 0, comm[R1 + R2][R1 + R2]{};
    int64_t u3size = 0, covered = 0, per_gen[R2]{};
    int64_t cov_by_w[MAXW + 1]{}, tgt_by_w[MAXW + 1]{};
    std::vector<uint32_t> leftover;
    int64_t all_total = 0, orig_det = 0, orig_und = 0, new_det = 0, still_und = 0, det11 = 0;
};
static S2Verify s2_verify(const Sol2& S) {
    S2Verify V;
    for (int a = 0; a < R2; ++a) { V.weights[a] = pwt(S.pau[a]); if (V.weights[a] != W2) V.wt_ok = false; }
    uint32_t all[R1 + R2];
    for (int i = 0; i < R1; ++i) all[i] = G1[i];
    for (int a = 0; a < R2; ++a) all[R1 + a] = S.pau[a];
    for (int i = 0; i < R1 + R2; ++i) for (int j = 0; j < R1 + R2; ++j) {
        V.comm[i][j] = symp(all[i], all[j]);
        if (i != j && V.comm[i][j]) { if (i < R1 || j < R1) V.comm_g_ok = false; else V.comm_h_ok = false; } }
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
    // rebuild U_3 from scratch and test every target directly against the eight new generators
    for_each_error([&](uint32_t E, int w) {
        int s0 = 0; for (int i = 0; i < R1; ++i) s0 |= symp(G1[i], E) << i;
        int s1 = 0; for (int a = 0; a < R2; ++a) s1 |= symp(S.pau[a], E) << a;
        V.all_total++;
        if (s0) V.orig_det++; else V.orig_und++;
        if (s0 || s1) V.det11++;
        if (!s0) {
            V.u3size++; V.tgt_by_w[w]++;
            if (s1) { V.covered++; V.cov_by_w[w]++; V.new_det++;
                for (int a = 0; a < R2; ++a) if ((s1 >> a) & 1) V.per_gen[a]++; }
            else { V.leftover.push_back(E); V.still_und++; }
        }
    });
    V.u3_ok = (V.u3size == (int64_t)U3.size());
    return V;
}

static void s2_report(const Sol2& S, const S2Verify& V, bool proven, int64_t ubound) {
    std::vector<uint64_t> uni(WORDS, 0), tm(WORDS);
    for (int a = 0; a < R2; ++a) { class_mask(S.cls[a], tm.data());
        for (int k = 0; k < WORDS; ++k) uni[k] |= tm[k]; }
    int64_t unionpc = 0; for (int k = 0; k < WORDS; ++k) unionpc += std::popcount(uni[k]);
    RB.clear();
    ap("============================================================\n");
    ap("STAGE 2 - EIGHT ADDITIONAL WEIGHT-10 GENERATORS\n");
    ap("============================================================\n\n");
    ap("Status: %s\n\n", proven ? "PROVEN GLOBAL OPTIMUM" : "BEST KNOWN - OPTIMUM NOT PROVEN");
    ap("Stage-1 target errors: %lld\n\n", (long long)V.u3size);
    ap("Detected by Stage 2:   %lld\n", (long long)S.coverage);
    ap("Still undetected:      %lld\n\n", (long long)(V.u3size - S.coverage));
    ap("Coverage percentage:   %.8f %%\n\n", 100.0 * double(S.coverage) / double(V.u3size));
    if (S.coverage == V.u3size) ap("ALL STAGE-1 UNDETECTED ERRORS DETECTED\n\n");
    for (int w = 1; w <= MAXW; ++w)
        ap("Weight %d:\nDetected / Stage-1 target : %lld / %lld\nRemaining                 : %lld\n\n",
           w, (long long)V.cov_by_w[w], (long long)V.tgt_by_w[w],
           (long long)(V.tgt_by_w[w] - V.cov_by_w[w]));
    ap("Stage-1 generators (frozen):\n");
    for (int i = 0; i < R1; ++i) ap("g%d = %s   (weight %d)\n", i + 1, pstr(G1[i]).c_str(), pwt(G1[i]));
    ap("\nAdditional generators:\n");
    for (int a = 0; a < R2; ++a) ap("h%d = %s   (weight %d)\n", a + 1, pstr(S.pau[a]).c_str(), V.weights[a]);
    ap("\nWeights:\n");
    for (int a = 0; a < R2; ++a) ap("h%d: %d\n", a + 1, V.weights[a]);
    ap("\nRank of complete 11-generator stabilizer:\n%d  (target %d)\n", V.rank, R1 + R2);
    ap("Qubits whose 11-generator column spans < 2 dimensions: %d"
       "  (equals the number of undetected weight-1 targets)\n", deficient_columns(S.pau));
    ap("\nCommutation (rows/cols g1..g3,h1..h8; 0 = commute):\n");
    for (int i = 0; i < R1 + R2; ++i) { ap("  ");
        for (int j = 0; j < R1 + R2; ++j) ap("%d ", V.comm[i][j]); ap("\n"); }
    ap("\nTarget-error detection by each new generator:\n");
    for (int a = 0; a < R2; ++a)
        ap("h%d: %lld   (algebraic cov1 = %d %s)\n", a + 1, (long long)V.per_gen[a],
           (int)COV1[S.cls[a]], V.per_gen[a] == (int64_t)COV1[S.cls[a]] ? "OK" : "*** MISMATCH ***");
    ap("\nUnion coverage:\n");
    ap("  |U_3|-bit bitset popcount     : %lld\n", (long long)unionpc);
    ap("  brute-force recount           : %lld\n", (long long)V.covered);
    ap("  algebraic (1/128)*sum cov1    : %lld\n", (long long)S.coverage);
    ap("  %s\n", (unionpc == V.covered && V.covered == S.coverage) ? "all three agree" : "*** MISMATCH ***");
    ap("\nINDEPENDENT VERIFICATION\n");
    ap("  U_3 rebuilt from scratch      : %lld errors %s\n", (long long)V.u3size, V.u3_ok ? "OK" : "*** MISMATCH ***");
    ap("  all weights == %d              : %s\n", W2, V.wt_ok ? "OK" : "FAIL");
    ap("  [h_a , g_i] = 0               : %s\n", V.comm_g_ok ? "OK" : "FAIL");
    ap("  [h_a , h_b] = 0               : %s\n", V.comm_h_ok ? "OK" : "FAIL");
    ap("  rank = 11                     : %s\n", V.rank_ok ? "OK" : "FAIL");
    ap("\nProvable upper bound on M_8 (any 8 commuting independent generators):\n");
    ap("  coverage = (1/128) * sum of cov1 over the 255 non-zero elements of the span\n");
    ap("  <= min( |U_3| , 255*max_v cov1(v)/128 ) = %lld\n", (long long)ubound);
    ap("\nCandidate classes: %zu\n", CAND.size());
    ap("Search nodes: %llu\nPruned:       %llu\nElapsed:      %.3f s\n",
       (unsigned long long)NODES.load(), (unsigned long long)PRUNED.load(), elapsed());
    ap("============================================================\n");
    dump("stage2_best.txt");

    RB.clear();
    ap("Stage-2 remaining (still undetected) errors: %zu of %lld Stage-1 targets\n\n",
       V.leftover.size(), (long long)V.u3size);
    for (int w = 1; w <= MAXW; ++w) {
        ap("Weight %d:\n", w); int n = 0;
        for (uint32_t E : V.leftover) if (pwt(E) == w) ap("  E%-6d = %s\n", ++n, pstr(E).c_str());
        if (!n) ap("  (none)\n");
        ap("\n");
    }
    dump("stage2_remaining.txt");
}

static void final_summary(const Sol2& S, const S1Verify& V1, const S2Verify& V2,
                          bool s1proven, bool s2proven) {
    RB.clear();
    ap("============================================================\n");
    ap("FINAL SUMMARY\n");
    ap("============================================================\n\n");
    ap("STAGE 1\n3 generators, all weight 8   [%s]\n\n",
       s1proven ? "PROVEN GLOBAL OPTIMUM" : "BEST KNOWN - OPTIMUM NOT PROVEN");
    ap("Total errors considered: %lld\n", (long long)TOTAL_ERRORS);
    ap("Detected:                %lld\n", (long long)V1.det);
    ap("Undetected:              %lld\n\n", (long long)V1.und);
    ap("Detection percentage:    %.8f %%\n\n", 100.0 * double(V1.det) / double(TOTAL_ERRORS));
    ap("STAGE 2\n8 additional generators, all weight 10   [%s]\n\n",
       s2proven ? "PROVEN GLOBAL OPTIMUM" : "BEST KNOWN - OPTIMUM NOT PROVEN");
    ap("Stage-1-undetected errors: %lld\n", (long long)V1.und);
    ap("Detected by Stage 2:       %lld\n", (long long)S.coverage);
    ap("Still undetected:          %lld\n\n", (long long)(V1.und - S.coverage));
    ap("Stage-2 coverage percentage: %.8f %%\n\n", 100.0 * double(S.coverage) / double(V1.und));
    ap("FINAL 11-GENERATOR CODE   [[14,3,d]]\n\n");
    ap("Total weight <=4 errors:            %lld\n", (long long)V2.all_total);
    ap("Originally detected by first 3:     %lld\n", (long long)V2.orig_det);
    ap("Newly detected by additional 8:     %lld\n", (long long)V2.new_det);
    ap("Still undetected:                   %lld\n\n", (long long)V2.still_und);
    ap("Overall detected:                   %lld\n", (long long)V2.det11);
    ap("Overall detection percentage:       %.8f %%\n\n",
       100.0 * double(V2.det11) / double(V2.all_total));
    ap("THE THREE QUANTITIES, SEPARATELY\n");
    ap("  M_3     = %lld   (weight<=4 errors detected by the optimal 3-generator set)\n", (long long)V1.det);
    ap("  |U_3|   = %lld   = 91770 - M_3\n", (long long)V1.und);
    ap("  M_8     = %lld   (Stage-1-undetected errors detected by the 8-generator extension)\n", (long long)S.coverage);
    ap("  M_final = %lld   (detected by the complete 11-generator stabilizer)\n\n", (long long)V2.det11);
    ap("  identity check  M_final == M_3 + M_8 :  %lld == %lld + %lld  ->  %s\n",
       (long long)V2.det11, (long long)V1.det, (long long)S.coverage,
       (V2.det11 == V1.det + S.coverage) ? "VERIFIED" : "*** FAILED ***");
    ap("  identity check  |U_3| == 91770 - M_3 :  %s\n",
       (V1.und == TOTAL_ERRORS - V1.det) ? "VERIFIED" : "*** FAILED ***");
    ap("============================================================\n");
    dump("summary.txt");

    RB.clear();
    ap("Final 11-generator stabilizer code  [[14,3]]\n");
    ap("rows 1-3  : Stage-1 generators, Pauli weight 8\n");
    ap("rows 4-11 : Stage-2 generators, Pauli weight 10\n\n");
    ap("Pauli form:\n");
    for (int i = 0; i < R1; ++i) ap("g%d = %s\n", i + 1, pstr(G1[i]).c_str());
    for (int a = 0; a < R2; ++a) ap("h%d = %s\n", a + 1, pstr(S.pau[a]).c_str());
    ap("\nParity-check matrix  H = [ X | Z ]   (11 x 28):\n");
    { uint32_t all[R1 + R2];
      for (int i = 0; i < R1; ++i) all[i] = G1[i];
      for (int a = 0; a < R2; ++a) all[R1 + a] = S.pau[a];
      for (int i = 0; i < R1 + R2; ++i)
          ap("  %s | %s\n", bstr(px(all[i])).c_str(), bstr(pz(all[i])).c_str()); }
    ap("\nrank = %d,  all 55 generator pairs commute,  weights = ", V2.rank);
    for (int i = 0; i < R1; ++i) ap("%d ", pwt(G1[i]));
    for (int a = 0; a < R2; ++a) ap("%d ", V2.weights[a]);
    ap("\n\nWeight <=4 errors detected by the complete code: %lld / %lld  (%.8f %%)\n",
       (long long)V2.det11, (long long)V2.all_total,
       100.0 * double(V2.det11) / double(V2.all_total));
    ap("Weight <=4 errors still undetected:              %lld\n", (long long)V2.still_und);
    dump("final_matrix.txt");
}

// ===========================================================================================
//  MAIN
// ===========================================================================================

// ===========================================================================================
//  ORBIT DRIVER -- run Stage 2 once for every inequivalent OPTIMAL Stage-1 solution
//
//  Stage 1 is already proven optimal, so the question this program answers is
//        max over { S3 : M_3(S3) = 80584 }   of   max over H8   M_8(S3, H8),
//  i.e. whether the residual errors are an artefact of one particular optimal Stage-1
//  choice or a structural obstruction shared by the whole optimal family.
//
//  The Stage-1 search state is a MULTISET of column SUBSPACE types, so the enumeration is
//  already the quotient by qubit permutations (multiset) and by local Cliffords (subspaces);
//  the S_3 lex-min prefix test additionally quotients generator permutations.  All three
//  group actions provably preserve generator weight, mutual commutation, rank and the whole
//  weight-<=4 detection objective, so distinct canonical multisets are genuinely
//  inequivalent under that group and each one needs its own Stage-2 run.  The quotient space
//  C(S3)/S3 depends on S3, so every orbit rebuilds its own class map, cov1 table and
//  candidate set from scratch -- nothing objective-dependent is cached across orbits.
// ===========================================================================================
static void reset_stage2() {
    COV1.clear(); REPR10.clear(); TCLS.clear(); MASKBASE.clear(); CAND.clear(); POOL.clear();
    COVMAX = 0; N_W10 = 0; NCLS10 = 0; TGTN = 0; WORDS = 0;
    BEST = Sol2(); BESTCOV.store(-1); NODES.store(0); PRUNED.store(0); ABORT.store(false);
    U3.clear(); for (int w = 0; w <= MAXW; ++w) U3_BY_W[w] = 0;
}

// ---- structure of the residual (still-undetected) set -------------------------------------
struct Residual {
    int64_t n = 0;
    int64_t by_w[MAXW + 1] = {0};
    bool is_subgroup = false;          // R + identity closed under multiplication
    int  dim = 0;                      // dimension of span(R) over F_2
    std::vector<uint32_t> gens;        // a basis of span(R)
    bool all_x = false;                // every residual error is X-type
    std::vector<int> qubits;           // physical qubits touched
};
static Residual analyse_residual(const std::vector<uint32_t>& R) {
    Residual r; r.n = (int64_t)R.size();
    for (uint32_t E : R) r.by_w[pwt(E)]++;
    // basis of span(R) over F_2 by incremental Gaussian elimination on the 28-bit vectors
    uint32_t piv[2 * NQ] = {0};
    for (uint32_t E : R) {
        uint32_t v = px(E) | (pz(E) << NQ);
        while (v) {
            int b = 31 - std::countl_zero(v);
            if (piv[b]) v ^= piv[b];
            else { piv[b] = v; r.gens.push_back(pk(v & QM, (v >> NQ) & QM)); break; }
        }
    }
    r.dim = (int)r.gens.size();
    // closure test: R together with the identity must be closed under multiplication
    std::unordered_set<uint32_t> S; S.insert(0u);
    for (uint32_t E : R) S.insert(E);
    r.is_subgroup = true;
    for (uint32_t a : S) { for (uint32_t b : S) if (!S.count(a ^ b)) { r.is_subgroup = false; break; }
        if (!r.is_subgroup) break; }
    if (r.is_subgroup && (int64_t)S.size() != (1LL << r.dim)) r.is_subgroup = false;
    r.all_x = true;
    uint32_t touched = 0;
    for (uint32_t E : R) { if (pz(E)) r.all_x = false; touched |= px(E) | pz(E); }
    for (int j = 0; j < NQ; ++j) if ((touched >> j) & 1) r.qubits.push_back(j);
    return r;
}
static std::string support_notation(uint32_t E) {
    std::string s; bool first = true;
    for (int j = 0; j < NQ; ++j) {
        int x = (px(E) >> j) & 1, z = (pz(E) >> j) & 1;
        if (!x && !z) continue;
        if (!first) s += " ";
        s += PAULI_CH[x][z]; s += "_" + std::to_string(j);
        first = false;
    }
    return s.empty() ? std::string("I") : s;
}
static void append_residual_block(const std::vector<uint32_t>& R, const Residual& r) {
    ap("Remaining undetected errors: %lld\n\n", (long long)r.n);
    for (int w = 1; w <= MAXW; ++w) {
        ap("Weight %d:\n", w); int n = 0;
        for (uint32_t E : R) if (pwt(E) == w) {
            ap("  %s = %s\n", pstr(E).c_str(), support_notation(E).c_str()); ++n; }
        if (!n) ap("  (none)\n");
        ap("\n");
    }
    ap("Residual weight distribution: %lld / %lld / %lld / %lld  (weights 1/2/3/4)\n",
       (long long)r.by_w[1], (long long)r.by_w[2], (long long)r.by_w[3], (long long)r.by_w[4]);
    ap("Residual set forms a subgroup (with the identity): %s\n", r.is_subgroup ? "YES" : "NO");
    ap("Residual subgroup dimension: %d\n", r.dim);
    ap("Generators:\n");
    for (uint32_t g : r.gens) ap("  %s = %s\n", pstr(g).c_str(), support_notation(g).c_str());
    ap("All residual errors are X-only: %s\n", r.all_x ? "YES" : "NO");
    ap("Physical qubits involved:");
    for (int q : r.qubits) ap(" %d", q);
    ap("\n");
}

// ---- one row of the global comparison table ------------------------------------------------
struct OrbitResult {
    int index = 0;
    uint32_t g[R1]{};
    std::array<uint8_t, NQ> types{};
    int64_t M3 = 0, U3n = 0, M8 = 0, Mfinal = 0;
    int64_t residual = 0, res_w[MAXW + 1] = {0};
    int64_t w4_total = 0, w4_det = 0, w4_und = 0;
    bool w4_all = false, proven = false, full = false;
    int  res_dim = 0; bool res_subgroup = false, res_allx = false;
    Sol2 sol;
    std::vector<uint32_t> leftover;
};

static void print_full_coverage_banner(const OrbitResult& O, const S2Verify& V) {
    RB.clear();
    ap("============================================================\n");
    ap("FULL COVERAGE FOUND\n");
    ap("============================================================\n\n");
    ap("All %lld Stage-1-undetected errors are detected\n", (long long)O.U3n);
    ap("by the eight additional weight-10 generators.\n\nTherefore:\n\n");
    ap("%lld / %lld errors detected.\n\n", (long long)O.Mfinal, (long long)TOTAL_ERRORS);
    ap("The complete 11-generator code detects every\nPauli error of weight <= 4.\n\n");
    ap("============================================================\n\n");
    ap("Stage-1 orbit index: %d\n\nStage-1 generators:\n", O.index);
    for (int i = 0; i < R1; ++i) ap("g%d = %s   (weight %d)\n", i + 1, pstr(O.g[i]).c_str(), pwt(O.g[i]));
    ap("\nStage-2 generators:\n");
    for (int a = 0; a < R2; ++a) ap("h%d = %s   (weight %d)\n", a + 1, pstr(O.sol.pau[a]).c_str(), V.weights[a]);
    ap("\nComplete 11-generator set / parity-check matrix [ X | Z ]:\n");
    uint32_t all[R1 + R2];
    for (int i = 0; i < R1; ++i) all[i] = O.g[i];
    for (int a = 0; a < R2; ++a) all[R1 + a] = O.sol.pau[a];
    for (int i = 0; i < R1 + R2; ++i)
        ap("  %-14s   %s | %s\n", pstr(all[i]).c_str(), bstr(px(all[i])).c_str(), bstr(pz(all[i])).c_str());
    ap("\nRank: %d\n\nCommutation matrix (0 = commute):\n", V.rank);
    for (int i = 0; i < R1 + R2; ++i) { ap("  ");
        for (int j = 0; j < R1 + R2; ++j) ap("%d ", V.comm[i][j]); ap("\n"); }
    ap("\nWeight distribution of the 11 generators:");
    for (int i = 0; i < R1 + R2; ++i) ap(" %d", pwt(all[i]));
    ap("\n\nVERIFICATION\n");
    ap("  all Stage-2 weights == 10 : %s\n", V.wt_ok ? "OK" : "FAIL");
    ap("  [h_a , g_i] = 0           : %s\n", V.comm_g_ok ? "OK" : "FAIL");
    ap("  [h_a , h_b] = 0           : %s\n", V.comm_h_ok ? "OK" : "FAIL");
    ap("  rank = 11                 : %s\n", V.rank_ok ? "OK" : "FAIL");
    ap("  U_3 rebuilt from scratch  : %lld errors\n", (long long)V.u3size);
    ap("  brute-force coverage      : %lld  (optimiser said %lld)\n",
       (long long)V.covered, (long long)O.sol.coverage);
    ap("============================================================\n");
    std::string keep = OUTDIR; OUTDIR = ".";
    dump("full_coverage_solution.txt", true);
    RB.clear(); ap("NONE\n"); dump("full_coverage_remaining_errors.txt");
    RB.clear();
    for (int i = 0; i < R1 + R2; ++i) ap("%s %s\n", bstr(px(all[i])).c_str(), bstr(pz(all[i])).c_str());
    dump("full_coverage_matrix.txt");
    OUTDIR = keep;
}

static void usage() {
    printf("stab14_2stage -- Stage 1: 3 weight-8 generators (proven optimal, ALL optima\n"
           "enumerated); Stage 2: 8 weight-10 generators run once per inequivalent optimum.\n"
           "  --mode heuristic|exact|hybrid    (default hybrid; Stage 1 is always exhaustive)\n"
           "  --threads N   --seed S   --deterministic\n"
           "  --time-limit S            global wall-clock cap\n"
           "  --stage2-time-limit S     heuristic seconds per Stage-1 orbit (default 45)\n"
           "  --exact-time S            branch-and-bound seconds per orbit (default 0)\n"
           "  --max-stage1-solutions N  cap the number of orbits examined (default: all)\n"
           "  --stage1-index N          run only orbit N (1-based)\n"
           "  --continue-after-full     keep going after a full-coverage certificate\n"
           "  --beam-size K --beam-pool N   --restarts N   --selftest N   --save-all\n");
}

int main(int argc, char** argv) {
    std::string mode = "hybrid";
    int nthreads = (int)std::max(1u, std::thread::hardware_concurrency());
    int restarts = 0, selftestN = 200, beamK = 0, beamPool = 40000;
    int maxOrbits = 0, onlyIndex = 0;
    uint64_t seed = 12345ull;
    double stage2Time = 45.0, exactTime = 0.0, globalTL = 1e18;
    bool deterministic = false, continueAfterFull = false, saveAll = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto nxt = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : std::string(); };
        if      (a == "--mode")                   mode = nxt();
        else if (a == "--threads")                nthreads = std::max(1, atoi(nxt().c_str()));
        else if (a == "--seed")                   seed = strtoull(nxt().c_str(), nullptr, 10);
        else if (a == "--time-limit")             globalTL = atof(nxt().c_str());
        else if (a == "--stage2-time-limit")      stage2Time = atof(nxt().c_str());
        else if (a == "--heur-time")              stage2Time = atof(nxt().c_str());
        else if (a == "--exact-time")             exactTime = atof(nxt().c_str());
        else if (a == "--max-stage1-solutions")   maxOrbits = atoi(nxt().c_str());
        else if (a == "--stage1-index")           onlyIndex = atoi(nxt().c_str());
        else if (a == "--continue-after-full")    continueAfterFull = true;
        else if (a == "--beam-size")              beamK = atoi(nxt().c_str());
        else if (a == "--beam-pool")              beamPool = atoi(nxt().c_str());
        else if (a == "--restarts")               restarts = atoi(nxt().c_str());
        else if (a == "--deterministic")          deterministic = true;
        else if (a == "--full-columns")           REQ_FULLCOL = true;
        else if (a == "--selftest")               selftestN = atoi(nxt().c_str());
        else if (a == "--save-all")               saveAll = true;
        else { usage(); return (a == "--help" || a == "-h") ? 0 : 1; }
    }
    (void)saveAll;
    DETERMINISTIC = deterministic;
    if (deterministic && restarts == 0) restarts = 3000;

    T0 = std::chrono::steady_clock::now();
    init_tables();

    // ------------------------------------------------------------------ STAGE 1
    printf("=== STAGE 1: all optimal three-generator weight-8 stabilizers ===\n");
    build_types1(); build_s3(); init_bounds1();
    printf("column types over F_2^3: %d   generator-permutation group S_3: %d elements\n",
           NTYPE1, NS3);
    // pass A -- no S_3 quotient: counts the optimal column multisets themselves
    S1_USE_SYM = false;
    std::memset(s1_c, 0, sizeof(s1_c)); s1_plx = 0; s1_dfs(0, 0);
    const uint64_t rawOpt = S1_RAWOPT;
    const int64_t optSumP = S1_bestSumP;
    const uint64_t rawNodes = S1_nodes;
    // pass B -- with the S_3 quotient: the inequivalent orbits
    S1_ALL.clear(); S1_RAWOPT = 0; S1_have = false; S1_bestSumP = INFP; S1_bestIdle = NQ + 1;
    S1_nodes = S1_pruned = S1_leaves = 0;
    S1_USE_SYM = true;
    std::memset(s1_c, 0, sizeof(s1_c)); s1_plx = 0; s1_dfs(0, 0);
    if (!S1_have || S1_bestSumP != optSumP) { printf("Stage-1 enumeration inconsistent\n"); return 1; }

    const int64_t M3 = TOTAL_ERRORS - ((91778 + optSumP) / 8 - 1);
    printf("\nNumber of raw optimal Stage-1 solutions:\n"
           "  %llu   (distinct optimal multisets of 14 column subspace types, i.e. already\n"
           "         modulo the 14! qubit permutations and the (S_3)^14 local Cliffords)\n",
           (unsigned long long)rawOpt);
    printf("\nNumber of inequivalent optimal Stage-1 solutions:\n"
           "  %zu   (additionally modulo the S_3 permutations of the three generators)\n",
           S1_ALL.size());
    printf("\nStage-1 optimum: M_3 = %lld / %lld,  |U_3| = %lld   [GLOBAL OPTIMUM PROVEN,"
           " tree exhausted: %llu nodes with symmetry, %llu without]\n",
           (long long)M3, (long long)TOTAL_ERRORS, (long long)(TOTAL_ERRORS - M3),
           (unsigned long long)S1_nodes, (unsigned long long)rawNodes);
    fflush(stdout);

    std::vector<std::array<uint8_t, NQ>> orbits(S1_ALL.begin(), S1_ALL.end());
    if (onlyIndex > 0) {
        if (onlyIndex > (int)orbits.size()) { printf("--stage1-index out of range\n"); return 1; }
        std::array<uint8_t, NQ> keep = orbits[onlyIndex - 1];
        orbits.assign(1, keep);
    } else if (maxOrbits > 0 && (int)orbits.size() > maxOrbits) {
        orbits.resize(maxOrbits);
        printf("(limited to the first %d orbits by --max-stage1-solutions)\n", maxOrbits);
    }

    // ------------------------------------------------------------------ STAGE 2, per orbit
    std::vector<OrbitResult> results;
    bool foundFull = false;
    for (size_t oi = 0; oi < orbits.size(); ++oi) {
        if (elapsed() > globalTL) { printf("\nglobal time limit reached, stopping orbit loop\n"); break; }
        const int idx = (onlyIndex > 0) ? onlyIndex : (int)oi + 1;
        char dirbuf[64]; snprintf(dirbuf, sizeof(dirbuf), "stage1_orbit_%04d", idx);
        std::error_code ec; std::filesystem::create_directories(dirbuf, ec);
        OUTDIR = dirbuf;

        reset_stage2();
        OrbitResult O; O.index = idx; O.types = orbits[oi];
        s1_build_matrix(O.types.data(), G1);
        for (int i = 0; i < R1; ++i) O.g[i] = G1[i];
        S1Verify V1 = s1_verify(G1, &U3);
        for (uint32_t E : U3) U3_BY_W[pwt(E)]++;
        O.M3 = V1.det; O.U3n = V1.und;
        if (!V1.ok_w || !V1.ok_comm || !V1.ok_rank || V1.det != M3) {
            printf("orbit %d FAILED Stage-1 verification\n", idx); continue; }
        std::copy(O.types.begin(), O.types.end(), S1_bestTypes);
        S1_bestSumP = optSumP;
        s1_report(G1, V1, true);
        s1_save_undetected();

        printf("\n--- orbit %04d / %zu ---  g1=%s g2=%s g3=%s   |U_3| = %lld\n",
               idx, orbits.size(), pstr(O.g[0]).c_str(), pstr(O.g[1]).c_str(),
               pstr(O.g[2]).c_str(), (long long)O.U3n);
        fflush(stdout);

        build_symplectic_basis2();
        build_cov1();
        build_maskbase();
        build_repr10();
        CAND.reserve(NCLS10);
        for (uint32_t v = 1; v < (uint32_t)VSZ2; ++v) if (REPR10[v] != NONE) CAND.push_back(v);
        std::stable_sort(CAND.begin(), CAND.end(), [](uint32_t a, uint32_t b) {
            if (COV1[a] != COV1[b]) return COV1[a] > COV1[b]; return a < b; });
        const int64_t UB8 = std::min<int64_t>(O.U3n, (int64_t(255) * COVMAX) / 128);
        printf("    C(S3)/S3 dim %d, %zu weight-10 candidate classes, max cov1 = %d, UB(M_8) = %lld\n",
               VDIM2, CAND.size(), COVMAX, (long long)UB8);
        fflush(stdout);

        if (oi == 0 && selftestN > 0) {
            printf("    [selftest] %d checks (cov1 table, class map, (1/128) identity vs bitset union) ... ",
                   selftestN); fflush(stdout);
            if (!self_test(selftestN, seed)) { printf("FAILED\n"); return 3; }
            printf("all pass\n");
        }
        if (beamK > 0) beam_search(beamK, beamPool, nthreads);
        if (mode == "heuristic" || mode == "hybrid")
            heuristic(nthreads, seed + 1000ull * idx, restarts,
                      std::min(elapsed() + stage2Time, globalTL));
        bool proven = (BESTCOV.load() >= UB8);
        if ((mode == "exact" || mode == "hybrid") && !proven && exactTime > 0) {
            TIME_LIMIT = std::min(globalTL, elapsed() + exactTime);
            exact_search(nthreads);
            if (!ABORT.load()) proven = true;
        }
        if (!BEST.valid) { printf("    no feasible 8-generator extension found\n"); continue; }
        if (BESTCOV.load() >= UB8) proven = true;

        S2Verify V2 = s2_verify(BEST);
        O.sol = BEST; O.M8 = BEST.coverage; O.Mfinal = V2.det11; O.proven = proven;
        O.leftover = V2.leftover;
        Residual r = analyse_residual(V2.leftover);
        O.residual = r.n; for (int w = 1; w <= MAXW; ++w) O.res_w[w] = r.by_w[w];
        O.res_dim = r.dim; O.res_subgroup = r.is_subgroup; O.res_allx = r.all_x;
        // weight-4 statistics for the complete 11-generator code, computed independently
        O.w4_total = TOTW[MAXW];
        { int64_t det = 0;
          uint32_t all[R1 + R2];
          for (int i = 0; i < R1; ++i) all[i] = G1[i];
          for (int a = 0; a < R2; ++a) all[R1 + a] = BEST.pau[a];
          for_each_error([&](uint32_t E, int w) {
              if (w != MAXW) return;
              int s = 0;
              for (int i = 0; i < R1 + R2; ++i) if (symp(all[i], E)) { s = 1; break; }
              det += s; });
          O.w4_det = det; O.w4_und = O.w4_total - det; O.w4_all = (O.w4_und == 0); }
        O.full = (O.M8 == O.U3n);

        s2_report(BEST, V2, proven, UB8);
        RB.clear(); append_residual_block(V2.leftover, r); dump("stage2_remaining.txt");
        final_summary(BEST, V1, V2, true, proven);
        // per-orbit summary.txt in the requested shape
        RB.clear();
        ap("Stage-1 orbit:\n%04d of %zu\n\n", idx, orbits.size());
        ap("Stage-1 generators:\n");
        for (int i = 0; i < R1; ++i) ap("g%d = %s\n", i + 1, pstr(O.g[i]).c_str());
        ap("\nStage-1 detected:\n%lld / %lld\n\n", (long long)O.M3, (long long)TOTAL_ERRORS);
        ap("Stage-1 undetected:\n%lld\n\n", (long long)O.U3n);
        ap("Stage-2 best:\n");
        for (int a = 0; a < R2; ++a) ap("h%d = %s\n", a + 1, pstr(O.sol.pau[a]).c_str());
        ap("\nStage-2 coverage:\n%lld / %lld\n\n", (long long)O.M8, (long long)O.U3n);
        ap("Final coverage:\n%lld / %lld\n\n", (long long)O.Mfinal, (long long)TOTAL_ERRORS);
        append_residual_block(V2.leftover, r);
        ap("\nWeight-4 errors (recomputed programmatically, not hard-coded):\n");
        ap("Total:      %lld\nDetected:   %lld\nUndetected: %lld\n",
           (long long)O.w4_total, (long long)O.w4_det, (long long)O.w4_und);
        ap("All weight-4 errors detected:\n%s\n\n", O.w4_all ? "YES" : "NO");
        ap("Stage-2 optimum:\n%s\n", proven ? "PROVEN" : "BEST KNOWN");
        dump("summary.txt");

        printf("    M_8 = %lld / %lld   residual %lld (w %lld/%lld/%lld/%lld, dim %d%s%s)"
               "   M_final = %lld   w4-undetected %lld%s\n",
               (long long)O.M8, (long long)O.U3n, (long long)O.residual,
               (long long)O.res_w[1], (long long)O.res_w[2], (long long)O.res_w[3],
               (long long)O.res_w[4], O.res_dim, r.is_subgroup ? ", subgroup" : "",
               r.all_x ? ", X-only" : "", (long long)O.Mfinal, (long long)O.w4_und,
               O.w4_all ? "  <== ALL WEIGHT-4 DETECTED" : "");
        fflush(stdout);
        if (O.w4_all) {
            printf("============================================================\n"
                   "ALL WEIGHT-4 ERRORS DETECTED\n"
                   "============================================================\n"
                   "Weight-4 errors:\nDetected: %lld\nUndetected: 0\n"
                   "(orbit %04d)\n", (long long)O.w4_det, idx);
        }
        results.push_back(O);
        if (O.full) {
            foundFull = true;
            print_full_coverage_banner(O, V2);
            if (!continueAfterFull) { printf("\nstopping after the full-coverage certificate"
                                             " (use --continue-after-full to keep going)\n"); break; }
        }
    }
    OUTDIR = ".";

    // ------------------------------------------------------------------ global comparison
    std::vector<OrbitResult> sorted = results;
    std::sort(sorted.begin(), sorted.end(), [](const OrbitResult& a, const OrbitResult& b) {
        if (a.M8 != b.M8) return a.M8 > b.M8;
        if (a.residual != b.residual) return a.residual < b.residual;
        if (a.res_w[MAXW] != b.res_w[MAXW]) return a.res_w[MAXW] < b.res_w[MAXW];
        return a.index < b.index; });
    RB.clear();
    ap("================================================================\n");
    ap("COMPARISON OF INEQUIVALENT OPTIMAL STAGE-1 SOLUTIONS\n");
    ap("================================================================\n\n");
    ap("Orbit     M3       |U3|     M8       Final     Residual  W4-undetected  Stage2\n");
    ap("----------------------------------------------------------------------------\n");
    for (const auto& o : sorted)
        ap("%04d      %-8lld %-8lld %-8lld %-9lld %-9lld %-14lld %s\n", o.index,
           (long long)o.M3, (long long)o.U3n, (long long)o.M8, (long long)o.Mfinal,
           (long long)o.residual, (long long)o.w4_und, o.proven ? "PROVEN" : "best known");
    ap("================================================================\n\n");
    if (!results.empty()) {
        const OrbitResult& B = sorted.front();
        ap("GLOBAL BEST KNOWN over the tested Stage-1 orbits:\n");
        ap("  orbit %04d,  M_8 = %lld / %lld,  M_final = %lld / %lld  (%.8f %%)\n",
           B.index, (long long)B.M8, (long long)B.U3n, (long long)B.Mfinal,
           (long long)TOTAL_ERRORS, 100.0 * double(B.Mfinal) / double(TOTAL_ERRORS));
        ap("  Stage-1 generators:\n");
        for (int i = 0; i < R1; ++i) ap("    g%d = %s\n", i + 1, pstr(B.g[i]).c_str());
        ap("  Stage-2 generators:\n");
        for (int a = 0; a < R2; ++a) ap("    h%d = %s\n", a + 1, pstr(B.sol.pau[a]).c_str());
        ap("\n");
        Residual br = analyse_residual(B.leftover);
        append_residual_block(B.leftover, br);
        ap("\n");
        if (foundFull) {
            ap("FULL COVERAGE FOUND\n\n");
            ap("An optimal Stage-1 solution exists whose eight weight-10 extensions detect\n");
            ap("all %lld Stage-1-undetected errors.\n\nFinal coverage = %lld / %lld.\n",
               (long long)B.U3n, (long long)B.Mfinal, (long long)TOTAL_ERRORS);
        } else {
            ap("NO FULL COVERAGE FOUND\n\n");
            ap("Across all %zu tested inequivalent globally optimal Stage-1 solutions, the best\n"
               "Stage-2 coverage is %lld / %lld, leaving %lld residual errors (listed above).\n",
               results.size(), (long long)B.M8, (long long)B.U3n, (long long)B.residual);
            ap("\nStage 1 is PROVEN globally optimal and its optimal family was enumerated\n"
               "exhaustively.  Each Stage-2 search is heuristic unless marked PROVEN, so this\n"
               "is a BEST-KNOWN result, not a proof of impossibility.\n");
        }
        bool anyW4 = false; for (const auto& o : results) if (o.w4_all) anyW4 = true;
        ap("\nAll weight-4 errors detected by some orbit: %s\n", anyW4 ? "YES" : "NO");
        ap("Weight-4 error total (computed, not hard-coded): %lld\n", (long long)TOTW[MAXW]);
    }
    ap("\nElapsed: %.2f s\n", elapsed());
    dump("orbit_comparison.txt", true);
    return 0;
}
