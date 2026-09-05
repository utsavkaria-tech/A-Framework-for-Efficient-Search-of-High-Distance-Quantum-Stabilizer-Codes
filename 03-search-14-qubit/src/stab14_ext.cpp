// ===========================================================================================
//  stab14_ext.cpp -- seven additional weight-10 stabilizer generators that maximally cover
//                    the 5430 weight-<=4 Pauli errors left undetected by a FIXED [[14,*]]
//                    four-generator stabilizer S0.
//
//  ------------------------------------------------------------------------------------
//  WHAT THE PROBLEM ACTUALLY IS (the reduction the whole program is built on)
//  ------------------------------------------------------------------------------------
//  S0 = <g1,g2,g3,g4> is a fixed 4-dimensional isotropic subspace of F_2^28 (Pauli group
//  modulo phases, symplectic form <a,b> = xa.zb + za.xb).
//
//  (R1)  Every admissible new generator h and every target error E lies in the centralizer
//        C(S0) = S0^perp, dim 28-4 = 24.  <h,E> is unchanged when h or E is shifted by an
//        element of S0, so the whole problem lives in the 20-dimensional NON-DEGENERATE
//        symplectic quotient  V = C(S0)/S0.  Everything -- coverage, mutual commutation,
//        independence -- is a function of the CLASSES only.  We therefore work with 20-bit
//        integers, using a symplectic basis u_1..u_10, w_1..w_10 of V:
//              class(v) = ( <w_i,v> )_i  |  ( <u_i,v> )_i << 10 ,
//              <v,v'>   = parity( (a & b') ^ (b & a') ).
//
//  (R2)  The seven new generators must be independent and pairwise commuting, i.e. their
//        classes span a 7-dimensional TOTALLY ISOTROPIC subspace D <= V, and rank 11 is
//        exactly dim D = 7.  A target E is UNDETECTED by the extension iff its class lies
//        in D^perp.  So the objective depends only on the subspace D, not on the basis.
//
//  (R3)  THE KEY IDENTITY.  Let cov1(v) = #{ E in U : <v,E> = 1 } for a class v.  Poisson
//        summation on V (equivalently, the Walsh transform of the target-class multiset)
//        gives, for any 7-dimensional D,
//
//              #covered  =  (1/64) * SUM_{v in D \ {0}} cov1(v).                    (*)
//
//        Proof: #uncovered = sum_{p in D^perp} f(p) = 2^-7 sum_{v in sigma(D)} fhat(v),
//        and fhat(sigma(w)) = 5430 - 2 cov1(w).  Hence #uncovered = 5430 - (1/64) sum cov1.
//        Scoring a complete 7-generator set is therefore 127 lookups in a 2^20 table -- no
//        5430-bit union, no symplectic products.  (The bitset engine is still implemented
//        and used for INDEPENDENT verification; the two agree on every self-test.)
//
//  (R4)  Equivalent code-theoretic statement.  L = S0 + span(h) is an 11-dimensional
//        isotropic subspace, i.e. an [[14,3]] stabilizer code containing S0, and
//              #uncovered = (1/2048) [ 91770 + SUM_{h in L, h != 0} P(wt(h)) ],
//        with the same P(m) table as the first-stage program.  Full coverage (5430/5430)
//        happens iff N(L) contains no non-zero Pauli of weight <= 4, i.e. iff L is a PURE
//        [[14,3,5]] code containing S0.  (Sanity check of the identity: for L = S0 itself,
//        (1/16)[91770 + 15*(-326)] = 5430, the known target count.)
//
//  (R5)  cov1 is computed for ALL 2^20 classes at once with a Walsh-Hadamard transform
//        (20 * 2^20 additions), and the set of classes reachable by a weight-10 (or 12)
//        Pauli is found by a Gray-code sweep of all 2^24 elements of C(S0).
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

// ------------------------------------------------------------------ problem constants
constexpr int NQ      = 14;                 // physical qubits
constexpr int NG0     = 4;                  // fixed generators
constexpr int NEW     = 7;                  // new generators to find
constexpr int MAXW    = 4;                  // error weight considered
constexpr int VDIM    = 2 * NQ - 2 * NG0;   // dim C(S0)/S0 = 20
constexpr int VSZ     = 1 << VDIM;          // 2^20 classes
constexpr int CDIM    = 2 * NQ - NG0;       // dim C(S0) = 24
constexpr uint32_t QM  = (1u << NQ) - 1u;   // 14-bit qubit mask
constexpr uint32_t NONE = 0xFFFFFFFFu;
constexpr int TARGETS_EXPECT = 5430;

// A Pauli (mod phase) is packed as  x | (z << 16),  x,z in [0, 2^14).
static inline uint32_t pk(uint32_t x, uint32_t z) { return x | (z << 16); }
static inline uint32_t px(uint32_t p) { return p & QM; }
static inline uint32_t pz(uint32_t p) { return (p >> 16) & QM; }
static inline int pwt(uint32_t p) { return std::popcount((p | (p >> 16)) & QM); }
static inline int symp(uint32_t a, uint32_t b) {           // symplectic product mod 2
    return (std::popcount((px(a) & pz(b)) ^ (pz(a) & px(b)))) & 1;
}
static const char PAULI_CH[2][2] = { {'I','Z'}, {'X','Y'} };   // [xbit][zbit]
static std::string pstr(uint32_t p) {
    std::string s; s.reserve(NQ);
    for (int j = 0; j < NQ; ++j) s += PAULI_CH[(px(p) >> j) & 1][(pz(p) >> j) & 1];
    return s;
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

// ------------------------------------------------------------------ the FIXED stabilizer
static const char* G0_STR[NG0] = {
    "XIIIXXIXIXIXYY",
    "IXIIXIXIXIXXYY",
    "IIXIIXXIIXXXXX",
    "IIIXIIIXXXXXZZ"
};
static uint32_t G0[NG0];

// ------------------------------------------------------------------ P(m) table (see R4)
static int64_t Cbin[32][32], Ptab[NQ + 1], Ew[NQ + 1][MAXW + 1], TOTW[MAXW + 1], TOTAL_ERRORS;
static void init_P() {
    for (int n = 0; n < 32; ++n) { Cbin[n][0] = 1;
        for (int k = 1; k <= n; ++k) Cbin[n][k] = Cbin[n-1][k-1] + (k <= n-1 ? Cbin[n-1][k] : 0); }
    auto C = [](int n, int k) { return (k < 0 || k > n) ? (int64_t)0 : Cbin[n][k]; };
    int64_t p3[MAXW + 1]; p3[0] = 1; for (int w = 1; w <= MAXW; ++w) p3[w] = p3[w-1] * 3;
    TOTAL_ERRORS = 0;
    for (int w = 0; w <= MAXW; ++w) { TOTW[w] = C(NQ, w) * p3[w]; if (w) TOTAL_ERRORS += TOTW[w]; }
    for (int m = 0; m <= NQ; ++m) {
        Ptab[m] = 0;
        for (int w = 0; w <= MAXW; ++w) {
            int64_t s = 0;
            for (int k = 0; k <= w; ++k) { int64_t t = C(m,k) * C(NQ-m, w-k) * p3[w-k]; s += (k&1) ? -t : t; }
            Ew[m][w] = s; if (w) Ptab[m] += s;
        }
    }
}

// ===========================================================================================
//  STEP 1 -- reconstruct the target set U from the four fixed generators (never trusted in)
// ===========================================================================================
static std::vector<uint32_t> TGT;        // the 5430 target Paulis
static std::vector<uint8_t>  TGTW;       // their Pauli weights
static int64_t TGT_BY_W[MAXW + 1] = {0};

static void build_targets() {
    static const uint32_t PXb[3] = {1,1,0}, PZb[3] = {0,1,1};   // X, Y, Z
    int idx[MAXW];
    for (int w = 1; w <= MAXW; ++w) {
        for (int i = 0; i < w; ++i) idx[i] = i;
        for (;;) {
            int npat = 1; for (int i = 0; i < w; ++i) npat *= 3;
            for (int p = 0; p < npat; ++p) {
                uint32_t x = 0, z = 0; int q = p;
                for (int i = 0; i < w; ++i) { int a = q % 3; q /= 3;
                    x |= PXb[a] << idx[i]; z |= PZb[a] << idx[i]; }
                uint32_t E = pk(x, z);
                int s = 0; for (int g = 0; g < NG0; ++g) s |= symp(G0[g], E) << g;
                if (s == 0) { TGT.push_back(E); TGTW.push_back(uint8_t(w)); TGT_BY_W[w]++; }
            }
            int i = w - 1; while (i >= 0 && idx[i] == NQ - w + i) --i;
            if (i < 0) break;
            ++idx[i]; for (int k = i + 1; k < w; ++k) idx[k] = idx[k-1] + 1;
        }
    }
}

// ===========================================================================================
//  STEP 2 -- symplectic basis of V = C(S0)/S0  and the class map
// ===========================================================================================
static uint32_t UB[VDIM / 2], WB[VDIM / 2];        // u_1..u_10 , w_1..w_10
static uint32_t S0EL[1 << NG0];                    // the 16 elements of S0

static inline uint32_t cls_of(uint32_t v) {        // v must lie in C(S0)
    uint32_t c = 0;
    for (int i = 0; i < VDIM / 2; ++i) {
        c |= uint32_t(symp(WB[i], v)) << i;
        c |= uint32_t(symp(UB[i], v)) << (i + VDIM / 2);
    }
    return c;
}
static inline uint32_t lift_of(uint32_t c) {       // a canonical Pauli representative
    uint32_t v = 0;
    for (int i = 0; i < VDIM / 2; ++i) {
        if ((c >> i) & 1)                v ^= UB[i];
        if ((c >> (i + VDIM / 2)) & 1)   v ^= WB[i];
    }
    return v;
}
// symplectic form directly on 20-bit classes
static inline int csymp(uint32_t a, uint32_t b) {
    constexpr uint32_t LO = (1u << (VDIM / 2)) - 1u;
    uint32_t t = ((a & LO) & (b >> (VDIM/2))) ^ ((a >> (VDIM/2)) & (b & LO));
    return std::popcount(t) & 1;
}

static void build_symplectic_basis() {
    // basis of C(S0) = { v : <g_i,v> = 0 } by Gaussian elimination on the 4 x 28 system
    std::vector<uint32_t> cur;
    {
        // 28 coordinate directions: x_j (bit j) and z_j (bit j+16)
        std::vector<uint32_t> gen;
        for (int j = 0; j < NQ; ++j) { gen.push_back(pk(1u << j, 0)); gen.push_back(pk(0, 1u << j)); }
        // reduce the 4 syndrome functionals
        std::vector<uint32_t> basis = gen;         // start from all of F_2^28
        for (int g = 0; g < NG0; ++g) {
            int piv = -1;
            for (size_t i = 0; i < basis.size(); ++i) if (symp(G0[g], basis[i])) { piv = (int)i; break; }
            if (piv < 0) { fprintf(stderr, "g%d is dependent -- rank < 4\n", g + 1); exit(1); }
            uint32_t pv = basis[piv]; basis.erase(basis.begin() + piv);
            for (auto& b : basis) if (symp(G0[g], b)) b ^= pv;
        }
        cur = basis;                                // 24 vectors spanning C(S0)
    }
    if (cur.size() != CDIM) { fprintf(stderr, "centralizer dim %zu != 24\n", cur.size()); exit(1); }
    // symplectic Gram-Schmidt: peel off hyperbolic pairs; the radical left over must be S0
    int np = 0;
    while (np < VDIM / 2) {
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
    if (np != VDIM / 2) { fprintf(stderr, "symplectic rank %d != 10\n", np); exit(1); }
    if (cur.size() != NG0) { fprintf(stderr, "radical dim %zu != 4\n", cur.size()); exit(1); }
    // the radical of the form on C(S0) is exactly S0 -- verify that it really is
    for (uint32_t r : cur) {
        // r must be a combination of g1..g4 : brute force over the 16 elements
        bool found = false;
        for (int m = 0; m < 16; ++m) { uint32_t t = 0;
            for (int g = 0; g < NG0; ++g) if ((m >> g) & 1) t ^= G0[g];
            if (t == r) { found = true; break; } }
        if (!found) { fprintf(stderr, "radical is not S0 -- basis construction failed\n"); exit(1); }
    }
    for (int m = 0; m < 16; ++m) { uint32_t t = 0;
        for (int g = 0; g < NG0; ++g) if ((m >> g) & 1) t ^= G0[g];
        S0EL[m] = t; }
}

// ===========================================================================================
//  STEP 3 -- cov1[v] for ALL 2^20 classes via one Walsh-Hadamard transform, and the sets of
//            classes reachable by a weight-10 / weight-12 Pauli via a Gray-code sweep of the
//            2^24 elements of C(S0).
// ===========================================================================================
static std::vector<uint16_t> COV1;        // cov1[v] = # targets anticommuting with class v
static std::vector<uint32_t> REPR10, REPR12;   // a canonical Pauli of that weight, or NONE
static std::vector<uint32_t> TCLS;        // class of each target
static int COVMAX = 0;                    // max_v cov1(v)   (admissible-bound ingredient)

static void build_cov1() {
    TCLS.resize(TGT.size());
    std::vector<int32_t> f(VSZ, 0);
    constexpr uint32_t LO = (1u << (VDIM / 2)) - 1u;
    for (size_t j = 0; j < TGT.size(); ++j) {
        uint32_t c = cls_of(TGT[j]); TCLS[j] = c;
        uint32_t s = ((c >> (VDIM/2)) & LO) | ((c & LO) << (VDIM/2));   // sigma(c)
        f[s] += 1;
    }
    for (int b = 0; b < VDIM; ++b) {                    // in-place Walsh-Hadamard
        int step = 1 << b;
        for (int i = 0; i < VSZ; i += step << 1)
            for (int k = i; k < i + step; ++k) {
                int32_t a = f[k], d = f[k + step];
                f[k] = a + d; f[k + step] = a - d;
            }
    }
    const int32_t N = (int32_t)TGT.size();
    COV1.resize(VSZ);
    for (int v = 0; v < VSZ; ++v) { int32_t c = (N - f[v]) / 2; COV1[v] = uint16_t(c);
        if (v && c > COVMAX) COVMAX = c; }
}

static int64_t N_W10 = 0, N_W12 = 0, NCLS10 = 0, NCLS12 = 0;
static void build_reprs() {
    REPR10.assign(VSZ, NONE); REPR12.assign(VSZ, NONE);
    uint32_t basis[CDIM], bcls[CDIM];
    for (int g = 0; g < NG0; ++g) { basis[g] = G0[g]; bcls[g] = 0; }
    for (int i = 0; i < VDIM / 2; ++i) { basis[NG0 + i] = UB[i]; bcls[NG0 + i] = 1u << i; }
    for (int i = 0; i < VDIM / 2; ++i) { basis[NG0 + VDIM/2 + i] = WB[i]; bcls[NG0 + VDIM/2 + i] = 1u << (i + VDIM/2); }
    uint32_t cur = 0, ccl = 0;
    const uint64_t TOT = 1ull << CDIM;
    for (uint64_t i = 1; i < TOT; ++i) {
        int b = std::countr_zero(i);
        cur ^= basis[b]; ccl ^= bcls[b];
        int w = pwt(cur);
        if (w == 10) { ++N_W10; if (ccl && (REPR10[ccl] == NONE || cur < REPR10[ccl])) {
                           if (REPR10[ccl] == NONE) ++NCLS10; REPR10[ccl] = cur; } }
        else if (w == 12) { ++N_W12; if (ccl && (REPR12[ccl] == NONE || cur < REPR12[ccl])) {
                           if (REPR12[ccl] == NONE) ++NCLS12; REPR12[ccl] = cur; } }
    }
}

// ===========================================================================================
//  STEP 4 -- the bitset engine (used for INDEPENDENT verification and for listing the
//            errors that stay uncovered).  TARGETS bits, one word-array per generator.
// ===========================================================================================
static int TARGETS = 0, WORDS = 0;
static std::vector<uint64_t> MASKBASE;    // VDIM masks: coverage of the basis classes

static void build_maskbase() {
    TARGETS = (int)TGT.size(); WORDS = (TARGETS + 63) / 64;
    MASKBASE.assign(size_t(VDIM) * WORDS, 0);
    for (int b = 0; b < VDIM; ++b) {
        uint64_t* m = &MASKBASE[size_t(b) * WORDS];
        uint32_t v = 1u << b;
        for (int j = 0; j < TARGETS; ++j) if (csymp(v, TCLS[j])) m[j >> 6] |= 1ull << (j & 63);
    }
}
// coverage bitset of a single class, built by linearity: mask(v) = XOR of the basis masks
static void class_mask(uint32_t v, uint64_t* out) {
    std::memset(out, 0, size_t(WORDS) * 8);
    while (v) { int b = std::countr_zero(v); v &= v - 1;
        const uint64_t* m = &MASKBASE[size_t(b) * WORDS];
        for (int k = 0; k < WORDS; ++k) out[k] ^= m[k]; }
}

// ===========================================================================================
//  STEP 5 -- scoring.  score_span() is the fast algebraic route (identity (*) above);
//            verify_bruteforce() is a completely independent recomputation.
// ===========================================================================================
struct Sol {
    uint32_t cls[NEW]{};        // the 7 classes in V (they span the 7-dim isotropic D)
    uint32_t pau[NEW]{};        // the concrete Pauli chosen for each class
    int      wt[NEW]{};         // its Pauli weight (10 or 12)
    int64_t  coverage = -1;
    int      nw12 = 0;
    int      ncol1 = 0;        // physical qubits whose 11-generator column is not 2-dimensional
    int64_t  sumIndiv = 0;
    int      minIndiv = 0;
    bool     valid = false;
};

// coverage of the subspace spanned by the d classes, using identity (*) generalised to
// dimension d:   #covered = (1/2^{d-1}) * SUM_{v in span\0} cov1(v)
static inline int64_t span_sum(const uint32_t* cls, int d) {
    int64_t s = 0;
    for (int m = 1; m < (1 << d); ++m) {
        uint32_t v = 0; int mm = m;
        while (mm) { int b = std::countr_zero(unsigned(mm)); mm &= mm - 1; v ^= cls[b]; }
        s += COV1[v];
    }
    return s;
}
static inline int64_t coverage_of(const uint32_t* cls, int d) {
    return span_sum(cls, d) / (1LL << (d - 1));
}

// ---------------------------------------------------------------- independent verification
struct Verify {
    bool wt_ok = true, comm_g_ok = true, comm_h_ok = true, rank_ok = true;
    int  rank = 0;
    int  weights[NEW]{};
    int  comm[NG0 + NEW][NG0 + NEW]{};
    int64_t detected = 0, per_gen[NEW]{};
    int64_t det_by_w[MAXW + 1]{}, tot_by_w[MAXW + 1]{};
    std::vector<uint32_t> leftover;              // target errors still uncovered
    // full 11-generator statistics (spec item 29)
    int64_t all_total = 0, all_orig_det = 0, all_orig_und = 0, all_new_det = 0,
            all_still_und = 0, all_det11 = 0;
};

static Verify verify_bruteforce(const Sol& S, const int* allowedW, int nAllowedW) {
    Verify V;
    for (int a = 0; a < NEW; ++a) {
        V.weights[a] = pwt(S.pau[a]);
        bool ok = false; for (int i = 0; i < nAllowedW; ++i) if (V.weights[a] == allowedW[i]) ok = true;
        if (!ok) V.wt_ok = false;
    }
    uint32_t all[NG0 + NEW];
    for (int i = 0; i < NG0; ++i) all[i] = G0[i];
    for (int a = 0; a < NEW; ++a) all[NG0 + a] = S.pau[a];
    for (int i = 0; i < NG0 + NEW; ++i) for (int j = 0; j < NG0 + NEW; ++j) {
        V.comm[i][j] = symp(all[i], all[j]);
        if (i != j && V.comm[i][j]) { if (i < NG0 || j < NG0) V.comm_g_ok = false; else V.comm_h_ok = false; }
    }
    // rank of the 11 x 28 matrix over F_2
    { uint32_t rows[NG0 + NEW]; for (int i = 0; i < NG0 + NEW; ++i)
          rows[i] = (px(all[i])) | (pz(all[i]) << NQ);
      int r = 0;
      for (int b = 0; b < 2 * NQ && r < NG0 + NEW; ++b) {
          int p = -1; for (int i = r; i < NG0 + NEW; ++i) if ((rows[i] >> b) & 1) { p = i; break; }
          if (p < 0) continue;
          std::swap(rows[r], rows[p]);
          for (int i = 0; i < NG0 + NEW; ++i) if (i != r && ((rows[i] >> b) & 1)) rows[i] ^= rows[r];
          ++r; }
      V.rank = r; V.rank_ok = (r == NG0 + NEW); }
    // re-enumerate every weight 1..4 error from scratch and score directly
    static const uint32_t PXb[3] = {1,1,0}, PZb[3] = {0,1,1};
    int idx[MAXW];
    for (int w = 1; w <= MAXW; ++w) {
        for (int i = 0; i < w; ++i) idx[i] = i;
        for (;;) {
            int npat = 1; for (int i = 0; i < w; ++i) npat *= 3;
            for (int p = 0; p < npat; ++p) {
                uint32_t x = 0, z = 0; int q = p;
                for (int i = 0; i < w; ++i) { int a = q % 3; q /= 3;
                    x |= PXb[a] << idx[i]; z |= PZb[a] << idx[i]; }
                uint32_t E = pk(x, z);
                int s0 = 0; for (int g = 0; g < NG0; ++g) s0 |= symp(G0[g], E) << g;
                int s1 = 0; for (int a = 0; a < NEW; ++a) s1 |= symp(S.pau[a], E) << a;
                V.all_total++;
                if (s0) V.all_orig_det++; else V.all_orig_und++;
                if (s0 || s1) V.all_det11++; 
                if (!s0) {
                    V.tot_by_w[w]++;
                    if (s1) { V.detected++; V.det_by_w[w]++; V.all_new_det++;
                              for (int a = 0; a < NEW; ++a) if ((s1 >> a) & 1) V.per_gen[a]++; }
                    else { V.leftover.push_back(E); V.all_still_und++; }
                }
            }
            int i = w - 1; while (i >= 0 && idx[i] == NQ - w + i) --i;
            if (i < 0) break;
            ++idx[i]; for (int k = i + 1; k < w; ++k) idx[k] = idx[k-1] + 1;
        }
    }
    return V;
}

// ===========================================================================================
//  STEP 6 -- candidate set, incumbent bookkeeping and the total order used for tie-breaks
// ===========================================================================================
static std::vector<uint32_t> CAND;        // candidate classes, sorted by cov1 descending
static std::vector<uint8_t>  CANDW12;     // 1 if this class needs a weight-12 representative
static int  WMODE = 10;                   // 10  or  1012
static int  COVMAX_CAND = 0;
static std::chrono::steady_clock::time_point T0;
static double elapsed() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count(); }

static void build_candidates(int candTop) {
    CAND.clear(); CANDW12.clear();
    for (uint32_t v = 1; v < (uint32_t)VSZ; ++v) {
        bool h10 = REPR10[v] != NONE, h12 = (WMODE == 1012) && REPR12[v] != NONE;
        if (!h10 && !h12) continue;
        CAND.push_back(v); CANDW12.push_back(uint8_t(h10 ? 0 : 1));
    }
    std::vector<uint32_t> ord(CAND.size()); std::iota(ord.begin(), ord.end(), 0u);
    std::stable_sort(ord.begin(), ord.end(), [](uint32_t a, uint32_t b) {
        if (COV1[CAND[a]] != COV1[CAND[b]]) return COV1[CAND[a]] > COV1[CAND[b]];
        return CAND[a] < CAND[b]; });
    std::vector<uint32_t> c2(CAND.size()); std::vector<uint8_t> w2(CAND.size());
    for (size_t i = 0; i < ord.size(); ++i) { c2[i] = CAND[ord[i]]; w2[i] = CANDW12[ord[i]]; }
    CAND.swap(c2); CANDW12.swap(w2);
    if (candTop > 0 && (int)CAND.size() > candTop) { CAND.resize(candTop); CANDW12.resize(candTop); }
    COVMAX_CAND = CAND.empty() ? 0 : COV1[CAND[0]];
}

// Number of physical qubits whose column of the FULL 11-generator matrix is not
// 2-dimensional.  A weight-w Krawtchouk computation gives B_1 = 3*n0 + n1 exactly, i.e. this
// count IS the number of undetected weight-1 target errors, so driving it to zero is a
// necessary condition for full coverage.  Used only by --full-columns.
static bool REQ_FULLCOL = false;
static int deficient_columns(const uint32_t* pau) {
    uint32_t all[NG0 + NEW];
    for (int i = 0; i < NG0; ++i) all[i] = G0[i];
    for (int a = 0; a < NEW; ++a) all[NG0 + a] = pau[a];
    int bad = 0;
    for (int j = 0; j < NQ; ++j) {
        uint32_t xc = 0, zc = 0;
        for (int i = 0; i < NG0 + NEW; ++i) {
            xc |= ((px(all[i]) >> j) & 1u) << i;
            zc |= ((pz(all[i]) >> j) & 1u) << i;
        }
        if (xc == 0 || zc == 0 || xc == zc) ++bad;      // span{x,z} has dimension < 2
    }
    return bad;
}

static Sol finalize_sol(const uint32_t* cls) {
    Sol S; std::copy(cls, cls + NEW, S.cls);
    std::sort(S.cls, S.cls + NEW);
    S.coverage = coverage_of(S.cls, NEW);
    S.sumIndiv = 0; S.minIndiv = 1 << 30; S.nw12 = 0;
    for (int a = 0; a < NEW; ++a) {
        uint32_t c = S.cls[a];
        if (REPR10[c] != NONE) { S.pau[a] = REPR10[c]; S.wt[a] = 10; }
        else { S.pau[a] = REPR12[c]; S.wt[a] = 12; ++S.nw12; }
        S.sumIndiv += COV1[c]; S.minIndiv = std::min<int>(S.minIndiv, COV1[c]);
    }
    S.ncol1 = REQ_FULLCOL ? deficient_columns(S.pau) : 0;
    S.valid = true; return S;
}
// strict total order: coverage > fewer weight-12 > larger individual sum > better balance > lex
static bool better(const Sol& a, const Sol& b) {
    if (!b.valid) return a.valid;
    if (REQ_FULLCOL && a.ncol1 != b.ncol1) return a.ncol1 < b.ncol1;
    if (a.coverage != b.coverage) return a.coverage > b.coverage;
    if (a.nw12     != b.nw12)     return a.nw12     < b.nw12;
    if (a.sumIndiv != b.sumIndiv) return a.sumIndiv > b.sumIndiv;
    if (a.minIndiv != b.minIndiv) return a.minIndiv > b.minIndiv;
    return std::lexicographical_compare(a.cls, a.cls + NEW, b.cls, b.cls + NEW);
}
static std::mutex BESTMTX;
static Sol BEST;
static std::atomic<int64_t> BESTCOV{-1};
// hard validity gate: the seven classes must be pairwise orthogonal (mutual commutation)
// and linearly independent (rank 11 for the complete stabilizer).  Nothing enters the
// incumbent without passing this, whatever the search did.
static bool valid_set(const uint32_t* cls) {
    for (int a = 0; a < NEW; ++a) {
        if (cls[a] == 0) return false;
        for (int b = a + 1; b < NEW; ++b) if (csymp(cls[a], cls[b])) return false;
    }
    uint32_t piv[VDIM] = {0};
    for (int a = 0; a < NEW; ++a) {
        uint32_t v = cls[a];
        while (v) { int b = 31 - std::countl_zero(v);
            if (piv[b]) v ^= piv[b]; else { piv[b] = v; break; } }
        if (!v) return false;                 // linearly dependent on the previous ones
    }
    return true;
}
static bool offer(const uint32_t* cls) {
    if (!valid_set(cls)) return false;
    Sol S = finalize_sol(cls);
    if (S.coverage < BESTCOV.load(std::memory_order_relaxed)) return false;
    std::lock_guard<std::mutex> lk(BESTMTX);
    if (!better(S, BEST)) return false;
    BEST = S; BESTCOV.store(S.coverage); return true;
}

// ===========================================================================================
//  STEP 7 -- randomized greedy construction + steepest-descent local search
//  Greedy gain of adding class v to a d-dimensional span D:  G(v) = SUM_{u in D} cov1(v^u),
//  which is EXACTLY the increment of the span sum, so greedy is optimising the true objective
//  at every step (not a surrogate).
// ===========================================================================================
struct Builder {
    std::vector<uint32_t> allow, tmp;
    uint32_t span[1 << NEW]; int nspan = 1;
    uint32_t spanSorted[1 << NEW];
    void reset() { span[0] = 0; nspan = 1; }
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

static bool greedy_build(std::mt19937_64& rng, int rndTop, uint32_t* out, Builder& B) {
    B.reset(); std::copy(B.span, B.span + 1, B.spanSorted);
    B.allow.assign(CAND.begin(), CAND.end());
    for (int d = 0; d < NEW; ++d) {
        if (B.allow.empty()) return false;
        int64_t bestG = -1; std::vector<std::pair<int64_t,uint32_t>> top;
        top.reserve(size_t(rndTop) + 4);
        for (uint32_t v : B.allow) {
            int64_t g = B.gain(v);
            if ((int)top.size() < rndTop) { top.emplace_back(g, v); std::push_heap(top.begin(), top.end(),
                    [](const std::pair<int64_t,uint32_t>&a, const std::pair<int64_t,uint32_t>&b){ return a.first > b.first; }); }
            else if (g > top.front().first) {
                std::pop_heap(top.begin(), top.end(), [](const std::pair<int64_t,uint32_t>&a, const std::pair<int64_t,uint32_t>&b){ return a.first > b.first; });
                top.back() = {g, v};
                std::push_heap(top.begin(), top.end(), [](const std::pair<int64_t,uint32_t>&a, const std::pair<int64_t,uint32_t>&b){ return a.first > b.first; });
            }
            bestG = std::max(bestG, g);
        }
        uint32_t pick = top[rng() % top.size()].second;
        out[d] = pick;
        B.add_to_span(pick);
        B.tmp.clear();
        for (uint32_t v : B.allow) if (!csymp(v, pick) && !B.in_span(v)) B.tmp.push_back(v);
        B.allow.swap(B.tmp);
    }
    return true;
}

// Steepest-descent 1-opt: repeatedly replace a single generator by the best compatible
// candidate.  om[i] caches, for candidate i, which of the seven current generators it
// commutes with; only the column of the generator that actually changed is recomputed, so a
// sweep costs 7 linear passes over the candidate array plus 64 table lookups per survivor.
struct LocalWS { std::vector<uint8_t> om; };

static void om_column(std::vector<uint8_t>& om, int a, uint32_t c) {
    const uint8_t bit = uint8_t(1u << a), clr = uint8_t(~bit);
    for (size_t i = 0; i < CAND.size(); ++i)
        om[i] = uint8_t((om[i] & clr) | (csymp(CAND[i], c) ? 0 : bit));
}

static bool local_improve(uint32_t* cls, LocalWS& WS) {
    if (WS.om.size() != CAND.size()) WS.om.assign(CAND.size(), 0);
    for (int a = 0; a < NEW; ++a) om_column(WS.om, a, cls[a]);
    Sol cur = finalize_sol(cls);
    bool improved = false;
    for (;;) {
        bool any = false;
        for (int a = 0; a < NEW; ++a) {
            uint32_t sp[1 << (NEW - 1)]; int ns = 1; sp[0] = 0;
            for (int b = 0; b < NEW; ++b) if (b != a) {
                for (int i = 0; i < ns; ++i) sp[ns + i] = sp[i] ^ cls[b];
                ns <<= 1; }
            uint32_t spS[1 << (NEW - 1)]; std::copy(sp, sp + ns, spS); std::sort(spS, spS + ns);
            const uint8_t need = uint8_t(((1u << NEW) - 1u) & ~(1u << a));
            uint32_t trial[NEW]; std::copy(cls, cls + NEW, trial);
            uint32_t bestV = cls[a]; Sol bestS = cur;
            for (size_t i = 0; i < CAND.size(); ++i) {
                if ((WS.om[i] & need) != need) continue;
                uint32_t v = CAND[i];
                if (std::binary_search(spS, spS + ns, v)) continue;
                trial[a] = v;
                Sol T = finalize_sol(trial);
                if (better(T, bestS)) { bestS = T; bestV = v; }
            }
            if (bestV != cls[a]) {
                cls[a] = bestV; cur = bestS; any = true; improved = true;
                om_column(WS.om, a, bestV);          // only this column became stale
            }
        }
        if (!any) break;
    }
    return improved;
}

// ---- large-neighbourhood search: drop `nk` generators and greedily rebuild them ----------
static bool rebuild_from(const uint32_t* keep, int nk, std::mt19937_64& rng, int rndTop,
                         uint32_t* out, Builder& B) {
    B.reset(); B.spanSorted[0] = 0;
    for (int i = 0; i < nk; ++i) { out[i] = keep[i]; B.add_to_span(keep[i]); }
    B.allow.clear();
    for (size_t i = 0; i < CAND.size(); ++i) {
        uint32_t v = CAND[i];
        bool ok = true;
        for (int k = 0; k < nk && ok; ++k) if (csymp(v, keep[k])) ok = false;
        if (ok && !B.in_span(v)) B.allow.push_back(v);
    }
    for (int d = nk; d < NEW; ++d) {
        if (B.allow.empty()) return false;
        std::vector<std::pair<int64_t, uint32_t>> top; top.reserve(size_t(rndTop) + 4);
        auto cmp = [](const std::pair<int64_t,uint32_t>& a, const std::pair<int64_t,uint32_t>& b) { return a.first > b.first; };
        for (uint32_t v : B.allow) {
            int64_t g = B.gain(v);
            if ((int)top.size() < rndTop) { top.emplace_back(g, v); std::push_heap(top.begin(), top.end(), cmp); }
            else if (g > top.front().first) { std::pop_heap(top.begin(), top.end(), cmp);
                top.back() = {g, v}; std::push_heap(top.begin(), top.end(), cmp); }
        }
        uint32_t pick = top[rng() % top.size()].second;
        out[d] = pick; B.add_to_span(pick);
        B.tmp.clear();
        for (uint32_t v : B.allow) if (!csymp(v, pick) && !B.in_span(v)) B.tmp.push_back(v);
        B.allow.swap(B.tmp);
    }
    return true;
}

// ===========================================================================================
//  STEP 8 -- exact branch and bound over 7-element increasing index sequences.
//  Admissible bound.  With D_d already fixed (2^d elements, span sum S_d), adding the next
//  generator contributes exactly G(v) = SUM_{u in D_d} cov1(v^u) and the remaining
//  generators contribute cosets totalling 128 - 2^{d+1} further elements, each of cov1 at
//  most COVMAX = max_v cov1(v).  Hence
//        coverage  <=  ( S_d + max_{v allowed} G(v) + (128 - 2^{d+1}) * COVMAX ) / 64,
//  which relaxes isotropy, independence and weight-realisability and is therefore valid.
// ===========================================================================================
static std::atomic<bool> ABORT{false};
static std::atomic<uint64_t> NODES{0}, PRUNED{0};
static double TIME_LIMIT = 1e18;

struct Exact {
    std::vector<uint32_t> allow[NEW + 1];       // indices into CAND
    std::vector<int64_t>  gains;
    uint32_t chosen[NEW]{};
    uint32_t span[1 << NEW]{}; 
    uint64_t nodes = 0, pruned = 0;

    void dfs(int d, int nspan, int64_t S) {
        if (ABORT.load(std::memory_order_relaxed)) return;
        ++nodes;
        if ((nodes & 0xFFFF) == 0 && elapsed() > TIME_LIMIT) { ABORT.store(true); return; }
        if (d == NEW) { offer(chosen); return; }
        auto& A = allow[d];
        if (A.empty()) { ++pruned; return; }
        gains.resize(A.size());
        int64_t gbest = -1;
        for (size_t i = 0; i < A.size(); ++i) {
            uint32_t v = CAND[A[i]];
            int64_t g = 0; for (int k = 0; k < nspan; ++k) g += COV1[span[k] ^ v];
            gains[i] = g; if (g > gbest) gbest = g;
        }
        const int64_t rest = int64_t(128 - (1 << (d + 1))) * COVMAX;
        if ((S + gbest + rest) / 64 <= BESTCOV.load(std::memory_order_relaxed)) { ++pruned; return; }
        for (size_t i = 0; i < A.size(); ++i) {
            if (ABORT.load(std::memory_order_relaxed)) return;
            int64_t S2 = S + gains[i];
            if ((S2 + rest) / 64 <= BESTCOV.load(std::memory_order_relaxed)) { ++pruned; continue; }
            uint32_t v = CAND[A[i]];
            chosen[d] = v;
            int ns2 = nspan << 1;
            for (int k = 0; k < nspan; ++k) span[nspan + k] = span[k] ^ v;
            uint32_t sorted[1 << NEW]; std::copy(span, span + ns2, sorted); std::sort(sorted, sorted + ns2);
            auto& B = allow[d + 1]; B.clear();
            for (size_t j = i + 1; j < A.size(); ++j) {
                uint32_t u = CAND[A[j]];
                if (csymp(u, v)) continue;
                if (std::binary_search(sorted, sorted + ns2, u)) continue;
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
        Exact E; E.span[0] = 0;
        for (;;) {
            size_t i = next.fetch_add(1);
            if (i >= CAND.size() || ABORT.load()) break;
            uint32_t v = CAND[i];
            E.chosen[0] = v; E.span[0] = 0; E.span[1] = v;
            int64_t S = COV1[v];
            const int64_t rest = int64_t(128 - 4) * COVMAX;
            if ((S + rest + 2 * COVMAX) / 64 <= BESTCOV.load()) { ++E.pruned; continue; }
            auto& B = E.allow[1]; B.clear();
            for (size_t j = i + 1; j < CAND.size(); ++j) {
                uint32_t u = CAND[j];
                if (csymp(u, v) || u == v) continue;
                B.push_back((uint32_t)j);
            }
            E.dfs(1, 2, S);
        }
        NODES += E.nodes; PRUNED += E.pruned;
    });
    for (auto& x : th) x.join();
}

// ===========================================================================================
//  STEP 9 -- reporting, provable upper bounds, self-tests, main
// ===========================================================================================
static std::string RB;
static void ap(const char* fmt, ...) {
    char buf[4096]; va_list a; va_start(a, fmt); vsnprintf(buf, sizeof(buf), fmt, a); va_end(a);
    RB += buf;
}

static int64_t UB_COVMAX = 0, UB_TOP127 = 0, UB_BEST = 0;
static void compute_upper_bounds() {
    UB_COVMAX = (int64_t(127) * COVMAX) / 64;
    std::vector<uint16_t> v(COV1.begin() + 1, COV1.end());        // skip class 0
    std::nth_element(v.begin(), v.begin() + 127, v.end(), std::greater<uint16_t>());
    int64_t s = 0; for (int i = 0; i < 127; ++i) s += v[i];
    UB_TOP127 = s / 64;
    UB_BEST = std::min<int64_t>({ (int64_t)TARGETS, UB_COVMAX, UB_TOP127 });
}

static void report(const Sol& S, bool proven, const char* fname) {
    static const int W10[1] = {10}; static const int W1012[2] = {10, 12};
    Verify V = verify_bruteforce(S, WMODE == 10 ? W10 : W1012, WMODE == 10 ? 1 : 2);
    // independent bitset union, computed from the class masks
    std::vector<uint64_t> uni(WORDS, 0), tmpm(WORDS, 0);
    for (int a = 0; a < NEW; ++a) { class_mask(S.cls[a], tmpm.data());
        for (int k = 0; k < WORDS; ++k) uni[k] |= tmpm[k]; }
    int64_t unionpc = 0; for (int k = 0; k < WORDS; ++k) unionpc += std::popcount(uni[k]);

    RB.clear();
    ap("============================================================\n");
    ap("%s\n", proven ? "OPTIMAL ADDITIONAL STABILIZER SET" : "BEST ADDITIONAL STABILIZER SET FOUND");
    ap("============================================================\n\n");
    ap("Weight mode: %s\n", WMODE == 10 ? "all new generators weight 10"
                                        : "each new generator weight 10 or 12");
    ap("Status: %s\n\n", proven ? "GLOBAL OPTIMUM PROVEN" : "BEST KNOWN SOLUTION - OPTIMUM NOT PROVEN");
    ap("Target errors: %d\n\n", TARGETS);
    ap("Detected:              %lld\n", (long long)S.coverage);
    ap("Undetected remaining:  %lld\n", (long long)(TARGETS - S.coverage));
    ap("Coverage:              %lld / %d\n", (long long)S.coverage, TARGETS);
    ap("Coverage percentage:   %.6f %%\n\n", 100.0 * double(S.coverage) / double(TARGETS));
    if (S.coverage == TARGETS) ap("ALL %d TARGET ERRORS DETECTED\n\n", TARGETS);
    ap("Original generators:\n");
    for (int i = 0; i < NG0; ++i) ap("g%d = %s   (weight %d)\n", i + 1, pstr(G0[i]).c_str(), pwt(G0[i]));
    ap("\nAdditional generators:\n");
    for (int a = 0; a < NEW; ++a) ap("h%d = %s   (weight %d)\n", a + 1, pstr(S.pau[a]).c_str(), V.weights[a]);
    ap("\nWeights:\n");
    for (int a = 0; a < NEW; ++a) ap("h%d: %d\n", a + 1, V.weights[a]);
    { int c10 = 0, c12 = 0; for (int a = 0; a < NEW; ++a) (V.weights[a] == 10 ? c10 : c12)++;
      ap("\nWeight distribution:\nweight 10 : %d\nweight 12 : %d\n", c10, c12); }
    ap("\nRank of complete stabilizer: %d  (target %d)\n", V.rank, NG0 + NEW);
    ap("Qubits whose 11-generator column is not 2-dimensional: %d"
       "  (this number equals the count of undetected weight-1 targets)\n", deficient_columns(S.pau));
    ap("\nCommutation matrix (rows/cols g1..g4,h1..h7; 0 = commute):\n");
    for (int i = 0; i < NG0 + NEW; ++i) { ap("  ");
        for (int j = 0; j < NG0 + NEW; ++j) ap("%d ", V.comm[i][j]); ap("\n"); }
    ap("\nTarget-error detection by each new generator:\n");
    for (int a = 0; a < NEW; ++a)
        ap("h%d: %lld   (algebraic cov1 = %d %s)\n", a + 1, (long long)V.per_gen[a],
           (int)COV1[S.cls[a]], V.per_gen[a] == (int64_t)COV1[S.cls[a]] ? "OK" : "*** MISMATCH ***");
    ap("\nUnion coverage:\n  bitset popcount            : %lld\n", (long long)unionpc);
    ap("  brute-force recount        : %lld\n", (long long)V.detected);
    ap("  algebraic (1/64)*sum cov1  : %lld\n", (long long)S.coverage);
    ap("  %s\n", (unionpc == V.detected && V.detected == S.coverage) ? "all three agree" : "*** MISMATCH ***");
    ap("\nDetected target errors by weight:\n");
    for (int w = 1; w <= MAXW; ++w)
        ap("  weight %d : %lld / %lld\n", w, (long long)V.det_by_w[w], (long long)V.tot_by_w[w]);
    ap("\nComplete 11 x 28 parity-check matrix [X | Z] (rows g1..g4, h1..h7):\n");
    { uint32_t all[NG0 + NEW];
      for (int i = 0; i < NG0; ++i) all[i] = G0[i];
      for (int a = 0; a < NEW; ++a) all[NG0 + a] = S.pau[a];
      for (int i = 0; i < NG0 + NEW; ++i) {
          ap("  ");
          for (int j = 0; j < NQ; ++j) ap("%d", (px(all[i]) >> j) & 1);
          ap(" | ");
          for (int j = 0; j < NQ; ++j) ap("%d", (pz(all[i]) >> j) & 1);
          ap("\n"); } }
    if (!V.leftover.empty()) {
        ap("\nRemaining undetected target errors: %zu\n", V.leftover.size());
        for (int w = 1; w <= MAXW; ++w) {
            ap("\nWeight %d:\n", w);
            int n = 0;
            for (uint32_t E : V.leftover) if (pwt(E) == w) { ap("  E%-4d = %s\n", ++n, pstr(E).c_str()); }
            if (n == 0) ap("  (none)\n");
        }
    } else ap("\nRemaining undetected target errors: none\n");
    ap("\nFull 11-generator check over all %lld errors of weight 1..4:\n", (long long)V.all_total);
    ap("  total errors            : %lld\n", (long long)V.all_total);
    ap("  originally detected     : %lld\n", (long long)V.all_orig_det);
    ap("  originally undetected   : %lld\n", (long long)V.all_orig_und);
    ap("  newly detected          : %lld\n", (long long)V.all_new_det);
    ap("  still undetected        : %lld\n", (long long)V.all_still_und);
    ap("  detected by all 11 gens : %lld\n", (long long)V.all_det11);
    ap("\nINDEPENDENT VERIFICATION\n");
    ap("  new-generator weights allowed   : %s\n", V.wt_ok ? "OK" : "FAIL");
    ap("  [h_a , g_i] = 0                 : %s\n", V.comm_g_ok ? "OK" : "FAIL");
    ap("  [h_a , h_b] = 0                 : %s\n", V.comm_h_ok ? "OK" : "FAIL");
    ap("  rank = 11                       : %s\n", V.rank_ok ? "OK" : "FAIL");
    ap("  target set rebuilt from scratch : %lld errors\n", (long long)V.all_orig_und);
    ap("  recomputed coverage             : %lld %s\n", (long long)V.detected,
       V.detected == S.coverage ? "(matches optimiser)" : "*** MISMATCH ***");
    ap("\nProvable upper bounds on the objective (valid for ANY admissible 7-set):\n");
    ap("  127 * max_v cov1(v) / 64            = %lld\n", (long long)UB_COVMAX);
    ap("  (sum of 127 largest cov1) / 64      = %lld\n", (long long)UB_TOP127);
    ap("  best upper bound used               = %lld\n", (long long)UB_BEST);
    ap("\nCandidate count: %zu classes\n", CAND.size());
    ap("Search nodes: %llu\nPruned nodes: %llu\nElapsed time: %.3f s\n",
       (unsigned long long)NODES.load(), (unsigned long long)PRUNED.load(), elapsed());
    ap("============================================================\n");
    fputs(RB.c_str(), stdout); fflush(stdout);
    FILE* f = fopen(fname, "w"); if (f) { fputs(RB.c_str(), f); fclose(f); }
}


// ===========================================================================================
//  BEAM SEARCH over the subspace construction.  A state is a d-dimensional isotropic subspace
//  given by an increasing sequence of candidate indices; its score is the exact partial span
//  sum S_d = sum_{v in D_d\0} cov1(v), which is the same quantity the final objective is
//  built from, so ranking by S_d is ranking by the true objective restricted to what is
//  already decided.  States are de-duplicated by a hash of their sorted span, which removes
//  the (2^d-1)(2^d-2).../ d! different bases of the same subspace.
// ===========================================================================================
struct BeamState {
    uint32_t cls[NEW];
    int32_t  last;          // index into POOL of the last generator (increasing order)
    int64_t  S;
    uint64_t key;           // canonical hash of the span
};

static std::vector<uint32_t> POOL;          // candidate classes used by the beam

static uint64_t span_key(const uint32_t* cls, int d) {
    uint32_t sp[1 << NEW]; int ns = 1; sp[0] = 0;
    for (int b = 0; b < d; ++b) { for (int i = 0; i < ns; ++i) sp[ns + i] = sp[i] ^ cls[b]; ns <<= 1; }
    std::sort(sp, sp + ns);
    uint64_t h = 1469598103934665603ull;
    for (int i = 1; i < ns; ++i) { h ^= sp[i]; h *= 1099511628211ull; }
    return h;
}

static void beam_search(int K, int poolN, int nthreads) {
    POOL.assign(CAND.begin(), CAND.begin() + std::min<size_t>(CAND.size(), size_t(poolN)));
    std::vector<BeamState> cur, nxt;
    for (int i = 0; i < (int)POOL.size() && (int)cur.size() < K; ++i) {
        BeamState s{}; s.cls[0] = POOL[i]; s.last = i; s.S = COV1[POOL[i]];
        s.key = span_key(s.cls, 1); cur.push_back(s);
    }
    for (int d = 1; d < NEW; ++d) {
        std::vector<std::vector<BeamState>> perThread(nthreads);
        std::atomic<size_t> next{0};
        std::vector<std::thread> th;
        for (int t = 0; t < nthreads; ++t) th.emplace_back([&, t]() {
            auto& out = perThread[t];
            uint32_t sp[1 << NEW];
            for (;;) {
                size_t si = next.fetch_add(1);
                if (si >= cur.size()) break;
                const BeamState& s = cur[si];
                int ns = 1; sp[0] = 0;
                for (int b = 0; b < d; ++b) { for (int i = 0; i < ns; ++i) sp[ns + i] = sp[i] ^ s.cls[b]; ns <<= 1; }
                uint32_t spS[1 << NEW]; std::copy(sp, sp + ns, spS); std::sort(spS, spS + ns);
                // keep this state's best few extensions
                constexpr int KPER = 24;
                std::vector<std::pair<int64_t, int>> best; best.reserve(KPER + 1);
                auto cmp = [](const std::pair<int64_t,int>& a, const std::pair<int64_t,int>& b) { return a.first > b.first; };
                for (int j = s.last + 1; j < (int)POOL.size(); ++j) {
                    uint32_t v = POOL[j];
                    bool ok = true;
                    for (int b = 0; b < d && ok; ++b) if (csymp(v, s.cls[b])) ok = false;
                    if (!ok) continue;
                    if (std::binary_search(spS, spS + ns, v)) continue;
                    int64_t g = 0; for (int i = 0; i < ns; ++i) g += COV1[sp[i] ^ v];
                    if ((int)best.size() < KPER) { best.emplace_back(g, j); std::push_heap(best.begin(), best.end(), cmp); }
                    else if (g > best.front().first) { std::pop_heap(best.begin(), best.end(), cmp);
                        best.back() = {g, j}; std::push_heap(best.begin(), best.end(), cmp); }
                }
                for (auto& b : best) {
                    BeamState n = s; n.cls[d] = POOL[b.second]; n.last = b.second;
                    n.S = s.S + b.first; n.key = span_key(n.cls, d + 1);
                    out.push_back(n);
                }
            }
        });
        for (auto& x : th) x.join();
        nxt.clear();
        for (auto& v : perThread) nxt.insert(nxt.end(), v.begin(), v.end());
        if (nxt.empty()) return;
        std::sort(nxt.begin(), nxt.end(), [](const BeamState& a, const BeamState& b) {
            if (a.S != b.S) return a.S > b.S;
            return std::lexicographical_compare(a.cls, a.cls + NEW, b.cls, b.cls + NEW); });
        std::vector<BeamState> ded; ded.reserve(std::min<size_t>(nxt.size(), size_t(K)));
        std::unordered_set<uint64_t> seen; seen.reserve(nxt.size() * 2);
        for (const auto& s : nxt) {
            if (!seen.insert(s.key).second) continue;
            ded.push_back(s);
            if ((int)ded.size() >= K) break;
        }
        cur.swap(ded);
    }
    for (const auto& s : cur) offer(s.cls);
}
// ---------------------------------------------------------------- self-tests
static uint64_t splitmix(uint64_t x) {
    x += 0x9E3779B97F4A7C15ull; x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBull; return x ^ (x >> 31);
}
static bool self_test(int n, uint64_t seed) {
    std::mt19937_64 rng(seed); bool ok = true;
    // (a) cov1 table vs direct symplectic counting
    for (int t = 0; t < n; ++t) {
        uint32_t v = uint32_t(rng() % VSZ);
        int64_t c = 0; for (int j = 0; j < TARGETS; ++j) c += csymp(v, TCLS[j]);
        if (c != (int64_t)COV1[v]) { printf("SELFTEST FAIL cov1[%u]: %lld vs %d\n", v, (long long)c, (int)COV1[v]); ok = false; break; }
    }
    // (b) class map: lift_of(cls_of(v)) must differ from v by an element of S0
    for (int t = 0; t < n && ok; ++t) {
        uint32_t v = 0;
        for (int k = 0; k < CDIM; ++k) if (rng() & 1)
            v ^= (k < NG0 ? G0[k] : (k < NG0 + VDIM / 2 ? UB[k - NG0] : WB[k - NG0 - VDIM / 2]));
        uint32_t r = v ^ lift_of(cls_of(v));
        bool in = false; for (int m = 0; m < 16; ++m) if (S0EL[m] == r) in = true;
        if (!in) { printf("SELFTEST FAIL: class map is not a section of C(S0) -> V\n"); ok = false; }
    }
    // (c) the (1/64)*sum-cov1 identity vs the explicit 5430-bit union, on real 7-sets
    Builder B; std::vector<uint64_t> uni(WORDS), tm(WORDS);
    for (int t = 0; t < std::max(1, n / 20) && ok; ++t) {
        uint32_t cls[NEW];
        if (!greedy_build(rng, 64, cls, B)) continue;
        std::fill(uni.begin(), uni.end(), 0ull);
        for (int a = 0; a < NEW; ++a) { class_mask(cls[a], tm.data());
            for (int k = 0; k < WORDS; ++k) uni[k] |= tm[k]; }
        int64_t pc = 0; for (int k = 0; k < WORDS; ++k) pc += std::popcount(uni[k]);
        int64_t alg = coverage_of(cls, NEW);
        if (pc != alg) { printf("SELFTEST FAIL: union %lld vs identity %lld\n", (long long)pc, (long long)alg); ok = false; }
    }
    return ok;
}

// ---------------------------------------------------------------- heuristic driver
static bool DETERMINISTIC = false;
static void heuristic(int nthreads, uint64_t seed, int restarts, double tlimit) {
    std::atomic<int> next{0};
    std::vector<std::thread> th;
    for (int t = 0; t < nthreads; ++t) th.emplace_back([&]() {
        Builder B; LocalWS WS; uint32_t cls[NEW], keep[NEW];
        Sol local;                                    // this worker's own incumbent
        for (;;) {
            int r = next.fetch_add(1);
            if (restarts > 0 && r >= restarts) break;
            if (elapsed() > tlimit) break;
            std::mt19937_64 rng(splitmix(seed ^ (uint64_t(r) * 0x1000193ull)));
            bool ok;
            if (DETERMINISTIC || !local.valid || (r % 4) == 0) {
                ok = greedy_build(rng, (r == 0) ? 1 : 1 + int(r % 12), cls, B);
            } else {
                // ruin and recreate: keep a random subset of the incumbent, rebuild the rest
                int nk = NEW - (1 + int(rng() % 3));
                int idx[NEW]; std::iota(idx, idx + NEW, 0);
                std::shuffle(idx, idx + NEW, rng);
                for (int i = 0; i < nk; ++i) keep[i] = local.cls[idx[i]];
                ok = rebuild_from(keep, nk, rng, 1 + int(rng() % 8), cls, B);
            }
            if (!ok) continue;
            local_improve(cls, WS);
            Sol S = finalize_sol(cls);
            if (valid_set(S.cls) && (!local.valid || better(S, local))) local = S;
            offer(cls);
            if (!DETERMINISTIC && (r % 32) == 31) {
                std::lock_guard<std::mutex> lk(BESTMTX);
                if (BEST.valid && (!local.valid || better(BEST, local))) local = BEST;
            }
        }
    });
    for (auto& x : th) x.join();
}

static void usage() {
    printf("stab14_ext -- seven extra weight-10 (or 10/12) stabilizer generators\n"
           "  --mode heuristic|exact|hybrid     (default hybrid)\n"
           "  --weight-mode 10|10or12|auto      (default auto: run 10 first, then 10or12)\n"
           "  --threads N        --seed S       --time-limit SECONDS\n"
           "  --heur-time S      seconds of randomized greedy + local search per problem\n"
           "  --restarts N       fixed number of randomized greedy restarts (0 = time based)\n"
           "  --cand-top N       restrict the candidate classes to the N with largest cov1\n"
           "  --deterministic    fixed restart count, reproducible result\n"
           "  --selftest N       (default 200)      --out PREFIX\n");
}

int main(int argc, char** argv) {
    std::string mode = "hybrid", wmodeArg = "auto", outPrefix = "extension";
    int nthreads = (int)std::max(1u, std::thread::hardware_concurrency());
    int restarts = 0, selftestN = 200, candTop = 0, beamK = 4000, beamPool = 40000;
    uint64_t seed = 12345ull; double heurTime = 20.0, exactTime = 60.0, globalTL = 1e18;
    bool deterministic = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto nxt = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : std::string(); };
        if      (a == "--mode")          mode = nxt();
        else if (a == "--weight-mode")   wmodeArg = nxt();
        else if (a == "--threads")       nthreads = std::max(1, atoi(nxt().c_str()));
        else if (a == "--seed")          seed = strtoull(nxt().c_str(), nullptr, 10);
        else if (a == "--time-limit")    globalTL = atof(nxt().c_str());
        else if (a == "--exact-time")    exactTime = atof(nxt().c_str());
        else if (a == "--heur-time")     heurTime = atof(nxt().c_str());
        else if (a == "--restarts")      restarts = atoi(nxt().c_str());
        else if (a == "--cand-top")      candTop = atoi(nxt().c_str());
        else if (a == "--beam")          beamK = atoi(nxt().c_str());
        else if (a == "--beam-pool")     beamPool = atoi(nxt().c_str());
        else if (a == "--deterministic") deterministic = true;
        else if (a == "--full-columns")  REQ_FULLCOL = true;
        else if (a == "--selftest")      selftestN = atoi(nxt().c_str());
        else if (a == "--out")           outPrefix = nxt();
        else { usage(); return (a == "--help" || a == "-h") ? 0 : 1; }
    }
    TIME_LIMIT = globalTL;
    if (deterministic && restarts == 0) restarts = 3000;
    DETERMINISTIC = deterministic;

    T0 = std::chrono::steady_clock::now();
    init_P();
    for (int i = 0; i < NG0; ++i) G0[i] = parse_pauli(G0_STR[i]);

    printf("=== extension of a fixed 14-qubit / 4-generator stabilizer ===\n");
    for (int i = 0; i < NG0; ++i) printf("g%d = %s  weight %d\n", i + 1, pstr(G0[i]).c_str(), pwt(G0[i]));
    for (int i = 0; i < NG0; ++i) if (pwt(G0[i]) != 8) { printf("FATAL: g%d does not have weight 8\n", i + 1); return 2; }
    for (int i = 0; i < NG0; ++i) for (int j = i + 1; j < NG0; ++j)
        if (symp(G0[i], G0[j])) { printf("FATAL: g%d and g%d do not commute\n", i + 1, j + 1); return 2; }

    build_targets();
    printf("\nreconstructed target set U = weight 1..4 errors with zero S0-syndrome\n");
    printf("  |U| = %zu   (weights 1/2/3/4 = %lld/%lld/%lld/%lld)\n", TGT.size(),
           (long long)TGT_BY_W[1], (long long)TGT_BY_W[2], (long long)TGT_BY_W[3], (long long)TGT_BY_W[4]);
    if ((int)TGT.size() != TARGETS_EXPECT) { printf("FATAL: |U| != %d -- refusing to optimise the wrong problem\n", TARGETS_EXPECT); return 3; }
    if (TGT_BY_W[1] != 12 || TGT_BY_W[2] != 69 || TGT_BY_W[3] != 528 || TGT_BY_W[4] != 4821) {
        printf("FATAL: target weight distribution does not match 12/69/528/4821\n"); return 3; }

    build_symplectic_basis();
    build_cov1();
    build_maskbase();
    printf("\ncentralizer C(S0): dim %d,  quotient V = C(S0)/S0: dim %d (%d classes)\n", CDIM, VDIM, VSZ);
    printf("max_v cov1(v) = %d   (the mean over all classes is %.1f)\n", COVMAX, TGT.size() / 2.0);
    build_reprs();
    printf("weight-10 Paulis in C(S0): %lld  reaching %lld of the %d non-zero classes\n",
           (long long)N_W10, (long long)NCLS10, VSZ - 1);
    printf("weight-12 Paulis in C(S0): %lld  reaching %lld of the %d non-zero classes\n",
           (long long)N_W12, (long long)NCLS12, VSZ - 1);
    compute_upper_bounds();
    printf("\nPROVABLE UPPER BOUNDS on the coverage of ANY 7 commuting independent generators:\n");
    printf("  coverage = (1/64) * sum of cov1 over the 127 non-zero elements of their span\n");
    printf("  <= 127 * %d / 64                = %lld\n", COVMAX, (long long)UB_COVMAX);
    printf("  <= (sum of 127 largest cov1)/64 = %lld\n", (long long)UB_TOP127);
    printf("  ==> M <= %lld   (the target set has %d errors)\n", (long long)UB_BEST, TARGETS);
    if (UB_BEST < TARGETS)
        printf("  ==> FULL COVERAGE OF ALL %d TARGET ERRORS IS PROVABLY IMPOSSIBLE\n", TARGETS);

    if (selftestN > 0) {
        printf("\n[selftest] %d checks (cov1 table, class map, identity vs 5430-bit union) ... ", selftestN);
        fflush(stdout);
        WMODE = 10; build_candidates(candTop);
        if (!self_test(selftestN, seed)) { printf("FAILED\n"); return 4; }
        printf("all pass\n");
    }

    std::vector<int> modes;
    if (wmodeArg == "10") modes.push_back(10);
    else if (wmodeArg == "10or12") modes.push_back(1012);
    else { modes.push_back(10); modes.push_back(1012); }

    int64_t M10 = -1, M1012 = -1;
    for (size_t mi = 0; mi < modes.size(); ++mi) {
        WMODE = modes[mi];
        if (WMODE == 1012 && M10 == TARGETS) {
            printf("\nthe weight-10 search already covers everything; relaxation not needed\n"); break; }
        BEST = Sol(); BESTCOV.store(-1); NODES.store(0); PRUNED.store(0); ABORT.store(false);
        build_candidates(candTop);
        printf("\n---------------------------------------------------------------\n");
        printf("PROBLEM %s : new generator weights %s\n", WMODE == 10 ? "A" : "B",
               WMODE == 10 ? "all = 10" : "each in {10,12}");
        printf("candidate classes: %zu   (max cov1 among them = %d)\n", CAND.size(), COVMAX_CAND);
        fflush(stdout);
        if (beamK > 0) {
            beam_search(beamK, beamPool, nthreads);
            printf("[beam] K=%d pool=%d -> coverage %lld / %d  (%.2f s)\n", beamK, beamPool,
                   (long long)BESTCOV.load(), TARGETS, elapsed());
            fflush(stdout);
        }
        double hEnd = elapsed() + heurTime;
        if (mode == "heuristic" || mode == "hybrid") {
            heuristic(nthreads, seed, restarts, std::min(hEnd, globalTL));
            printf("[heuristic] best coverage = %lld / %d after %.2f s\n",
                   (long long)BESTCOV.load(), TARGETS, elapsed());
            fflush(stdout);
        }
        bool proven = (BESTCOV.load() >= UB_BEST);
        if ((mode == "exact" || mode == "hybrid") && !proven) {
            TIME_LIMIT = std::min(globalTL, elapsed() + exactTime);
            printf("[exact] branch and bound over %zu candidate classes (budget %.0f s) ...\n",
                   CAND.size(), exactTime);
            // Honest statement about what the bound can and cannot do.  The coset relaxation
            // prunes the root only when 127*COVMAX/64 <= incumbent, i.e. COVMAX <= 64*inc/127.
            printf("[exact] note: this relaxation can prune only while max_v cov1(v) <= %.1f;"
                   " the actual value is %d, so it cannot prune near the root and an exhaustive\n"
                   "        proof of optimality is out of reach with it -- expect the time limit.\n",
                   64.0 * double(BESTCOV.load()) / 127.0, COVMAX);
            fflush(stdout);
            exact_search(nthreads);
            if (!ABORT.load()) proven = true;
            printf("[exact] %s  nodes %llu  pruned %llu  (%.2f s)\n",
                   ABORT.load() ? "TIME LIMIT / INCOMPLETE" : "SEARCH TREE EXHAUSTED",
                   (unsigned long long)NODES.load(), (unsigned long long)PRUNED.load(), elapsed());
        }
        if (!BEST.valid) { printf("no feasible 7-generator set found\n"); continue; }
        if (BESTCOV.load() >= UB_BEST) proven = true;      // meets a provable upper bound
        std::string fn = outPrefix + (WMODE == 10 ? "_w10.txt" : "_w10or12.txt");
        report(BEST, proven, fn.c_str());
        if (WMODE == 10) M10 = BEST.coverage; else M1012 = BEST.coverage;
    }
    printf("\n=== SUMMARY ===\n");
    if (M10   >= 0) printf("M_10    = %lld / %d\n", (long long)M10, TARGETS);
    if (M1012 >= 0) printf("M_10/12 = %lld / %d\n", (long long)M1012, TARGETS);
    printf("provable upper bound for both = %lld\n", (long long)UB_BEST);
    return 0;
}
