// ===========================================================================================
//  stab14.cpp -- Exact optimizer for a 4-generator, 14-qubit stabilizer check matrix H=[X|Z]
//
//  OBJECTIVE (exact, never a proxy):
//      N_det = #{ Pauli errors e : 1 <= wt(e) <= 4 , syndrome s(e) != 0 },  wt(Y)=1.
//      Total error universe = sum_{w=1..4} C(14,w) 3^w = 91770.
//
//  CONSTRAINTS: every generator has Pauli weight exactly GENERATOR_WEIGHT (default 8),
//      the 4 generators pairwise commute, and rank_F2(H) = 4.
//
//  -----------------------------------------------------------------------------------------
//  MATHEMATICAL CORE  (see the long derivation in README_MATH.txt)
//  -----------------------------------------------------------------------------------------
//  Column view.  Qubit j contributes a column (x_j, z_j) in F_2^4 x F_2^4, where x_j[i] / z_j[i]
//  are the X/Z bits of generator i on qubit j.
//
//  (1) Per-qubit local Clifford acts as GL(2,2) on (x_j,z_j):  (x,z) -> (ax+bz, cx+dz), ad+bc=1.
//      It preserves (a) the per-row Pauli weight, (b) the commutation form, (c) the whole
//      objective.  The complete GL(2,2)-invariant of a column is the SUBSPACE
//                 W_j = span{x_j, z_j}  <= F_2^4 ,     dim W_j in {0,1,2}.
//      Hence the search space is the multiset { W_1..W_14 } of 51 subspaces (1+15+35), not
//      256 column types and certainly not 2^112 matrices.
//
//  (2) Group element g_chi = prod_i g_i^{chi_i}, chi in F_2^4.  g_chi acts non-trivially on
//      qubit j  <=>  chi is NOT orthogonal to W_j.  Therefore
//                 m(chi) := wt(g_chi) = 14 - c(chi),      c(chi) := #{ j : W_j <= chi^perp }.
//      Writing U_j := W_j^perp,  c(chi) = #{ j : chi in U_j }: a point-degree in PG(3,2).
//
//  (3) Walsh/Hadamard over the 16-element syndrome space gives the EXACT objective as a
//      function of the stabilizer-group WEIGHT ENUMERATOR only:
//            #{ e : wt(e)=w, s(e)=0 } = (1/16) * sum_{chi} e_w( F_1(chi),...,F_14(chi) )
//      where F_j(chi) = 3 if g_chi acts trivially on qubit j, else -1.  So for chi != 0 the
//      multiset {F_j} is (3 with multiplicity 14-m, -1 with multiplicity m) and
//            e_w = sum_k C(m,k)(-1)^k C(14-m, w-k) 3^{w-k}   =: Ew[m][w].
//      Scoring a candidate is 15 table lookups instead of 91770 symplectic products.
//
//  (4) Constraints in column language:
//        weight:        c(e_i) = 14 - GENERATOR_WEIGHT   for the four basis points e_i,
//        commutation:   XOR_j plucker(W_j) = 0 in F_2^6   (plucker(W)=x z^T + z x^T is a
//                       basis-independent invariant of the 2-dim subspace W; 0 if dim W<2),
//        rank 4:        c(chi) <= 13 for all chi != 0.
//
//  Residual symmetry after all of the above: S_4 permuting the four generators (qubit
//  permutations are already quotiented out by using a MULTISET, local Cliffords by using
//  subspaces).  S_4 is used as a lex-min prefix pruning rule inside the DFS.
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

// ------------------------------------------------------------------ compile-time parameters
constexpr int NQ            = 14;      // physical qubits
constexpr int NGEN          = 4;       // stabilizer generators (rank of H)
constexpr int MAXW          = 4;       // maximum error weight considered
constexpr int GENERATOR_WEIGHT_DEFAULT = 8;   // <-- change to 4 for the alternative spec
constexpr int NCHI          = 16;      // |F_2^4|
constexpr int NTYPE_MAX     = 51;      // subspaces of F_2^4 of dim <= 2

static int GENERATOR_WEIGHT = GENERATOR_WEIGHT_DEFAULT;
static int TARGET_C         = NQ - GENERATOR_WEIGHT_DEFAULT;   // required c(e_i)

// ------------------------------------------------------------------ small combinatorics
static int64_t Cbin[64][64];
static void init_binom() {
    for (int n = 0; n < 64; ++n) { Cbin[n][0] = 1; for (int k = 1; k <= n; ++k)
        Cbin[n][k] = Cbin[n-1][k-1] + (k <= n-1 ? Cbin[n-1][k] : 0); }
}
static inline int64_t Cb(int n, int k) { return (k < 0 || k > n || n < 0) ? 0 : Cbin[n][k]; }

// Ew[m][w] = e_w of the multiset (3 x (14-m), -1 x m) -- see (3) above.
static int64_t Ew[NQ + 1][MAXW + 1];
static int64_t Ptab[NQ + 1];          // P(m) = sum_{w=1..4} Ew[m][w]
static int64_t PbyC[NQ + 1];          // PbyC[c] = P(14-c)   (indexed by the degree c)
static int64_t TOTAL_ERRORS = 0;      // 91770
static int64_t TOTW[MAXW + 1];        // C(14,w) 3^w

static void init_tables() {
    init_binom();
    int64_t p3[MAXW + 1]; p3[0] = 1; for (int w = 1; w <= MAXW; ++w) p3[w] = p3[w-1] * 3;
    for (int w = 0; w <= MAXW; ++w) TOTW[w] = Cb(NQ, w) * p3[w];
    TOTAL_ERRORS = 0; for (int w = 1; w <= MAXW; ++w) TOTAL_ERRORS += TOTW[w];
    for (int m = 0; m <= NQ; ++m) {
        for (int w = 0; w <= MAXW; ++w) {
            int64_t s = 0;
            for (int k = 0; k <= w; ++k) {
                int64_t t = Cb(m, k) * Cb(NQ - m, w - k) * p3[w - k];
                s += (k & 1) ? -t : t;
            }
            Ew[m][w] = s;
        }
        Ptab[m] = 0; for (int w = 1; w <= MAXW; ++w) Ptab[m] += Ew[m][w];
    }
    for (int c = 0; c <= NQ; ++c) PbyC[c] = Ptab[NQ - c];
}

// ===========================================================================================
//  SUBSPACE (COLUMN-TYPE) TABLES
//  A column type is a subspace W <= F_2^4 with dim W in {0,1,2}: 1 + 15 + 35 = 51 types.
//  Everything the objective and the constraints need is a function of W alone.
// ===========================================================================================
struct TypeInfo {
    uint16_t wmask;      // bitmask over F_2^4 of the elements of W
    uint16_t umask;      // bitmask over F_2^4 of the elements of U = W^perp
    uint8_t  dimW;       // 0,1,2
    uint8_t  sig;        // bit i set  <=>  e_i in U  <=>  generator i is IDENTITY on this qubit
    uint8_t  pl;         // 6-bit Plucker invariant  x z^T + z x^T   (0 when dim W < 2)
    uint8_t  degU;       // |U| - 1 = number of chi != 0 with chi in U  (15, 7 or 3)
    uint8_t  bx, bz;     // a basis (x,z) of W, used only when reconstructing H
    uint16_t chi[15];    // the nonzero elements of U (degU of them)
};
static TypeInfo TY[NTYPE_MAX];
static int NTYPE = 0;

static inline int par(unsigned a) { return std::popcount(a) & 1; }

// pair index for the 6 off-diagonal entries (i<k) of the 4x4 symplectic form
static const int PI0[6] = {0,0,0,1,1,2};
static const int PI1[6] = {1,2,3,2,3,3};

static uint8_t plucker_of(int x, int z) {
    uint8_t p = 0;
    for (int b = 0; b < 6; ++b) {
        int i = PI0[b], k = PI1[b];
        int v = (((x >> i) & 1) & ((z >> k) & 1)) ^ (((z >> i) & 1) & ((x >> k) & 1));
        if (v) p |= uint8_t(1u << b);
    }
    return p;
}

static void build_types() {
    std::vector<uint16_t> masks;
    masks.push_back(1u);                                   // dim 0 : {0}
    for (int v = 1; v < 16; ++v) masks.push_back(uint16_t(1u | (1u << v)));   // dim 1
    for (int v = 1; v < 16; ++v) for (int w = v + 1; w < 16; ++w) {           // dim 2
        uint16_t m = uint16_t(1u | (1u << v) | (1u << w) | (1u << (v ^ w)));
        if (std::find(masks.begin(), masks.end(), m) == masks.end()) masks.push_back(m);
    }
    std::vector<TypeInfo> tmp;
    for (uint16_t wm : masks) {
        TypeInfo t{}; t.wmask = wm;
        t.dimW = uint8_t(std::countr_zero(unsigned(std::popcount(wm))));   // |W| = 2^dim
        // U = W^perp
        uint16_t um = 0;
        for (int chi = 0; chi < 16; ++chi) {
            bool ok = true;
            for (int w = 0; w < 16; ++w) if ((wm >> w) & 1) if (par(unsigned(chi & w))) { ok = false; break; }
            if (ok) um = uint16_t(um | (1u << chi));
        }
        t.umask = um;
        t.degU = uint8_t(std::popcount(um) - 1);
        t.sig = 0;
        for (int i = 0; i < 4; ++i) if ((um >> (1 << i)) & 1) t.sig |= uint8_t(1u << i);
        int n = 0; for (int chi = 1; chi < 16; ++chi) if ((um >> chi) & 1) t.chi[n++] = uint16_t(chi);
        // basis of W and its Plucker invariant (independent of the chosen basis)
        int b[2] = {0,0}, nb = 0;
        for (int v = 1; v < 16 && nb < 2; ++v) if ((wm >> v) & 1) b[nb++] = v;
        t.bx = uint8_t(b[0]); t.bz = uint8_t(b[1]);
        t.pl = (t.dimW == 2) ? plucker_of(b[0], b[1]) : uint8_t(0);
        tmp.push_back(t);
    }
    // Ordering: dim W ascending (so that once the DFS enters the 35 "line" types every
    // remaining column adds exactly 3 to sum_chi c(chi) -- an EXACT budget for the bound),
    // then |sig| descending (types that can still satisfy the four weight equations first).
    std::sort(tmp.begin(), tmp.end(), [](const TypeInfo& a, const TypeInfo& b) {
        if (a.dimW != b.dimW) return a.dimW < b.dimW;
        int pa = std::popcount(a.sig), pb = std::popcount(b.sig);
        if (pa != pb) return pa > pb;
        return a.wmask < b.wmask;
    });
    NTYPE = (int)tmp.size();
    for (int i = 0; i < NTYPE; ++i) TY[i] = tmp[i];
}

// ---- S_4 : permutations of the four generators.  Acts on F_2^4 by permuting bits, hence on
// ---- subspaces, hence on column types.  This is the ONLY symmetry left after the multiset
// ---- (qubit permutations) and subspace (local Clifford) reductions.
static int S4act[24][NTYPE_MAX];
static int NS4 = 0;
static void build_s4() {
    int p[4] = {0,1,2,3}; NS4 = 0;
    std::sort(p, p + 4);
    do {
        int map16[16];
        for (int v = 0; v < 16; ++v) { int r = 0;
            for (int i = 0; i < 4; ++i) if ((v >> i) & 1) r |= 1 << p[i];
            map16[v] = r; }
        for (int t = 0; t < NTYPE; ++t) {
            uint16_t nm = 0;
            for (int v = 0; v < 16; ++v) if ((TY[t].wmask >> v) & 1) nm = uint16_t(nm | (1u << map16[v]));
            int found = -1;
            for (int u = 0; u < NTYPE; ++u) if (TY[u].wmask == nm) { found = u; break; }
            S4act[NS4][t] = found;
        }
        ++NS4;
    } while (std::next_permutation(p, p + 4));
}

// ===========================================================================================
//  SCORING FROM THE DEGREE VECTOR c(.)   -- 15 table lookups, no error enumeration
// ===========================================================================================
static const int EPT[4] = {1, 2, 4, 8};   // the four generator points e_1..e_4 inside F_2^4

struct Score {
    int64_t sumP;                 // sum_{chi != 0} P(m(chi))   (the quantity we MINIMISE)
    int64_t undetected;           // non-identity errors with zero syndrome
    int64_t detected;             // = TOTAL_ERRORS - undetected   (the quantity we MAXIMISE)
    int64_t det_w[MAXW + 1];      // detected per weight
    int64_t und_w[MAXW + 1];      // undetected per weight
};

static inline int64_t sumP_of(const uint8_t* c) {          // c[1..15]
    int64_t s = 0; for (int chi = 1; chi < 16; ++chi) s += PbyC[c[chi]];
    return s;
}

static Score score_of(const uint8_t* c) {
    Score S{}; int cnt[NQ + 1] = {0};
    for (int chi = 1; chi < 16; ++chi) cnt[NQ - c[chi]]++;      // histogram of m(chi)
    S.sumP = sumP_of(c);
    S.undetected = 0; S.detected = 0;
    for (int w = 1; w <= MAXW; ++w) {
        int64_t acc = TOTW[w];
        for (int m = 0; m <= NQ; ++m) if (cnt[m]) acc += int64_t(cnt[m]) * Ew[m][w];
        int64_t n0 = acc / 16;                                   // exact: acc is divisible by 16
        S.und_w[w] = n0; S.det_w[w] = TOTW[w] - n0;
        S.undetected += n0; S.detected += TOTW[w] - n0;
    }
    return S;
}

// ===========================================================================================
//  BRANCH-AND-BOUND MACHINERY (admissible lower bounds on the final sumP)
// ===========================================================================================
constexpr int64_t RSCALE = 1 << 20;
constexpr int64_t INFP   = (int64_t)1 << 40;
static int64_t Pmin[NQ + 2][NQ + 1];        // min_{0<=d<=r, c+d<=13} P(14-c-d)
static int64_t Rmax[NQ + 2][NQ + 1];        // max gain per unit of budget, scaled, rounded up
static uint8_t sigAvail[NTYPE_MAX + 1];     // OR of sig over all types with index >= t
static uint8_t sigPopMax[NTYPE_MAX + 1];    // max popcount(sig) over types with index >= t
static uint8_t degMax[NTYPE_MAX + 1];       // max degU over types with index >= t

static void init_bound_tables() {
    for (int c = 0; c <= NQ; ++c) for (int r = 0; r <= NQ; ++r) {
        int64_t best = (c <= NQ - 1) ? PbyC[c] : INFP;          // c == 14 is rank-deficient
        int64_t rat  = 0;
        int dmax = std::min(r, (NQ - 1) - c);
        for (int d = 1; d <= dmax; ++d) {
            best = std::min(best, PbyC[c + d]);
            int64_t gain = PbyC[c] - PbyC[c + d];               // >0 means improvement
            if (gain > 0) { int64_t rr = (gain * RSCALE + d - 1) / d; rat = std::max(rat, rr); }
        }
        Pmin[c][r] = best; Rmax[c][r] = rat;
    }
    sigAvail[NTYPE] = 0; sigPopMax[NTYPE] = 0; degMax[NTYPE] = 0;
    for (int t = NTYPE - 1; t >= 0; --t) {
        sigAvail[t]  = uint8_t(sigAvail[t + 1] | TY[t].sig);
        sigPopMax[t] = uint8_t(std::max<int>(sigPopMax[t + 1], std::popcount(TY[t].sig)));
        degMax[t]    = uint8_t(std::max<int>(degMax[t + 1], TY[t].degU));
    }
}

// Admissible lower bound on the final sumP for a node with degree vector c, r columns still
// to place, all of them of type index >= t0.  Returns INFP when the node is infeasible.
static inline int64_t node_bound(const uint8_t* c, int r, int t0) {
    // --- hard feasibility on the four generator-weight equations c(e_i) = TARGET_C
    int need = 0;
    for (int i = 0; i < 4; ++i) {
        int d = TARGET_C - c[EPT[i]];
        if (d < 0) return INFP;                       // overshoot: generator too light
        if (d > r) return INFP;                       // cannot be repaired: at most +1 per column
        if (d > 0 && !((sigAvail[t0] >> i) & 1)) return INFP;   // no remaining type touches e_i
        need += d;
    }
    // (r == 0 falls through: the loop below returns the exact sumP of the leaf)
    if (r > 0 && need > r * (int)sigPopMax[t0]) return INFP;
    // --- objective bounds
    const int64_t base4 = 4 * PbyC[TARGET_C];
    int64_t s1 = 0, s2 = 0, rat = 0;
    for (int chi = 1; chi < 16; ++chi) {
        if (chi == 1 || chi == 2 || chi == 4 || chi == 8) continue;
        int cc = c[chi];
        if (cc > NQ - 1) return INFP;                 // rank would drop below 4
        int64_t pm = Pmin[cc][r];
        if (pm >= INFP) return INFP;
        s1 += pm; s2 += PbyC[cc];
        rat = std::max(rat, Rmax[cc][r]);
    }
    int64_t Bmax = int64_t(degMax[t0]) * r - need;     // budget usable by the 11 other points
    if (Bmax < 0) return INFP;
    int64_t lb2 = s2 - (Bmax * rat + RSCALE - 1) / RSCALE;
    return base4 + std::max(s1, lb2);
}

// ===========================================================================================
//  RECONSTRUCTION  (column multiset  ->  concrete 4 x 28 check matrix)  and
//  INDEPENDENT VERIFICATION  (explicit enumeration of all 91770 errors)
// ===========================================================================================
struct Matrix { uint16_t Hx[NGEN]; uint16_t Hz[NGEN]; };

static Matrix build_matrix(const uint8_t* types) {
    Matrix M{}; for (int i = 0; i < NGEN; ++i) { M.Hx[i] = 0; M.Hz[i] = 0; }
    for (int j = 0; j < NQ; ++j) {
        const TypeInfo& t = TY[types[j]];
        for (int i = 0; i < NGEN; ++i) {
            if ((t.bx >> i) & 1) M.Hx[i] = uint16_t(M.Hx[i] | (1u << j));
            if ((t.bz >> i) & 1) M.Hz[i] = uint16_t(M.Hz[i] | (1u << j));
        }
    }
    return M;
}

struct Verdict {
    bool ok_weight, ok_commute, ok_rank, ok_nonidentity;
    int  weights[NGEN];
    int  rank;
    int64_t synd[16];       // syndrome histogram over the 91770 errors
    int64_t synd_w[MAXW + 1][16];
    int64_t detected, undetected;
    int64_t det_w[MAXW + 1];
};

// --- completely independent scorer: brute force over every weight<=4 Pauli error ----------
static Verdict verify_bruteforce(const Matrix& M) {
    Verdict V{}; std::memset(V.synd, 0, sizeof(V.synd)); std::memset(V.synd_w, 0, sizeof(V.synd_w));
    for (int i = 0; i < NGEN; ++i) {
        V.weights[i] = std::popcount(unsigned(M.Hx[i] | M.Hz[i]));
    }
    V.ok_weight = true;
    for (int i = 0; i < NGEN; ++i) if (V.weights[i] != GENERATOR_WEIGHT) V.ok_weight = false;
    V.ok_nonidentity = true;
    for (int i = 0; i < NGEN; ++i) if ((M.Hx[i] | M.Hz[i]) == 0) V.ok_nonidentity = false;
    V.ok_commute = true;
    for (int i = 0; i < NGEN; ++i) for (int k = i + 1; k < NGEN; ++k) {
        int s = par(unsigned(M.Hx[i] & M.Hz[k])) ^ par(unsigned(M.Hz[i] & M.Hx[k]));
        if (s) V.ok_commute = false;
    }
    // rank over F_2 of the 4 x 28 matrix
    uint32_t rows[NGEN]; for (int i = 0; i < NGEN; ++i) rows[i] = uint32_t(M.Hx[i]) | (uint32_t(M.Hz[i]) << 14);
    int rank = 0;
    for (int b = 0; b < 28 && rank < NGEN; ++b) {
        int piv = -1; for (int i = rank; i < NGEN; ++i) if ((rows[i] >> b) & 1) { piv = i; break; }
        if (piv < 0) continue;
        std::swap(rows[rank], rows[piv]);
        for (int i = 0; i < NGEN; ++i) if (i != rank && ((rows[i] >> b) & 1)) rows[i] ^= rows[rank];
        ++rank;
    }
    V.rank = rank; V.ok_rank = (rank == NGEN);
    // exhaustive enumeration: choose the support, then one of {X,Y,Z} per support qubit
    static const uint16_t PX[3] = {1, 1, 0};    // x-bit of X, Y, Z
    static const uint16_t PZ[3] = {0, 1, 1};    // z-bit of X, Y, Z
    int idx[MAXW];
    for (int w = 1; w <= MAXW; ++w) {
        for (int i = 0; i < w; ++i) idx[i] = i;
        while (true) {
            int npat = 1; for (int i = 0; i < w; ++i) npat *= 3;
            for (int p = 0; p < npat; ++p) {
                uint16_t ex = 0, ez = 0; int q = p;
                for (int i = 0; i < w; ++i) { int a = q % 3; q /= 3;
                    ex = uint16_t(ex | (uint16_t(PX[a]) << idx[i]));
                    ez = uint16_t(ez | (uint16_t(PZ[a]) << idx[i])); }
                int s = 0;
                for (int g = 0; g < NGEN; ++g)
                    s |= (par(unsigned(M.Hx[g] & ez)) ^ par(unsigned(M.Hz[g] & ex))) << g;
                V.synd[s]++; V.synd_w[w][s]++;
            }
            int i = w - 1; while (i >= 0 && idx[i] == NQ - w + i) --i;
            if (i < 0) break;
            ++idx[i]; for (int k = i + 1; k < w; ++k) idx[k] = idx[k-1] + 1;
        }
    }
    V.undetected = V.synd[0];
    V.detected = 0; for (int s = 1; s < 16; ++s) V.detected += V.synd[s];
    for (int w = 1; w <= MAXW; ++w) { int64_t d = 0; for (int s = 1; s < 16; ++s) d += V.synd_w[w][s];
        V.det_w[w] = d; }
    return V;
}

static const char PAULI[2][2] = { {'I','Z'}, {'X','Y'} };   // [xbit][zbit]

static std::string row_string(const Matrix& M, int i) {
    std::string s; s.reserve(NQ);
    for (int j = 0; j < NQ; ++j) s += PAULI[(M.Hx[i] >> j) & 1][(M.Hz[i] >> j) & 1];
    return s;
}
static std::string bits_string(uint16_t v) {
    std::string s; for (int j = 0; j < NQ; ++j) s += char('0' + ((v >> j) & 1));
    return s;
}

// ===========================================================================================
//  SHARED STATE HELPERS
// ===========================================================================================
static inline void compute_c(const uint8_t* types, int n, uint8_t* c, uint8_t& plx) {
    std::memset(c, 0, 16); plx = 0;
    for (int j = 0; j < n; ++j) {
        const TypeInfo& t = TY[types[j]];
        for (int k = 0; k < t.degU; ++k) c[t.chi[k]]++;
        plx ^= t.pl;
    }
}
static inline bool feasible_full(const uint8_t* c, uint8_t plx) {
    if (plx) return false;                                   // generators must commute
    for (int i = 0; i < 4; ++i) if (c[EPT[i]] != TARGET_C) return false;
    for (int chi = 1; chi < 16; ++chi) if (c[chi] > NQ - 1) return false;   // rank 4
    return true;
}

// S_4 canonical form of a (sorted, ascending) multiset of column types.
static void canon_form(const uint8_t* types, int n, uint8_t* out) {
    uint8_t buf[NQ], best[NQ]; bool first = true;
    for (int p = 0; p < NS4; ++p) {
        const int* A = S4act[p];
        for (int j = 0; j < n; ++j) buf[j] = uint8_t(A[types[j]]);
        std::sort(buf, buf + n);
        if (first || std::lexicographical_compare(buf, buf + n, best, best + n)) {
            std::copy(buf, buf + n, best); first = false; }
    }
    std::copy(best, best + n, out);
}
// Valid prefix pruning: if some permutation maps the ascending prefix to a lexicographically
// smaller ascending sequence, then the whole multiset cannot be the S_4 canonical one.
static inline bool is_lexmin_prefix(const uint8_t* types, int n) {
    uint8_t buf[NQ];
    for (int p = 1; p < NS4; ++p) {
        const int* A = S4act[p];
        for (int j = 0; j < n; ++j) buf[j] = uint8_t(A[types[j]]);
        std::sort(buf, buf + n);
        if (std::lexicographical_compare(buf, buf + n, types, types + n)) return false;
    }
    return true;
}

// ===========================================================================================
//  GLOBAL INCUMBENT + REPORTING
// ===========================================================================================
static std::mutex           g_mtx;
static std::atomic<int64_t> g_bestSumP{INFP};
static uint8_t              g_bestTypes[NQ];
static bool                 g_have = false;
static std::atomic<uint64_t> g_nodes{0}, g_pruned{0}, g_leaves{0};
static std::chrono::steady_clock::time_point g_t0;
static std::string          g_outfile = "best_stabilizer_14q.txt";

static double elapsed() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - g_t0).count();
}

static char g_rbuf[65536]; static int g_rn = 0;
static void ap(const char* fmt, ...) {
    va_list a; va_start(a, fmt);
    int room = (int)sizeof(g_rbuf) - g_rn - 1;
    if (room > 0) { int w = vsnprintf(g_rbuf + g_rn, (size_t)room, fmt, a); if (w > 0) g_rn += (w < room ? w : room); }
    va_end(a);
}

static void report(const uint8_t* types, const char* tag) {
    uint8_t c[16], plx; compute_c(types, NQ, c, plx);
    Score S = score_of(c);
    Matrix M = build_matrix(types);
    Verdict V = verify_bruteforce(M);
    g_rn = 0;

    ap("==================================================\n%s\n", tag);
    ap("Detected errors (wt <= 4): %lld / %lld\n", (long long)S.detected, (long long)TOTAL_ERRORS);
    ap("Undetected errors: %lld\n", (long long)S.undetected);
    ap("Detection percentage: %.8f %%\n\n", 100.0 * double(S.detected) / double(TOTAL_ERRORS));
    for (int w = 1; w <= MAXW; ++w)
        ap("Weight-%d detected: %lld / %lld\n", w, (long long)S.det_w[w], (long long)TOTW[w]);
    ap("\nGenerator matrix:\n");
    for (int i = 0; i < NGEN; ++i) ap("g%d = %s   (weight %d)\n", i + 1, row_string(M, i).c_str(), V.weights[i]);
    ap("\nX part:\n");
    for (int i = 0; i < NGEN; ++i) ap("%s\n", bits_string(M.Hx[i]).c_str());
    ap("\nZ part:\n");
    for (int i = 0; i < NGEN; ++i) ap("%s\n", bits_string(M.Hz[i]).c_str());
    ap("\nStabilizer group weight enumerator  (m(chi) = weight of prod g_i^chi_i):\n ");
    for (int chi = 1; chi < 16; ++chi) ap(" %d", NQ - c[chi]);
    ap("\n\nSyndrome distribution:\n");
    for (int s = 0; s < 16; ++s)
        ap("%d%d%d%d %lld\n", (s>>3)&1, (s>>2)&1, (s>>1)&1, s&1, (long long)V.synd[s]);
    ap("\nVERIFICATION\n");
    ap("  dimensions 4 x 28              : OK\n");
    ap("  generator weights == %-2d        : %s\n", GENERATOR_WEIGHT, V.ok_weight ? "OK" : "FAIL");
    ap("  pairwise commutation (6 pairs)  : %s\n", V.ok_commute ? "OK" : "FAIL");
    ap("  rank(H) = 4                     : %s (rank %d)\n", V.ok_rank ? "OK" : "FAIL", V.rank);
    ap("  no identity generator           : %s\n", V.ok_nonidentity ? "OK" : "FAIL");
    ap("  brute-force detected            : %lld  %s\n", (long long)V.detected,
       V.detected == S.detected ? "(matches algebraic score)" : "*** MISMATCH ***");
    for (int w = 1; w <= MAXW; ++w)
        ap("  brute-force weight-%d detected   : %lld %s\n", w, (long long)V.det_w[w],
           V.det_w[w] == S.det_w[w] ? "OK" : "*** MISMATCH ***");
    ap("\nSearch nodes: %llu\nPruned nodes: %llu\nLeaves: %llu\nElapsed time: %.3f s\n",
       (unsigned long long)g_nodes.load(), (unsigned long long)g_pruned.load(),
       (unsigned long long)g_leaves.load(), elapsed());
    ap("==================================================\n");
    fputs(g_rbuf, stdout); fflush(stdout);
    FILE* f = fopen(g_outfile.c_str(), "w"); if (f) { fputs(g_rbuf, f); fclose(f); }
}

// Tie-breaking among configurations with the identical (optimal) objective.  Purely
// cosmetic: first prefer solutions that leave no physical qubit completely untouched
// (dim W_j = 0 means generators act as identity there), then take the S_4-lex-minimum so
// that the reported representative does not depend on thread scheduling.
static int idle_qubits(const uint8_t* types) {
    int n = 0; for (int j = 0; j < NQ; ++j) if (TY[types[j]].dimW == 0) ++n; return n;
}
static int g_bestIdle = NQ + 1;

// returns true if the candidate became the new incumbent
static bool offer(const uint8_t* types, bool quiet = false) {
    uint8_t c[16], plx; compute_c(types, NQ, c, plx);
    if (!feasible_full(c, plx)) return false;
    int64_t sp = sumP_of(c);
    if (sp > g_bestSumP.load(std::memory_order_relaxed)) return false;
    std::lock_guard<std::mutex> lk(g_mtx);
    uint8_t cn[NQ]; canon_form(types, NQ, cn);
    int idle = idle_qubits(cn);
    if (g_have) {
        int64_t bs = g_bestSumP.load();
        if (sp > bs) return false;
        if (sp == bs) {
            if (idle > g_bestIdle) return false;
            if (idle == g_bestIdle &&
                !std::lexicographical_compare(cn, cn + NQ, g_bestTypes, g_bestTypes + NQ)) return false;
        }
    }
    bool improved = (!g_have) || sp < g_bestSumP.load();
    std::copy(cn, cn + NQ, g_bestTypes); g_bestSumP.store(sp); g_bestIdle = idle; g_have = true;
    if (!quiet && improved) report(g_bestTypes, "NEW BEST");
    return true;
}

// ===========================================================================================
//  HEURISTIC PHASE 1 : simulated annealing over the 14-column multiset
//  Moves change one column's subspace type; the objective is recomputed with 15 lookups and
//  constraint violations are handled with a penalty, so the walk may cross infeasible ground.
// ===========================================================================================
static int64_t penalty_of(const uint8_t* c, uint8_t plx) {
    int64_t p = 0;
    for (int i = 0; i < 4; ++i) p += std::abs(int(c[EPT[i]]) - TARGET_C);
    p += 3 * std::popcount(unsigned(plx));
    for (int chi = 1; chi < 16; ++chi) if (c[chi] > NQ - 1) p += c[chi] - (NQ - 1);
    return p;
}
static void anneal(uint64_t seed, double seconds, int restarts) {
    std::mt19937_64 rng(seed);
    const int64_t LAMBDA = 4000;
    uint8_t types[NQ], c[16], plx;
    for (int rs = 0; rs < restarts; ++rs) {
        if (elapsed() > seconds) break;
        for (int j = 0; j < NQ; ++j) types[j] = uint8_t(rng() % unsigned(NTYPE));
        compute_c(types, NQ, c, plx);
        int64_t E = sumP_of(c) + LAMBDA * penalty_of(c, plx);
        const int ITERS = 400000;
        double T0 = 900.0, T1 = 0.35;
        for (int it = 0; it < ITERS; ++it) {
            double T = T0 * std::pow(T1 / T0, double(it) / ITERS);
            int j = int(rng() % NQ);
            uint8_t oldt = types[j], newt = uint8_t(rng() % unsigned(NTYPE));
            if (newt == oldt) continue;
            const TypeInfo& a = TY[oldt]; const TypeInfo& b = TY[newt];
            for (int k = 0; k < a.degU; ++k) c[a.chi[k]]--;
            for (int k = 0; k < b.degU; ++k) c[b.chi[k]]++;
            uint8_t np = uint8_t(plx ^ a.pl ^ b.pl);
            int64_t E2 = sumP_of(c) + LAMBDA * penalty_of(c, np);
            int64_t d = E2 - E;
            bool acc = (d <= 0) || (std::generate_canonical<double, 24>(rng) < std::exp(double(-d) / T));
            if (acc) { types[j] = newt; plx = np; E = E2;
                if (penalty_of(c, np) == 0) offer(types);
            } else {
                for (int k = 0; k < b.degU; ++k) c[b.chi[k]]--;
                for (int k = 0; k < a.degU; ++k) c[a.chi[k]]++;
            }
            if ((it & 8191) == 0 && elapsed() > seconds) break;
        }
    }
}

// ===========================================================================================
//  HEURISTIC PHASE 2 : beam search over the same DFS tree (states = ascending type prefixes)
// ===========================================================================================
struct BState {
    uint8_t types[NQ]; uint8_t c[16]; uint8_t plx; int64_t key;
    std::array<uint8_t, NQ> canon;
};
static void beam_search(int K) {
    std::vector<BState> cur, nxt;
    { BState s{}; std::memset(s.c, 0, 16); s.plx = 0; s.key = 0; cur.push_back(s); }
    for (int depth = 0; depth < NQ; ++depth) {
        nxt.clear();
        int r = NQ - depth - 1;
        for (const BState& s : cur) {
            int t0 = (depth == 0) ? 0 : s.types[depth - 1];
            for (int t = t0; t < NTYPE; ++t) {
                BState ch = s; ch.types[depth] = uint8_t(t);
                const TypeInfo& ti = TY[t];
                bool bad = false;
                for (int k = 0; k < ti.degU; ++k) if (++ch.c[ti.chi[k]] > NQ - 1) bad = true;
                ch.plx = uint8_t(ch.plx ^ ti.pl);
                if (bad) continue;
                int64_t b = node_bound(ch.c, r, t);
                if (b >= INFP) continue;
                ch.key = b;
                canon_form(ch.types, depth + 1, ch.canon.data());
                for (int q = depth + 1; q < NQ; ++q) ch.canon[q] = 0;
                nxt.push_back(ch);
            }
        }
        if (nxt.empty()) return;
        std::sort(nxt.begin(), nxt.end(), [](const BState& a, const BState& b) {
            if (a.key != b.key) return a.key < b.key; return a.canon < b.canon; });
        // deduplicate S_4-equivalent partial configurations, then truncate to the beam width
        std::vector<BState> ded; ded.reserve(nxt.size());
        for (size_t i = 0; i < nxt.size(); ++i) {
            if (i && nxt[i].canon == nxt[i-1].canon && nxt[i].key == nxt[i-1].key) continue;
            ded.push_back(nxt[i]);
            if ((int)ded.size() >= K) break;
        }
        cur.swap(ded);
    }
    for (const BState& s : cur) offer(s.types);
}

// ===========================================================================================
//  EXACT BRANCH-AND-BOUND
//  The tree enumerates every multiset of 14 column subspaces exactly once (columns are added
//  in non-decreasing type index, which is precisely the quotient by the S_14 qubit
//  permutations).  Nodes die from: (i) the four generator-weight equations c(e_i)=TARGET_C
//  becoming unreachable, (ii) rank loss c(chi)=14, (iii) the admissible objective bound,
//  (iv) the S_4 lex-minimality test on the prefix.  Nothing else is discarded, so a completed
//  run proves global optimality.
// ===========================================================================================
static std::atomic<bool> g_abort{false};
static double g_timeLimit = 1e18;
static int    g_symDepth  = 6;
static bool   g_tieMode   = false;
static int    g_split     = 2;

struct Searcher {
    uint8_t types[NQ]{}; uint8_t c[16]{}; uint8_t plx = 0;
    uint64_t nodes = 0, pruned = 0, leaves = 0;
    int splitAt = -1;
    std::vector<std::array<uint8_t, NQ>>* sink = nullptr;

    void dfs(int depth, int t0) {
        if (g_abort.load(std::memory_order_relaxed)) return;
        ++nodes;
        if ((nodes & 0x3FFFF) == 0 && elapsed() > g_timeLimit) { g_abort.store(true); return; }
        const int r = NQ - depth;
        const int64_t best = g_bestSumP.load(std::memory_order_relaxed);
        const int64_t b = node_bound(c, r, t0);
        if (b >= INFP || (g_tieMode ? (b > best) : (b >= best))) { ++pruned; return; }
        if (r == 0) { ++leaves; if (plx == 0) offer(types); return; }
        if (splitAt >= 0 && depth == splitAt) {
            std::array<uint8_t, NQ> a{}; std::copy(types, types + depth, a.begin());
            sink->push_back(a); return;
        }
        for (int t = t0; t < NTYPE; ++t) {
            const TypeInfo& ti = TY[t];
            bool bad = false;
            for (int k = 0; k < ti.degU; ++k) if (++c[ti.chi[k]] > NQ - 1) bad = true;
            plx = uint8_t(plx ^ ti.pl); types[depth] = uint8_t(t);
            if (!bad && (depth + 1 > g_symDepth || is_lexmin_prefix(types, depth + 1)))
                dfs(depth + 1, t);
            plx = uint8_t(plx ^ ti.pl);
            for (int k = 0; k < ti.degU; ++k) c[ti.chi[k]]--;
        }
    }
};

static void exact_search(int nthreads) {
    std::vector<std::array<uint8_t, NQ>> tasks;
    { Searcher S; S.splitAt = std::min(g_split, NQ); S.sink = &tasks; S.dfs(0, 0);
      g_nodes += S.nodes; g_pruned += S.pruned; }
    printf("[exact] %zu independent subtrees after depth-%d split; %d threads\n",
           tasks.size(), g_split, nthreads);
    fflush(stdout);
    std::atomic<size_t> next{0};
    std::vector<std::thread> th;
    for (int w = 0; w < nthreads; ++w) th.emplace_back([&]() {
        Searcher S;
        uint64_t n = 0, p = 0, l = 0;
        for (;;) {
            size_t i = next.fetch_add(1);
            if (i >= tasks.size() || g_abort.load()) break;
            const auto& pre = tasks[i];
            int d = std::min(g_split, NQ);
            std::copy(pre.begin(), pre.begin() + d, S.types);
            compute_c(S.types, d, S.c, S.plx);
            S.nodes = S.pruned = S.leaves = 0;
            S.dfs(d, S.types[d - 1]);
            n += S.nodes; p += S.pruned; l += S.leaves;
        }
        g_nodes += n; g_pruned += p; g_leaves += l;
    });
    for (auto& t : th) t.join();
}

// ===========================================================================================
//  SELF-TESTS AND BENCHMARK
// ===========================================================================================
static bool self_test(uint64_t seed, int trials) {
    std::mt19937_64 rng(seed); bool ok = true;
    for (int it = 0; it < trials; ++it) {
        uint8_t types[NQ], c[16], plx;
        for (int j = 0; j < NQ; ++j) types[j] = uint8_t(rng() % unsigned(NTYPE));
        compute_c(types, NQ, c, plx);
        Score S = score_of(c);
        Matrix M = build_matrix(types);
        Verdict V = verify_bruteforce(M);
        // (a) algebraic weight enumerator == brute-force generator weights
        for (int i = 0; i < NGEN; ++i) if (V.weights[i] != NQ - c[EPT[i]]) {
            printf("SELFTEST FAIL: generator weight mismatch\n"); ok = false; }
        // (b) algebraic Plucker commutation == brute-force commutation
        if ((plx == 0) != V.ok_commute) {
            printf("SELFTEST FAIL: commutation mismatch plx=%02x bruteforce_ok=%d\n  types:", plx, (int)V.ok_commute);
            for (int j = 0; j < NQ; ++j) printf(" t%d(d%d,pl%02x)", types[j], TY[types[j]].dimW, TY[types[j]].pl);
            printf("\n  bf pairs:");
            for (int i2 = 0; i2 < NGEN; ++i2) for (int k2 = i2 + 1; k2 < NGEN; ++k2)
                printf(" %d", par(unsigned(M.Hx[i2] & M.Hz[k2])) ^ par(unsigned(M.Hz[i2] & M.Hx[k2])));
            printf("\n");
            ok = false;
        }
        // (c) algebraic objective == brute-force objective, total and per weight
        if (S.detected != V.detected) { printf("SELFTEST FAIL: detected %lld vs %lld\n",
            (long long)S.detected, (long long)V.detected); ok = false; }
        for (int w = 1; w <= MAXW; ++w) if (S.det_w[w] != V.det_w[w]) {
            printf("SELFTEST FAIL: weight-%d %lld vs %lld\n", w,
                   (long long)S.det_w[w], (long long)V.det_w[w]); ok = false; }
        if (!ok) return false;
    }
    return ok;
}

static void bench(uint64_t seed) {
    std::mt19937_64 rng(seed);
    const int N = 20000;
    std::vector<std::array<uint8_t, NQ>> cfg(N);
    for (auto& a : cfg) for (int j = 0; j < NQ; ++j) a[j] = uint8_t(rng() % unsigned(NTYPE));
    auto t1 = std::chrono::steady_clock::now();
    int64_t acc = 0;
    for (auto& a : cfg) { uint8_t c[16], plx; compute_c(a.data(), NQ, c, plx); acc += sumP_of(c); }
    auto t2 = std::chrono::steady_clock::now();
    for (int i = 0; i < N / 200; ++i) { Matrix M = build_matrix(cfg[i].data());
        Verdict V = verify_bruteforce(M); acc += V.detected; }
    auto t3 = std::chrono::steady_clock::now();
    double a1 = std::chrono::duration<double>(t2 - t1).count() / N;
    double a2 = std::chrono::duration<double>(t3 - t2).count() / (N / 200);
    printf("[bench] column-type algebraic score : %10.1f ns / candidate\n", a1 * 1e9);
    printf("[bench] 91770-error brute force     : %10.1f ns / candidate  (speed-up %.0fx)\n",
           a2 * 1e9, a2 / a1);
    printf("[bench] checksum %lld\n", (long long)acc);
}

// ===========================================================================================
//  MAIN
// ===========================================================================================
static void usage() {
    printf("stab14 -- exact optimiser for a 4-generator, 14-qubit stabilizer check matrix\n"
           "  --mode heuristic|exact|hybrid   (default hybrid)\n"
           "  --threads N                     (default hardware concurrency)\n"
           "  --beam-size K                   (default 20000)\n"
           "  --time-limit SECONDS            (default none)\n"
           "  --seed S                        (default 12345)\n"
           "  --generator-weight W            (default %d)\n"
           "  --sym-depth D                   S_4 lex-min prefix pruning depth (default 6)\n"
           "  --split D                       depth at which the tree is cut for threads (2)\n"
           "  --deterministic                 keep ties, report the S_4-lex-min optimum\n"
           "  --selftest N                    formula-vs-brute-force trials (default 200)\n"
           "  --bench                         benchmark the two scoring routes and exit\n"
           "  --out FILE                      (default best_stabilizer_14q.txt)\n",
           GENERATOR_WEIGHT_DEFAULT);
}

int main(int argc, char** argv) {
    std::string mode = "hybrid";
    int nthreads = (int)std::max(1u, std::thread::hardware_concurrency());
    int beamK = 20000, selftestN = 200;
    uint64_t seed = 12345ULL;
    bool doBench = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto nxt = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : std::string(); };
        if      (a == "--mode")             mode = nxt();
        else if (a == "--threads")          nthreads = std::max(1, atoi(nxt().c_str()));
        else if (a == "--beam-size")        beamK = std::max(1, atoi(nxt().c_str()));
        else if (a == "--time-limit")       g_timeLimit = atof(nxt().c_str());
        else if (a == "--seed")             seed = strtoull(nxt().c_str(), nullptr, 10);
        else if (a == "--generator-weight") GENERATOR_WEIGHT = atoi(nxt().c_str());
        else if (a == "--sym-depth")        g_symDepth = atoi(nxt().c_str());
        else if (a == "--split")            g_split = std::max(1, atoi(nxt().c_str()));
        else if (a == "--deterministic")  { g_tieMode = true; }
        else if (a == "--selftest")         selftestN = atoi(nxt().c_str());
        else if (a == "--bench")            doBench = true;
        else if (a == "--out")              g_outfile = nxt();
        else { usage(); return (a == "--help" || a == "-h") ? 0 : 1; }
    }
    if (GENERATOR_WEIGHT < 1 || GENERATOR_WEIGHT > NQ) { printf("bad generator weight\n"); return 1; }
    TARGET_C = NQ - GENERATOR_WEIGHT;

    g_t0 = std::chrono::steady_clock::now();
    init_tables(); build_types(); build_s4(); init_bound_tables();

    printf("=== 14-qubit / 4-generator stabilizer detection optimiser ===\n");
    printf("GENERATOR_WEIGHT = %d   (required c(e_i) = 14 - w = %d)\n", GENERATOR_WEIGHT, TARGET_C);
    printf("error universe   = %lld errors of Pauli weight 1..4 (wt(Y)=1)\n", (long long)TOTAL_ERRORS);
    printf("column types     = %d subspaces of F_2^4 (dim 0/1/2 = 1/15/35)\n", NTYPE);
    printf("mode = %s, threads = %d, beam = %d, seed = %llu, sym-depth = %d%s\n",
           mode.c_str(), nthreads, beamK, (unsigned long long)seed, g_symDepth,
           g_tieMode ? ", deterministic" : "");
    printf("P(m) table (contribution of one group element of weight m):\n ");
    for (int m = 0; m <= NQ; ++m) printf(" %lld", (long long)Ptab[m]);
    printf("\n\n");

    if (doBench) { bench(seed); return 0; }

    if (selftestN > 0) {
        printf("[selftest] %d random configurations, algebraic score vs brute force ... ", selftestN);
        fflush(stdout);
        if (!self_test(seed, selftestN)) { printf("FAILED\n"); return 2; }
        printf("all identical, OK\n\n");
    }

    if (mode == "heuristic" || mode == "hybrid") {
        printf("[heuristic] beam search (K = %d) ...\n", beamK); fflush(stdout);
        beam_search(beamK);
        printf("[heuristic] simulated annealing ...\n"); fflush(stdout);
        anneal(seed, std::min(g_timeLimit, elapsed() + 20.0), 24);
        printf("[heuristic] best sumP = %lld after %.2f s\n\n",
               (long long)g_bestSumP.load(), elapsed());
    }
    if (mode == "exact" || mode == "hybrid") {
        printf("[exact] branch and bound over all multisets of 14 column subspaces ...\n");
        fflush(stdout);
        exact_search(nthreads);
        if (g_abort.load()) printf("\n*** TIME LIMIT HIT -- result is the best found, NOT proven optimal ***\n");
        else                printf("\n*** SEARCH TREE EXHAUSTED -- the solution below is the GLOBAL OPTIMUM ***\n");
    }
    if (!g_have) { printf("no feasible configuration found\n"); return 3; }
    report(g_bestTypes, g_abort.load() ? "FINAL BEST (not proven optimal)" : "FINAL BEST (proven global optimum)");
    return 0;
}
