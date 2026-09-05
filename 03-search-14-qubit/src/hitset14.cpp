// =========================================================================================
//  hitset14.cpp -- residual-targeted COLUMN replacement search for a [[14,3,5]] stabilizer
//
//  Independent of every earlier program in this project.  Nothing outside residual_columns/
//  is read or written.
//
// -----------------------------------------------------------------------------------------
//  0.  THE REFORMULATION THIS PROGRAM IS BUILT ON
//
//  The check matrix is H = [X|Z], 11 x 28.  Read it COLUMN-wise: qubit j contributes
//  x_j, z_j in F_2^11 (the j-th columns of X and of Z).  Put
//
//        W_j  =  span{ x_j , z_j }  <=  F_2^11 .
//
//  Syndromes.  An error E with Pauli X^{a_j} Z^{b_j} on qubit j has syndrome
//  s(E) = SUM_j ( a_j z_j + b_j x_j ).  So the syndrome contributed by qubit j is
//        X -> z_j ,   Z -> x_j ,   Y -> x_j + z_j ,
//  i.e. the three NON-ZERO elements of W_j, each exactly once.  Therefore
//
//     d >= 5   <=>   no 1 <= k <= 4 distinct qubits j_1..j_k and non-zero v_i in W_{j_i}
//                    have  v_1 + ... + v_k = 0 .                                    (*)
//
//  Everything follows from (*).  Three immediate consequences:
//
//    * dim W_j = 2 for every j.  If dim W_j <= 1 then one of the three contributions is 0,
//      which is a weight-1 undetected error.  So each column is a LINE of PG(10,2).
//    * k=2 forces W_i ∩ W_j = 0 for i != j.
//    * the condition depends only on the 14 subspaces, NOT on the basis (x_j,z_j) chosen
//      inside each.  Changing that basis is exactly a local Clifford on qubit j
//      (GL(2,2) = S_3 permuting the labels X,Y,Z), so local Clifford is quotiented out
//      for free by working with subspaces.  That is a 6^14 = 7.8e10 reduction.
//
//  Commutation.  <g_a,g_b> = SUM_j (x_{aj} z_{bj} + z_{aj} x_{bj}), so with the PLUCKER
//  form  pl(W) = x z^T + z x^T  (symmetric, zero diagonal, hence a point of the 55-dimensional
//  Lambda^2(F_2^11)), the whole commutation requirement is the single linear equation
//
//        SUM_{j=1..14}  pl(W_j)  =  0    in  Lambda^2(F_2^11) = F_2^55 .            (**)
//
//  pl(W) is basis-independent: a change of basis multiplies it by det = 1.  Under the
//  Pluecker embedding, LINES of PG(10,2) correspond bijectively to the DECOMPOSABLE
//  non-zero 2-vectors x^z, equivalently to the alternating forms of rank exactly 2.
//
//  Rank.  rank H = 11  <=>  W_1 + ... + W_14 = F_2^11.
//
//  Generator basis change.  A in GL(11,2) sends W_j -> A W_j and pl -> A pl A^T, so it
//  preserves (*) and (**).  It is a symmetry of the problem but it moves ALL columns at
//  once, so it is broken as soon as we freeze some columns; it is not used for pruning.
//
// -----------------------------------------------------------------------------------------
//  1.  WHY THE HITTING-SET CONDITION IS EXACTLY RIGHT (AND WHERE IT COMES FROM)
//
//  Freeze the columns outside a set S and let the r = |S| columns in S vary.  A weight-<=4
//  error supported entirely OUTSIDE S has a syndrome that does not involve any free column,
//  so it is unchanged by the search.  Hence
//
//        S must meet the support of every currently-undetected weight-<=4 error,
//
//  which for the five given residuals is precisely a hitting-set condition.  It is
//  NECESSARY, and it is not an extra heuristic: it is the f = 0 case of the decomposition
//  below.  Split an error's support into its free part (f qubits) and its fixed part, and
//  let C_m = { syndromes of errors of weight <= m supported on FIXED qubits } (with the
//  empty error contributing 0).  Then (*) becomes, for the free points v_i,
//
//        f = 0 :  0 not in C_4 \ {0}          <- feasibility of S: the hitting condition
//        f = 1 :  v_1 not in C_3
//        f = 2 :  v_1 + v_2 not in C_2
//        f = 3 :  v_1 + v_2 + v_3 not in C_1
//        f = 4 :  v_1 + v_2 + v_3 + v_4 != 0
//
//  C_3, C_2, C_1 are computed once per S as 2048-bit tables.  The f=1 row is the powerful
//  one: EVERY free column must be a line all three of whose points avoid C_3, and that set
//  is the same for all free qubits.  Call the surviving lines the ALLOWED SET A.
//
// -----------------------------------------------------------------------------------------
//  2.  THE STRONGEST FORMULATION I FOUND: DECOMPOSING AN ALTERNATING FORM
//
//  By (**) the free columns must satisfy   SUM_{j in S} pl(W_j) = T,   where
//  T = SUM_{j not in S} pl(W_j) is a fixed alternating form (equivalently T is the Pluecker
//  sum of the ORIGINAL columns of S, since the original code satisfies (**)).
//
//  So the search is: WRITE T AS A SUM OF r DECOMPOSABLE 2-VECTORS, each drawn from A,
//  subject to the f>=2 cross conditions above.  This is a classical object, and it gives an
//  exact prune that costs one 11x11 F_2 rank computation per node:
//
//        a sum of t decomposables has rank <= 2t  (each has rank 2),
//        and conversely every form of rank 2m IS a sum of m decomposables (Darboux),
//        so  min #decomposables summing to T  =  rank(T)/2 .
//
//     PRUNE (exact):  after choosing k of the r columns, the remaining target T_k must obey
//                     rank(T_k) <= 2 (r - k).
//
//  Nothing that could complete is discarded, because any completion writes T_k as a sum of
//  r-k decomposables and therefore has rank <= 2(r-k).  Two consequences make the search
//  collapse:
//
//    * at depth r-1 exactly one column is left, so rank(T_{r-1}) must be EXACTLY 2 and the
//      last line is FORCED: it is the row space of T_{r-1}.  No enumeration at the last level.
//    * subtracting a decomposable changes the rank by -2, 0 or +2, and it drops by 2 exactly
//      when the line is a HYPERBOLIC plane of T_k (T_k(x,z) = 1).  When rank(T) = 2r every
//      single step must drop, which is extremely restrictive.
//
//  This is why the search is exhaustive and fast rather than a 698027^r enumeration.
//
// -----------------------------------------------------------------------------------------
//  3.  DEGENERACY
//
//  (*) forbids ANY weight-<=4 Pauli in the centralizer N(S), which is the PURE (non-
//  degenerate) condition.  A degenerate code is allowed to contain weight-<=4 elements of
//  the stabilizer group S itself, since those act trivially.  The element of S with
//  coefficient vector a in F_2^11 acts non-trivially on qubit j iff a is NOT orthogonal to
//  W_j, so
//        wt(E_a) = 14 - #{ j : a _|_ W_j } ,
//  and a weight-<=4 stabilizer element needs a non-zero a orthogonal to at least 10 columns.
//  At most r of those can be free, so a must be orthogonal to at least 10 - r of the FIXED
//  columns.  The program checks this for every a in F_2^11 for each S: when the maximum
//  achievable orthogonality count is < 10, NO degenerate solution can exist for that S and
//  the pure search is therefore exhaustive over ALL solutions, degenerate ones included.
//  That closure is reported per subset, so the exhaustiveness claim is never assumed.
// =========================================================================================

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cstdarg>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <algorithm>
#include <numeric>
#include <atomic>
#include <mutex>
#include <thread>
#include <chrono>
#include <bit>
#include <filesystem>
#include <random>
#include <ctime>

constexpr int NQ   = 14;          // physical qubits
constexpr int RG   = 11;          // stabilizer generators
constexpr int MAXW = 4;           // error weight window
constexpr int VN   = 1 << RG;     // 2048 syndrome vectors
constexpr uint32_t VMASK = VN - 1;

// ---------------------------------------------------------------- timing / logging
static long long now_ns() {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(
               std::chrono::steady_clock::now().time_since_epoch()).count();
}
static std::atomic<long long> T0_NS{0};
static void set_start_time() { T0_NS.store(now_ns(), std::memory_order_relaxed); }
static double elapsed() {
    long long t0 = T0_NS.load(std::memory_order_relaxed);
    if (t0 == 0) { long long n = now_ns(), e = 0;
        t0 = T0_NS.compare_exchange_strong(e, n) ? n : e; }
    return double(now_ns() - t0) * 1e-9;
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
static std::string RB;
static void ap(const char* fmt, ...) {
    char buf[16384]; va_list a; va_start(a, fmt); vsnprintf(buf, sizeof(buf), fmt, a); va_end(a);
    RB += buf;
}
static void dump(const std::string& path) {
    FILE* f = fopen(path.c_str(), "w"); if (f) { fputs(RB.c_str(), f); fclose(f); }
}

// ---------------------------------------------------------------- Pauli rows (for I/O only)
// A generator is stored as two 14-bit masks (X part, Z part).
struct Gen { uint32_t x = 0, z = 0; };
static const char PCH[2][2] = { {'I','Z'}, {'X','Y'} };
static std::string gstr(const Gen& g) {
    std::string s;
    for (int j = 0; j < NQ; ++j) s += PCH[(g.x >> j) & 1][(g.z >> j) & 1];
    return s;
}
static Gen parse_gen(const char* s) {
    Gen g;
    for (int j = 0; j < NQ; ++j) {
        switch (s[j]) {
            case 'I': break;
            case 'X': g.x |= 1u << j; break;
            case 'Z': g.z |= 1u << j; break;
            case 'Y': g.x |= 1u << j; g.z |= 1u << j; break;
            default: fprintf(stderr, "bad Pauli char in %s\n", s); exit(1);
        }
    }
    return g;
}
static int gwt(const Gen& g) { return std::popcount(g.x | g.z); }
static int par(uint32_t a) { return std::popcount(a) & 1; }
static int commute(const Gen& a, const Gen& b) { return par((a.x & b.z) ^ (a.z & b.x)); }

// ---------------------------------------------------------------- the starting code
static const char* START_G[RG] = {
    "IIIIXXXXIIZZXX",   // g1  weight 8
    "IIIIXXIIXXXXZZ",   // g2  weight 8
    "IIIIIIZZZZZZZZ",   // g3  weight 8
    "IXIXXXIZYZIXXX",   // h1  weight 10
    "IZXIXXYIXYIYXX",   // h2
    "IIIZYXZXYXZIZY",   // h3
    "IIIYXZYIYZZYZX",   // h4
    "IIZIXZXYZXZYIZ",   // h5
    "XIIIXYXZXYIXYX",   // h6
    "IIYIZXZXYXXIYX",   // h7
    "ZIIIZXIZZYZXYX"    // h8
};

// ---------------------------------------------------------------- columns as lines of PG(10,2)
// Column j is the pair (x_j, z_j) in F_2^11 x F_2^11; its three syndrome points are
//    X -> z_j,  Z -> x_j,  Y -> x_j ^ z_j.
struct Col { uint32_t x = 0, z = 0; };
static inline void col_points(const Col& c, uint32_t p[3]) {
    p[0] = c.z; p[1] = c.x; p[2] = c.x ^ c.z;
}
static inline bool col_ok(const Col& c) {          // dim W_j must be 2
    return c.x && c.z && c.x != c.z;
}

// ---------------------------------------------------------------- alternating forms
// An 11x11 symmetric zero-diagonal F_2 matrix, stored as 11 rows of 11 bits.
struct Alt {
    uint16_t r[RG] = {0,0,0,0,0,0,0,0,0,0,0};
    bool operator==(const Alt& o) const {
        for (int i = 0; i < RG; ++i) if (r[i] != o.r[i]) return false; return true; }
};
static inline void alt_xor(Alt& a, const Alt& b) { for (int i = 0; i < RG; ++i) a.r[i] ^= b.r[i]; }
static inline bool alt_zero(const Alt& a) {
    for (int i = 0; i < RG; ++i) if (a.r[i]) return false; return true; }
// pl(W) = x z^T + z x^T
static Alt plucker(uint32_t x, uint32_t z) {
    Alt a;
    for (int i = 0; i < RG; ++i) {
        uint16_t row = 0;
        if ((x >> i) & 1) row ^= (uint16_t)z;
        if ((z >> i) & 1) row ^= (uint16_t)x;
        row &= (uint16_t)~(1u << i);              // zero diagonal
        a.r[i] = row;
    }
    return a;
}
static int alt_rank(const Alt& a) {
    uint16_t m[RG]; for (int i = 0; i < RG; ++i) m[i] = a.r[i];
    int rank = 0;
    for (int col = 0; col < RG; ++col) {
        int piv = -1;
        for (int i = rank; i < RG; ++i) if ((m[i] >> col) & 1) { piv = i; break; }
        if (piv < 0) continue;
        std::swap(m[rank], m[piv]);
        for (int i = 0; i < RG; ++i)
            if (i != rank && ((m[i] >> col) & 1)) m[i] ^= m[rank];
        ++rank;
    }
    return rank;
}
// Row space of an alternating form = its support.  Returns the dimension and a basis.
static int alt_support(const Alt& a, uint32_t* basis) {
    uint16_t m[RG]; for (int i = 0; i < RG; ++i) m[i] = a.r[i];
    int n = 0;
    for (int col = 0; col < RG; ++col) {
        int piv = -1;
        for (int i = n; i < RG; ++i) if ((m[i] >> col) & 1) { piv = i; break; }
        if (piv < 0) continue;
        std::swap(m[n], m[piv]);
        for (int i = 0; i < RG; ++i)
            if (i != n && ((m[i] >> col) & 1)) m[i] ^= m[n];
        basis[n] = m[n];
        ++n;
    }
    return n;
}
// For a rank-2 form the row space IS the line; return its two basis vectors.
static bool alt_support2(const Alt& a, uint32_t& u, uint32_t& v) {
    uint16_t m[RG]; for (int i = 0; i < RG; ++i) m[i] = a.r[i];
    uint32_t b[2]; int n = 0;
    for (int col = 0; col < RG && n < 3; ++col) {
        int piv = -1;
        for (int i = n; i < RG; ++i) if ((m[i] >> col) & 1) { piv = i; break; }
        if (piv < 0) continue;
        std::swap(m[n], m[piv]);
        for (int i = 0; i < RG; ++i)
            if (i != n && ((m[i] >> col) & 1)) m[i] ^= m[n];
        if (n < 2) b[n] = m[n];
        ++n;
    }
    if (n != 2) return false;
    u = b[0]; v = b[1];
    return true;
}

// =========================================================================================
//  THE STARTING CODE, ITS COLUMNS, AND FULLY INDEPENDENT VERIFICATION
// =========================================================================================
static Gen G0[RG];
static Col COL0[NQ];

static void gens_to_cols(const Gen* g, Col* c) {
    for (int j = 0; j < NQ; ++j) {
        uint32_t X = 0, Z = 0;
        for (int i = 0; i < RG; ++i) {
            if ((g[i].x >> j) & 1) X |= 1u << i;
            if ((g[i].z >> j) & 1) Z |= 1u << i;
        }
        c[j].x = X; c[j].z = Z;
    }
}
static void cols_to_gens(const Col* c, Gen* g) {
    for (int i = 0; i < RG; ++i) { g[i].x = 0; g[i].z = 0; }
    for (int j = 0; j < NQ; ++j)
        for (int i = 0; i < RG; ++i) {
            if ((c[j].x >> i) & 1) g[i].x |= 1u << j;
            if ((c[j].z >> i) & 1) g[i].z |= 1u << j;
        }
}
// rank of the 11 x 28 check matrix over F_2
static int gen_rank(const Gen* g) {
    uint64_t m[RG];
    for (int i = 0; i < RG; ++i) m[i] = (uint64_t)g[i].x | ((uint64_t)g[i].z << 14);
    int rank = 0;
    for (int b = 0; b < 28; ++b) {
        int piv = -1;
        for (int i = rank; i < RG; ++i) if ((m[i] >> b) & 1) { piv = i; break; }
        if (piv < 0) continue;
        std::swap(m[rank], m[piv]);
        for (int i = 0; i < RG; ++i) if (i != rank && ((m[i] >> b) & 1)) m[i] ^= m[rank];
        ++rank;
    }
    return rank;
}
// is the Pauli (ex,ez) in the stabilizer GROUP generated by g?
static bool in_stabilizer(const Gen* g, uint32_t ex, uint32_t ez) {
    uint64_t m[RG];
    for (int i = 0; i < RG; ++i) m[i] = (uint64_t)g[i].x | ((uint64_t)g[i].z << 14);
    uint64_t t = (uint64_t)ex | ((uint64_t)ez << 14);
    int rank = 0;
    for (int b = 0; b < 28; ++b) {
        int piv = -1;
        for (int i = rank; i < RG; ++i) if ((m[i] >> b) & 1) { piv = i; break; }
        if (piv < 0) continue;
        std::swap(m[rank], m[piv]);
        for (int i = rank + 1; i < RG; ++i) if ((m[i] >> b) & 1) m[i] ^= m[rank];
        if ((t >> b) & 1) t ^= m[rank];
        ++rank;
    }
    return t == 0;
}
// minimum non-zero weight of the stabilizer group (2048 elements)
static int stab_min_weight(const Gen* g) {
    int best = NQ + 1;
    for (int a = 1; a < VN; ++a) {
        uint32_t ex = 0, ez = 0;
        for (int i = 0; i < RG; ++i) if ((a >> i) & 1) { ex ^= g[i].x; ez ^= g[i].z; }
        int w = std::popcount(ex | ez);
        if (w && w < best) best = w;
    }
    return best;
}

// ---------------------------------------------------------------- error enumeration
template <class F> static void for_each_error(F&& f) {
    static const uint32_t PX3[3] = {1,1,0}, PZ3[3] = {0,1,1};   // X, Y, Z
    int idx[MAXW];
    for (int w = 1; w <= MAXW; ++w) {
        for (int i = 0; i < w; ++i) idx[i] = i;
        for (;;) {
            int npat = 1; for (int i = 0; i < w; ++i) npat *= 3;
            for (int p = 0; p < npat; ++p) {
                uint32_t ex = 0, ez = 0; int q = p;
                for (int i = 0; i < w; ++i) { int a = q % 3; q /= 3;
                    ex |= PX3[a] << idx[i]; ez |= PZ3[a] << idx[i]; }
                f(ex, ez, w);
            }
            int i = w - 1; while (i >= 0 && idx[i] == NQ - w + i) --i;
            if (i < 0) break;
            ++idx[i]; for (int k = i + 1; k < w; ++k) idx[k] = idx[k-1] + 1;
        }
    }
}
static uint32_t syndrome(const Gen* g, uint32_t ex, uint32_t ez) {
    uint32_t s = 0;
    for (int i = 0; i < RG; ++i) if (par((g[i].x & ez) ^ (g[i].z & ex))) s |= 1u << i;
    return s;
}
static std::string estr(uint32_t ex, uint32_t ez) {
    std::string s; bool first = true;
    for (int j = 0; j < NQ; ++j) {
        int a = (ex >> j) & 1, b = (ez >> j) & 1;
        if (!a && !b) continue;
        if (!first) s += " ";
        s += PCH[a][b]; s += "_" + std::to_string(j); first = false;
    }
    return s.empty() ? std::string("I") : s;
}

struct CodeReport {
    bool weights_ok = false, commute_ok = false, rank_ok = false, cols_ok = false;
    bool plucker_ok = false;
    int rank = 0, stabmin = 0;
    int64_t total_err = 0;
    std::vector<std::pair<uint32_t,uint32_t>> residual;      // zero syndrome, weight<=4
    std::vector<std::pair<uint32_t,uint32_t>> residual_nondeg; // and NOT in the stabilizer
    int64_t by_w[MAXW + 1] = {0,0,0,0,0};
    int gw[RG] = {0};
};
// Independent: rebuilds every error from scratch and uses direct symplectic products.
static CodeReport analyse_code(const Gen* g) {
    CodeReport R;
    R.weights_ok = true;
    for (int i = 0; i < RG; ++i) R.gw[i] = gwt(g[i]);
    R.commute_ok = true;
    for (int i = 0; i < RG; ++i) for (int k = i + 1; k < RG; ++k)
        if (commute(g[i], g[k])) R.commute_ok = false;
    R.rank = gen_rank(g); R.rank_ok = (R.rank == RG);
    Col c[NQ]; gens_to_cols(g, c);
    R.cols_ok = true;
    for (int j = 0; j < NQ; ++j) if (!col_ok(c[j])) R.cols_ok = false;
    Alt sum;
    for (int j = 0; j < NQ; ++j) { Alt p = plucker(c[j].x, c[j].z); alt_xor(sum, p); }
    R.plucker_ok = alt_zero(sum);
    for_each_error([&](uint32_t ex, uint32_t ez, int w) {
        ++R.total_err;
        if (syndrome(g, ex, ez) == 0) {
            R.residual.push_back({ex, ez});
            R.by_w[w]++;
            if (!in_stabilizer(g, ex, ez)) R.residual_nondeg.push_back({ex, ez});
        }
    });
    R.stabmin = stab_min_weight(g);
    return R;
}
// Second, completely different distance routine: works only from the column subspaces and
// checks condition (*) directly.  Returns the number of violating (k <= 4)-tuples.
static int64_t column_violations(const Col* c) {
    uint32_t P[NQ][3];
    for (int j = 0; j < NQ; ++j) col_points(c[j], P[j]);
    int64_t bad = 0;
    for (int a = 0; a < NQ; ++a) for (int ia = 0; ia < 3; ++ia) {
        if (P[a][ia] == 0) { ++bad; continue; }
        for (int b = a + 1; b < NQ; ++b) for (int ib = 0; ib < 3; ++ib) {
            uint32_t s2 = P[a][ia] ^ P[b][ib];
            if (s2 == 0) ++bad;
            for (int d = b + 1; d < NQ; ++d) for (int id = 0; id < 3; ++id) {
                uint32_t s3 = s2 ^ P[d][id];
                if (s3 == 0) ++bad;
                for (int e = d + 1; e < NQ; ++e) for (int ie = 0; ie < 3; ++ie)
                    if ((s3 ^ P[e][ie]) == 0) ++bad;
            }
        }
    }
    return bad;
}

// =========================================================================================
//  PER-SUBSET SETUP:  C_1, C_2, C_3  AND THE ALLOWED LINE SET
// =========================================================================================
struct Bits {
    uint64_t w[VN / 64];
    void clear() { memset(w, 0, sizeof(w)); }
    inline void set(uint32_t v) { w[v >> 6] |= 1ull << (v & 63); }
    inline bool test(uint32_t v) const { return (w[v >> 6] >> (v & 63)) & 1; }
    int count() const { int c = 0; for (uint64_t q : w) c += std::popcount(q); return c; }
};
struct Line {
    uint32_t p[3];       // the three non-zero points, p[0] < p[1] < p[2]
    Alt pl;
    uint32_t key;        // p[0] * 2048 + p[1]
};
static inline Line make_line(uint32_t a, uint32_t b) {
    Line L; uint32_t t[3] = { a, b, a ^ b };
    std::sort(t, t + 3);
    L.p[0] = t[0]; L.p[1] = t[1]; L.p[2] = t[2];
    L.pl = plucker(t[0], t[1]);
    L.key = t[0] * VN + t[1];
    return L;
}

// ---- the running forbidden-point set --------------------------------------------------
// F[w] = syndromes of the weight-w errors supported on the FIXED qubits together with the
// free columns already chosen.  A further column is admissible exactly when none of its
// three points lies in F[0] u F[1] u F[2] u F[3]: adding that point would complete a
// weight-<=4 error of syndrome zero.  This is the same set of conditions as the f = 1..4
// rows of the decomposition, but collapsed into ONE table, so testing a line costs three
// bit lookups instead of re-deriving every pair, triple and quadruple.
struct FSets {
    Bits f[4];
    void forbidden(Bits& D) const {
        for (size_t i = 0; i < VN / 64; ++i)
            D.w[i] = f[0].w[i] | f[1].w[i] | f[2].w[i] | f[3].w[i];
    }
    // extend by one new column with the three points p[0..2]
    void extend(const uint32_t* p, FSets& out) const {
        out = *this;
        for (int w = 3; w >= 1; --w) {
            const Bits& src = f[w - 1];
            for (uint32_t base = 0; base < VN; base += 64) {
                uint64_t word = src.w[base >> 6];
                while (word) {
                    uint32_t v = base + (uint32_t)std::countr_zero(word);
                    word &= word - 1;
                    out.f[w].set(v ^ p[0]);
                    out.f[w].set(v ^ p[1]);
                    out.f[w].set(v ^ p[2]);
                }
            }
        }
    }
};
static inline bool line_ok(const Bits& D, const Line& L) {
    return !D.test(L.p[0]) && !D.test(L.p[1]) && !D.test(L.p[2]);
}
struct SubsetCtx {
    int r = 0;
    std::vector<int> freeQ, fixedQ;
    Bits C1, C2, C3;                  // syndromes of weight<=1, <=2, <=3 errors on FIXED qubits
    FSets F0;                         // the same, split by exact weight (see FSets)
    bool feasible = false;            // no weight<=4 all-fixed error with zero syndrome
    Alt T;                            // required Pluecker sum of the free columns
    int rankT = 0;
    std::vector<Line> A;              // allowed lines, sorted by key
    std::vector<uint32_t> Akey;
    std::vector<std::vector<int>> through;   // through[p] = allowed lines containing point p
    int maxOrth = 0;                  // max over a!=0 of #{fixed j : a _|_ W_j}
    bool degen_impossible = false;    // maxOrth + r < 10  =>  pure search is fully exhaustive
};

static void build_subset(const Col* c, const std::vector<int>& S, SubsetCtx& X) {
    X.r = (int)S.size();
    X.freeQ = S;
    X.fixedQ.clear();
    for (int j = 0; j < NQ; ++j)
        if (std::find(S.begin(), S.end(), j) == S.end()) X.fixedQ.push_back(j);
    const int nf = (int)X.fixedQ.size();
    uint32_t P[NQ][3];
    for (int j = 0; j < NQ; ++j) col_points(c[j], P[j]);

    X.C1.clear(); X.C2.clear(); X.C3.clear();
    for (int w = 0; w < 4; ++w) X.F0.f[w].clear();
    X.F0.f[0].set(0);
    X.C1.set(0); X.C2.set(0); X.C3.set(0);          // the empty (weight-0) fixed error
    X.feasible = true;
    for (int a = 0; a < nf; ++a) {
        int qa = X.fixedQ[a];
        for (int ia = 0; ia < 3; ++ia) {
            uint32_t s1 = P[qa][ia];
            X.C1.set(s1); X.C2.set(s1); X.C3.set(s1); X.F0.f[1].set(s1);
            if (s1 == 0) X.feasible = false;                       // weight-1 all-fixed
            for (int b = a + 1; b < nf; ++b) {
                int qb = X.fixedQ[b];
                for (int ib = 0; ib < 3; ++ib) {
                    uint32_t s2 = s1 ^ P[qb][ib];
                    X.C2.set(s2); X.C3.set(s2); X.F0.f[2].set(s2);
                    if (s2 == 0) X.feasible = false;               // weight-2 all-fixed
                    for (int d = b + 1; d < nf; ++d) {
                        int qd = X.fixedQ[d];
                        for (int id = 0; id < 3; ++id) {
                            uint32_t s3 = s2 ^ P[qd][id];
                            X.C3.set(s3); X.F0.f[3].set(s3);
                            if (s3 == 0) X.feasible = false;       // weight-3 all-fixed
                            for (int e = d + 1; e < nf; ++e) {
                                int qe = X.fixedQ[e];
                                for (int ie = 0; ie < 3; ++ie)
                                    if ((s3 ^ P[qe][ie]) == 0) X.feasible = false;  // weight-4
                            }
                        }
                    }
                }
            }
        }
    }
    // required Pluecker sum for the free columns
    Alt T;
    for (int j : X.fixedQ) { Alt p = plucker(c[j].x, c[j].z); alt_xor(T, p); }
    X.T = T; X.rankT = alt_rank(T);

    // allowed lines: all three points must avoid C_3
    std::vector<uint32_t> pts;
    for (uint32_t v = 1; v < VN; ++v) if (!X.C3.test(v)) pts.push_back(v);
    X.A.clear();
    for (size_t i = 0; i < pts.size(); ++i)
        for (size_t k = i + 1; k < pts.size(); ++k) {
            uint32_t a = pts[i], b = pts[k], d = a ^ b;
            if (d < b) continue;                       // count each line once: a<b<d
            if (X.C3.test(d)) continue;
            X.A.push_back(make_line(a, b));
        }
    std::sort(X.A.begin(), X.A.end(), [](const Line& u, const Line& v) { return u.key < v.key; });
    X.Akey.clear(); X.Akey.reserve(X.A.size());
    for (const Line& L : X.A) X.Akey.push_back(L.key);
    X.through.assign(VN, {});
    for (int i = 0; i < (int)X.A.size(); ++i)
        for (int t = 0; t < 3; ++t) X.through[X.A[i].p[t]].push_back(i);

    // degeneracy closure: a weight-<=4 stabilizer element needs a _|_ to >= 10 columns
    X.maxOrth = 0;
    for (uint32_t a = 1; a < VN; ++a) {
        int cnt = 0;
        for (int j : X.fixedQ) if (!par(a & c[j].x) && !par(a & c[j].z)) ++cnt;
        if (cnt > X.maxOrth) X.maxOrth = cnt;
    }
    X.degen_impossible = (X.maxOrth + X.r < NQ - MAXW);   // < 10
}
static int find_line(const SubsetCtx& X, uint32_t key) {
    auto it = std::lower_bound(X.Akey.begin(), X.Akey.end(), key);
    if (it == X.Akey.end() || *it != key) return -1;
    return (int)(it - X.Akey.begin());
}

// =========================================================================================
//  THE SEARCH
// =========================================================================================
struct Stats {
    uint64_t subsets = 0, hitting = 0, feasible = 0;
    uint64_t nodes = 0, prune_rank = 0, prune_cross = 0;
    uint64_t forced_rank = 0, forced_absent = 0, forced_order = 0, forced_cross = 0;
    uint64_t candidates = 0, verified = 0, solutions = 0, degenerate_open = 0;
    double seconds = 0;
    void add(const Stats& o) {
        subsets += o.subsets; hitting += o.hitting; feasible += o.feasible;
        nodes += o.nodes; prune_rank += o.prune_rank; prune_cross += o.prune_cross;
        forced_rank += o.forced_rank; forced_absent += o.forced_absent;
        forced_order += o.forced_order; forced_cross += o.forced_cross;
        candidates += o.candidates; verified += o.verified; solutions += o.solutions;
        degenerate_open += o.degenerate_open;
    }
};
// all four cross conditions that involve the NEW column L together with the k already chosen
static bool DFS_NOCROSS = false;   // selftest only: isolates the rank/Pluecker core
static inline bool cross_ok(const SubsetCtx& X, const Line* ch, int k, const Line& L) {
    if (DFS_NOCROSS) return true;
    for (int i = 0; i < k; ++i)
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
            if (X.C2.test(L.p[a] ^ ch[i].p[b])) return false;                    // f = 2
    for (int i = 0; i < k; ++i) for (int j = i + 1; j < k; ++j)
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) for (int d = 0; d < 3; ++d)
            if (X.C1.test(L.p[a] ^ ch[i].p[b] ^ ch[j].p[d])) return false;       // f = 3
    for (int i = 0; i < k; ++i) for (int j = i + 1; j < k; ++j) for (int m = j + 1; m < k; ++m)
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
            for (int d = 0; d < 3; ++d) for (int e = 0; e < 3; ++e)
                if ((L.p[a] ^ ch[i].p[b] ^ ch[j].p[d] ^ ch[m].p[e]) == 0) return false; // f = 4
    return true;
}
static std::string OUTDIR = ".";
static std::atomic<uint64_t> SOLCOUNT{0};

// full, independent acceptance test of one completed assignment
static bool accept_candidate(const Col* base, const SubsetCtx& X, const Line* ch,
                             Stats& st, const std::string& tag) {
    Col c[NQ];
    for (int j = 0; j < NQ; ++j) c[j] = base[j];
    for (int i = 0; i < X.r; ++i) { c[X.freeQ[i]].x = ch[i].p[0]; c[X.freeQ[i]].z = ch[i].p[1]; }
    Gen g[RG]; cols_to_gens(c, g);
    ++st.verified;
    CodeReport R = analyse_code(g);
    int64_t colv = column_violations(c);
    bool pure = R.residual.empty();
    bool dfive = R.residual_nondeg.empty();
    if (!dfive) return false;
    ++st.solutions;
    uint64_t id = ++SOLCOUNT;
    RB.clear();
    ap("SOLUTION %llu   (%s)\n\n", (unsigned long long)id, tag.c_str());
    ap("free qubits :");
    for (int q : X.freeQ) ap(" %d", q);
    ap("\n\ngenerators:\n");
    for (int i = 0; i < RG; ++i) ap("  %s   (weight %d)\n", gstr(g[i]).c_str(), R.gw[i]);
    ap("\nH = [ X | Z ]:\n");
    for (int i = 0; i < RG; ++i) {
        std::string xs, zs;
        for (int j = 0; j < NQ; ++j) { xs += char('0' + ((g[i].x >> j) & 1));
                                       zs += char('0' + ((g[i].z >> j) & 1)); }
        ap("  %s | %s\n", xs.c_str(), zs.c_str());
    }
    ap("\nVERIFICATION (independent: all errors rebuilt, syndromes by direct symplectic products)\n");
    ap("  pairwise commutation                 : %s\n", R.commute_ok ? "OK" : "FAIL");
    ap("  rank                                 : %d %s\n", R.rank, R.rank_ok ? "OK" : "FAIL");
    ap("  every column is a 2-dim subspace     : %s\n", R.cols_ok ? "OK" : "FAIL");
    ap("  Pluecker sum vanishes                : %s\n", R.plucker_ok ? "OK" : "FAIL");
    ap("  weight-<=4 errors examined           : %lld\n", (long long)R.total_err);
    ap("  zero-syndrome weight-<=4 errors      : %zu\n", R.residual.size());
    ap("  of those NOT in the stabilizer group : %zu\n", R.residual_nondeg.size());
    ap("  second, column-only distance routine : %lld violating tuples\n", (long long)colv);
    ap("  minimum non-zero stabilizer weight   : %d\n", R.stabmin);
    ap("  code is                              : %s\n",
       pure ? "PURE (non-degenerate), d >= 5" : "DEGENERATE, d >= 5");
    if (!R.residual.empty()) {
        ap("\n  degenerate zero-syndrome errors (all lie in the stabilizer group):\n");
        for (auto& e : R.residual) ap("    %s\n", estr(e.first, e.second).c_str());
    }
    char fn[256];
    snprintf(fn, sizeof(fn), "%s/results/solution_%04llu.txt", OUTDIR.c_str(),
             (unsigned long long)id);
    dump(fn);
    logline("\n*** SOLUTION %llu FOUND  (%s)  %s ***\n", (unsigned long long)id, tag.c_str(),
            pure ? "PURE d>=5" : "DEGENERATE d>=5");
    for (int i = 0; i < RG; ++i) logline("    %s\n", gstr(g[i]).c_str());
    logline("    written to %s\n\n", fn);
    return true;
}

// depth-first over the free columns; the last one is forced by the rank-2 condition
static void dfs(const Col* base, const SubsetCtx& X, int k, Alt T, int lastIdx,
                Line* ch, Stats& st) {
    ++st.nodes;
    int rk = alt_rank(T);
    if (rk > 2 * (X.r - k)) { ++st.prune_rank; return; }      // EXACT: sum of t decomposables
    if (k == X.r - 1) {                                        // last column is determined
        if (rk != 2) { ++st.forced_rank; return; }
        uint32_t u, v;
        if (!alt_support2(T, u, v)) { ++st.forced_rank; return; }
        Line L = make_line(u, v);
        int idx = find_line(X, L.key);
        if (idx < 0) { ++st.forced_absent; return; }
        if (idx <= lastIdx) { ++st.forced_order; return; }     // canonical: strictly increasing
        if (!cross_ok(X, ch, k, X.A[idx])) { ++st.forced_cross; return; }
        ch[k] = X.A[idx];
        ++st.candidates;
        accept_candidate(base, X, ch, st, "search");
        return;
    }
    const int n = (int)X.A.size();
    for (int idx = lastIdx + 1; idx < n; ++idx) {
        const Line& L = X.A[idx];
        if (!cross_ok(X, ch, k, L)) { ++st.prune_cross; continue; }
        ch[k] = L;
        Alt T2 = T; alt_xor(T2, L.pl);
        dfs(base, X, k + 1, T2, idx, ch, st);
    }
}

// ---------------------------------------------------------------------------------------
// The production search.  Same mathematics as dfs() above, but it carries
//   (a) the running forbidden-point table F, so testing a line costs three bit lookups, and
//   (b) the surviving candidate list, so a line eliminated at depth 2 is never re-tested at
//       depths 3,4,5 of that branch.
// dfs() is retained because the self-test uses it to validate this one.
static void dfs2(const Col* base, const SubsetCtx& X, int k, Alt T,
                 const std::vector<int>& cand, const FSets& F, Line* ch, Stats& st) {
    ++st.nodes;
    const int rem = X.r - k;
    int rk = alt_rank(T);
    if (rk > 2 * rem) { ++st.prune_rank; return; }
    Bits D; F.forbidden(D);

    if (rem == 2) {                       // last two columns: solved, never enumerated
        if (rk == 0) { ++st.forced_rank; return; }
        auto try_pair = [&](int idx1) {
            const Line& L1 = X.A[idx1];
            if (!line_ok(D, L1)) return;
            Alt T2 = T; alt_xor(T2, L1.pl);
            if (alt_rank(T2) != 2) return;
            uint32_t a2, b2;
            if (!alt_support2(T2, a2, b2)) return;
            int idx2 = find_line(X, make_line(a2, b2).key);
            if (idx2 < 0) { ++st.forced_absent; return; }
            if (idx2 <= idx1) return;                    // each unordered pair once
            FSets F2; F.extend(L1.p, F2);
            Bits D2; F2.forbidden(D2);
            if (!line_ok(D2, X.A[idx2])) { ++st.forced_cross; return; }
            ch[k] = L1; ch[k + 1] = X.A[idx2];
            ++st.candidates;
            accept_candidate(base, X, ch, st, "search");
        };
        if (rk == 2) {
            uint32_t u, v;
            if (!alt_support2(T, u, v)) { ++st.forced_rank; return; }
            uint32_t pts[3] = { u, v, u ^ v };
            for (uint32_t p : pts)
                for (int idx1 : X.through[p])
                    if (std::binary_search(cand.begin(), cand.end(), idx1)) try_pair(idx1);
        } else {                                          // rk == 4
            uint32_t bas[RG]; int dim = alt_support(T, bas);
            if (dim != 4) { ++st.forced_rank; return; }
            for (int m = 1; m < 16; ++m) {
                uint32_t p1 = 0;
                for (int i = 0; i < 4; ++i) if ((m >> i) & 1) p1 ^= bas[i];
                for (int m2 = m + 1; m2 < 16; ++m2) {
                    uint32_t p2 = 0;
                    for (int i = 0; i < 4; ++i) if ((m2 >> i) & 1) p2 ^= bas[i];
                    if (p2 == p1) continue;
                    Line L = make_line(p1, p2);
                    if ((L.p[0] != p1 && L.p[0] != p2) || (L.p[1] != p1 && L.p[1] != p2)) continue;
                    int idx1 = find_line(X, L.key);
                    if (idx1 >= 0 && std::binary_search(cand.begin(), cand.end(), idx1))
                        try_pair(idx1);
                }
            }
        }
        return;
    }
    std::vector<int> next;
    next.reserve(cand.size());
    for (size_t pos = 0; pos < cand.size(); ++pos) {
        int idx = cand[pos];
        const Line& L = X.A[idx];
        Alt T2 = T; alt_xor(T2, L.pl);
        if (alt_rank(T2) > 2 * (rem - 1)) { ++st.prune_rank; continue; }
        FSets F2; F.extend(L.p, F2);
        Bits D2; F2.forbidden(D2);
        next.clear();
        for (size_t q = pos + 1; q < cand.size(); ++q)
            if (line_ok(D2, X.A[cand[q]])) next.push_back(cand[q]);
            else ++st.prune_cross;
        if ((int)next.size() < rem - 1) { ++st.prune_cross; continue; }
        ch[k] = L;
        dfs2(base, X, k + 1, T2, next, F2, ch, st);
    }
}
static void run_subset(const Col* base, const SubsetCtx& X, Line* ch, Stats& st) {
    Bits D; X.F0.forbidden(D);
    std::vector<int> cand;
    for (int i = 0; i < (int)X.A.size(); ++i) if (line_ok(D, X.A[i])) cand.push_back(i);
    dfs2(base, X, 0, X.T, cand, X.F0, ch, st);
}

// =========================================================================================
//  SUBSET ENUMERATION
// =========================================================================================
static std::vector<std::pair<uint32_t,uint32_t>> RESIDUALS;   // the starting code's residuals
static bool hits_all_residuals(const std::vector<int>& S) {
    uint32_t m = 0; for (int j : S) m |= 1u << j;
    for (auto& e : RESIDUALS) if (((e.first | e.second) & m) == 0) return false;
    return true;
}
template <class F> static void for_each_subset(int r, F&& f) {
    std::vector<int> S(r);
    for (int i = 0; i < r; ++i) S[i] = i;
    for (;;) {
        f(S);
        int i = r - 1;
        while (i >= 0 && S[i] == NQ - r + i) --i;
        if (i < 0) break;
        ++S[i];
        for (int k = i + 1; k < r; ++k) S[k] = S[k-1] + 1;
    }
}

// =========================================================================================
//  VERIFICATION OF THE STARTING CODE AND THE SELF-TEST SUITE
// =========================================================================================
static void print_report(const CodeReport& R, const char* what) {
    logline("---- %s ----\n", what);
    logline("  generator weights            :");
    for (int i = 0; i < RG; ++i) logline(" %d", R.gw[i]);
    logline("\n  all %d generator pairs commute: %s\n", RG * (RG - 1) / 2,
            R.commute_ok ? "OK" : "FAIL");
    logline("  rank of the check matrix     : %d %s\n", R.rank, R.rank_ok ? "OK" : "FAIL");
    logline("  every column spans a 2-dim W : %s\n", R.cols_ok ? "OK" : "FAIL");
    logline("  Pluecker sum over columns = 0: %s\n", R.plucker_ok ? "OK" : "FAIL");
    logline("  weight-<=4 errors enumerated : %lld\n", (long long)R.total_err);
    logline("  undetected (zero syndrome)   : %zu   [w1=%lld w2=%lld w3=%lld w4=%lld]\n",
            R.residual.size(), (long long)R.by_w[1], (long long)R.by_w[2],
            (long long)R.by_w[3], (long long)R.by_w[4]);
    logline("  of those outside stabilizer  : %zu  (these are the real distance failures)\n",
            R.residual_nondeg.size());
    logline("  min non-zero stabilizer weight: %d\n", R.stabmin);
    for (auto& e : R.residual)
        logline("      %s = %s%s\n", gstr(Gen{e.first, e.second}).c_str(),
                estr(e.first, e.second).c_str(),
                std::find(R.residual_nondeg.begin(), R.residual_nondeg.end(), e)
                    == R.residual_nondeg.end() ? "   (in the stabilizer: harmless)" : "");
}
static bool verify_start(bool verbose) {
    CodeReport R = analyse_code(G0);
    if (verbose) print_report(R, "STARTING 11-GENERATOR CODE");
    bool ok = R.commute_ok && R.rank_ok && R.cols_ok && R.plucker_ok;
    int64_t colv = column_violations(COL0);
    if (verbose)
        logline("  independent column-only distance routine: %lld violating tuples "
                "(must equal %zu) %s\n", (long long)colv, R.residual.size(),
                colv == (int64_t)R.residual.size() ? "OK" : "MISMATCH");
    if (colv != (int64_t)R.residual.size()) ok = false;
    RESIDUALS = R.residual;
    return ok;
}
struct SubsetInfo {
    std::vector<int> S;
    int nA = 0, rankT = 0, maxOrth = 0;
    bool degen_impossible = false;
};
static void collect_subsets(int r, std::vector<SubsetInfo>& out, Stats& st, bool build);
static bool selftest(int nrand) {
    bool pass = true;
    logline("=============================== SELF TEST ===============================\n");
    if (!verify_start(true)) { logline("  starting code FAILED verification\n"); pass = false; }

    // the five residuals the user reported, checked literally
    static const char* WANT[5] = {
        "IIIXIIIIIIIIZZ", "IIIIIIXZIIIXII",
        "YIIIIZIIIYIXII", "IXIIIIZYXIIIII", "IIZIXIIIIIXXII" };
    std::set<std::string> got, want;
    for (auto& e : RESIDUALS) got.insert(gstr(Gen{e.first, e.second}));
    for (int i = 0; i < 5; ++i) want.insert(WANT[i]);
    logline("  residual set matches the five reported errors : %s (%zu found)\n",
            got == want ? "OK" : "MISMATCH", got.size());
    if (got != want) { pass = false;
        for (auto& s : got) if (!want.count(s)) logline("      extra:   %s\n", s.c_str());
        for (auto& s : want) if (!got.count(s)) logline("      missing: %s\n", s.c_str()); }

    // Pluecker / rank algebra
    std::mt19937_64 rng(20260831);
    int bad = 0;
    for (int t = 0; t < nrand; ++t) {
        uint32_t a = 1 + (uint32_t)(rng() % (VN - 1)), b = 1 + (uint32_t)(rng() % (VN - 1));
        if (a == b) continue;
        Line L = make_line(a, b);
        if (alt_rank(L.pl) != 2) { ++bad; continue; }
        uint32_t u, v;
        if (!alt_support2(L.pl, u, v)) { ++bad; continue; }
        Line L2 = make_line(u, v);
        if (L2.key != L.key) ++bad;
    }
    logline("  rank(pl(line)) == 2 and support recovers the line : %d failures  %s\n",
            bad, bad == 0 ? "OK" : "FAIL");
    if (bad) pass = false;
    // a sum of t decomposables has rank <= 2t   (the basis of the exact prune)
    bad = 0;
    for (int t = 0; t < nrand; ++t) {
        int m = 1 + (int)(rng() % 5);
        Alt S;
        for (int i = 0; i < m; ++i) {
            uint32_t a = 1 + (uint32_t)(rng() % (VN - 1)), b = 1 + (uint32_t)(rng() % (VN - 1));
            if (a == b) b ^= 1;
            Line L = make_line(a, b); alt_xor(S, L.pl);
        }
        if (alt_rank(S) > 2 * m) ++bad;
    }
    logline("  rank(sum of t decomposables) <= 2t                : %d failures  %s\n",
            bad, bad == 0 ? "OK" : "FAIL");
    if (bad) pass = false;

    // deliberate breakage: duplicating a column must create a weight-2 undetected error
    { Col c[NQ]; memcpy(c, COL0, sizeof(c)); c[0] = c[1];
      Gen g[RG]; cols_to_gens(c, g);
      CodeReport R2 = analyse_code(g);
      bool caught = R2.by_w[2] > 0;
      logline("  deliberate break (column 0 := column 1) caught    : %s (w2 residuals = %lld)\n",
              caught ? "OK" : "NOT CAUGHT", (long long)R2.by_w[2]);
      if (!caught) pass = false; }
    // deliberate breakage: make a column degenerate (dim W = 1)
    { Col c[NQ]; memcpy(c, COL0, sizeof(c)); c[5].z = c[5].x;
      Gen g[RG]; cols_to_gens(c, g);
      CodeReport R2 = analyse_code(g);
      bool caught = (R2.by_w[1] > 0) || !R2.cols_ok;
      logline("  deliberate break (dim W_5 = 1) caught             : %s (w1 residuals = %lld)\n",
              caught ? "OK" : "NOT CAUGHT", (long long)R2.by_w[1]);
      if (!caught) pass = false; }
    // deliberate breakage: a column change that breaks commutation must be visible
    { Col c[NQ]; memcpy(c, COL0, sizeof(c)); c[3].x ^= 1u;
      Gen g[RG]; cols_to_gens(c, g);
      CodeReport R2 = analyse_code(g);
      logline("  deliberate break (commutation) caught            : %s\n",
              !R2.commute_ok || !R2.plucker_ok ? "OK" : "NOT CAUGHT");
      if (R2.commute_ok && R2.plucker_ok) pass = false; }

    // KNOWN CONFIGURATION: the original free columns must satisfy the Pluecker equation
    // exactly, and must be rejected by the allowed-set filter (they carry the residuals).
    {
        std::vector<int> S = {3, 6, 11};
        SubsetCtx X; build_subset(COL0, S, X);
        Alt sum;
        for (int j : S) { Alt p = plucker(COL0[j].x, COL0[j].z); alt_xor(sum, p); }
        alt_xor(sum, X.T);
        logline("  original free columns reproduce the target T     : %s\n",
                alt_zero(sum) ? "OK" : "FAIL");
        if (!alt_zero(sum)) pass = false;
        int inA = 0;
        for (int j : S) {
            Line L = make_line(COL0[j].x, COL0[j].z);
            if (find_line(X, L.key) >= 0) ++inA;
        }
        logline("  original free columns excluded by the C_3 filter : %d of 3 still allowed %s\n",
                inA, inA < 3 ? "OK" : "SUSPICIOUS");
    }
    // END-TO-END: put the original lines back into A and confirm the search machinery
    // reaches exactly that assignment and then correctly REJECTS it (5 residuals).
    {
        std::vector<int> S = {3, 6, 11};
        SubsetCtx X; build_subset(COL0, S, X);
        for (int j : S) X.A.push_back(make_line(COL0[j].x, COL0[j].z));
        std::sort(X.A.begin(), X.A.end(),
                  [](const Line& u, const Line& v) { return u.key < v.key; });
        X.Akey.clear(); for (const Line& L : X.A) X.Akey.push_back(L.key);
        // walk the original assignment by hand through the same acceptance routine
        Line ch[3];
        for (int i = 0; i < 3; ++i) ch[i] = make_line(COL0[S[i]].x, COL0[S[i]].z);
        Stats st;
        bool acc = accept_candidate(COL0, X, ch, st, "selftest-known-config");
        logline("  known (original) configuration is rebuilt and REJECTED : %s\n",
                !acc ? "OK" : "WRONGLY ACCEPTED");
        if (acc) pass = false;
    }
    // THE CRITICAL TEST.  A depth-first search that wrongly prunes everything would report
    // "no solution" for every r and look exactly like a proof.  So: plant the original three
    // columns back into A, switch off ONLY the distance cross-conditions (which the original
    // columns must fail, since they carry the residuals), and check that the rank prune, the
    // forced last column, the line lookup and the increasing-index rule together RECOVER that
    // known assignment.  The search core is validated independently of the filters.
    for (std::vector<int> S : { std::vector<int>{3,6,11}, std::vector<int>{1,3,11},
                                std::vector<int>{7,11,12} }) {
        SubsetCtx X; build_subset(COL0, S, X);
        for (int j : S) X.A.push_back(make_line(COL0[j].x, COL0[j].z));
        std::sort(X.A.begin(), X.A.end(),
                  [](const Line& u, const Line& v) { return u.key < v.key; });
        X.A.erase(std::unique(X.A.begin(), X.A.end(),
                  [](const Line& u, const Line& v) { return u.key == v.key; }), X.A.end());
        X.Akey.clear(); for (const Line& L : X.A) X.Akey.push_back(L.key);
        Line ch[8]; Stats st;
        DFS_NOCROSS = true;
        dfs(COL0, X, 0, X.T, -1, ch, st);          // reference engine
        DFS_NOCROSS = false;
        logline("  DFS core recovers the planted assignment on {%d,%d,%d} : %s "
                "(%llu candidates, %llu nodes)\n", S[0], S[1], S[2],
                st.candidates >= 1 ? "OK" : "FAIL -- SEARCH CORE IS BROKEN",
                (unsigned long long)st.candidates, (unsigned long long)st.nodes);
        if (st.candidates < 1) pass = false;
    }
    // CROSS-VALIDATION OF THE TWO SEARCH ENGINES.  dfs() enumerates to depth r-1 and forces
    // only the last column, testing every pair/triple/quadruple condition explicitly.
    // dfs2() carries a running forbidden-point table and solves the last TWO columns
    // algebraically.  They are independent implementations of the same predicate, so on
    // every subset they must agree on the number of candidates and of accepted codes.
    {
        int mism = 0, tested = 0;
        for (int r = 3; r <= 5; ++r) {
            Stats stc; std::vector<SubsetInfo> subs;
            collect_subsets(r, subs, stc, true);
            int step = std::max<int>(1, (int)subs.size() / 6);
            for (size_t i = 0; i < subs.size(); i += step) {
                SubsetCtx X; build_subset(COL0, subs[i].S, X);
                Line ch1[8], ch2[8]; Stats s1, s2;
                dfs(COL0, X, 0, X.T, -1, ch1, s1);
                run_subset(COL0, X, ch2, s2);
                ++tested;
                if (s1.candidates != s2.candidates || s1.solutions != s2.solutions) {
                    ++mism;
                    logline("      [MISMATCH] r=%d subset %zu : reference %llu/%llu vs "
                            "optimized %llu/%llu (candidates/solutions)\n", r, i,
                            (unsigned long long)s1.candidates, (unsigned long long)s1.solutions,
                            (unsigned long long)s2.candidates, (unsigned long long)s2.solutions);
                }
            }
        }
        logline("  reference engine vs optimized engine, %d subsets : %d mismatches  %s\n",
                tested, mism, mism == 0 ? "OK" : "FAIL");
        if (mism) pass = false;
    }
    logline("========================= SELF TEST %s =========================\n\n",
            pass ? "PASSED" : "FAILED");
    return pass;
}

// =========================================================================================
//  ANALYSIS / CALIBRATION AND THE SEARCH DRIVER
// =========================================================================================
static void collect_subsets(int r, std::vector<SubsetInfo>& out, Stats& st, bool build) {
    out.clear();
    for_each_subset(r, [&](const std::vector<int>& S) {
        ++st.subsets;
        if (!hits_all_residuals(S)) return;
        ++st.hitting;
        SubsetCtx X; build_subset(COL0, S, X);
        if (!X.feasible) return;
        ++st.feasible;
        SubsetInfo I; I.S = S; I.nA = (int)X.A.size(); I.rankT = X.rankT;
        I.maxOrth = X.maxOrth; I.degen_impossible = X.degen_impossible;
        if (build) out.push_back(I);
    });
}
static void analyze(int rmin, int rmax) {
    logline("=========================== SEARCH-SPACE ANALYSIS ===========================\n");
    logline("residual supports of the starting code:\n");
    for (auto& e : RESIDUALS)
        logline("   %-16s support {", estr(e.first, e.second).c_str());
    logline("\n");
    for (auto& e : RESIDUALS) {
        uint32_t m = e.first | e.second;
        std::string s;
        for (int j = 0; j < NQ; ++j) if ((m >> j) & 1) s += std::to_string(j) + " ";
        logline("   %-22s -> { %s}\n", estr(e.first, e.second).c_str(), s.c_str());
    }
    logline("\n r   C(14,r)  hitting  feasible   |A| min/med/max   rank(T) min/max  "
            "degeneracy ruled out\n");
    logline(" --------------------------------------------------------------------------"
            "----------------\n");
    for (int r = rmin; r <= rmax; ++r) {
        Stats st; std::vector<SubsetInfo> subs;
        double t0 = elapsed();
        collect_subsets(r, subs, st, true);
        double t1 = elapsed();
        if (subs.empty()) {
            logline(" %2d  %7llu  %7llu  %8llu   (none feasible)\n", r,
                    (unsigned long long)st.subsets, (unsigned long long)st.hitting,
                    (unsigned long long)st.feasible);
            continue;
        }
        std::vector<int> nas, rks; int dimp = 0;
        for (auto& I : subs) { nas.push_back(I.nA); rks.push_back(I.rankT);
                               if (I.degen_impossible) ++dimp; }
        std::sort(nas.begin(), nas.end()); std::sort(rks.begin(), rks.end());
        logline(" %2d  %7llu  %7llu  %8llu   %6d/%6d/%6d   %4d / %-4d      %d of %zu"
                "   (setup %.1f s)\n",
                r, (unsigned long long)st.subsets, (unsigned long long)st.hitting,
                (unsigned long long)st.feasible,
                nas.front(), nas[nas.size()/2], nas.back(), rks.front(), rks.back(),
                dimp, subs.size(), t1 - t0);
    }
    logline("\n|A| is the number of admissible lines for EVERY free column (the f=1 condition).\n");
    logline("rank(T) is the rank of the required Pluecker sum; a sum of r decomposables has\n");
    logline("rank <= 2r, so rank(T) = 2r forces every DFS step to drop the rank by 2.\n");
    logline("'degeneracy ruled out' counts subsets where no weight-<=4 stabilizer element can\n");
    logline("exist at all, so for those the pure search is exhaustive over ALL codes.\n");
}
// time a few subsets to project the cost of a full run at this r
static void calibrate(int r, int sample, int threads) {
    Stats st0; std::vector<SubsetInfo> subs;
    collect_subsets(r, subs, st0, true);
    logline("---- calibration, r = %d : %zu feasible subsets ----\n", r, subs.size());
    if (subs.empty()) { logline("   nothing to do\n"); return; }
    int n = std::min((int)subs.size(), sample);
    double t0 = elapsed();
    Stats acc;
    // sample EVENLY across the subset list: |A| varies by a factor of several, so the first
    // n subsets are not representative and would give a misleading projection
    for (int t = 0; t < n; ++t) {
        int i = (int)((size_t)t * subs.size() / n);
        SubsetCtx X; build_subset(COL0, subs[i].S, X);
        Line ch[8]; Stats st;
        double a = elapsed();
        run_subset(COL0, X, ch, st);
        double b = elapsed();
        acc.add(st);
        logline("   subset %4d {", i);
        for (int q : subs[i].S) logline(" %d", q);
        logline(" }  |C3| = %4d  free pts = %4d  |A| = %6zu  rank(T) = %2d  "
                "nodes %-12llu  %.3f s\n",
                X.C3.count(), VN - X.C3.count(), X.A.size(), X.rankT,
                (unsigned long long)st.nodes, b - a);
    }
    double t1 = elapsed();
    double per = (t1 - t0) / n;
    logline("   sampled %d of %zu subsets in %.2f s  (%.4f s each)\n",
            n, subs.size(), t1 - t0, per);
    logline("   nodes %llu, rank-pruned %llu, cross-pruned %llu, candidates %llu, "
            "solutions %llu\n", (unsigned long long)acc.nodes,
            (unsigned long long)acc.prune_rank, (unsigned long long)acc.prune_cross,
            (unsigned long long)acc.candidates, (unsigned long long)acc.solutions);
    logline("   PROJECTED FULL r=%d RUN : %s on 1 thread, %s on %d threads\n\n", r,
            hms(per * subs.size()).c_str(), hms(per * subs.size() / threads).c_str(), threads);
}

// =========================================================================================
//  IMPROVEMENT SCAN
//
//  The main search only accepts d >= 5, so it is blind to a replacement that merely reduces
//  the residual count (5 -> 4 -> ...).  This mode drops the distance filter entirely and
//  enumerates EVERY valid commuting rank-11 replacement of the r columns, reporting the
//  residual count of each.  It is exhaustive because of one exact fact:
//
//     if rank(T_k) = 2 (r-k), then supp(T_k) has dimension 2(r-k) and every remaining line
//     must lie INSIDE supp(T_k).
//
//  (T_k is a sum of r-k decomposables, so supp(T_k) is contained in their span, which has
//  dimension at most 2(r-k); equality forces the span to BE supp(T_k).)  So the candidate
//  set at each level is the set of lines of a known 2(r-k)-dimensional space -- 651 lines
//  when that is 6-dimensional, 35 when it is 4-dimensional.  When rank(T_k) < 2(r-k) there
//  is slack and this shortcut is not available; such branches are counted and reported as
//  NOT scanned, so the exhaustiveness claim stays honest.
// =========================================================================================
static void lines_in_subspace(const uint32_t* bas, int dim, std::vector<Line>& out) {
    out.clear();
    const int n = 1 << dim;
    std::vector<uint32_t> pts;
    for (int m = 1; m < n; ++m) {
        uint32_t p = 0;
        for (int i = 0; i < dim; ++i) if ((m >> i) & 1) p ^= bas[i];
        pts.push_back(p);
    }
    for (size_t i = 0; i < pts.size(); ++i)
        for (size_t k = i + 1; k < pts.size(); ++k) {
            uint32_t a = pts[i], b = pts[k], d = a ^ b;
            if (a > b || d < b) continue;              // keep a < b < a^b : each line once
            out.push_back(make_line(a, b));
        }
}
struct ScanState { int64_t best = 1 << 30; uint64_t leaves = 0, skipped = 0, invalid = 0; };
static void scan_dfs(const Col* base, const SubsetCtx& X, int k, Alt T, Line* ch,
                     ScanState& S) {
    const int rem = X.r - k;
    int rk = alt_rank(T);
    if (rk > 2 * rem) return;
    if (rk < 2 * rem) { ++S.skipped; return; }         // slack: shortcut not valid here
    if (rem == 1) {
        uint32_t u, v;
        if (!alt_support2(T, u, v)) return;
        ch[k] = make_line(u, v);
        ++S.leaves;
        Col c[NQ];
        for (int j = 0; j < NQ; ++j) c[j] = base[j];
        for (int i = 0; i < X.r; ++i) { c[X.freeQ[i]].x = ch[i].p[0]; c[X.freeQ[i]].z = ch[i].p[1]; }
        for (int j = 0; j < NQ; ++j) if (!col_ok(c[j])) { ++S.invalid; return; }
        Gen g[RG]; cols_to_gens(c, g);
        if (gen_rank(g) != RG) { ++S.invalid; return; }
        CodeReport R = analyse_code(g);
        int64_t nres = (int64_t)R.residual_nondeg.size();
        if (nres < S.best) {
            S.best = nres;
            logline("   improvement scan: residual count %lld  (free qubits", (long long)nres);
            for (int q : X.freeQ) logline(" %d", q);
            logline(")\n");
            if (nres < 5) {
                for (int i = 0; i < RG; ++i) logline("      %s\n", gstr(g[i]).c_str());
                RB.clear();
                ap("IMPROVED CODE: %lld undetected weight-<=4 errors (was 5)\n\n", (long long)nres);
                for (int i = 0; i < RG; ++i)
                    ap("  %s   (weight %d)\n", gstr(g[i]).c_str(), R.gw[i]);
                ap("\nfree qubits :");
                for (int q : X.freeQ) ap(" %d", q);
                ap("\nremaining undetected errors:\n");
                for (auto& e : R.residual_nondeg) ap("  %s\n", estr(e.first, e.second).c_str());
                ap("\nrank %d, all pairs commute %s, min stabilizer weight %d\n",
                   R.rank, R.commute_ok ? "yes" : "NO", R.stabmin);
                char fn[256];
                snprintf(fn, sizeof(fn), "%s/results/improvement_%lld_%d.txt",
                         OUTDIR.c_str(), (long long)nres, (int)S.leaves);
                dump(fn);
            }
        }
        return;
    }
    uint32_t bas[RG]; int dim = alt_support(T, bas);
    std::vector<Line> cand; lines_in_subspace(bas, dim, cand);
    for (const Line& L : cand) {
        Alt T2 = T; alt_xor(T2, L.pl);
        if (alt_rank(T2) != 2 * (rem - 1)) continue;
        ch[k] = L;
        scan_dfs(base, X, k + 1, T2, ch, S);
    }
}
static void scan_improvements(int r) {
    Stats stc; std::vector<SubsetInfo> subs;
    collect_subsets(r, subs, stc, true);
    logline("=========== IMPROVEMENT SCAN, r = %d : %zu feasible subsets ===========\n",
            r, subs.size());
    ScanState S;
    double t0 = elapsed();
    for (size_t i = 0; i < subs.size(); ++i) {
        SubsetCtx X; build_subset(COL0, subs[i].S, X);
        Line ch[8];
        scan_dfs(COL0, X, 0, X.T, ch, S);
    }
    logline("  complete replacements examined : %llu\n", (unsigned long long)S.leaves);
    logline("  rejected (bad column or rank)  : %llu\n", (unsigned long long)S.invalid);
    logline("  branches with rank slack (NOT scanned) : %llu\n", (unsigned long long)S.skipped);
    logline("  BEST residual count found      : %lld   (the starting code has 5)\n",
            (long long)S.best);
    logline("  runtime                        : %s\n", hms(elapsed() - t0).c_str());
    if (S.skipped == 0)
        logline("  STATUS: EXHAUSTIVE over all valid %d-column replacements.\n", r);
    else
        logline("  STATUS: complete except for %llu slack branches.\n",
                (unsigned long long)S.skipped);
    logline("=====================================================================\n\n");
}

// ---------------------------------------------------------------- checkpoint
static void save_checkpoint(int r, size_t done, const Stats& st) {
    RB.clear();
    ap("# hitset14 checkpoint\nversion 1\ntimestamp %s\nelapsed %.2f\n",
       now_stamp().c_str(), elapsed());
    ap("r %d\nsubsets_done %zu\n", r, done);
    ap("subsets %llu\nhitting %llu\nfeasible %llu\n", (unsigned long long)st.subsets,
       (unsigned long long)st.hitting, (unsigned long long)st.feasible);
    ap("nodes %llu\nprune_rank %llu\nprune_cross %llu\n", (unsigned long long)st.nodes,
       (unsigned long long)st.prune_rank, (unsigned long long)st.prune_cross);
    ap("candidates %llu\nverified %llu\nsolutions %llu\n", (unsigned long long)st.candidates,
       (unsigned long long)st.verified, (unsigned long long)st.solutions);
    char fn[256]; snprintf(fn, sizeof(fn), "%s/checkpoint_r%d.txt", OUTDIR.c_str(), r);
    dump(fn);
}
static size_t load_checkpoint(int r, Stats& st) {
    char fn[256]; snprintf(fn, sizeof(fn), "%s/checkpoint_r%d.txt", OUTDIR.c_str(), r);
    FILE* f = fopen(fn, "r"); if (!f) return 0;
    char line[512]; size_t done = 0;
    while (fgets(line, sizeof(line), f)) {
        if (!strncmp(line, "subsets_done ", 13)) done = strtoull(line + 13, nullptr, 10);
        else if (!strncmp(line, "nodes ", 6)) st.nodes = strtoull(line + 6, nullptr, 10);
        else if (!strncmp(line, "prune_rank ", 11)) st.prune_rank = strtoull(line + 11, nullptr, 10);
        else if (!strncmp(line, "prune_cross ", 12)) st.prune_cross = strtoull(line + 12, nullptr, 10);
        else if (!strncmp(line, "candidates ", 11)) st.candidates = strtoull(line + 11, nullptr, 10);
        else if (!strncmp(line, "verified ", 9)) st.verified = strtoull(line + 9, nullptr, 10);
        else if (!strncmp(line, "solutions ", 10)) st.solutions = strtoull(line + 10, nullptr, 10);
    }
    fclose(f);
    return done;
}

static void run_search(int r, int threads, bool resume, double reportEvery) {
    Stats st;
    std::vector<SubsetInfo> subs;
    Stats stc; collect_subsets(r, subs, stc, true);
    st.subsets = stc.subsets; st.hitting = stc.hitting; st.feasible = stc.feasible;
    size_t start = resume ? load_checkpoint(r, st) : 0;
    logline("=========================================================================\n");
    logline("SEARCH  r = %d   free columns per subset\n", r);
    logline("  C(14,%d) = %llu subsets, %llu hit all five residual supports, %llu feasible\n",
            r, (unsigned long long)stc.subsets, (unsigned long long)stc.hitting,
            (unsigned long long)stc.feasible);
    logline("  starting at subset %zu of %zu, %d threads\n", start, subs.size(), threads);
    logline("=========================================================================\n");
    if (subs.empty()) {
        logline("No feasible subset at r = %d.  PROVEN: no d>=5 code can be reached by\n", r);
        logline("replacing any %d columns of this stabilizer.\n\n", r);
        return;
    }
    std::atomic<size_t> next{start};
    std::vector<char> done(subs.size(), 0);
    for (size_t i = 0; i < start; ++i) done[i] = 1;
    std::mutex mtx;
    size_t completed = start;
    double lastReport = elapsed();
    double t0 = elapsed();
    int allDegenClosed = 1;
    auto worker = [&]() {
        for (;;) {
            size_t i = next.fetch_add(1);
            if (i >= subs.size()) break;
            SubsetCtx X; build_subset(COL0, subs[i].S, X);
            Line ch[8]; Stats local;
            run_subset(COL0, X, ch, local);
            std::lock_guard<std::mutex> lk(mtx);
            st.add(local);
            if (!X.degen_impossible) { ++st.degenerate_open; allDegenClosed = 0; }
            done[i] = 1;
            while (completed < subs.size() && done[completed]) ++completed;
            double now = elapsed();
            if (now - lastReport > reportEvery) {
                lastReport = now;
                double frac = double(completed - start) / double(subs.size() - start + 1e-9);
                double rem = frac > 1e-9 ? (now - t0) * (1.0 / frac - 1.0) : 0;
                logline("\n---- STATUS %s ----\n", now_stamp().c_str());
                logline("  r = %d   subsets %zu / %zu = %.2f %%\n", r, completed, subs.size(),
                        100.0 * completed / subs.size());
                logline("  elapsed %s   estimated remaining %s\n",
                        hms(now - t0).c_str(), hms(rem).c_str());
                logline("  DFS nodes %llu   rank-pruned %llu   cross-pruned %llu\n",
                        (unsigned long long)st.nodes, (unsigned long long)st.prune_rank,
                        (unsigned long long)st.prune_cross);
                logline("  candidates %llu   fully verified %llu   SOLUTIONS %llu\n",
                        (unsigned long long)st.candidates, (unsigned long long)st.verified,
                        (unsigned long long)st.solutions);
                logline("--------------------------------------------------\n\n");
                save_checkpoint(r, completed, st);
            }
        }
    };
    std::vector<std::thread> th;
    for (int i = 0; i < std::max(1, threads); ++i) th.emplace_back(worker);
    for (auto& t : th) t.join();
    st.seconds = elapsed() - t0;
    save_checkpoint(r, subs.size(), st);

    logline("\n=================== RESULT FOR r = %d ===================\n", r);
    logline("  column subsets examined                 : %llu\n", (unsigned long long)st.subsets);
    logline("  eliminated by the residual-hitting test : %llu\n",
            (unsigned long long)(st.subsets - st.hitting));
    logline("  eliminated by full fixed-part feasibility: %llu\n",
            (unsigned long long)(st.hitting - st.feasible));
    logline("  subsets actually searched               : %zu\n", subs.size());
    logline("  DFS nodes (candidate column assignments): %llu\n", (unsigned long long)st.nodes);
    logline("  eliminated by the exact rank/commutation bound : %llu\n",
            (unsigned long long)st.prune_rank);
    logline("  eliminated by the cross distance conditions    : %llu\n",
            (unsigned long long)st.prune_cross);
    logline("  last column forced but wrong rank       : %llu\n", (unsigned long long)st.forced_rank);
    logline("  last column forced but not admissible   : %llu\n", (unsigned long long)st.forced_absent);
    logline("  reached full distance verification      : %llu\n", (unsigned long long)st.verified);
    logline("  codes with d >= 5                       : %llu\n", (unsigned long long)st.solutions);
    logline("  runtime                                 : %s\n", hms(st.seconds).c_str());
    logline("  subsets where a degenerate solution was not ruled out : %llu\n",
            (unsigned long long)st.degenerate_open);
    if (st.solutions == 0) {
        logline("\n  STATUS: PROVEN / EXHAUSTIVE for r = %d.\n", r);
        logline("  No replacement of any %d columns of this stabilizer yields a code with no\n", r);
        logline("  weight-<=4 undetected error%s.\n",
                allDegenClosed ? " (degenerate solutions included: ruled out algebraically)"
                               : " that is PURE; see the degeneracy note above");
    } else {
        logline("\n  STATUS: %llu solution(s) found and independently verified; see results/.\n",
                (unsigned long long)st.solutions);
    }
    logline("========================================================\n\n");
}

// =========================================================================================
static void usage() {
    printf(
    "hitset14 -- residual-targeted column replacement search for a [[14,3,5]] code\n\n"
    "  --verify            verify the starting code and its five residuals, then exit\n"
    "  --selftest [N]      full correctness suite (N random algebra trials, default 2000)\n"
    "  --analyze [a b]     search-space analysis for r = a..b (default 3..7), then exit\n"
    "  --calibrate R [n]   time n subsets at r = R and project the full runtime\n"
    "  --search R          run the exhaustive search for r = R\n"
    "  --threads N         default: all hardware threads\n"
    "  --resume            continue from checkpoint_rR.txt\n"
    "  --report-every S    status interval in seconds (default 900)\n"
    "  --out DIR           output directory (default: the current directory)\n"
    "  --help\n");
}
int main(int argc, char** argv) {
    set_start_time();
    int threads = (int)std::thread::hardware_concurrency(); if (threads <= 0) threads = 4;
    int mode = 0, rA = 3, rB = 7, searchR = 0, calR = 0, calN = 4, nrand = 2000;
    bool resume = false;
    double reportEvery = 900.0;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto has = [&](int d = 1) { return i + d < argc && argv[i+d][0] != '-'; };
        if (a == "--help") { usage(); return 0; }
        else if (a == "--verify") mode = 1;
        else if (a == "--selftest") { mode = 2; if (has()) nrand = atoi(argv[++i]); }
        else if (a == "--analyze") { mode = 3;
            if (has() && has(2)) { rA = atoi(argv[i+1]); rB = atoi(argv[i+2]); i += 2; } }
        else if (a == "--calibrate") { mode = 4;
            if (has()) calR = atoi(argv[++i]);
            if (has()) calN = atoi(argv[++i]); }
        else if (a == "--search") { mode = 5; if (has()) searchR = atoi(argv[++i]); }
        else if (a == "--scan") { mode = 6; if (has()) searchR = atoi(argv[++i]); }
        else if (a == "--threads") { if (has()) threads = std::max(1, atoi(argv[++i])); }
        else if (a == "--resume") resume = true;
        else if (a == "--report-every") { if (has()) reportEvery = atof(argv[++i]); }
        else if (a == "--out") { if (has()) OUTDIR = argv[++i]; }
        else { fprintf(stderr, "unknown option %s\n", a.c_str()); usage(); return 1; }
    }
    std::error_code ec;
    std::filesystem::create_directories(OUTDIR, ec);
    std::filesystem::create_directories(OUTDIR + "/results", ec);
    LOGF = fopen((OUTDIR + "/search.log").c_str(), resume ? "a" : "w");

    for (int i = 0; i < RG; ++i) G0[i] = parse_gen(START_G[i]);
    gens_to_cols(G0, COL0);

    logline("=========================================================================\n");
    logline("hitset14 -- residual-targeted column replacement, [[14,3,5]] search\n");
    logline("started %s   threads %d\n", now_stamp().c_str(), threads);
    logline("=========================================================================\n");

    if (mode == 0 || mode == 1) { verify_start(true); if (mode == 1) return 0; }
    if (mode != 1) verify_start(mode == 2 ? false : true);

    switch (mode) {
        case 2: { bool ok = selftest(nrand); return ok ? 0 : 1; }
        case 3: analyze(rA, rB); return 0;
        case 4: if (calR < 3) { logline("give a rank: --calibrate R\n"); return 1; }
                calibrate(calR, calN, threads); return 0;
        case 6: if (searchR < 3) { logline("give a rank: --scan R\n"); return 1; }
                scan_improvements(searchR); return 0;
        case 5: if (searchR < 3) { logline("give a rank: --search R\n"); return 1; }
                run_search(searchR, threads, resume, reportEvery); return 0;
        default: usage(); return 0;
    }
}
