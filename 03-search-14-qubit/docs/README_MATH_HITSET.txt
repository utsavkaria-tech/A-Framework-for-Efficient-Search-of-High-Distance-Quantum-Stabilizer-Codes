=========================================================================================
RESIDUAL-TARGETED COLUMN REPLACEMENT -- searching for a [[14,3,5]] stabilizer code
=========================================================================================
SELF-CONTAINED.  Everything this experiment produces lives in residual_columns/.  No file,
binary, log, checkpoint or result of any earlier search in this project is read or written.
The starting 11 generators are compiled into the source and re-verified at every startup.

-----------------------------------------------------------------------------------------
1.  THE PROBLEM, AND THE REFORMULATION EVERYTHING RESTS ON
-----------------------------------------------------------------------------------------
We start from 11 commuting generators on 14 qubits (three of weight 8, eight of weight 10)
with rank 11, for which exactly five non-identity Pauli errors of weight <= 4 are
undetected.  Those five are the whole problem; the code has distance 3 instead of 5.

Read the check matrix H = [X|Z] COLUMN-wise.  Qubit j contributes x_j, z_j in F_2^11 (the
j-th columns of X and of Z); set

        W_j  =  span{ x_j , z_j }   <=   F_2^11 .

An error E acting as X^{a_j} Z^{b_j} on qubit j has syndrome s(E) = SUM_j (a_j z_j + b_j x_j).
So the syndrome contributed by qubit j is

        X -> z_j ,      Z -> x_j ,      Y -> x_j + z_j ,

which are precisely the three NON-ZERO elements of W_j, each occurring once.  Therefore

  d >= 5  <=>  no 1 <= k <= 4 distinct qubits j_1..j_k and non-zero v_i in W_{j_i} satisfy
               v_1 + ... + v_k = 0 .                                                   (*)

Three consequences are immediate and are used throughout:

  * dim W_j = 2 for every j.  If dim W_j <= 1 one of the three contributions is 0, i.e. a
    weight-1 undetected error.  So every column is a LINE of PG(10,2); there are 698027.
  * k = 2 forces W_i ∩ W_j = 0 for i != j; in particular no two columns may be equal.
  * (*) depends only on the 14 SUBSPACES, never on the basis (x_j, z_j) chosen inside one.
    Changing that basis is exactly a local Clifford on qubit j (GL(2,2) = S_3 permuting the
    labels X, Y, Z).  Working with subspaces therefore quotients out the entire local
    Clifford group for free -- a factor 6^14 = 7.8e10 -- with no risk of over-quotienting,
    because the objective is literally a function of the subspaces.

Commutation.  <g_a, g_b> = SUM_j (x_{aj} z_{bj} + z_{aj} x_{bj}), so with the PLUCKER form
pl(W) = x z^T + z x^T (symmetric, zero diagonal, i.e. a point of Lambda^2(F_2^11) = F_2^55)
the whole commutation requirement collapses to ONE linear equation

        SUM_{j=1..14}  pl(W_j)  =  0 .                                                (**)

pl(W) is basis-independent (a change of basis multiplies it by det = 1), and under the
Pluecker embedding the lines of PG(10,2) correspond bijectively to the DECOMPOSABLE non-zero
2-vectors x^z, equivalently to the alternating forms of rank exactly 2.

Rank.  rank H = 11  <=>  W_1 + ... + W_14 = F_2^11.

Generator basis change.  A in GL(11,2) sends W_j -> A W_j and pl -> A pl A^T, so it
preserves both (*) and (**).  It is a genuine symmetry, but it moves ALL columns at once and
is therefore broken the moment any column is frozen; it is deliberately NOT used for pruning.

-----------------------------------------------------------------------------------------
2.  WHY THE HITTING-SET IDEA IS EXACTLY RIGHT -- AND WHAT IT REALLY IS
-----------------------------------------------------------------------------------------
Freeze the columns outside a set S and let the r = |S| columns in S vary.  An undetected
weight-<=4 error supported entirely outside S has a syndrome built only from frozen columns,
so it survives untouched.  Hence

        S must meet the support of every currently undetected weight-<=4 error.

For the five residuals that is a hitting-set condition, and it is NECESSARY, not heuristic.
But it is not an extra ingredient either: it is the f = 0 case of the following exact
decomposition.  Split an error's support into its free part (f qubits) and its frozen part,
and let C_m be the set of syndromes of errors of weight <= m supported on FROZEN qubits
(the empty error contributing 0).  Then (*) is equivalent to, for free points v_i,

        f = 0 :  0 not in C_4 \ {0}      <-  feasibility of S (the hitting condition)
        f = 1 :  v_1 not in C_3
        f = 2 :  v_1 + v_2 not in C_2
        f = 3 :  v_1 + v_2 + v_3 not in C_1
        f = 4 :  v_1 + v_2 + v_3 + v_4 != 0

The f = 1 row is the powerful one.  It says every free column must be a line ALL THREE of
whose points avoid C_3, and that set does not depend on which free qubit we are filling.
Call the surviving lines the ADMISSIBLE SET A.  Measured (see section 6), |A| collapses from
698027 to between 3 and 10 at r = 3, and to a few tens of thousands at r = 6.

The whole f = 1..4 family is then collapsed into a single running table.  Let F[w] be the
set of syndromes of weight-w errors supported on the frozen qubits TOGETHER WITH the free
columns already chosen.  A further column is admissible exactly when none of its three
points lies in F[0] u F[1] u F[2] u F[3] -- adding such a point would complete a weight-<=4
error of syndrome zero.  F extends in one pass when a column is fixed:

        F'[w] = F[w]  u  ( F[w-1] + p_i ),   i = 1,2,3 .

Testing a line therefore costs three bit lookups instead of re-deriving every pair, triple
and quadruple, and a line eliminated at depth 2 is never re-tested deeper in that branch.
Measured effect at r = 5: 36.1 million search nodes fall to 421 thousand.

-----------------------------------------------------------------------------------------
3.  THE STRONGEST FORMULATION FOUND: DECOMPOSING AN ALTERNATING FORM
-----------------------------------------------------------------------------------------
By (**) the free columns must satisfy SUM_{j in S} pl(W_j) = T, where T is the Pluecker sum
of the frozen columns -- equivalently, of the ORIGINAL columns of S, since the starting code
satisfies (**).  So the search is exactly:

        WRITE T AS A SUM OF r DECOMPOSABLE 2-VECTORS, EACH DRAWN FROM A,
        subject to the f >= 2 conditions above.

That is a classical object and it yields an exact prune costing one 11x11 F_2 rank
computation per node.  A sum of t decomposables has rank <= 2t because each summand has
rank 2; conversely (Darboux / symplectic normal form) every alternating form of rank 2m IS a
sum of m decomposables.  Hence

        min #decomposables summing to T   =   rank(T) / 2 ,

    PRUNE (exact):  after choosing k of the r columns the remaining target T_k must satisfy
                    rank(T_k) <= 2 (r - k).

Nothing that could complete is discarded: any completion writes T_k as a sum of r-k
decomposables and so has rank at most 2(r-k).  Two consequences make the search collapse:

  * THE LAST TWO COLUMNS ARE NEVER ENUMERATED.  With two left, rank(T) <= 4 and
      rank 0 : the two lines would have to be equal, which (*) forbids -- dead branch;
      rank 4 : dim(W_1 + W_2) <= 4 and supp(T) <= W_1 + W_2, so W_1 (+) W_2 = supp(T)
               exactly; both lines lie in a KNOWN 4-space, only 35 lines to test;
      rank 2 : W_1 ∩ W_2 is 1-dimensional, spanned by some u, and T = u ^ (a+b), so u lies
               in supp(T); enumerating the admissible lines through the three points of
               supp(T) covers every possibility.
    In each case the partner line is then forced as the support of T - pl(W_1).
  * subtracting a decomposable changes the rank by -2, 0 or +2, and it drops by 2 exactly
    when the line is a HYPERBOLIC plane of T_k.  When rank(T) = 2r EVERY step must drop,
    which is extremely restrictive.  Measured: rank(T) = 2r for all subsets at r = 3, 4 and
    for nearly all at r = 5, which is why those searches finish in seconds.  At r = 6,
    2r = 12 exceeds the maximum possible rank 10 of an alternating form on F_2^11, so there
    is always slack of at least one step -- and that, not the size of A, is what makes r = 6
    three orders of magnitude more expensive than r = 5.

This is why the search is exhaustive and fast rather than a 698027^r enumeration.

-----------------------------------------------------------------------------------------
4.  SYMMETRY THAT IS USED, AND SYMMETRY THAT IS DELIBERATELY NOT
-----------------------------------------------------------------------------------------
USED -- local Clifford (S_3)^14: absorbed exactly, by searching subspaces (section 1).
USED -- permutations of the FREE qubits among themselves.  Two assignments differing by such
        a permutation are carried into each other by a qubit permutation fixing every frozen
        column, so they are equivalent codes with the same distance.  The search therefore
        enumerates the chosen lines in strictly increasing index order, an r! reduction.
        This is valid for the objective "does a d >= 5 code exist"; solutions are reported in
        canonical form and each represents an orbit of size up to r!.
NOT USED -- GL(11,2) generator basis change: a symmetry of the code, but it moves the frozen
        columns too, so it cannot be applied once part of the matrix is fixed.
NOT USED -- permutations mixing free and frozen qubits, for the same reason.
NOT USED -- any quotient by the residual errors' own symmetries: the five residuals are
        input data, not a group, and identifying configurations through them would not be a
        valid quotient for the objective.

-----------------------------------------------------------------------------------------
5.  DEGENERACY -- WHAT IS AND IS NOT PROVEN
-----------------------------------------------------------------------------------------
Condition (*) forbids ANY weight-<=4 Pauli in the centralizer N(S).  That is the PURE
(non-degenerate) requirement, and it is exactly the acceptance condition as stated: "no
non-identity undetected Pauli error of weight <= 4".  A degenerate code is additionally
allowed to contain weight-<=4 elements of the stabilizer group S itself, since those act
trivially on the code space.

The program handles this honestly on both sides:

  * ACCEPTANCE is the true distance.  A candidate is accepted when no weight-<=4 error lies
    in N(S) \ S; the report states whether the code is pure or degenerate.  So any solution
    found is a genuine d >= 5 code, degenerate or not.
  * EXHAUSTIVENESS is claimed only over PURE codes, because the enumeration prunes on (*).
    A degenerate solution would need a non-zero a in F_2^11 with wt(E_a) <= 4, i.e.
    wt(E_a) = 14 - #{ j : a _|_ W_j } <= 4, so a must be orthogonal to at least 10 of the 14
    columns.  At most r of those are free, so a must be orthogonal to at least 10 - r of the
    FROZEN columns.  The program computes, for every a, how many frozen columns it is
    orthogonal to, and reports for each subset whether that closes the gap.  Measured: it
    closes for 2 of the 12 subsets at r = 3 and for none at r >= 4.  So the r = 3, 4, 5
    results below are proofs about PURE codes, and the degenerate case is left open -- this
    is stated in the program's own output rather than glossed over.

For reference the STARTING code has minimum non-zero stabilizer weight 6, so all five of its
residuals are genuine distance failures and none is a degenerate freebie.

-----------------------------------------------------------------------------------------
6.  WHAT THE PROGRAM MEASURED
-----------------------------------------------------------------------------------------
Minimum hitting set.  Residual supports are {3,12,13}, {6,7,11}, {0,5,9,11}, {1,6,7,8},
{2,4,10,11}.  The first and fourth are disjoint, so two columns can never hit all five, and
the third column must lie in {0,5,9,11} ∩ {2,4,10,11} = {11}.  Hence the minimum is 3 and
every r = 3 subset contains qubit 11 -- which is exactly what the enumeration produces.

  r   C(14,r)  hit all 5  feasible   |A| min/median/max      rank(T)   2r
  --------------------------------------------------------------------------
  3      364        12        12        3 /      9 /     10   6         6
  4     1001       156       156      134 /    200 /    544   6 - 8     8
  5     2002       691       691     1785 /   2881 /   5017   8 - 10   10
  6     3003      1644      1644    15211 /  20572 /  34499   6 - 10   12

-----------------------------------------------------------------------------------------
7.  VERIFICATION -- AND WHY A "PROOF" HERE IS NOT SELF-CONFIRMING
-----------------------------------------------------------------------------------------
A depth-first search that wrongly prunes everything reports "no solution" for every r and
looks exactly like a proof.  The self-test is built to make that failure mode impossible to
miss:

  1. the starting code is re-verified from scratch: weights, all 55 commuting pairs, rank 11,
     every column 2-dimensional, the Pluecker identity (**), and all 91770 weight-<=4 errors
     enumerated and their syndromes recomputed by direct symplectic products;
  2. the five residuals produced must equal the five reported errors, literally;
  3. a SECOND, structurally different distance routine works only from the column subspaces
     and checks (*) directly; its count of violating tuples must equal the number of
     zero-syndrome errors found by the first;
  4. rank(pl(line)) = 2 and the support recovers the line, on thousands of random lines;
  5. rank(sum of t decomposables) <= 2t on thousands of random sums -- the prune's premise;
  6. three DELIBERATE breakages (duplicate a column, collapse a column to dim 1, break
     commutation) must each be caught;
  7. a KNOWN configuration -- the original columns of a subset -- is planted back into A and
     the search core must RECOVER it, with the distance filters switched off so that the
     rank prune, the forced last column, the line lookup and the ordering rule are validated
     independently of the filters.  Without this test the null results would be worthless;
  8. TWO INDEPENDENT SEARCH ENGINES.  dfs() enumerates to depth r-1, forces only the last
     column and tests every pair/triple/quadruple condition explicitly.  dfs2() carries the
     running forbidden-point table and solves the last two columns algebraically.  They are
     compared on 19 subsets spanning r = 3, 4, 5 and must agree on the number of candidates
     and of accepted codes.

Every accepted candidate is then re-verified from scratch by the independent routine before
being written to results/.

-----------------------------------------------------------------------------------------
8.  RESULTS SO FAR AND HOW TO RUN THE REST
-----------------------------------------------------------------------------------------
    r = 3   PROVEN EXHAUSTIVE (pure)   0 solutions      12 nodes,      < 1 s
    r = 4   PROVEN EXHAUSTIVE (pure)   0 solutions     601 nodes,      < 1 s
    r = 5   PROVEN EXHAUSTIVE (pure)   0 solutions  421211 nodes,        2 s

No candidate at r = 3, 4 or 5 ever reached the distance verifier: the commutation equation
(**) and the admissibility filter are jointly unsatisfiable there.  Replacing three, four or
five columns of this stabilizer CANNOT produce a pure d >= 5 code.

IMPROVEMENT SCAN (--scan), which drops the distance filter and enumerates every valid
commuting rank-11 replacement, so that it can see improvements short of d >= 5:

    r = 3   EXHAUSTIVE, 21435 complete replacements, 0 slack branches, 100 s
            BEST residual count reachable = 5 -- the starting code is already optimal
            among ALL 3-column replacements.  Not merely "no d >= 5": no reduction at all.

r >= 4 is out of reach for the unfiltered scan: at r = 4 the candidate set at depth 0 is the
5440 hyperbolic lines of an 8-dimensional space, and the leaf count per subset is about
3.6e7, so a full scan would be ~1e10 complete replacements.  The d >= 5 search at r = 4 is
unaffected -- it is the admissibility filter that makes that one cheap.

    cd "...\stabilizer_14q\residual_columns"
    .\hitset14.exe --search 6 --threads 12          <- the next real run, ~5h 35m
    .\hitset14.exe --search 6 --threads 12 --resume <- after any interruption

r = 6 calibration (5 evenly spread subsets, 129-169 s each, mean 146.8 s):
    1644 feasible subsets, 67h 02m on 1 thread, 5h 35m on 12 threads.
The cost jump from r = 5 (2 s) to r = 6 (5.6 h) is structural, not incidental: an
alternating form on F_2^11 has rank at most 10, so at r = 6 the requirement
rank(T) <= 2r = 12 is slack by at least one step and the "every step must drop the rank"
regime that collapses r <= 5 no longer applies.  Expect r = 7 to be far worse again.

Other modes: --verify, --selftest [N], --analyze a b, --calibrate R n, --scan R,
--report-every S, --out DIR.  checkpoint_rR.txt is rewritten as subsets complete and
--resume skips exactly the finished ones.

STATUS LABELS the program uses and means literally:
    PROVEN / EXHAUSTIVE  the tree for that r was completed; no pure d >= 5 code exists there
    BEST KNOWN           never used for a null result, only for a found code before it has
                         passed independent verification
A null result at one r says nothing about any other r, and nothing about degenerate codes
except on the subsets where section 5's closure applies.
