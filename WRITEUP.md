# Independent verification and rank bounds for Gaetz’s counterexample to Brenti’s log-concavity conjecture

**Note.** This is an independent computational follow-up to Christian Gaetz, *A counterexample to a log-concavity conjecture of Brenti*, [arXiv:2608.13836v1](https://arxiv.org/abs/2608.13836). It treats finite symmetric groups (type \(A\)) and the small exceptional/Weyl types \(H_3\), \(H_4\), \(F_4\), and \(B_2\)–\(B_5\). The lead computational result is a complete census of \(H_4\).

## 1. The claim being checked

For \(u\le v\) in a Coxeter group, the Kazhdan–Lusztig \(R\)-polynomials determine non-negative polynomials \(\widetilde{R}_{u,v}(q)\) and
\[
Q_{u,v}(q)=\sum_{i}a_i q^i
\]
by grouping even or odd powers of \(\widetilde{R}_{u,v}\) according to the parity of \(\ell(v)-\ell(u)\). Brenti conjectured that \((a_i)\) is log-concave. Gaetz disproves this in type \(A\) by exhibiting an explicit pair in \(S_{14}\), reports that an exhaustive search of \(S_8\) found no counterexample, and gives a separate example in type \(H_4\).

## 2. Complete classification in \(H_4\) (and \(H_3\), \(F_4\), \(B_2\)–\(B_5\))

Exceptional and Weyl types are built from the geometric representation over \(\mathbb{Z}[\varphi]\) (\(\varphi^2=\varphi+1\)) for \(H_3/H_4\), and from the integer Cartan matrix (asymmetric double bond) for \(F_4\) and \(B_n\). Bruhat order is the recursive right-descent criterion, tabulated by dynamic programming in increasing length. \(\widetilde{R}\) uses the same right-descent recurrence as type \(A\), including the shared `uint64_t` overflow abort (`add_shifted`).

**Counting convention.** “Intervals” means all comparable pairs \(u\le v\) (including short ones). “Reduced failures” means comparable, \(\ell(v)-\ell(u)\ge 4\), and **no common right descent**. A shared right descent has the same \(Q\) as the strictly shorter pair \((us,vs)\), which is counted when the reduced upper is processed (the same skip as type \(A\)). So 20,163 is a count of reduced failures, not of all intervals.

Checks independent of the census: the length histogram of the 14,400 enumerated elements matches the Poincaré polynomial of \(H_4\). Exponents are 1, 11, 19, 29, so \(P(q)=[2]_q[12]_q[20]_q[30]_q\) with \([n]_q=1+q+\cdots+q^{n-1}\), and \(\lvert W\rvert=P(1)=14400\), \(\ell(w_0)=60\). The same check holds for \(H_3\) (\([2]_q[6]_q[10]_q\)), \(F_4\), and \(B_2\)–\(B_3\). Tabulated Bruhat agrees with `le_recursive` on all of \(H_3\) and on 8,000 random \(H_4\) pairs, and agrees with the transitive closure of covering relations (length-\(+1\) comparable pairs) on both \(H_3\) and \(H_4\). First versus last right descent produce identical \(\widetilde{R}\) on all of \(H_3\).

| type | \(\lvert W\rvert\) | \(\ell(w_0)\) | intervals \(u\le v\) | \(Q\) log-concave fails | unimodal fails | internal zeros | min fail \(\ell(v)-\ell(u)\) |
|------|-------------------:|--------------:|---------------------:|------------------------:|---------------:|---------------:|-----------------------------:|
| \(H_3\) | 120 | 15 | 5,491 | 0 | 0 | 0 | — |
| \(B_2\) | 8 | 4 | 33 | 0 | 0 | 0 | — |
| \(B_3\) | 48 | 9 | 847 | 0 | 0 | 0 | — |
| \(B_4\) | 384 | 16 | 40,249 | 0 | 0 | 0 | — |
| \(B_5\) | 3,840 | 25 | 3,089,459 | 0 | 0 | 0 | — |
| \(F_4\) | 1,152 | 24 | 396,809 | 0 | 0 | 0 | — |
| \(H_4\) | 14,400 | 60 | 75,539,433 | **20,163** | 0 | 0 | **16** |

So among these groups, **log-concavity of \(Q\) fails only in \(H_4\)**. It fails often: 20,163 reduced intervals, 1,591 distinct \(Q\)-vectors. Every failure has `violation_at=2` (the same \(a_2^2<a_1 a_3\) pattern as the \(S_{14}\) example, where \(11^2<1\cdot 123\)). All of those \(Q\) remain unimodal and have no internal zeros.

### 2.1 Gaetz’s \(H_4\) example

With \(m(s_1,s_2)=5\), \(m(s_2,s_3)=m(s_3,s_4)=3\), the paper’s reduced words
\[
\begin{align*}
u&=s_2s_3s_2s_1s_2s_1s_4s_3s_2s_1s_2s_1s_3s_2s_1s_2s_3s_4s_3s_2s_1,\\
v&=s_1s_2s_1s_2s_3s_2s_1s_2s_1s_3s_2s_4s_3s_2s_1s_2s_1s_3s_2s_1s_2s_3s_4s_3s_2s_1s_2s_1s_3s_2s_1s_2s_4s_3s_2s_1s_2s_3s_4
\end{align*}
\]
are comparable of relative length 18, and
\[
Q_{u,v}(q)=q^9+14q^8+78q^7+220q^6+326q^5+234q^4+67q^3+8q^2+q,
\]
i.e. ascending coefficients \([0,1,8,67,234,326,220,78,14,1]\). Then \(8^2=64<67=1\cdot 67\).

That pair **shares the right descent** \(s_4\), so it is not itself one of the 20,163 reduced failures; the enumerator records the equivalent pair with common right descents peeled, same \(Q\). Certificate lines store a greedy first-right-descent peeling, which is a reduced word for the inverse element, so the published words appear verbatim in `H4.cert` as the encoding of the inverse pair \((u^{-1},v^{-1})\), which *is* reduced. This \(Q\) occurs for exactly 10 reduced intervals. Both the right-reduced pair and the inverse pair are in the certificate.

### 2.2 Smaller \(H_4\) failures

The shortest failures have relative length **16** (56 intervals, three \(Q\)-vectors):

| copies | \(Q\) (ascending) | inequality |
|-------:|-------------------|------------|
| 24 | \([0,1,7,51,121,119,55,12,1]\) | \(7^2=49<51\) |
| 16 | \([0,1,7,51,123,120,55,12,1]\) | \(7^2=49<51\) |
| 16 | \([0,1,7,52,124,120,55,12,1]\) | \(7^2=49<52\) |

One explicit length-16 pair, as reduced words in \(s_1,\ldots,s_4\) (re-verified; this pair is already right-reduced, and both it and its inverse appear in the certificate):
\[
\begin{align*}
u&=s_3s_2s_1s_2s_1s_3s_2s_4s_3s_2s_1s_2s_1s_3s_2s_4,\\
v&=s_2s_1s_2s_3s_2s_1s_2s_1s_3s_2s_1s_4s_3s_2s_1s_2s_1s_3s_2s_1s_2s_4s_3s_2s_1s_2s_1s_3s_2s_1s_4s_3,
\end{align*}
\]
with \(Q=[0,1,7,52,124,120,55,12,1]\).

The 56 length-16 failures partition into 6 inversion orbits of size 2 and 44 of size 1 on the right-reduced set (the size-1 mates share a right descent after inversion). Relative-length histogram of the 20,163 failures starts at 16 and reaches 58; even lengths dominate because the typical violation still sits at the \(q^1,q^2,q^3\) end of \(Q\). See `H4_SUMMARY.md` for the full length table, unique min-length \(Q\)s, and orbit counts.

### 2.3 Inversion orbits

In \(H_4\), \(w_0\) is central, so conjugation by \(w_0\) is trivial; the remaining involution on pairs is \((u,v)\mapsto(u^{-1},v^{-1})\), which preserves \(Q\). It does not preserve the right-reduced subset: 13,213 of the 20,163 inverses share a right descent, so they are counted only after peeling. On the reduced set this gives 13,307 orbits of size 1 (94 truly self-inverse, 13,213 unpaired) and 3,428 of size 2, accounting for all 20,163. The odd total comes from the odd number of size-1 orbits. Every reduced failure’s right-reduced inverse is a reduced failure with the same \(Q\).

This is a data point on why the first type-\(A\) failure is so large, but it does **not** answer that question: it **sharpens** it. In \(H_4\) the first failure is already at relative length 16, two less than Gaetz’s published example, and far below the type-\(A\) length 29. The gap between 16 (rank 4, non-crystallographic) and the type-\(A\) minimum in \(\{10,\ldots,14\}\) is now wider and stranger.

## 3. Independent check of the \(S_{14}\) pair

Let
\[
\begin{align*}
u&=[1,2,5,7,9,3,11,4,6,12,13,8,10,14],\\
v&=[8,4,6,12,10,13,1,14,2,11,3,7,5,9].
\end{align*}
\]
An independent C++17 implementation of the right-descent recurrence for \(\widetilde{R}\) (with \(\widetilde{R}_{x,y}=0\) off Bruhat order, and type-\(A\) northwest rank comparison) yields

\[
Q_{u,v}(q)=q^{14}+16q^{13}+101q^{12}+333q^{11}+630q^{10}+695q^9+425q^8+123q^7+11q^6+q^5,
\]

in agreement with the paper. In particular \(11^2=121<123=1\cdot 123\), so log-concavity fails.

Two small remarks:

- \(u\le v\) in Bruhat order, as required.
- The relative length is \(\ell(v)-\ell(u)=46-17=\mathbf{29}\), not 28. Degree 14 of \(Q\) is compatible with either parity of \(\ell(v)-\ell(u)\); the odd case is the one that occurs. This does not affect the counterexample (possible typo in the note).

The pair is already reduced (no common left or right descent) and is pattern-minimal in the following sense: every nonempty proper value-pattern and position-pattern of \((u,v)\) in \(S_k\) for \(5\le k\le 13\) is either Bruhat-incomparable or has log-concave \(Q\). The Bruhat-cover neighbourhood of the pair consists of four combinatorially equivalent copies (the orbit under inversion and diagram reversal), all with the same \(Q\) and all still in \(S_{14}\).

## 4. Exhaustive type-\(A\) search

The same program enumerates, for each \(v\in S_n\), the principal order ideal \(\{u:u\le v\}\) by covering transpositions, then evaluates \(Q_{u,v}\) exactly. Pairs that share a right descent are not evaluated: the recurrence gives them the same \(\widetilde{R}\) as the strictly shorter pair \((us,vs)\), which is visited when \(vs\) is processed. Coefficient arithmetic is exact `uint64_t` with overflow abort. Log-concavity tests use 128-bit products. While enumerating we also record unimodality failures and internal zeros of \(Q\); none occurred through \(S_9\).

| \(n\) | comparable pairs \(u\le v\) | log-concavity failures | wall time (16 threads, i5-12600K) |
|------:|----------------------------:|-----------------------:|----------------------------------:|
| 5 | 3,781 | 0 | \(<0.1\) s |
| 6 | 98,407 | 0 | \(0.15\) s |
| 7 | 3,550,919 | 0 | \(1.8\) s |
| 8 | 170,288,585 | 0 | \(22\) s |
| 9 | 10,501,351,657 | 0 | \(27.3\) min |

The \(S_8\) result confirms Gaetz’s exhaustive search. The \(S_9\) result is new: **there is no type-\(A\) counterexample in \(S_n\) for \(n\le 9\)**. Combined with the verified \(S_{14}\) example, the smallest possible rank of a type-\(A\) counterexample is therefore some \(n\) with
\[
10\le n\le 14.
\]

`--exhaustive N --symmetry` further restricts to canonical upper elements under the maps \(v\mapsto v^{-1}\) and conjugation by \(w_0\) (diagram reversal), a \(\le 4\)-element orbit of automorphisms that preserve \(\widetilde{R}\). That quotient is used for the \(S_{10}\) run below.

## 5. Directed search in \(S_{10}\)–\(S_{13}\)

Directed searches in \(S_{10}\) through \(S_{13}\) found no failure:

- all \(2^{14}-2\) proper value and position patterns of the \(S_{14}\) pair;
- inversion-table scalings of that pair, with small coordinate perturbations;
- \(6\times 10^5\) random reduced intervals per rank \(n=10,11,12,13\), biased toward long intervals.

Typical sampled \(Q\)-vectors in these ranks remain comfortably log-concave (the closest ratio \(a_i^2/(a_{i-1}a_{i+1})\) seen was about \(1.24\) in \(S_{13}\), versus \(121/123\approx 0.984\) for the \(S_{14}\) example). This is negative evidence only; it does not rule out a counterexample in \(S_{10}\)–\(S_{13}\).

## 6. Exhaustive \(S_{10}\) with symmetry quotient

`--exhaustive 10 --symmetry` enumerates only canonical upper elements under inversion and \(w_0\)-conjugation. A full run is hours rather than days; the certificate and checkpoint are `S10.cert` / `S10.done`. Outcome will be either a type-\(A\) counterexample of rank 10, or the bound \(11\le n\le 14\). This note does not wait for that run.

## 7. Code

The implementation is C++17/OpenMP (`src/main.cpp`, `src/coxeter.hpp`):

[https://github.com/guillochon/2608.13836v1-num](https://github.com/guillochon/2608.13836v1-num)

```sh
make
./brenti_search --self-test
./brenti_search --verify-s14
./brenti_search --verify-h4
./brenti_search --coxeter-orbits H4.cert --summary H4_SUMMARY.md
OMP_NUM_THREADS=4 ./brenti_search --coxeter H4 --out H4.cert
OMP_NUM_THREADS=$(nproc) ./brenti_search --exhaustive 10 --symmetry --out S10.cert --checkpoint S10.done
```

`--self-test` checks group orders, Poincaré polynomials, covering-closure Bruhat, the \(H_4\) examples, \(H_3\) Bruhat (tabulated versus recursive), first-versus-last descent \(\widetilde{R}\) on \(H_3\), overflow abort of `add_shifted`, left/right/inverse \(\widetilde{R}\) on random \(S_7\) pairs, and the \(S_{14}\) polynomial. `--verify-h4` also confirms that Gaetz’s pair and the length-16 example sit in `H4.cert` when that file is present. The 20,163-line certificate is gitignored; regenerate it as above, or see `H4_SUMMARY.md`.

## 8. What remains open

- The exact smallest \(n\) with a type-\(A\) counterexample (some integer in \(\{10,11,12,13,14\}\); \(S_{10}\) is running).
- Whether the \(S_{14}\) interval is isomorphic to an interval already present in some smaller \(S_n\) (pattern deletion does not produce one; a poset embedding could still exist).
- A conceptual reason the first type-\(A\) failure is so large. That \(H_4\) already fails at relative length 16 does not answer this — it sharpens the question, because the gap between rank-4 non-crystallographic length 16 and the type-\(A\) minimum in \(\{10,\ldots,14\}\) is now wider.
- Whether some inflation of the \(S_{14}\) pair makes the log-concavity ratio arbitrarily small.
