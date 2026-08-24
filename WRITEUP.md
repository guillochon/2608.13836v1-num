# Independent verification and rank bounds for Gaetz’s counterexample to Brenti’s log-concavity conjecture

**Note.** This is an independent computational follow-up to Christian Gaetz, *A counterexample to a log-concavity conjecture of Brenti*, [arXiv:2608.13836v1](https://arxiv.org/abs/2608.13836). It concerns only finite symmetric groups (type \(A\)). The paper’s separate type \(H_4\) example is not treated.

## 1. The claim being checked

For \(u\le v\) in a Coxeter group, the Kazhdan–Lusztig \(R\)-polynomials determine non-negative polynomials \(\widetilde{R}_{u,v}(q)\) and
\[
Q_{u,v}(q)=\sum_{i}a_i q^i
\]
by grouping even or odd powers of \(\widetilde{R}_{u,v}\) according to the parity of \(\ell(v)-\ell(u)\). Brenti conjectured that \((a_i)\) is log-concave. Gaetz disproves this in type \(A\) by exhibiting an explicit pair in \(S_{14}\), and reports that an exhaustive search of \(S_8\) found no counterexample (the largest symmetric group he searched exhaustively).

## 2. Independent check of the \(S_{14}\) pair

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
- The relative length is \(\ell(v)-\ell(u)=46-17=\mathbf{29}\), not 28. Degree 14 of \(Q\) is compatible with either parity of \(\ell(v)-\ell(u)\); the odd case is the one that occurs. This does not affect the counterexample.

The pair is already reduced (no common left or right descent) and is pattern-minimal in the following sense: every nonempty proper value-pattern and position-pattern of \((u,v)\) in \(S_k\) for \(5\le k\le 13\) is either Bruhat-incomparable or has log-concave \(Q\). The Bruhat-cover neighbourhood of the pair consists of four combinatorially equivalent copies (the orbit under inversion and diagram reversal), all with the same \(Q\) and all still in \(S_{14}\).

## 3. Exhaustive type-\(A\) search

The same program enumerates, for each \(v\in S_n\), the principal order ideal \(\{u:u\le v\}\) by covering transpositions, then evaluates \(Q_{u,v}\) exactly. Pairs that share a right descent are not evaluated: the recurrence gives them the same \(\widetilde{R}\) as the strictly shorter pair \((us,vs)\), which is visited when \(vs\) is processed. Coefficient arithmetic is exact `uint64_t` with overflow abort. Log-concavity tests use 128-bit products.

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

## 4. Directed search in \(S_{10}\)–\(S_{13}\)

A full enumeration of \(S_{10}\) is still out of reach at this scale (the interval count appears to grow by a factor of roughly \(50\)–\(60\) per rank). Directed searches in \(S_{10}\) through \(S_{13}\) found no failure:

- all \(2^{14}-2\) proper value and position patterns of the \(S_{14}\) pair;
- inversion-table scalings of that pair, with small coordinate perturbations;
- \(6\times 10^5\) random reduced intervals per rank \(n=10,11,12,13\), biased toward long intervals.

Typical sampled \(Q\)-vectors in these ranks remain comfortably log-concave (the closest ratio \(a_i^2/(a_{i-1}a_{i+1})\) seen was about \(1.24\) in \(S_{13}\), versus \(121/123\approx 0.984\) for the \(S_{14}\) example). This is negative evidence only; it does not rule out a counterexample in \(S_{10}\)–\(S_{13}\).

## 5. Code

The implementation is a single C++17/OpenMP file:

[https://github.com/guillochon/2608.13836v1-num](https://github.com/guillochon/2608.13836v1-num)

```sh
make
./brenti_search --verify-s14
OMP_NUM_THREADS=$(nproc) ./brenti_search --exhaustive 9 --out S9.cert --checkpoint S9.done
```

It is source-only (no Sage/coxeter3 dependency), intended so that the \(S_{14}\) polynomial and the \(S_{\le 9}\) enumerations can be reproduced on any machine with a C++17 OpenMP compiler.

## 6. What remains open

- The exact smallest \(n\) with a type-\(A\) counterexample (some integer in \(\{10,11,12,13,14\}\)).
- Whether the \(S_{14}\) interval is isomorphic to an interval already present in some smaller \(S_n\) (pattern deletion does not produce one; a poset embedding could still exist).
- A conceptual reason the first type-\(A\) failure is so large, given that type \(H_4\) already fails.
