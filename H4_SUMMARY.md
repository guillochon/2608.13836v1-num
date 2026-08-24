# H4 census summary

Regenerate the 20,163-line certificate with

```sh
OMP_NUM_THREADS=4 ./brenti_search --coxeter H4 --out H4.cert
```

and re-check orbits with `./brenti_search --coxeter-orbits H4.cert --summary H4_SUMMARY.md`.

## Counting convention

- **intervals** = all comparable pairs \(u\le v\) (including short ones): **75,539,433** in \(H_4\).
- **reduced failures** = comparable, \(\ell(v)-\ell(u)\ge 4\), **no common right descent**. A shared right descent has the same \(Q\) as the strictly shorter pair \((us,vs)\), which is counted when the reduced upper element is processed (same skip as type \(A\)).
- So **20,163** is the number of reduced failures, not of all intervals.

## Totals

| quantity | value |
|----------|------:|
| \(|W|\) | 14,400 |
| \(\ell(w_0)\) | 60 |
| intervals | 75,539,433 |
| reduced log-concavity failures | 20,163 |
| distinct \(Q\)-vectors | 1591 |
| unimodal / internal-zero failures | 0 / 0 |
| violation index | 2 in every failure |
| min relative length | 16 |
| Gaetz Example 4 \(Q\) copies | 10 |
| self-inverse reduced pairs | 94 |
| reduced pairs whose inverse is also reduced | 6950 |
| inversion orbits of size 1 on the reduced set | 13307 |
| inversion orbits of size 2 on the reduced set | 3428 |

Every reduced failure's right-reduced inverse is a reduced failure with the same \(Q\). In \(H_4\), \(w_0\) is central, so conjugation by \(w_0\) is trivial. The involution \((u,v)\mapsto(u^{-1},v^{-1})\) preserves \(Q\) but does not preserve the right-reduced subset: 13213 inverses share a right descent, so they are counted in reduced form rather than as a second copy. On the reduced set this gives 13307 orbits of size 1 (94 truly self-inverse, 13213 unpaired) and 3428 of size 2, accounting for all 20,163. The odd total comes from the odd number of size-1 orbits.

## Minimum length (relative length 16)

56 reduced failures. The raw inversion map closes this slice (12 of 56 inverses are already right-reduced and still length 16), with 0 self-inverse pairs, 44 orbits of size 1 and 6 of size 2 (44 + 2×6 = 56).

| copies | \(Q\) (ascending) |
|-------:|-------------------|
| 24 | `[0,1,7,51,121,119,55,12,1]` |
| 16 | `[0,1,7,51,123,120,55,12,1]` |
| 16 | `[0,1,7,52,124,120,55,12,1]` |

One explicit length-16 pair, as reduced words (left-to-right product; this pair is reduced and both it and its inverse appear in `H4.cert`):

```
u = 3,2,1,2,1,3,2,4,3,2,1,2,1,3,2,4
v = 2,1,2,3,2,1,2,1,3,2,1,4,3,2,1,2,1,3,2,1,2,4,3,2,1,2,1,3,2,1,4,3
Q = [0,1,7,52,124,120,55,12,1]   (7^2=49<52)
```

## Gaetz Example 4

Multiplying the published reduced words (left to right, \(s_1,\ldots,s_4\) with \(m(s_1,s_2)=5\)) recovers a comparable pair of relative length 18 with the published \(Q\). That pair **shares the right descent** \(s_4\), so it is not itself one of the 20,163 reduced failures; the enumerator records the equivalent pair with common right descents peeled, which has the same \(Q\).

Certificate words are the greedy first-right-descent peeling, which is a reduced word for the **inverse** element. Consequently the published words appear verbatim in `H4.cert` as the encoding of the inverse pair \((u^{-1},v^{-1})\), which *is* reduced and is one of the 10 intervals with this \(Q\):

```
u = 2,3,2,1,2,1,4,3,2,1,2,1,3,2,1,2,3,4,3,2,1
v = 1,2,1,2,3,2,1,2,1,3,2,4,3,2,1,2,1,3,2,1,2,3,4,3,2,1,2,1,3,2,1,2,4,3,2,1,2,3,4
length = 18
Q = [0,1,8,67,234,326,220,78,14,1]   (8^2=64<67)
```

Gaetz inverse pair in cert: yes; right-reduced pair in cert: yes. That \(Q\) occurs for 10 reduced intervals.

## Failures by relative length

| length | reduced failures |
|-------:|-----------------:|
| 16 | 56 |
| 18 | 362 |
| 20 | 669 |
| 22 | 1140 |
| 24 | 1267 |
| 26 | 1736 |
| 27 | 44 |
| 28 | 1905 |
| 29 | 124 |
| 30 | 2019 |
| 31 | 82 |
| 32 | 1992 |
| 33 | 159 |
| 34 | 1521 |
| 35 | 159 |
| 36 | 1347 |
| 37 | 65 |
| 38 | 1511 |
| 39 | 49 |
| 40 | 1108 |
| 41 | 41 |
| 42 | 1001 |
| 43 | 13 |
| 44 | 569 |
| 45 | 1 |
| 46 | 406 |
| 48 | 402 |
| 50 | 196 |
| 52 | 145 |
| 54 | 50 |
| 56 | 21 |
| 58 | 3 |

## Poincaré check

Exponents of \(H_4\) are 1, 11, 19, 29, so the length generating function is \(P(q)=[2]_q[12]_q[20]_q[30]_q\) with \([n]_q=1+q+\cdots+q^{n-1}\). The histogram of the 14,400 enumerated lengths matches this polynomial (PASS).
