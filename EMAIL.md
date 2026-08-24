# Draft email to Christian Gaetz

**To:** gaetz@berkeley.edu  
**Subject:** Independent check of arXiv:2608.13836 — \(H_4\) classified (20,163 reduced failures; first at length 16), \(S_9\) clean, \(S_{10}\) running

---

Dear Professor Gaetz,

I recently read your note *A counterexample to a log-concavity conjecture of Brenti* (arXiv:2608.13836v1) and wrote a small independent C++/OpenMP implementation of the \(\widetilde{R}\) recurrence.

The item I wanted to share first is a complete census of \(H_4\). Your Example 4 reproduces, with
\[
Q=q^9+14q^8+78q^7+220q^6+326q^5+234q^4+67q^3+8q^2+q
\]
(relative length 18, violation \(8^2<67\) at the second nonzero coefficient). A full pass over the group is small enough to finish: **75,539,433** Bruhat intervals \(u\le v\), of which **20,163 reduced failures** of log-concavity (1,591 distinct \(Q\)-vectors). A pair is counted as a reduced failure only if it is comparable, \(\ell(v)-\ell(u)\ge 4\), and has **no common right descent** — the same skip as in type \(A\), since a shared right descent has the same \(Q\) as a strictly shorter pair that is counted when the reduced upper is processed. So 20,163 is not a count of all intervals.

Every one of those 20,163 failures violates at index 2, and all remain unimodal with no internal zeros. That is the same \(a_2^2<a_1 a_3\) pattern as your \(S_{14}\) example (\(11^2<1\cdot 123\)). By contrast \(H_3\), \(F_4\), and \(B_2\)–\(B_5\) are entirely clean.

Your published pair shares the right descent \(s_4\), so it is identified with a right-reduced pair of the same \(Q\); both that reduced pair and the inverse pair \((u^{-1},v^{-1})\) appear among the 20,163. The published reduced words appear verbatim in the certificate as the encoding of the inverse pair. (Certificate words are a greedy right-descent peeling, which is a reduced word for the inverse element.) Your \(Q\) occurs for exactly 10 reduced intervals.

The shortest failures have relative length **16** (56 intervals, three \(Q\)-vectors). One explicit length-16 pair, as reduced words in \(s_1,\ldots,s_4\) with \(m(s_1,s_2)=5\):
\[
\begin{align*}
u&=s_3s_2s_1s_2s_1s_3s_2s_4s_3s_2s_1s_2s_1s_3s_2s_4,\\
v&=s_2s_1s_2s_3s_2s_1s_2s_1s_3s_2s_1s_4s_3s_2s_1s_2s_1s_3s_2s_1s_2s_4s_3s_2s_1s_2s_1s_3s_2s_1s_4s_3,
\end{align*}
\]
with \(Q=[0,1,7,52,124,120,55,12,1]\) and \(7^2=49<52\). I re-verified this pair. The 56 length-16 failures partition into 6 inversion orbits of size 2 and 44 of size 1 on the right-reduced set (the size-1 mates share a right descent after inversion, so they are counted only in reduced form).

That \(H_4\) already fails at length 16 does not explain why the first type-\(A\) failure is so large — it sharpens the question. The gap between 16 (rank 4, non-crystallographic) and the type-\(A\) minimum in \(\{10,\ldots,14\}\) is now wider and stranger.

On type \(A\): your \(S_{14}\) pair checks out, including \(11^2<1\cdot 123\). One minor numerical remark, which does not affect the counterexample: the relative length is \(\ell(v)-\ell(u)=46-17=29\) rather than 28 (degree 14 of \(Q\) is compatible with either parity). I also repeated the exhaustive type-\(A\) search. There are no log-concavity failures among the \(170{,}288{,}585\) intervals in \(S_8\), in agreement with your report, nor among the \(10{,}501{,}351{,}657\) intervals in \(S_9\). So any type-\(A\) counterexample has rank at least 10, and with your example the smallest possible rank is some \(n\) with \(10\le n\le 14\).

I have an exhaustive \(S_{10}\) search running with the inverse/\(w_0\)-conjugation symmetry quotient. Are you or anyone nearby already enumerating \(S_{10}\)? I would rather not duplicate a multi-hour run. I will follow up either way (a rank-10 example, or the bound \(11\le n\le 14\)).

The 20,163-line \(H_4\) certificate can be regenerated with `./brenti_search --coxeter H4 --out H4.cert`; I am happy to send the file. A compact summary (unique min-length \(Q\)s, orbit counts, Poincaré check) is in `H4_SUMMARY.md` in the repository:

https://github.com/guillochon/2608.13836v1-num

A short write-up is in `WRITEUP.md` there. Please feel free to ignore or reuse any of this.

Best regards,  
James Guillochon  
guillochon@gmail.com
