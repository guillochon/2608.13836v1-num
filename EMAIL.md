# Draft email to Christian Gaetz

**To:** gaetz@berkeley.edu  
**Subject:** Independent check of arXiv:2608.13836 — \(S_9\) clean, \(H_4\) classified, \(S_{10}\) running

---

Dear Professor Gaetz,

I recently read your note *A counterexample to a log-concavity conjecture of Brenti* (arXiv:2608.13836v1) and wrote a small independent C++/OpenMP implementation of the \(\widetilde{R}\) recurrence, first in type \(A\) and then for the small exceptional/Weyl groups.

Your \(S_{14}\) pair checks out: I recover the same \(Q\)-polynomial
\[
q^{14}+16q^{13}+101q^{12}+333q^{11}+630q^{10}+695q^9+425q^8+123q^7+11q^6+q^5,
\]
and the same failure \(11^2<1\cdot 123\). One minor numerical remark: the relative length is \(\ell(v)-\ell(u)=46-17=29\) rather than 28. That does not affect the counterexample (degree 14 of \(Q\) is compatible with either parity).

I also repeated the exhaustive type-\(A\) search. There are no log-concavity failures among the \(170{,}288{,}585\) Bruhat intervals in \(S_8\), in agreement with your report. Pushing one rank further, there are none among the \(10{,}501{,}351{,}657\) intervals in \(S_9\) either (about 27 minutes on a 16-thread desktop). So any type-\(A\) counterexample has rank at least 10, and with your example the smallest possible rank is some \(n\) with \(10\le n\le 14\).

Your \(H_4\) Example 4 likewise reproduces, with
\[
Q=q^9+14q^8+78q^7+220q^6+326q^5+234q^4+67q^3+8q^2+q
\]
(relative length 18). A complete census of \(H_4\) is small enough to finish: 75,539,433 intervals, of which **20,163** fail log-concavity (1,591 distinct \(Q\)-vectors). All failures violate at index 2, and all remain unimodal with no internal zeros. The shortest failures have relative length **16** (56 intervals, three \(Q\)-vectors, e.g. \([0,1,7,52,124,120,55,12,1]\) with \(7^2<52\)). Your published \(Q\) occurs for 10 intervals. By contrast \(H_3\), \(F_4\), and \(B_2\)–\(B_5\) are entirely clean.

I have an exhaustive \(S_{10}\) search running with the inverse/\(w_0\)-conjugation symmetry quotient on the same machine. I would be glad to know if you or anyone in your orbit is already enumerating \(S_{10}\), so as not to duplicate a multi-hour run; either outcome (a rank-10 example, or the bound \(11\le n\le 14\)) seems publishable.

The code is here if it is of any use:

https://github.com/guillochon/2608.13836v1-num

A short write-up of the counts, the length-29 check, and the \(H_4\) census is in `WRITEUP.md` in that repository. Please feel free to ignore or reuse any of this.

Best regards,  
James Guillochon  
guillochon@gmail.com
