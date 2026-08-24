# Draft email to Christian Gaetz

**To:** gaetz@berkeley.edu  
**Subject:** Independent check of arXiv:2608.13836 — no type-A counterexample through \(S_9\)

---

Dear Professor Gaetz,

I recently read your note *A counterexample to a log-concavity conjecture of Brenti* (arXiv:2608.13836v1) and wrote a small independent C++/OpenMP implementation of the type-\(A\) \(\widetilde{R}\) recurrence, with the aim of reproducing the \(S_{14}\) example and seeing how far an exhaustive search would go.

Your \(S_{14}\) pair checks out: I recover the same \(Q\)-polynomial
\[
q^{14}+16q^{13}+101q^{12}+333q^{11}+630q^{10}+695q^9+425q^8+123q^7+11q^6+q^5,
\]
and the same failure \(11^2<1\cdot 123\). One minor numerical remark: the relative length is \(\ell(v)-\ell(u)=46-17=29\) rather than 28. That does not affect the counterexample (degree 14 of \(Q\) is compatible with either parity).

I also repeated the exhaustive type-\(A\) search. There are no log-concavity failures among the \(170{,}288{,}585\) Bruhat intervals in \(S_8\), in agreement with your report. Pushing one rank further, there are none among the \(10{,}501{,}351{,}657\) intervals in \(S_9\) either (about 27 minutes on a 16-thread desktop). So any type-\(A\) counterexample has rank at least 10, and with your example the smallest possible rank is some \(n\) with \(10\le n\le 14\).

Directed searches in \(S_{10}\)–\(S_{13}\) (all proper patterns of your pair, inversion-table mimics, and several hundred thousand random reduced intervals per rank) did not turn up another example; your \(S_{14}\) pair also appears to be pattern-minimal. I did not attempt a full enumeration of \(S_{10}\).

The code is here if it is of any use:

https://github.com/guillochon/2608.13836v1-num

A short write-up of the counts and the length-29 check is in `WRITEUP.md` in that repository. Please feel free to ignore or reuse any of this; I mainly wanted to let you know that the \(S_{14}\) polynomial independently reproduces and that \(S_9\) is clean.

Best regards,  
James Guillochon  
guillochon@gmail.com
