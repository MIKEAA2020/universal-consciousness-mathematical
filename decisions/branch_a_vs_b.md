# Branch A vs. Branch B — Decision Memo

**Prepared:** 2026-10-09, after the Sen equation-level verification (`entropy/sen_verification.md`), the real-contraction MERA program (`mera/MERA_REPORT.md`), and the full VD notebook (`vd/vd_cancellation.ipynb`). The user asked: *"we decide Branch A vs B."* This memo lays the evidence on the table with the metaphysical commitments visible, as the documents themselves prescribe (O:2679–2683: *"You must decide which axiom is non-negotiable: is pure inducement part of the metaphysics, or was it an extra assumption?"*).

---

## 1. What the branches are

| | **Branch A** | **Branch B** |
|---|---|---|
| Posture | R0 only; $G_{\rm bare} > 0$ allowed as a counterterm inside App(A) | Pure Sakharov inducement: $G_{\rm bare} = 0$ elevated to a non-negotiable postulate |
| New fields | none | $\Delta N_{\rm Weyl} \ge 24$ (12 vector-like Dirac), light but dark |
| Matter-only $N_{\rm eff}$ | $-23/12$ (sign irrelevant: $G_{\rm bare}$ absorbs it) | must be $> 0$; $+1/12$ at the minimum $\Delta N_{\rm Weyl} = 24$ |
| c_log (Sen conventions, non-extremal) | $C_{\rm local}^{\rm SM} = +3.083\,\ln a$ | $C_{\rm local}^{\rm SM+24W} = +4.017\,\ln a$ |
| R0 status | $G_{\rm bare}=0$ is *not* a consequence of R0 (the R0 Clarification Memo, O:2696) | overrides the memo; R0 is re-scoped |

## 2. What this session's three computations change

**(a) The Sen audit removes Branch B's advertised numerical payoff.** The opus's motivation for Branch B was a distinct falsifiable number: "c_log^B = −12.53" (later "−5.50") vs Branch A's "−11.06" (later "−5.03"). Under Sen's actual conventions both numbers are wrong in sign and magnitude; the corrected pair is **+3.08 vs +4.02 ln a** — the A/B gap survives (+0.93 ln a ≈ +0.47 ln A_H) but is produced by 12 new fermions' heat-kernel traces, and both values now sit on the *same* side (positive) of every comparison the files made (LQG −2 ln a, "string 0"). Branch B's selling point — "a sharper, more falsifiable number" — is intact in structure but every published candidate number is dead.

**(b) The VD notebook confirms the matter-sector exclusion but withdraws its "locked" status.** $N_{\rm eff}^{\rm matter} = -23/12$ is now actually derived line-by-line (per-field weights re-derived, VD $\alpha$-independence verified at $\alpha=1$ and algebraically; the opus's own $\alpha^{1-d/2}$ exponent is inconsistent with its own $D(\alpha)$ — corrected to $\alpha^{d/2-1}$ in the notebook, which does not affect the $\alpha=1$ value). But the graviton sector — omitted in both files, in every variant — contributes to the same channel: the TT part derives to exactly zero, while the non-TT + ghost + Nielsen–Kallosh parts do not vanish. Until the full spin-2 heat kernel is done (Christensen–Duff cross-check), the honest statement is *"SM matter + G_bare=0 is excluded unless the graviton sector overturns the sign"* — not "gauge invariant, locked".

**(c) The graviton sector is mandatory in the entropy channel.** Sen's Eq. (1.2) carries +424/90 from the graviton loop — the largest single term. Omitting it (as both files do) flips the sign of the SM c_log. Whichever branch is chosen, its published entropy numbers must include it.

## 3. The decision axes (as the documents frame them)

1. **Is pure inducement part of the intuition?** The R0 Clarification Memo (O:2696–2700, one of the documents' genuinely strong artifacts) answers: R0 forbids $A$ having $R, G, x, t$; it does not forbid App(A) having a counterterm $G_{\rm bare}$. Making $G_{\rm bare}=0$ non-negotiable is an *additional aesthetic postulate* — Lakatos-wise, an extra belt element, not hard core.
2. **Generativity vs. minimality.** Branch B is more generative (predicts 12 dark vector-like Dirac fermions — in principle testable) but pays: anomaly-freedom constraints, asymptotic-freedom protection (they must be heavier than ~TeV to not spoil AF, hence decoupled from low-energy inducement — the documents' own trap, O:2657), and now the graviton-sector caveat on the very sign that motivated the branch.
3. **Track record.** The correction cascade's one "locked-then-survived" quantitative result is the matter-sector sign — now verified but with its lock withdrawn pending the graviton completion. Every other locked number died (see `analysis/prediction_ledger.md`). Prior odds favor the branch that requires fewer new commitments per falsifiable number.

## 4. Recommendation

**Hold Branch A as the working frame; do not lock Branch B.** Reasons:

- Branch A survives all three audits structurally: R0 is untouched, no new fields are forced, and its (corrected) falsifiable content — c_log = +3.08 ln a (SM, non-extremal, Sen conventions) with the graviton included, r conditional, the MERA poset→metric layer with its *real* numbers — remains testable in principle.
- Branch B's distinctive commitments (12 dark Dirac fermions *and* the sign-flip requirement *and* the AF/decoupling tightrope) now rest on a channel whose graviton part is uncomputed. Adopting it now means locking a sign before its own arithmetic is complete — precisely the failure mode the documents' history punishes (four locked numbers, four withdrawals).
- The decision the documents deferred ("which axiom is non-negotiable") is, on the evidence, not forced: nothing derived this session *requires* $G_{\rm bare}=0$; the R0 memo stands.

**Condition to revisit:** Branch B becomes the rational choice if (i) the completed spin-2 heat kernel leaves $N_{\rm eff}^{\rm SM+graviton} \le 0$ (pure inducement then *fails* — note this would strengthen, not weaken, the exclusion claim while killing Branch B), or (ii) the user regards pure inducement as intrinsic to the metaphysical intuition (a frame-level commitment that no computation can settle — but then it should be stated as a postulate in the axioms, not smuggled in as a consequence of R0).

**Either way, the next irreversible steps are the same:** finish the graviton sector of the VD notebook (Christensen–Duff cross-check), replace every c_log in the documents with the Sen-convention values including +424/90, and re-run the falsification table in `analysis/prediction_ledger.md` before any number is called a prediction.
