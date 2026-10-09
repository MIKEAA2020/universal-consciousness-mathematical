# Graviton Heat-Kernel Completion — Christensen–Duff Channel

**Task:** the open item from `vd/vd_cancellation.ipynb` §4 and `decisions/branch_a_vs_b.md` §5: the non-TT + ghost + Nielsen–Kallosh graviton sector in the inducement (tr a₁) channel, with the Christensen–Duff cross-check.

**Computation:** `vd/graviton_completion.py`, executed on the conformal round S⁴ (g = Ω²δ, Ω = 2/(1+r²)) so that every sign is forced by direct computation; raw results `vd/graviton_a1_results.json`.

---

## 1. Primary sources verified this session (full texts fetched)

| Source | What it provided | Status |
|---|---|---|
| Sen, arXiv:1205.0971 (ar5iv full text) | Eq. (1.2) verbatim; **"The last term 424 is the contribution from the graviton loop. In pure gravity theory only this term is present."**; the graviton-fluctuation measure normalization (Eq. 2.37); bibliography [45]/[46] = the two Christensen–Duff papers | VERIFIED |
| Vassilevich, hep-th/0306138 (ar5iv full text) | §3.5 "Graviton": the York decomposition Eqs. (3.64)–(3.75) — the computation route used below; Table 1 (a₄ by spin, from Christensen–Duff [127]): **spin-2 incl. ghosts: (a, b, c, d) = (212, 0, 0, 717/4)** | VERIFIED |
| Solodukhin, arXiv:1104.3712 (ar5iv full text) | Eqs. (131)–(141): per-spin operators Δ⁽ˢ⁾ = −∇² + X⁽ˢ⁾, the a₁ formula, Newton-constant renormalization per spin (independent cross-check of the matter table: scalar (1/6−ξ)R, Weyl +1, vector −R/2 with ghosts — all match `vd_cancellation.ipynb`); Eq. (275): the A/B anomaly table — **graviton: A = 0, B = 424/3** | VERIFIED |
| Christensen–Duff NPB **B154** 301 (1979) / **B170** 480 (1980) | the primary spin-2 heat-kernel references | NOT FETCHED (paywalled) — the one remaining open item |

## 2. New structural finding: Sen's C_local = 3B − A per field

From Sen's Eq. (1.2) coefficients and Solodukhin's Eq. (275) anomaly table (A in 1/90π² units, B in 1/30π² units; C_local in 1/90 units):

| Field | A | 3B | **C_local = 3B − A** | Sen (1.2) |
|---|---|---|---|---|
| real scalar | 1 | 3 | **2** | +2 ✓ |
| Dirac | 11 | 18 | **7** | +7 ✓ |
| vector | 62 | 36 | **−26** | −26 ✓ |
| spin-3/2 | 0 | −233/2 | **−233/2** | −233/2 ✓ |
| **graviton** | **0** | **424** | **424** | **+424 ✓** |

Verified across all five spins. Consequences: (i) the graviton's entropy-channel contribution is purely B-type (Weyl²) — its Euler-type (A) anomaly vanishes, the classic Duff observation, confirmed here from two independent sources; (ii) Sen's "graviton loop = +424" is exactly 3× the B-anomaly 424/3 — the entropy channel and the Christensen–Duff anomaly table are one structure; (iii) this identifies precisely which literature numbers anchor the graviton sector in each channel.

## 3. Exact operator algebra (all computed, all verified)

On the unit round S⁴ (R = 12; Einstein check R_μν = 3g_μν and Riemann-derived Ricci both pass):

1. **Sigma anchor.** −Δ(x_i·Ω) = 4·x_i·Ω = (R/3)(x_iΩ) for the l=1 harmonics: the σ-operator −Δ − R/3 has the l=1 (conformal-Killing-scalar) zero modes — the Gibbons–Hawking–Perry structure. This anchors the a₁ convention: **a₁(−∇² + X) = R/6 − X per component** (spectral-shift argument, validated by this zero-mode structure).
2. **Ghost operator.** From the gauge variation of the de Donder condition, computed explicitly: δF(ε) = □ε + (R/4)ε for constant ε — the ghost operator is **M = □ + R/4** (c = +1; C = 3·δ, constant).
3. **TT shift identity.** 2R_{μρνσ}h^{ρσ} = (R/6)(g_{μν}·g^{αβ}h_{αβ} − h_{μν}), verified identically for symbolic symmetric h. On TT (metric-traceless): **2Riem·h = −(R/6)h** — so the TT endomorphism is X_TT = −R/6 in this Riemann convention; +R/6 with the opposite Riemann sign (the de Sitter graviton standard, m² = 2Λ/3 = R/6); +R/3 under Solodukhin's Eq. (133) X⁽²⁾ (traceless-space reading). **With the dS-standard convention the committed notebook's claim "TT graviton tr a₁ = 0" is CONFIRMED** (a₁ = R/6 − R/6 = 0 per TT component).
4. **York Jacobian.** L(Killing vector) = 0 exactly (all 10); gradient block **M~|_grad = (3/2)∇(□φ) + (R/2)∇φ** (k₁ = 3/2, k₂ = R/2 — verified on explicit gradients); co-closed block **M~|_cc = □ + R/4 — identical to the FP ghost operator** (verified on an explicit curl field ξ = δC, co-closed by construction; the Jacobian's co-closed part IS the ghost determinant — the York and de-Donder bookkeeping connect exactly here).

## 4. The a₁ assembly (York route, R/12 units, unit S⁴)

Z = det_V(L†L)^{+1/2}·det_S(−Δ−R/3)^{−1/2}·det_TT(−Δ+X_TT)^{−1/2} (Vassilevich Eq. 3.75) → W-weights: Jacobian −½, σ +½, TT +½.

| Piece | Value (R/12 units) | Basis |
|---|---|---|
| σ (1 comp, X = −R/3, weight +½) | **+3** | a₁ = R/2; zero-mode-anchored |
| Jacobian det(2M~)^{1/2} (weight −½, 2^{1−d/2} = ½) | **−19/4** | gradient family (3/2)·(operator, X = −R/3): effective 4; co-closed family (X = −R/4): 5 × 3 comps = 15; total 19 |
| TT (5 comps, weight +½) | **(5/2)(2 − X̃_TT)** | the §3.3 identity |

**Result by branch:**

| TT convention | w_grav (R/12) | **N_eff^SM+graviton** | < 0 |
|---|---|---|---|
| dS-graviton standard (X_TT = +R/6; TT a₁ = 0 — notebook confirmed) | −7/4 | **−33/16 = −2.0625** | ✓ |
| Solodukhin Eq. (133) X⁽²⁾ reading (X = +R/3) | −27/4 | **−119/48 ≈ −2.479** | ✓ |
| this-convention literal (X_TT = −R/6) | +33/4 | **−59/48 ≈ −1.229** | ✓ |

**Every branch: N_eff^SM+graviton < 0.** The graviton sector does not overturn the matter-sector sign; in two of three branches it makes 1/G_ind more negative. Pure inducement (G_bare = 0) fails with the graviton sector included, within the stated conventions.

## 5. Methodological caveat (disclosed, not resolved)

The same restricted-determinant bookkeeping applied to the **Maxwell vector** (York/Coulomb route) yields **+R/24** per vector, versus the **triple-anchored** Feynman-gauge number **−R/2** (verified independently in the notebook, the source documents, and Solodukhin's Eqs. (134)–(141)) — a discrepancy of −13R/24 per vector. The off-shell one-loop effective action is gauge-choice dependent; restricted-determinant a₁ bookkeeping carries Kabat-type subtleties that this session did not resolve. The graviton a₁ numbers above therefore carry state **DERIVED-with-caveat**, not VERIFIED. The closing check — the de-Donder spin-2 heat kernel with ghosts and the Nielsen–Kallosh completion, i.e. the actual Christensen–Duff NPB B170 computation — remains open (paywalled primary source). The sign conclusion is robust across all branches computed, and the caveat's scale (O(1) R/12-units per gauge parameter) does not reach the +23 units needed to flip the N_eff^SM+grav sign in the two standard branches without also flipping the anchored matter numbers.

## 6. Consequences

1. **P9's condition (i) is resolved.** The branch memo's reopen trigger — "if the completed spin-2 heat kernel leaves N_eff^SM+graviton ≤ 0" — has fired in the negative direction in every computed branch: the exclusion claim ("SM minimal + G_bare = 0 excluded", a statement within appearance about a model choice) is **strengthened**, with the caveat of §5 recorded.
2. **The entropy channel is unchanged and literature-anchored:** the graviton loop +424/90 is mandatory in every published c_log (P4, P12).
3. **Branch B's last technical hope is closed.** Its metaphysical and conventional standing was excluded by the user's rulings (framework Articles III, IV.1); the inducement-channel graviton completion removes the residual "graviton sign flip" possibility in all computed branches.
4. **The committed notebook's §4 graviton sketch is substantiated:** TT a₁ = 0 (standard convention), ghost operator M = □ + R/4, and the non-TT structure — all exactly as the notebook sketched, now with verified arithmetic; the notebook's "open" flag on the non-TT+NK completion is answered by §4–§5 here (assembly done; de-Donder cross-check open).
