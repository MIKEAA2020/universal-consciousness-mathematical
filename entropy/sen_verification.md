# Sen Equation-Level Verification — c_log, Field Content, and the Graviton Sector

**Verifiers:** Super Z (GLM), session 2026-10-09. Online verification against primary sources, per the open item left in the opus files ("Exact value pending verification against Sen primary sources", O:2869/2935) and the user's instruction to include the graviton sector, which both files omit.

**Sources actually consulted (fetched full text unless noted):**

| Source | What it actually is | Fetched |
|---|---|---|
| arXiv:1205.0971 (A. Sen) | *Logarithmic Corrections to Schwarzschild and Other Non-extremal Black Hole Entropy in Different Dimensions* | full text (ar5iv) |
| arXiv:1005.3044 (S. Banerjee, R.K. Gupta, A. Sen) | *Logarithmic Corrections to Extremal Black Hole Entropy from Quantum Entropy Function* (N=4 matter multiplets) | full text (ar5iv) |
| arXiv:0809.3304 (A. Sen) | *Quantum Entropy Function from AdS(2)/CFT(1) Correspondence* | abstract |
| arXiv:1109.3706 (R.K. Gupta, A. Sen) | *Logarithmic Corrections to Rotating Extremal Black Hole Entropy in Four and Five Dimensions* | abstract |
| arXiv:0905.0932 (T. Nishioka, S. Ryu, T. Takayanagi) | *Holographic Entanglement Entropy: An Overview* — **not a Sen paper at all** | arXiv page |
| A. Sen, ICTP Spring School 2012 slides, "Extremal Black Hole Entropy" | summary table of extremal log-correction results | full PDF |

---

## 1. Citation verdicts

**V1 — "Sen 0905.0932" does not exist as cited.** arXiv:0905.0932 is *Holographic Entanglement Entropy: An Overview* by **Nishioka, Ryu and Takayanagi**. It is not Sen's quantum entropy function paper, and it contains none of the equations the opus attributes to it ("Eq. 2.10–2.15, 4.1–4.5, Table 1/2, 5.10–5.15"). Sen's quantum entropy function paper is **arXiv:0809.3304**. Every "Sen 0905.0932" reference in both files (O:2417, O:2759, O:2875, O:2938, D:5098, D:5112 ff.) is a wrong citation. This is the same error class as the "Washburn 2023" hallucination (evaluation Finding 22).

**V2 — "Sen 1205.0971" is cited correctly.** Title as above; it is the non-extremal paper, and it is indeed the relevant one for S = A/4G + c_log log A — but the equation numbers the opus lists ("Eq. 3.1–3.8", "Sec. 4 Dirac vs Weyl examples") do not contain the formula the opus claims. The formula actually appears as **Eq. (1.2)** (and equivalently (2.34)-derived form), quoted next.

**V3 — The extremal companion literature is a different set of papers.** The extremal log-correction results the opus's checklist gestures at live in 1005.3044 (matter multiplets), 1109.3706 (rotating), and further papers cited therein — not in 0905.0932.

---

## 2. What Sen's papers actually say (equation level)

### 2.1 Non-extremal black holes (arXiv:1205.0971)

**Eq. (1.1) (microcanonical):**
S_mc(M, J, Q) = S_BH(M, J, Q) + ln a · ( C_local − ½(D−4) − ½(D−2) N_C − ½(D−4) n_V )

**Eq. (1.2) (field content, D = 4, uncharged):**

> C_local = (1/90) · ( 2 n_S − 26 n_V + 7 n_F − (233/2) n_{3/2} + 424 )

with the paper's own words: *"if the theory contains, besides gravity, n_S massless scalar fields, n_V massless vector fields, n_F massless **Dirac** fields and n_{3/2} massless spin 3/2 fields, all minimally coupled to gravity without any other interactions. **The last term 424 is the contribution from the graviton loop.** In pure gravity theory only this term is present."*

**Eq. (1.5) (singlet ensemble):** ΔS_singlet = ln a · ( C_local − ½(D−4) − ½(D−2) N_R − ½(D−4) n_V ), with N_R = (D−1)(D−2)/2 rotation generators, N_C = [(D−1)/2] Cartan generators, and 'a' the black-hole size parameter, A_H ~ a^{D−2} (so ln a = ½ ln A_H up to constants).

**Pure gravity (D=4), from the paper:** C_local = 212/45, hence ΔS_singlet = (212/45 − 3) ln a = (77/45) ln a. The LQG prediction after the same ensemble conversion is −2 ln a. Sen's abstract: *"For Schwarzschild black holes in four space-time dimensions the macroscopic result seems to disagree with the existing result in loop quantum gravity."*

### 2.2 Extremal black holes (arXiv:1005.3044 and Sen's 2012 summary table)

**Eq. (2.22) of 1005.3044 (a single massless scalar on AdS₂×S²):**

> ΔS_BH = −(1/180) · ln(a²/ε)

This is the true origin of the "1/180" prefactor the opus uses — it is the *extremal, per-scalar* coefficient. The paper then shows the N=4 matter-multiplet contributions **cancel** exactly, and explicitly defers the gravitational sector: *"While a similar calculation is possible in principle for the fields in the gravitational sector, the computation is technically involved, and we have not [carried it out]"*.

**Sen's own summary table (1205.0971 Eq.-table and ICTP 2012 slides), extremal, D=4, A_H ~ Λ²:**

| Theory | log correction |
|---|---|
| N=4 with n_v matter multiplets | 0 |
| N=8 / type II on T⁶ | −8 ln Λ |
| N=2 with n_V vector, n_H hyper multiplets | **(1/6)(23 + n_H − n_V) ln Λ** |
| N=6 | −4 ln Λ |
| N=5 | −2 ln Λ |
| N=3 | +2 ln Λ |
| BMPV (5D, J=0) | −(1/4)(n_V + 3) ln Λ |
| **Extremal Kerr (pure gravity)** | **16/45 ln A_H** (from 1109.3706) |

The "23" in the N=2 formula and the 16/45 for extremal Kerr are gravitational-sector contributions — i.e., in the extremal program the graviton sector is **included** in the tested results, and it is not zero.

---

## 3. Verdict on the opus's c_log formula

The opus's final pending state (O:2921): S_log = −(1/180)(N_s + α N_f + 62 N_v) log A, α = 7 (Dirac) or 7/2 (Weyl), ensemble ±½, candidate values −5.03 / −5.91.

**Finding — this formula appears in no Sen paper.** It is a hybrid that does not match either the extremal or the non-extremal literature results:

| Field | opus weight (per field, in 1/180 units, sign as used) | Sen non-extremal (per field, in 1/90 units) | Sen extremal (per field) |
|---|---|---|---|
| real scalar | −1/180 | **+2/90 = +4/180** | −1/180 per ln(a²) (1005.3044 Eq. 2.22) |
| Dirac fermion | −7/180 | **+7/90 = +14/180** | (not quoted per-Dirac in fetched sources; constrained by N=4 cancellation) |
| vector | **+62/180** | **−26/90 = −52/180** | (N=4 multiplet cancels; per-vector not isolated here) |
| spin-3/2 | absent | −(233/2)/90 | — |
| **graviton + ghosts** | **ABSENT** | **+424/90 = +848/180** | 23-term / 16/45 (extremal Kerr) |
| ensemble terms | "±½" | −N_C, −N_R, −n_V/2-type, exactly enumerated | (different structure) |

Three independent structural errors: (i) the vector coefficient has the **wrong sign and wrong magnitude** (+62/180 vs Sen's −52/180 non-extremal; the opus's 62 is the flat-space trace-anomaly a₂ coefficient mis-imported into a black-hole heat-kernel computation); (ii) the scalar coefficient is the extremal one (−1/180) pasted into a formula otherwise resembling nothing in the literature; (iii) the graviton sector — the single largest term in the non-extremal result (+424/90 = +4.71, larger than the entire SM contribution in magnitude) — is missing entirely, in both files, in every variant. The sign of the opus's total is negative; Sen's SM total is positive.

---

## 4. Corrected numbers (Sen conventions, minimally coupled, uncharged D=4)

Using Sen Eq. (1.2) with the opus's own field counts (SM: n_S = 4 real Higgs d.o.f., n_F = 22.5 Dirac = 45 Weyl, n_V = 12, n_{3/2} = 0):

**Branch A (SM, non-extremal, Sen Eq. 1.2):**
C_local = (2·4 + 7·22.5 − 26·12 + 424)/90 = (8 + 157.5 − 312 + 424)/90 = **277.5/90 = +3.0833**

- microcanonical (D=4, uncharged): ΔS = (C_local − N_C) ln a = 2.0833 ln a ≈ **+1.042 ln A_H**
- singlet: (C_local − 3) ln a = 0.0833 ln a ≈ **+0.042 ln A_H**
- graviton-only part: +424/90 = +4.711 (in C_local units)

**Branch B (SM + 24 Weyl = +12 Dirac vector-like, n_F = 34.5):**
C_local = (8 + 241.5 − 312 + 424)/90 = 361.5/90 = **+4.0167**; microcanonical ≈ +1.508 ln A_H.

**Comparison with the opus's pending values:**

| Quantity | opus (pending) | Sen-verified |
|---|---|---|
| c_log (SM) | −5.03 (or −5.91) | **+3.08 ln a = +1.54 ln A_H (C_local); +1.04 ln A_H (microcanonical); +0.042 ln A_H (singlet)** |
| c_log (Branch B) | −5.50 (or −6.38) | **+4.02 ln a** |
| ensemble shift | "±½" | −N_C = −1 (mc) / −N_R = −3 (singlet) in ln a units, D=4 |
| sign vs LQG | "opposite to LQG's −0.5/−3" | Sen: SM+gravity **+3.08 ln a** vs LQG −2 ln a — still opposite, but the opus's own number is on the wrong side of zero |
| graviton sector | omitted | +424/90 non-extremal; the "23" (N=2) / 16/45 (extremal Kerr) extremal |

Caveats stated by Sen himself: Eq. (1.2) assumes minimal coupling and no other interactions; the charged (Kerr–Newman) variant is Eq. (2.34); non-Abelian gauge multiplets and charged fermions (the actual SM) require care beyond the strictly minimal case. These caveats shift details, not the three structural verdicts above.

---

## 5. The graviton sector (the user's flagged omission) — consequences

1. **Non-extremal:** the graviton loop contributes +424/90 = +4.711 to C_local — the largest single coefficient in Eq. (1.2). Omitting it (as both files do in every variant) is not a small correction; it exceeds the entire SM matter contribution in magnitude (SM matter net: 277.5 − 424 = −146.5 → −1.628/90... i.e., SM matter alone gives C_local = −146.5/90 = −1.628, and the graviton flips the total to +3.083). The opus's negative sign is exactly what one gets by dropping the graviton and using wrong vector signs — the sign of its "load-bearing" result is an artifact of the omission.
2. **Extremal:** the tested extremal results (N=2's "23", extremal Kerr's 16/45, N=8's −8) all include the gravitational multiplet. Any MCE/quantum-entropy-function-style claim about extremal BH log corrections must state its graviton sector; the opus never does.
3. **For the Branch A/B decision:** the induced-gravity/VD analysis (separate notebook) and this entropy analysis are the two places the graviton sector enters. In the entropy channel, including it changes the SM c_log by +4.71 (ln a units) and Branch B by the same amount — the A-vs-B gap (3.08 vs 4.02) is unchanged by the graviton (it adds equally), but the *absolute* number and its *sign* versus LQG/string comparisons change qualitatively.

---

## 6. Resolution of the opus's open checklist (O:2873–2895, D:5095–5181)

| Checklist item | Resolution |
|---|---|
| Sen 1205.0971 Eq. 3.1–3.8 "C = 1/180(N_s+7N_f+62N_v)" | **Not present.** The actual formula is Eq. (1.2): (1/90)(2n_S − 26n_V + 7n_F − (233/2)n_{3/2} + 424). |
| Dirac vs Weyl convention | **Dirac** (n_F counts Dirac fields — verbatim from 1205.0971). The opus guessed this correctly. |
| Prefactor 1/180 vs 1/360 | Neither, for the non-extremal total: **1/90** with the coefficients above. The 1/180 is the *extremal per-scalar* coefficient (1005.3044 Eq. 2.22). |
| Ensemble ±½ (Sen "Eq. 5.10–5.15") | The cited equations are in the wrong paper. Sen's actual ensemble terms: microcanonical −N_C ln a, singlet −N_R ln a (D=4: −1 and −3 in ln a units), plus −½(D−4)n_V (zero in D=4). |
| "−5.03 vs −5.91" | **Neither.** Sen-convention SM values are positive: +3.08 ln a (C_local), +2.08 ln a (mc), +0.083 ln a (singlet). |
| Solodukhin cross-check (a,c vs Euler/Weyl basis) | Not fetched this session. Structural note: the trace-anomaly/entanglement route (Solodukhin's program) and Sen's partition-function route differ by the temperature-derivative and zero-mode terms; the opus's c_log = −4a(4π)² and c = 283/120 mixing of the two channels is not a valid identity in either. |

**Bottom line:** the "c_log ≈ −5.03, still distinct from LQG/string" line that the opus left pending is unsalvageable as stated. Under Sen's actual conventions the SM (+gravity) non-extremal log coefficient is **positive** (+3.08 ln a ≈ +1.54 ln A_H before ensemble terms), the graviton sector is mandatory and dominant, and the LQG comparison remains a disagreement but with both signs and magnitudes different from what the files claim. The falsifiability structure survives (the coefficient is still a sharp, low-energy-determined number), but every number and the sign must be replaced before any ET/CE note or "lock".

---

## Addendum (2026-10-09, later session): the C_local = 3B − A identification and the graviton-sector anomaly structure

Cross-checking Sen's Eq. (1.2) coefficients against the anomaly table in Solodukhin's review (arXiv:1104.3712, Eq. (275), from Duff/Birrell–Davies; A in 1/90π² units, B in 1/30π² units) yields an exact identification, verified across all five spins:

**C_local = 3B − A per field** (in 1/90 units).

| Field | A | 3B | C_local = 3B − A | Sen (1.2) |
|---|---|---|---|---|
| real scalar | 1 | 3 | 2 | +2 ✓ |
| Dirac | 11 | 18 | 7 | +7 ✓ |
| vector | 62 | 36 | −26 | −26 ✓ |
| spin-3/2 | 0 | −233/2 | −233/2 | −233/2 ✓ |
| graviton | 0 | 424 | 424 | +424 ✓ |

Consequences: (1) the graviton's entropy-channel contribution is purely B-type (Weyl²) — its A-type (Euler) anomaly vanishes (Duff's observation, confirmed from two independent sources: Solodukhin Eq. (275) and, at the a₄-operator level, Vassilevich's Table 1 spin-2-with-ghosts row (212, 0, 0, 717/4) derived from Christensen–Duff NPB B154); (2) Sen's +424 = 3 × (424/3): the entropy channel and the Christensen–Duff anomaly table are one structure; (3) the graviton sector is now literature-anchored in both channels — see `vd/graviton_completion.md` for the inducement-channel computation. Sources fetched this session: Sen 1205.0971 full text (Eq. 1.2 and the "graviton loop" attribution verbatim; CD papers located in his bibliography as [45]/[46]); Vassilevich hep-th/0306138 full text (§3.5 York decomposition; Table 1); Solodukhin 1104.3712 full text (Eqs. 131–141, 275). Christensen–Duff NPB B154/B170 themselves remain unfetched (paywalled).
