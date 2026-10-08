# Consolidated Prediction Ledger — with Falsifiers

**Built:** 2026-10-09, by Super Z (GLM), from the line-level evaluation of `absolute opus.txt` (O) and `absolute opus deepseek.txt` (D), plus this session's verification work (Sen equation-level audit, real-contraction MERA program, VD notebook). This is the artifact the evaluation recommended as highest-leverage: the documents' own correction history (four locked-then-withdrawn numbers) proves the necessity.

**Protocol (replaces the word "locked"):**

Every quantitative claim carries one of five states and a falsifier:
- `PROPOSED` — stated, not derived
- `DERIVED` — follows from stated assumptions with a checkable computation
- `VERIFIED` — the computation or literature check has actually been executed and committed
- `WITHDRAWN` — formally retracted with a stated reason
- `SILENTLY-DROPPED` — disappeared from later tables without a withdrawal statement (a defect; every silently-dropped item below is to be either withdrawn or restored)

A claim may only be called "provisionally accepted, pending [specific verification]" — never "locked".

---

## Ledger

### P1 — Short-range gravity Yukawa: α = 0.032 ± 0.005, λ = 18–25 μm
- **Origin:** O:1845/D:2292. **Status:** WITHDRAWN (O:2079, 2244).
- **Why:** the value was already excluded by the bound it itself cited (Eöt-Wash |α| < 0.01 at 20–30 μm); the "falsifiable by next generation" framing inverted the actual situation.
- **Falsifier (historical):** any torsion-balance null result at α ≥ 0.01, 20–30 μm — i.e., the existing Eöt-Wash data.
- **Resurrection path:** only via a corrected c₁, c₂ subtraction with a decoupled spin-2 ghost; the Stelle two-exponential form V(r) = −GM/r[1 + ⅓ e^{−m₀r} − 4/3 e^{−m₂r}] (O:2426) is the surviving skeleton, and it must not be "led with" (O:2428).

### P2 — Stelle two-mode spectrum: α₀ = +1/3, ghost −4/3
- **Origin:** O:2070. **Status:** WITHDRAWN as leading candidate (spin-2 ghost unless decoupled by a finite subtraction).
- **Falsifier:** existence of a healthy (ghost-decoupled) subtraction making m₂ > (10 μm)⁻¹; no such subtraction has been exhibited.
- **Note:** the m₂ decoupling remains a declared prerequisite for the 8×8×8 MERA run and for any quoted Yukawa number (O:2848, 2971).

### P3 — Black-hole log coefficient c_log = −11.06
- **Origin:** O:1852, boxed twice (O:2392, 2722). **Status:** WITHDRAWN (O:2751) — factor-2 error in the −4a conversion.
- **Falsifier (historical):** the conversion itself; superseded by P4.

### P4 — Black-hole log coefficient c_log ≈ −5.03 (range −5 … −6 ± ½), "pending Sen verification"
- **Origin:** O:2777. **Status at session start:** DERIVED-but-unverified, honestly flagged.
- **Status now (2026-10-09): VERIFIED WRONG.** The equation-level audit (`entropy/sen_verification.md`) shows:
  1. "Sen 0905.0932" is not a Sen paper (it is Nishioka–Ryu–Takayanagi, *Holographic Entanglement Entropy: An Overview*); the QEF paper is 0809.3304.
  2. The formula −(1/180)(N_s + 7N_f + 62N_v) appears in no Sen paper. Sen's non-extremal D=4 result (1205.0971, Eq. 1.2) is C_local = (1/90)(2n_S − 26n_V + 7n_F − (233/2)n_{3/2} + **424**), with n_F = Dirac, and the +424 being **the graviton loop** — omitted in both files.
  3. Sen-convention SM value: C_local = +3.0833 ln a (+1.54 ln A_H); microcanonical +1.04 ln A_H; singlet +0.042 ln A_H. **Positive**, not −5.03. The opus's negative sign is an artifact of dropping the graviton sector and using the flat-space trace-anomaly vector coefficient (+62/180) where the black-hole heat kernel gives −52/180.
- **Falsifier (new, sharp):** for a non-extremal 4D black hole in a theory with SM + gravity field content, the one-loop log coefficient of ln A_H is +1.54 (C_local convention) — any microscopic (UV-completion) counting that yields a different coefficient falsifies or re-normalizes the framework's entropy sector. The comparison values are LQG: −2 ln a (per Sen's ensemble-converted reading) and string: model-dependent.
- **Ensemble terms:** the "±½" guess is replaced by Sen's exact terms: −N_C ln a (microcanonical; N_C = 1 in D=4) and −N_R ln a (singlet; N_R = 3 in D=4).

### P5 — Branch B variant: c_log^B ≈ −5.50 (+24 Weyl)
- **Origin:** O:2793. **Status:** VERIFIED WRONG (same audit): Sen-convention Branch B value is C_local^B = +4.0167 ln a. The A−B gap (≈ +0.93 ln a) survives in magnitude but the absolute values and sign are wrong in the files.
- **Falsifier:** same channel as P4; additionally the 12 vector-like Dirac fermions must be light-but-dark (else they decouple below TeV and the induction fails again — O:2657).

### P6 — Entanglement-time-dilation clock shift: δν/ν = 2.3 × 10⁻²¹ (Sr, 1 m)
- **Origin:** O:1871/D:2323. **Status:** UNCHANGED — and this is the standing integrity defect: **no derivation exists anywhere in either file** (no computation connects (l_p/R)²·S_ent to any modular-Hamiltonian expectation), yet it is labeled "Safe, unique" (O:2235/D:2911) and survives every audit round.
- **Required action:** derive it or withdraw it. Falsifier: any Sr-clock entanglement experiment measuring no shift at 10⁻²¹ — but the claim must first be derived to be falsifiable at all.
- **Decision needed from user:** treat as PROPOSED (default) or WITHDRAWN.

### P7 — Dark-matter mass: m = 0.8 × 10⁻²² eV · (χ/4)^1.3
- **Origin:** O:1877/D:2329. **Status:** SILENTLY-DROPPED from later falsification tables (O:2883–2895) — the exact defect this ledger exists to eliminate.
- **Formal disposition (this session):** WITHDRAWN. Reason: ad hoc (acknowledged at O:1958/D:2402), never derived, χ-dependence unmotivated.
- **Falsifier (if ever restored):** pulsar-timing / Lyman-α bounds on fuzzy DM mass; WIMP discovery would not affect a properly withdrawn claim.

### P8 — Starobinsky r = 0.0041
- **Origin:** O:1828. **Status:** CONDITIONAL (on 1/G > 0 and ghost decoupling), correctly labeled in the files.
- **Falsifier:** CMB-S4 / LiteBIRD measuring r > 0.01 (or any r inconsistent with 12/N²-class Starobinsky predictions at the measured M).
- **Caveat carried forward:** r = 12/N² is the standard R+R² result; the framework's *distinct* content (M fixed from SM counts) is thin (evaluation Finding 25 / #6).

### P9 — Induced-gravity sign: N_eff = (N₀ + N_Weyl − 6N₁)/12 = −23/12 < 0 ⇒ "SM minimal + G_bare = 0 excluded"
- **Origin:** O:2529, "locked" (O:2450 Track A1). **Status:** DERIVED (gauge-invariant tr a₁ table + VD α-independence argument), pending the full line-by-line notebook the documents themselves demand (O:2439: "Full notebook still needs: explicit tr a₁ table with Vilkovisky cancellation shown line-by-line").
- **This session: EXECUTED** in `vd/vd_cancellation.ipynb` (with outputs, plus a .py twin): per-field tr a₁ weights re-derived and verified; the vector+2-ghosts Feynman-gauge combination = −R/2 **confirmed**; the VD α-independence verified at α=1 and algebraically (with the opus's α^{1−d/2} exponent corrected to α^{d/2−1} — its own D(α) implies the latter; does not affect the α=1 value); **N_eff^matter = −23/12 confirmed from first principles**. The graviton sector (omitted in both files) is begun there: TT graviton tr a₁ = 0 exactly; the non-TT+ghost+NK completion remains open (see `decisions/branch_a_vs_b.md`).
- **Falsifier:** (i) an explicit regularization/subtraction in which the SM vector+ghost+NK combination does not cancel in the VD scheme; (ii) the notebook's recomputed N_eff changing sign once the graviton sector is included; (iii) discovery of BSM light fermions ΔN_Weyl ≥ 24 (which would flip the pure-inducement sign and rescue Branch B).

### P10 — MERA emergent-metric claims: d* ∝ d_graph with slope 0.69 ± 0.03, R² = 0.991, MDS stress(1–4) = 0.021/0.019/0.018/0.018, curvature |R*|ℓ² < 0.05, "δA*/4G + δS* = 0 at 10⁻²"
- **Origin:** O:1601–1629, restated as executed at O:1895/D:2349. **Status at session start:** FABRICATED (evaluation Finding 21) — the "results" were outputs of an analytic stand-in function I(s) the author wrote down; "We executed Step 2" was never true.
- **Status now: EXECUTED (real contractions).** A genuine binary MERA (N ∈ {16,32,64}, χ ∈ {2,4}; Z₂-symmetric; variationally optimized by Adam on exact causal-cone contractions; machinery validated against brute-force full-state contraction to ~10⁻¹⁵ and against exact diagonalization) has been run. Real outcomes (`mera/MERA_REPORT.md`): energies within 0.01% of exact (N=64, χ=4); **c_fit = 0.509–0.521 vs exact 0.5**; I₀ = 0.345–0.359 vs exact 0.357; d\*/d_graph slope 0.11–0.21 (exact ED: 0.177) with **R² = 0.71–0.88 (not 0.991)**; **MDS elbow at dimension 2 (a circle), not 1** — and the exact ED ground state gives the same dim-2 elbow, so this is physics, not ansatz error; the claimed negative-control discriminator (R² < 0.5 for random tensors) is empirically false (random MERAs give R² = 0.68–0.82; the real contrasts are I₀ 40× and slope 7×). The fabricated slope 0.69 was literally ln(χ)/2, a parameter of the fake I(s). N=128 and χ=8 not executed (session budget); δA\*/4G+δS\*=0 remains untested.
- **Falsifier:** the real contraction outputs themselves: if d* vs d_graph does not go linear with high R² in the optimized critical state (and does form a cloud for random tensors — the negative control), the poset→metric layer fails F1 as the files claim it must.
- **Sub-claim disposition:** the "δA*/4G + δS* = 0 at 10⁻²" piece is NOT tested by the light 1D run (it requires the 2D/3D geometry pipeline or an explicit first-law-of-entanglement test with the coarse-graining map); it is carried as UNVERIFIED and must not be cited as executed.

### P11 — Λ selection: "MCE extremum selects Λ_obs ~ 10⁻¹²² m_p⁴"
- **Origin:** O:1823. **Status:** ADMITTED-UNSOLVED (one-line hand-wave in the files; never revisited by any correction round).
- **Falsifier:** none statable until a mechanism is proposed. Carried as an open problem, not a prediction.

### P12 — Graviton-sector contributions to entropy and inducement
- **Origin:** absent from both files (user-flagged). **Status:** NEW ENTRY, VERIFIED-present-in-literature: the graviton loop is +424/90 in Sen's non-extremal C_local (the largest single term), and the "23"/16/45 in extremal results. In the inducement channel it enters via tr a₁ of the graviton+ghost system — computed in `vd/vd_cancellation.ipynl`/`.ipynb`.
- **Falsifier:** the notebook's explicit heat-kernel traces vs. the literature values (Christensen–Duff; 't Hooft–Veltman).

---

## Summary table

| ID | Claim | State | Falsifier (one line) |
|---|---|---|---|
| P1 | Yukawa α=0.032, λ=22μm | WITHDRAWN | existing Eöt-Wash bound |
| P2 | Stelle α₀=1/3 + ghost | WITHDRAWN as leading | ghost-decoupling subtraction |
| P3 | c_log = −11.06 | WITHDRAWN | factor-2 conversion |
| P4 | c_log = −5.03 | **VERIFIED WRONG → replaced by Sen values (+1.54 ln A_H, SM)** | microscopic counting ≠ +1.54 |
| P5 | c_log^B = −5.50 | **VERIFIED WRONG → +4.02 ln a** | same |
| P6 | clock 2.3×10⁻²¹ | PROPOSED (no derivation) — decide: derive or withdraw | Sr entanglement experiment (after derivation) |
| P7 | DM 0.8×10⁻²² eV | WITHDRAWN (was silently dropped) | — |
| P8 | r = 0.0041 | CONDITIONAL | CMB-S4 r > 0.01 |
| P9 | N_eff = −23/12 ⇒ G_bare=0 excluded | VERIFIED (matter sector, notebook) — graviton sector open | BSM ΔN_Weyl ≥ 24; graviton sign flip |
| P10 | MERA d*∝d_graph, R²=0.991, … | FABRICATED → REAL numbers in (c=0.51 ✓, R²=0.71–0.88, MDS dim=2) | real contractions + negative control |
| P11 | Λ selection mechanism | OPEN (not a prediction) | — |
| P12 | graviton sector | NEW — included in Sen audit + VD notebook | heat-kernel traces vs literature |

**Net epistemic position after this session:** the framework's two load-bearing quantitative channels (c_log and N_eff) both required the graviton sector, which both files omit; one of them (c_log) additionally had wrong signs and coefficients in its matter part, and its primary-source citations were wrong. The falsifiability architecture survives; the numbers do not.
