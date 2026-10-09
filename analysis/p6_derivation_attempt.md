# P6 — Pre-Registered Derivation Attempt

**Protocol source:** user ruling, 2026-10-09, verbatim:

> "P6: attempt derivation under pre-specified assumptions, or withdraw.
> Rules:
> 1. Specify the model first: Lagrangian, cutoff, coupling, clock type. Commit to these before computing.
> 2. Compute forward. Do not tune any parameter to reach 2.3×10⁻²¹.
> 3. If the number comes out, report the model in full. Reword the ledger entry as a model prediction, not a framework prediction.
> 4. If the number does not come out, withdraw P6. Do not adjust the model to rescue the number.
> 5. If no model can be specified without reverse-engineering, withdraw P6 immediately. That is the hallucination case.
> Report which of the three outcomes occurred: derived, not derived, or not specifiable."

**The claim under test** (ledger P6; origin O:1871/D:2323): an entanglement-time-dilation clock shift δν/ν = 2.3 × 10⁻²¹ for strontium clocks at 1 m separation, labeled "Safe, unique" in the source documents. No derivation exists in either source file: no computation connects (l_p/R)²·S_ent to any modular-Hamiltonian expectation.

**Integrity mechanism:** this Part I is committed to git **before** the computation in Part II is run. Any post-hoc change to the models below would be visible in the commit history. No parameter specified in Part I may be adjusted after computation begins.

---

## Part I — Pre-specified models (committed before computation)

Common experimental setting (fixed for all models): two identical ⁸⁷Sr optical-lattice clocks, transition ¹S₀–³P₀, frequency ν_Sr = 429.228 THz; vertical separation h = 1 m; Earth's surface, standard gravity g = 9.80665 m/s². Fixed constants: c = 2.99792458×10⁸ m/s, ħ = 1.054571817×10⁻³⁴ J·s, G = 6.67430×10⁻¹¹ m³kg⁻¹s⁻², m(⁸⁷Sr) = 86.9089 u = 1.44316×10⁻²⁵ kg (mc² = 81.04 GeV), E_clock = h ν_Sr = 1.776 eV, l_p = 1.616255×10⁻³⁵ m, E_p = 1.2209×10¹⁹ GeV.

Model selection criterion (fixed before computation): the five classes below are the canonical mechanisms in the existing literature for clock-frequency shifts or clock-state modulations — classical gravitational redshift (the mandatory baseline), quantum time dilation of composite clocks, gravitationally mediated entanglement, Planck-suppressed modified dispersion, and the source documents' own proposed mechanism. Selection is by canonicity, not by proximity to any target value.

### M1 — Classical gravitational redshift (baseline; general relativity)
- **Lagrangian:** Einstein–Hilbert + matter; clocks modeled as two-level systems at rest at heights differing by h.
- **Cutoff:** none (classical GR; the effect is infrared).
- **Coupling:** gravitational potential only.
- **Clock type:** ⁸⁷Sr optical transition.
- **Prediction to compute:** δν/ν = g h / c².
- **Role:** the reference magnitude; any distinct mechanism must be a correction relative to this, with its own independent scale.

### M2 — Quantum time dilation of an entangled composite clock (Smith–Ahmadi formalism)
- **Lagrangian:** relativistic quantum mechanics of the clock atom: Dirac Hamiltonian + internal two-level Hamiltonian H_int = E_clock |e⟩⟨e|.
- **Coupling:** special-relativistic coupling of internal energy to the proper-time operator; corrections bounded by internal-energy over rest-energy ratios for the occupied internal states.
- **Cutoff:** none (effect is of relative order ≤ E_int/mc²; no UV sensitivity).
- **Clock type:** ⁸⁷Sr, internal superposition/entanglement between ground and excited states.
- **Quantity to compute:** the maximal possible entanglement-induced fractional shift, bounded by E_clock/mc² for the two-level system; and the general state-dependence structure (corrections vanish for product/classical states, are bounded by the internal-energy scale otherwise).

### M3 — Gravitationally mediated entanglement between the two clocks (Bose et al. 2017 / Marletto–Vedral 2017 class)
- **Lagrangian:** non-relativistic quantum mechanics of two Sr atoms with Newtonian pair potential V = −G m²/d, d = 1 m.
- **Coupling:** Newton's G.
- **Cutoff:** none (non-relativistic; no UV sensitivity).
- **Clock type:** the atoms' internal transitions as clocks, entangled via the gravitational interaction.
- **Quantity to compute:** the dimensionless gravitational coupling λ_G = G m²/(ħc); the phase accumulation rate G m²/(ħ d) at d = 1 m; the induced fractional clock modulation implied by either.

### M4 — Planck-suppressed modified dispersion relation
- **Lagrangian/model:** modified dispersion E² = m²c⁴ + p²c²[1 + ξ (E/E_p)ⁿ⁻²], n ∈ {1,2}, with the dimensionless coefficient ξ fixed a priori at ξ = 1 (no tuning permitted; ξ is the naturalness parameter and the pre-registered value is 1).
- **Cutoff:** E_p (the Planck energy is the suppression scale).
- **Coupling:** ξ = 1.
- **Clock type:** ⁸⁷Sr transition frequency as affected by the dispersion modification; two energies to compute: the rest energy mc² (generous bound) and the transition energy E_clock (the relevant one).
- **Quantity to compute:** δν/ν = (E/E_p)ⁿ for the two energies and both n.

### M5 — The source documents' own mechanism: δν/ν = (l_p/R)² · S_ent
- **Model:** as stated at O:1871: the shift is (l_p/R)² times an entanglement entropy, R = 1 m.
- **Required specification the source never provides:** the system whose entanglement entropy S_ent is. Pre-registration fixes the defensible range: S_ent = O(1) (few-qubit clock entanglement); S_ent ≤ 10³⁰ (Avogadro-scale many-body entropy of a 1 m laboratory system); S_ent = A/(4 l_p²) (holographic area bound for a 1 m² region).
- **Quantity to compute:** (l_p/R)²·S_ent across that range; and the value of S_ent that would be required to produce 2.3×10⁻²¹.

### Decision rule (fixed before computation)
- If any model M1–M5 produces 2.3×10⁻²¹ (within the significant digits stated in the claim) with its pre-registered parameters, outcome = **derived**; the ledger entry is reworded as a prediction of that model, stated in full.
- If the models produce well-defined numbers different from 2.3×10⁻²¹, outcome = **not derived**; P6 is withdrawn (rule 4).
- If producing 2.3×10⁻²¹ requires introducing or adjusting any parameter not fixed above, outcome = **not specifiable**; P6 is withdrawn immediately (rule 5).

*Part II (computation and report) follows in a separate commit after this one. Part I is frozen from this commit onward.*
