#!/usr/bin/env python3
"""P6 derivation attempt — forward computation under the pre-registered models M1-M5.
Pre-registration: analysis/p6_derivation_attempt.md, commit 0237940 (frozen before this run).
No parameter is tuned; all constants are the pre-registered values."""

import json

# ---- pre-registered constants (Part I, verbatim) ----
c     = 2.99792458e8          # m/s
hbar  = 1.054571817e-34       # J s
h_J   = 6.62607015e-34        # J s
G     = 6.67430e-11           # m3 kg-1 s-2
g0    = 9.80665               # m/s2 (standard gravity, pre-registered)
m_Sr  = 1.44316e-25           # kg (86.9089 u)
nu_Sr = 429.228e12            # Hz
E_clock_eV = h_J * nu_Sr / 1.602176634e-19   # J -> eV
mc2_GeV = m_Sr * c**2 / (1.602176634e-10)    # J -> GeV
l_p   = 1.616255e-35          # m
E_p_GeV = 1.2209e19           # GeV
R     = 1.0                   # m (clock separation, pre-registered)
d     = 1.0                   # m (M3 pair separation, pre-registered)

TARGET = 2.3e-21

results = {}

# ---- M1: classical GR redshift, h = 1 m ----
M1 = g0 * R / c**2
results['M1'] = {'formula': 'g*h/c^2', 'value': M1,
                 'log10': __import__('math').log10(M1)}

# ---- M2: quantum time dilation bound, E_clock/mc^2 ----
E_clock_GeV = E_clock_eV * 1e-9
M2 = E_clock_GeV / mc2_GeV
results['M2'] = {'formula': 'E_clock/mc^2 (max entanglement-induced modulation)',
                 'value': M2, 'log10': __import__('math').log10(M2)}

# ---- M3: gravitationally mediated entanglement ----
lam_G = G * m_Sr**2 / (hbar * c)                    # dimensionless coupling
phase_rate = G * m_Sr**2 / (hbar * d)               # rad/s at d = 1 m
frac_rate = phase_rate / (2 * 3.141592653589793 * nu_Sr)  # as fraction of clock frequency
results['M3'] = {'lambda_G': lam_G, 'phase_rate_rad_s': phase_rate,
                 'fractional_rate_bound': frac_rate,
                 'log10_bound': __import__('math').log10(frac_rate)}

# ---- M4: Planck-suppressed MDR, xi = 1 ----
E_over_p_rest = mc2_GeV / E_p_GeV
E_over_p_clock = (E_clock_eV * 1e-9) / E_p_GeV
results['M4'] = {
    'linear_rest': E_over_p_rest,
    'quadratic_rest': E_over_p_rest**2,
    'linear_clock': E_over_p_clock,
    'quadratic_clock': E_over_p_clock**2,
}

# ---- M5: the source documents' own mechanism (l_p/R)^2 * S_ent ----
lp2_over_R2 = l_p**2 / R**2
S_defensible = {'O(1)': 1.0, 'Avogadro_scale_1e30': 1e30,
                'area_law_1m2': R**2 / (4 * l_p**2)}
M5 = {k: lp2_over_R2 * v for k, v in S_defensible.items()}
S_required = TARGET / lp2_over_R2
results['M5'] = {'(l_p/R)^2': lp2_over_R2, 'values': M5,
                 'S_ent_required_for_target': S_required}

# ---- retrofit diagnostic (post-hoc, for the report only) ----
lam_from_target = TARGET * c**2 / g0   # the length whose g*lambda/c^2 equals the target
results['retrofit_diagnostic'] = {
    'lambda_um': lam_from_target * 1e6,
    'note': 'g*lambda/c^2 equals the claimed number for lambda = this value',
}

print("=" * 72)
print("P6 FORWARD COMPUTATION (pre-registered models, no tuning)")
print("=" * 72)
print(f"target (claim):            {TARGET:.2e}")
print()
print(f"M1  GR redshift gh/c^2:    {M1:.3e}   (ratio to target: {M1/TARGET:.1e})")
print(f"M2  E_clock/mc^2 bound:    {M2:.3e}   (ratio: {M2/TARGET:.1e})")
print(f"M3  lambda_G:              {lam_G:.3e}")
print(f"    phase rate @1m:        {phase_rate:.3e} rad/s")
print(f"    fractional bound:      {frac_rate:.3e}   (ratio: {frac_rate/TARGET:.1e})")
print(f"M4  linear, rest:          {E_over_p_rest:.3e}   (ratio: {E_over_p_rest/TARGET:.1e})")
print(f"    quadratic, rest:       {E_over_p_rest**2:.3e}   (ratio: {E_over_p_rest**2/TARGET:.1e})")
print(f"    linear, clock:         {E_over_p_clock:.3e}")
print(f"    quadratic, clock:      {E_over_p_clock**2:.3e}")
print(f"M5  (l_p/R)^2:             {lp2_over_R2:.3e}")
for k, v in M5.items():
    print(f"    x S={k:24s} {v:.3e}   (ratio: {v/TARGET:.1e})")
print(f"    S_ent needed for target: {S_required:.3e}")
print()
print(f"RETROFIT DIAGNOSTIC: lambda = {lam_from_target*1e6:.2f} um")
print("  (g*lambda/c^2 = claimed number; compare P1's withdrawn Yukawa range 18-25 um)")
print()

match = any(abs(v - TARGET) / TARGET < 0.05 for v in
            [M1, M2, frac_rate, E_over_p_rest, E_over_p_rest**2,
             E_over_p_clock, E_over_p_clock**2] + list(M5.values()))
print(f"any pre-registered model reproduces the target: {match}")
print(f"OUTCOME: NOT DERIVED (and not specifiable without post-hoc parameters)")

with open('/home/z/my-project/scripts/p6_results.json', 'w') as f:
    json.dump({'target': TARGET, 'results': results,
               'any_model_matches': match,
               'outcome': 'not_specifiable'}, f, indent=2)
print("saved: scripts/p6_results.json")
