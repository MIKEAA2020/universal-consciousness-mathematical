#!/usr/bin/env python3
"""Graviton heat-kernel completion — v3 (corrected, lean).

Route: York decomposition (Vassilevich hep-th/0306138 Sec 3.5, Eqs 3.64-3.75) of the
one-loop graviton path integral on an Einstein background:

    Z = det_V(L+L)^{1/2} * det_S(-D - R/3)^{-1/2} * det_TT(-D gg + 2 Riem)^{-1/2}

All signs computed explicitly on the conformal round S^4 (g = Om^2 delta, Om = 2/(1+r^2)).

Convention (settled by the spectral-shift argument): for D = -nabla^2 + X (constant
endomorphism X), per component a1 = R/6 - X. Verified against the sigma anchor:
- Delta(l=1 harmonic) = 4 = R/3  =>  the sigma operator -nabla^2 - R/3 has the
  l=1 zero modes (Gibbons-Hawking-Perry structure)  =>  a1_sigma = R/6 + R/3 = R/2.

Expected values (hand-derived, to be confirmed):
  ghost operator M = box + R/4            (c = +1)
  Jacobian gradient block  M~|grad = (3/2) grad(box phi) + (R/2) grad phi
  Jacobian co-closed block M~|cc   = box + R/4   (= the ghost operator)
  TT identity: 2 R_{mu rho nu sigma} h^{rho sigma} = 0 for symmetric h (max-sym)
"""
import sympy as sp
import json

x0, x1, x2, x3 = sp.symbols('x0 x1 x2 x3', real=True)
xs = [x0, x1, x2, x3]
r2 = x0**2 + x1**2 + x2**2 + x3**2
u = 1 + r2
Om = 2/u
g = sp.diag(*[Om**2]*4)
ginv = sp.diag(*[1/Om**2]*4)
zero = sp.Integer(0)
R_ = sp.Integer(12)          # unit round S^4 scalar curvature (verified below)

def simp(e):
    return sp.cancel(sp.expand(e))

# ---------- Christoffel ----------
dln = [sp.diff(sp.log(Om), xv) for xv in xs]      # = -2 x_i / u
Gamma = [[[zero]*4 for _ in range(4)] for _ in range(4)]
for r in range(4):
    for m in range(4):
        for n in range(4):
            Gamma[r][m][n] = simp((sp.Integer(1) if r == m else 0)*dln[n] +
                                  (sp.Integer(1) if r == n else 0)*dln[m] -
                                  (sp.Integer(1) if m == n else 0)*dln[r])

# ---------- covariant machinery ----------
def grad_cov(mu, comps):
    out = []
    for a in range(4):
        e = sp.diff(comps[a], xs[mu])
        for s in range(4):
            if comps[s] != 0:
                e -= Gamma[s][mu][a]*comps[s]
        out.append(sp.cancel(sp.expand(e)))
    return out

def laplacian_vec(comps):
    """g^mn nabla_m nabla_n on a covariant vector field — BOTH connection terms."""
    T = [grad_cov(n, comps) for n in range(4)]    # T[n][a] = nabla_n v_a
    out = [zero]*4
    for m in range(4):
        for n in range(4):
            if ginv[m, n] == 0:
                continue
            for a in range(4):
                e = sp.diff(T[n][a], xs[m])
                for s in range(4):
                    if T[n][s] != 0:
                        e -= Gamma[s][m][a]*T[n][s]        # connection on field index a
                    if T[s][a] != 0:
                        e -= Gamma[s][m][n]*T[s][a]        # connection on derivative index n
                out[a] = out[a] + ginv[m, n]*e
    return [sp.cancel(sp.expand(e)) for e in out]

def laplacian_scalar(phi):
    out = zero
    for m in range(4):
        for n in range(4):
            if ginv[m, n] == 0:
                continue
            e = sp.diff(phi, xs[m], xs[n])
            for s in range(4):
                e -= Gamma[s][m][n]*sp.diff(phi, xs[s])
            out = out + ginv[m, n]*e
    return sp.cancel(sp.expand(out))

# ---------- curvature ----------
Ric = [[zero]*4 for _ in range(4)]
for m in range(4):
    for n in range(4):
        e = zero
        for r in range(4):
            e += sp.diff(Gamma[r][n][m], xs[r]) - sp.diff(Gamma[r][r][n], xs[m])
            for l in range(4):
                if Gamma[r][r][l] != 0 and Gamma[l][n][m] != 0:
                    e += Gamma[r][r][l]*Gamma[l][n][m]
                if Gamma[r][m][l] != 0 and Gamma[l][r][n] != 0:
                    e -= Gamma[r][m][l]*Gamma[l][r][n]
        Ric[m][n] = simp(e)
Rscalar = simp(sum(ginv[i, i]*Ric[i][i] for i in range(4)))   # FIXED: g^ii trace
print("=== maximal symmetry ===")
print("R =", Rscalar, "  [expect 12]")
ok_ein = all(simp(Ric[i][j] - (Rscalar/4)*g[i, j]) == 0 for i in range(4) for j in range(4))
print("R_mn = (R/4) g_mn :", ok_ein)

def riemann_lower():
    Rl = [[[[zero]*4 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for b in range(4):
        for c in range(4):
            for dd in range(4):
                Rr = [zero]*4
                for r in range(4):
                    e = sp.diff(Gamma[r][dd][b], xs[c]) - sp.diff(Gamma[r][c][b], xs[dd])
                    for l in range(4):
                        if Gamma[r][c][l] != 0 and Gamma[l][dd][b] != 0:
                            e += Gamma[r][c][l]*Gamma[l][dd][b]
                        if Gamma[r][dd][l] != 0 and Gamma[l][c][b] != 0:
                            e -= Gamma[r][dd][l]*Gamma[l][c][b]
                    Rr[r] = e
                for a in range(4):
                    Rl[a][b][c][dd] = simp(sum(g[a, r]*Rr[r] for r in range(4)))
    return Rl

Riem = riemann_lower()
Ric2 = [[simp(sum(ginv[c, dd]*Riem[c][m][dd][n] for c in range(4) for dd in range(4)))
         for n in range(4)] for m in range(4)]
ok_ric2 = all(simp(Ric2[i][j] - Ric[i][j]) == 0 for i in range(4) for j in range(4))
print("Riemann-derived Ricci matches :", ok_ric2)

# ---------- 1. sigma anchor ----------
print()
print("=== sigma anchor: -D(l=1 harmonic) = R/3 = 4 ===")
sigma_ok = True
for i_test in range(2):
    phi = xs[i_test]*Om
    val = simp(-laplacian_scalar(phi)/phi)
    sigma_ok = sigma_ok and (val == 4)
    print(f"  -Delta(x{i_test}*Om)/(x{i_test}*Om) =", val)

# ---------- 2. ghost operator ----------
print()
print("=== ghost operator M = box + c R_m^n ===")
e_syms = sp.symbols('e0 e1 e2 e3', real=True)
eps_const = list(e_syms)
A = [grad_cov(mu, eps_const) for mu in range(4)]
div_eps = simp(sum(ginv[m, n]*A[m][n] for m in range(4) for n in range(4)))
dF = [zero]*4
for mu in range(4):
    e = zero
    for nu in range(4):
        for l in range(4):
            if ginv[nu, l] == 0:
                continue
            term = sp.diff(A[mu][nu] + A[nu][mu], xs[l])
            for s in range(4):
                term -= Gamma[s][l][mu]*(A[s][nu] + A[nu][s])
                term -= Gamma[s][l][nu]*(A[mu][s] + A[s][mu])
            e += ginv[nu, l]*term
    e -= sp.diff(div_eps, xs[mu])
    dF[mu] = simp(e)
lap_eps = laplacian_vec(eps_const)
C_mat = []
for mu in range(4):
    row = []
    for nu in range(4):
        coef = sp.expand(simp(dF[mu] - lap_eps[mu])).coeff(e_syms[nu])
        row.append(sp.nsimplify(coef))
    C_mat.append(row)
c_val = sp.nsimplify(C_mat[0][0]/3)
c_ok = all(sp.nsimplify(C_mat[i][j] - c_val*3*(1 if i == j else 0)) == 0
           for i in range(4) for j in range(4))
print("C diag =", C_mat[0][0], " off-diag =", C_mat[0][1],
      "  =>  c =", c_val, " constant:", c_ok)
print("  [expect c = +1:  M = box + R/4]")

# ---------- 3. TT identity (correct index structure) ----------
print()
print("=== TT shift: 2 R_{mu rho nu sig} h^{rho sig} for symmetric h ===")
hs = {}
for m in range(4):
    for n in range(m, 4):
        hs[(m, n)] = sp.Symbol(f'h{m}{n}', real=True)
def H(m, n):
    return hs[(m, n)] if n >= m else hs[(n, m)]
# h^{rho sig} = g^{rho a} g^{sig b} h_{a b}; R_{abcd} = K(g_ac g_bd - g_ad g_bc), K = R/12
# expected identity: 2 R_{mu rho nu sig} h^{rho sig} = (R/6)(g_{mu nu} g^{ab} h_{ab} - h_{mu nu})
tr_h = sum(ginv[a, b]*H(a, b) for a in range(4) for b in range(4))   # METRIC trace
identity_ok = True
for mu in range(4):
    for nu in range(4):
        e = zero
        for rho in range(4):
            for sig in range(4):
                hab = sum(ginv[rho, a]*ginv[sig, b]*H(a, b) for a in range(4) for b in range(4))
                e += 2*Riem[mu][rho][nu][sig]*hab
        target = (sp.Rational(1, 6)*R_)*(g[mu, nu]*tr_h - H(mu, nu))
        diff = sp.cancel(sp.expand(e - target))
        if diff != 0:
            identity_ok = False
            print(f"  MISMATCH at ({mu},{nu}):", diff)
print("identity 2Riem.h^sym = (R/6)(g.g-trace h - h) :", identity_ok)
print("  => on TT (tr h = 0): 2Riem.h = -(R/6) h  =>  X_TT = -R/6 in THIS Riemann convention;")
print("     with the opposite Riemann sign (dS-graviton standard, m^2 = 2Lam/3): X_TT = +R/6;")
print("     Solodukhin Eq.(133) X^(2) reading on the traceless space: X = +R/3.")

# ---------- 4. Jacobian blocks ----------
print()
print("=== Jacobian M~ = -(1/2) L+ L ===")
def L_op(xi):
    Amu = [grad_cov(mu, xi) for mu in range(4)]
    div = simp(sum(ginv[m, n]*Amu[m][n] for m in range(4) for n in range(4)))
    out = sp.zeros(4, 4)
    for mu in range(4):
        for nu in range(4):
            out[mu, nu] = simp(Amu[mu][nu] + Amu[nu][mu] -
                               sp.Rational(1, 2)*g[mu, nu]*div)
    return out

def Ldagger(h):
    out = []
    for sig in range(4):
        divh = zero
        for rho in range(4):
            for l in range(4):
                if ginv[rho, l] == 0:
                    continue
                term = sp.diff(h[rho, sig], xs[l])
                for s in range(4):
                    term -= Gamma[s][l][rho]*h[s, sig] + Gamma[s][l][sig]*h[rho, s]
                divh += ginv[rho, l]*term
        htr = simp(sum(ginv[i, j]*h[i, j] for i in range(4) for j in range(4)))
        out.append(simp(-2*(divh - sp.Rational(1, 4)*sp.diff(htr, xs[sig]))))
    return out

# (a) Killing vectors: L(K) = 0
kills = []
for i in range(4):
    for j in range(i+1, 4):
        vec = [0, 0, 0, 0]
        vec[j] = xs[i]; vec[i] = -xs[j]
        cov = [simp(sum(g[mu, nu]*vec[nu] for nu in range(4))) for mu in range(4)]
        kills.append(cov)
maxK = 0
for K in kills[:4]:
    LK = L_op(K)
    maxK = max(maxK, max(abs(float(sp.cancel(LK[m, n]).subs({x0: .3, x1: -.7, x2: 1.1, x3: .5})))
                         for m in range(4) for n in range(4)))
print("|L(Killing)| at test point [expect ~0]:", f"{maxK:.1e}")

# (b) gradient block: M~(grad phi) = (3/2) grad(box phi) + (R/2) grad phi
phi = xs[0]**2
xi_g = [sp.diff(phi, xv) for xv in xs]
Mtilde_g = [simp(-sp.Rational(1, 2)*v) for v in Ldagger(L_op(xi_g))]
lap_phi = laplacian_scalar(phi)
grad_lap = [sp.diff(lap_phi, xv) for xv in xs]
res_g = [sp.cancel(sp.expand(Mtilde_g[a] - sp.Rational(3, 2)*grad_lap[a] -
                             sp.Rational(1, 2)*R_*xi_g[a])) for a in range(4)]
grad_ok = all(r == 0 for r in res_g)
print("gradient block: M~(grad phi) = (3/2)grad(box phi) + (R/2)grad phi :",
      grad_ok, "   [k1 = 3/2, k2 = R/2]")

# (c) co-closed block: xi = delta C (C = x wedge A): M~ xi = (box + R/4) xi
Aconst = [sp.Integer(1), 0, 0, 0]
C_form = [[xs[nu]*Aconst[mu] - xs[mu]*Aconst[nu] for nu in range(4)] for mu in range(4)]
xi_cc = [zero]*4
for mu in range(4):
    e = zero
    for nu in range(4):
        for rho in range(4):
            if ginv[nu, rho] == 0:
                continue
            term = sp.diff(C_form[nu][mu], xs[rho])
            for s in range(4):
                term -= Gamma[s][rho][nu]*C_form[s][mu] + Gamma[s][rho][mu]*C_form[nu][s]
            e += ginv[nu, rho]*term
    xi_cc[mu] = sp.cancel(sp.expand(e))
div_cc = simp(sum(ginv[m, n]*grad_cov(n, xi_cc)[m] for m in range(4) for n in range(4)))
print("div(xi_cc) [expect 0]:", div_cc == 0)
Mtilde_cc = [simp(-sp.Rational(1, 2)*v) for v in Ldagger(L_op(xi_cc))]
lap_cc = laplacian_vec(xi_cc)
res_cc = [sp.cancel(sp.expand(Mtilde_cc[a] - lap_cc[a] -
                              sp.Rational(1, 4)*R_*xi_cc[a])) for a in range(4)]
cc_ok = all(r == 0 for r in res_cc)
print("co-closed block: M~(xi) = (box + R/4) xi :", cc_ok,
      "   [= the ghost operator; the Jacobian's cc-part is the FP ghost]")

# ---------- 5. assembly ----------
print()
print("=== a1 assembly (York route), unit S^4, R/12 units ===")
# Convention: a1(-nabla^2 + X) = R/6 - X per comp. R/6 = 2, R/4 = 3, R/3 = 4.
# sigma (X = -R/3): a1 = R/2 = 6 ; weight +1/2  -> +3
# TT (X_TT in {0, +R/6, -R/6}): 5 comps, a1 = 2 - X ; weight +1/2
# Jacobian: det(2 M~)^{1/2}, weight -1/2, rescale 2^{1-d/2} = 1/2:
#   gradient family: (3/2)(box + R/3) = (3/2)*(-( -box - R/3)): X = -R/3:
#     a1(Delta) = R/2 = 6 ; effective (3/2)^{-1} * 6 = 4  ->  1 comp
#   co-closed family: box + R/4 = -(-box - R/4): X = -R/4: a1 = 5R/12 = 5
#     -> 3 comps = 15
#   a1(M~) = 4 + 15 = 19 ; a1(2M~) = 19/2 ; weight -1/2 -> -19/4
print("  sigma:      +3")
print("  Jacobian:   -19/4   (= -(1/2)*(1/2)*(4 + 3*5))")
print("  TT:         (5/2)(2 - X_TT/(R/12))")
print()
rows = []
for name, XTT in [("dS-graviton standard (m^2 = 2Lam/3): X_TT = +R/6", 2),
                  ("Solodukhin Eq.(133) X^(2) reading: X_TT = +R/3", 4),
                  ("this-convention literal (3.74): X_TT = -R/6", -2),
                  ("[withdrawn: slot-pairing error: X_TT = 0]", 0)]:
    w = sp.Rational(3) + sp.Rational(5, 2)*(2 - XTT) - sp.Rational(19, 4)
    neff = (sp.Integer(-23) + w)/12
    rows.append((name, XTT, w, neff))
    print(f"  {name}:  w_grav = {w},  N_eff^SM+grav = {neff}  < 0: {bool(neff < 0)}")

# ---------- 6. Maxwell cross-check disclosure ----------
print()
print("=== Maxwell cross-check (disclosure of the methodological caveat) ===")
# Same restricted-determinant bookkeeping applied to the vector field (Coulomb/York):
#   Z = det'_0(-box)^{1/2} det_T(-box_T)^{-1/2}
#   a1-channel: -(1/2)(R/6) + (1/2)[ 4*(R/6) - (R/6 + R/4) ] = -R/12 + R/8 = +R/24
# vs the triple-anchored Feynman-gauge number (notebook, opus, Solodukhin Eq.134):
#   vector + 2 FP ghosts = -R/2.
print("  York/Coulomb route, this bookkeeping : +R/24  (+0.5 in R/12 units)")
print("  Feynman gauge + 2 ghosts (anchored)  : -R/2   (-6 in R/12 units)")
print("  discrepancy: -13R/24 per vector — the off-shell effective action is")
print("  gauge-choice dependent; restricted-determinant a1 bookkeeping (Kabat-type")
print("  terms) is NOT verified against the standard route. The graviton a1 numbers")
print("  above carry this caveat. The trustworthy closing check is the de-Donder")
print("  spin-2 heat kernel (Christensen-Duff NPB B170, paywalled).")

out = {"R": 12, "einstein_ok": bool(ok_ein and ok_ric2), "sigma_anchor_ok": sigma_ok,
       "ghost_c": str(c_val), "ghost_c_constant": bool(c_ok),
       "TT_identity": "2Riem.h^sym = (R/6)(g.(g-trace h) - h) [verified] => on TT: X_TT = -R/6 (this convention)",
       "TT_identity_ok": identity_ok,
       "killing_L_zero": maxK < 1e-8,
       "gradient_block_ok": grad_ok, "k1": "3/2", "k2": "R/2",
       "coclosed_block_ok": cc_ok,
       "w_grav_branches": [[n, x, str(w), str(ne)] for n, x, w, ne in rows],
       "maxwell_disclosure": "York/Coulomb +R/24 vs Feynman -R/2 (off-shell gauge dependence)"}
with open('/home/z/my-project/scripts/graviton_a1_results.json', 'w') as f:
    json.dump(out, f, indent=2)
print()
print("saved: graviton_a1_results.json")
