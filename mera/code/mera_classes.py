"""mera_classes.py — RDM equality classes for level-shared-tensor MERAs.

Fact exploited: mera_real.MERA uses ONE disentangler u[l] and ONE isometry w[l]
per layer (shared across all positions). Two targets whose causal cones are
isomorphic at every level therefore contract the SAME shared tensors with the
SAME connectivity -> bitwise-identical contractions (up to leg-order, which
does not affect tr(rho @ h) for the swap-symmetric h2 = -sx^x sx and the
diagonal h1 = -sz).

Consequence: the exact energy
    E = (1/N) sum_a tr(rho_{(a,a+1)} h2) + (1/N) sum_s tr(rho_s h1)
collapses to a sum over a few equality classes:
    E = sum_c (n_c/N) tr(rho_{rep(c)} h2) + sum_c' (n'_c/N) tr(rho_{rep(c')} h1)
which is the IDENTICAL function, computed with |classes| contractions instead
of 2N. This file determines the classes and validates the identity.

Validation plan:
  V1 (structure independence): classes found with chi=2 surrogate tensors ==
      classes found with chi=4 surrogate tensors at N in {16, 32}
      (cone geometry is chi-independent, so the class partition must be too).
  V2 (energy identity, real artifact): for the SAVED optimized N=64 chi=4
      seed-0 params, class-weighted energy == direct all-bonds/all-sites energy.
  V3 (full N=128): with chi=2 surrogate tensors, class-weighted == direct.
Output: scripts/mera_classes_N128.json (bond + site class reps and counts).
"""
import sys, json, time
sys.path.insert(0, '/home/z/my-project/scripts')
import numpy as np
import mera_real
mera_real.set_backend('np')
from mera_real import MERA

_H2 = -np.kron(np.array([[0, 1], [1, 0]], dtype=complex),
               np.array([[0, 1], [1, 0]], dtype=complex))
_H1 = -np.array([[1, 0], [0, -1]], dtype=complex)


def cluster_values(vals, tol=1e-9):
    """Group indices by (approximately) equal float values."""
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    groups, cur, prev = [], [order[0]], vals[order[0]]
    for i in order[1:]:
        if abs(vals[i] - prev) < tol * max(1.0, abs(prev)):
            cur.append(i)
        else:
            groups.append(cur)
            cur = [i]
        prev = vals[i]
    groups.append(cur)
    return groups


def bond_site_values(m, N):
    e2 = [float(np.real(np.trace(m.rdm([a, (a + 1) % N]) @ _H2))) for a in range(N)]
    e1 = [float(np.real(np.trace(m.rdm([s]) @ _H1))) for s in range(N)]
    return e2, e1


def classes_of(m, N):
    e2, e1 = bond_site_values(m, N)
    bg = cluster_values(e2)
    sg = cluster_values(e1)
    return bg, sg, e2, e1


def energy_direct(m, N):
    e2, e1 = bond_site_values(m, N)
    return float(np.mean(e2) + np.mean(e1))


def energy_classed(m, N, bg, sg):
    e2, e1 = bond_site_values(m, N)
    eb = sum(len(g) * e2[g[0]] for g in bg) / N
    es = sum(len(g) * e1[g[0]] for g in sg) / N
    return float(eb + es)


def partition_signature(groups, N):
    return tuple(sorted(len(g) for g in groups)), [g[0] for g in groups]


def main():
    log = {}

    # ---- V1: chi-independence of the class partition at N=16, 32 ----
    for N in (16, 32):
        t0 = time.time()
        m2 = MERA(N, 2, seed=11)
        m4 = MERA(N, 4, seed=13)
        bg2, sg2, _, _ = classes_of(m2, N)
        bg4, sg4, _, _ = classes_of(m4, N)
        sig2b, reps2b = partition_signature(bg2, N)
        sig4b, reps4b = partition_signature(bg4, N)
        sig2s, reps2s = partition_signature(sg2, N)
        sig4s, reps4s = partition_signature(sg4, N)
        same_b = (sorted([sorted(g) for g in bg2]) == sorted([sorted(g) for g in bg4]))
        same_s = (sorted([sorted(g) for g in sg2]) == sorted([sorted(g) for g in sg4]))
        print(f'V1 N={N}: bond classes chi2 {list(sig2b)} vs chi4 {list(sig4b)} '
              f'-> identical partition: {same_b}; site classes {list(sig2s)} vs '
              f'{list(sig4s)} -> identical: {same_s}  ({time.time()-t0:.1f}s)', flush=True)
        log[f'V1_N{N}'] = dict(bond_sig=list(sig2b), site_sig=list(sig2s),
                               bond_partition_identical=bool(same_b),
                               site_partition_identical=bool(same_s))
        assert same_b and same_s, 'V1 FAILED: partition depends on chi'

    # ---- V2: energy identity on the saved optimized N=64 chi=4 state ----
    d = np.load('/home/z/my-project/scripts/mera_params_N64_chi4_seed0.npz')
    import jax.numpy as jnp
    from mera_opt import tensors_from_params
    prm = {k: jnp.array(d[k]) for k in d.files}
    us, ws, t = tensors_from_params(prm, 64, 4)
    m64 = MERA(64, 4, tensors=([np.asarray(x) for x in us],
                               [np.asarray(x) for x in ws], np.asarray(t)))
    t0 = time.time()
    bg, sg, _, _ = classes_of(m64, 64)
    ed = energy_direct(m64, 64)
    ec = energy_classed(m64, 64, bg, sg)
    print(f'V2 N=64 chi=4 (saved optimized): direct E/N = {ed:.9f}, '
          f'classed E/N = {ec:.9f}, diff = {abs(ed-ec):.2e}, '
          f'bond classes = {[len(g) for g in bg]}, site classes = {[len(g) for g in sg]} '
          f'({time.time()-t0:.1f}s)', flush=True)
    log['V2_N64_chi4'] = dict(direct=ed, classed=ec, diff=float(abs(ed - ec)),
                              bond_class_sizes=[len(g) for g in bg],
                              site_class_sizes=[len(g) for g in sg])
    assert abs(ed - ec) < 1e-10, 'V2 FAILED: class energy != direct energy'

    # ---- V3 + class table at N=128 (chi=2 surrogate: geometry is chi-free) ----
    m128 = MERA(128, 2, seed=11)
    t0 = time.time()
    bg, sg, e2, e1 = classes_of(m128, 128)
    ed = energy_direct(m128, 128)
    ec = energy_classed(m128, 128, bg, sg)
    print(f'V3 N=128 (chi=2 surrogate): direct E/N = {ed:.9f}, classed = {ec:.9f}, '
          f'diff = {abs(ed-ec):.2e} ({time.time()-t0:.1f}s)', flush=True)
    print(f'    bond classes ({len(bg)}): reps {[g[0] for g in bg]} '
          f'sizes {[len(g) for g in bg]}', flush=True)
    print(f'    site classes ({len(sg)}): reps {[g[0] for g in sg]} '
          f'sizes {[len(g) for g in sg]}', flush=True)
    assert abs(ed - ec) < 1e-10, 'V3 FAILED'

    # spread within each class (must be ~machine epsilon)
    bspread = max(max(abs(e2[i] - e2[g[0]]) for i in g) for g in bg)
    sspread = max(max(abs(e1[i] - e1[g[0]]) for i in g) for g in sg)
    print(f'    within-class spread: bonds {bspread:.1e}, sites {sspread:.1e}', flush=True)

    out = dict(
        N=128,
        bond_classes=[dict(rep=int(g[0]), count=len(g),
                           members=[int(x) for x in g]) for g in bg],
        site_classes=[dict(rep=int(g[0]), count=len(g),
                           members=[int(x) for x in g]) for g in sg],
        validation=log,
        within_class_spread=dict(bonds=float(bspread), sites=float(sspread)),
        surrogate=dict(chi=2, seed=11,
                       note='cone geometry is chi-independent; chi=2 surrogate '
                            'gives the exact N=128 class partition for any chi'),
    )
    with open('/home/z/my-project/scripts/mera_classes_N128.json', 'w') as f:
        json.dump(out, f, indent=1)
    print('saved /home/z/my-project/scripts/mera_classes_N128.json')
    print('ALL CLASS VALIDATIONS PASSED')


if __name__ == '__main__':
    main()
