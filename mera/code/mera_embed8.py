"""mera_embed8.py — exact warm-start embedding: N=128 chi=4 params -> chi=8 params.

Both chi values share the SAME tree (N=128, L=7 levels, n_l = 128>>l); only
site dims at levels >= 2 differ (4 vs 8). Site-dim embedding (preserving the
Z2 even/odd split, even sector first):

    old even {0,1} -> new {0,1};   old odd {2,3} -> new {4,5}
    complement: even {2,3}, odd {6,7}

With this, every tensor embeds block-diagonally in (old | complement):
  u_l (l>=2):  H8_even = [[H4_even, R], [R^dag, C]] in the parity-sorted pair
               basis (R = coupling, C = complement generator) -> expm gives
               U8 = old-block = U4 exactly when R = 0, C = 0.
  w_l (l>=2):  X8_even = [[X4_even (8x2), 0], [0, V (24x2)]]  (V random
               orthonormal) -> isometry_of(X8) = X8; old upper-site channels
               receive exactly the X4 columns.
  w1:          output side only: X8e = [X4e | Ve] with Ve orthonormal
               completion orthogonal to X4e's columns.
  u0, w0, u1:  dims identical -> params copied verbatim.
  t:           t8_even = [t4_even, 0, 0] (pad with zeros).

With R = 0, C = 0 the chi=8 state is EXACTLY the chi=4 state (every layer
preserves the old/complement split and t selects the old channel), so the
chi=8 energy equals the chi=4 energy to numerical precision — validated below.
With small random R (coupling) the state is perturbed O(R^2) and all
complement parameters acquire gradients: that is the Phase B initialization.

Usage: python3 mera_embed8.py            (validates exact embedding, saves both)
"""
import sys, time, json
sys.path.insert(0, '/home/z/my-project/scripts')
import numpy as np
import jax.numpy as jnp
from mera_opt import pair_parity_order, tensors_from_params
import mera_real
mera_real.set_backend('np')
from mera_real import MERA

N = 128
_H2 = -np.kron(np.array([[0, 1], [1, 0]]), np.array([[0, 1], [1, 0]])).astype(complex)
_H1 = -np.array([[1, 0], [0, -1]], dtype=complex)

# ---------------- position maps (parity-sorted pair bases) ----------------
EV4, OD4 = pair_parity_order(4, 2)[:(4 * 4 // 2)], None
_order4 = pair_parity_order(4, 2)
_ev4 = [x for x in _order4[:8]]
_od4 = [x for x in _order4[8:]]
_order8 = pair_parity_order(8, 4)
_ev8 = [x for x in _order8[:32]]
_od8 = [x for x in _order8[32:]]

M_SITE = {0: 0, 1: 1, 2: 4, 3: 5}   # old site index (d=4) -> new (d=8)

def _map_positions(old_pairs, new_pairs):
    out = []
    for p4 in old_pairs:
        i1, i2 = divmod(p4, 4)
        p8 = M_SITE[i1] * 8 + M_SITE[i2]
        out.append(new_pairs.index(p8))
    return out

MAP_EV = _map_positions(_ev4, _ev8)     # len 8, values in 0..31
MAP_OD = _map_positions(_od4, _od8)
COMP_EV = [a for a in range(32) if a not in MAP_EV]
COMP_OD = [a for a in range(32) if a not in MAP_OD]

# ---------------- block extraction / packing helpers ----------------

def unpack_herm(params, m):
    H = (params[:m * m] + 1j * params[m * m:2 * m * m]).reshape(m, m)
    return (H + H.conj().T) / 2.0

def pack_herm(H):
    m = H.shape[0]
    return np.concatenate([H.real.reshape(-1), H.imag.reshape(-1)])

def unpack_X(params, big, small):
    return (params[:big * small] + 1j * params[big * small:2 * big * small]).reshape(big, small)

def pack_X(X):
    big, small = X.shape
    return np.concatenate([X.real.reshape(-1), X.imag.reshape(-1)])

def orthonormal_completion(A, k, rng):
    """k random orthonormal columns orthogonal to A's columns (A: (m, n))."""
    m = A.shape[0]
    G = rng.standard_normal((m, k)) + 1j * rng.standard_normal((m, k))
    G = G - A @ (A.conj().T @ G)
    Q, _ = np.linalg.qr(G)
    return Q[:, :k]

# ---------------- the embedding ----------------

def embed_params(prm4, us4, ws4, coupling=0.0, seed=7):
    """us4/ws4: ASSEMBLED chi=4 tensors (from tensors_from_params) — the w
    embedding must use the assembled isometries (isometry_of re-orthonormalizes
    raw params), while the u embedding works exactly at the generator level."""
    rng = np.random.default_rng(seed)
    prm8 = {}
    prm8['u0'] = np.array(prm4['u0'])
    prm8['u1'] = np.array(prm4['u1'])
    prm8['w0'] = np.array(prm4['w0'])
    # assembled w-blocks in the parity-sorted basis: Xe4[a, c] = W4[ev4[a], c]
    def w_blocks(W4):
        W4 = np.asarray(W4)
        Xe = np.array([[W4[p, c] for c in (0, 1)] for p in _ev4])
        Xo = np.array([[W4[p, 2 + c] for c in (0, 1)] for p in _od4])
        return Xe, Xo
    # ---- w1: output side only (input pair basis is identical, d_in = 4) ----
    Xe4, Xo4 = w_blocks(ws4[1])
    Ve = orthonormal_completion(Xe4, 2, rng)
    Vo = orthonormal_completion(Xo4, 2, rng)
    Xe1 = np.hstack([Xe4, Ve])
    Xo1 = np.hstack([Xo4, Vo])
    if coupling > 0:
        # degeneracy guard: exactly-orthonormal columns make X^dag X = I with
        # exactly tied eigenvalues -> eigh backward is 0/0 (NaN). A 1e-4 random
        # perturbation breaks the tie; energy shift ~1e-4, well under the coupling.
        for X in (Xe1, Xo1):
            X += 1e-4 * (rng.standard_normal(X.shape) + 1j * rng.standard_normal(X.shape))
    prm8['w1'] = np.concatenate([pack_X(Xe1), pack_X(Xo1)])
    # ---- u_l, w_l for l = 2..6 ----
    for l in range(2, 7):
        # u: generators 8x8 -> 32x32 in (old | comp) blocks
        p = np.array(prm4['u%d' % l])
        H4e = unpack_herm(p[:128], 8)
        H4o = unpack_herm(p[128:], 8)
        H8e = np.zeros((32, 32), dtype=complex)
        H8o = np.zeros((32, 32), dtype=complex)
        H8e[np.ix_(MAP_EV, MAP_EV)] = H4e
        H8o[np.ix_(MAP_OD, MAP_OD)] = H4o
        if coupling > 0:
            # normalized: Frobenius norm of each coupling block = `coupling`
            # (per-entry scale would blow up the collective mixing amplitude)
            Re = rng.standard_normal((8, 24)) + 1j * rng.standard_normal((8, 24))
            Ro = rng.standard_normal((8, 24)) + 1j * rng.standard_normal((8, 24))
            Re = coupling * Re / np.linalg.norm(Re)
            Ro = coupling * Ro / np.linalg.norm(Ro)
            H8e[np.ix_(MAP_EV, COMP_EV)] = Re
            H8e[np.ix_(COMP_EV, MAP_EV)] = Re.conj().T
            H8o[np.ix_(MAP_OD, COMP_OD)] = Ro
            H8o[np.ix_(COMP_OD, MAP_OD)] = Ro.conj().T
        prm8['u%d' % l] = np.concatenate([pack_herm(H8e), pack_herm(H8o)])
        # w: assembled X 8x2 -> 32x4 with old rows/cols + complement completion
        Xe4, Xo4 = w_blocks(ws4[l])
        X8e = np.zeros((32, 4), dtype=complex)
        X8o = np.zeros((32, 4), dtype=complex)
        X8e[np.ix_(MAP_EV, [0, 1])] = Xe4
        X8o[np.ix_(MAP_OD, [0, 1])] = Xo4
        X8e[np.ix_(COMP_EV, [2, 3])] = orthonormal_completion(
            np.zeros((24, 0), dtype=complex), 2, rng)
        X8o[np.ix_(COMP_OD, [2, 3])] = orthonormal_completion(
            np.zeros((24, 0), dtype=complex), 2, rng)
        if coupling > 0:
            for X in (X8e, X8o):
                X += 1e-4 * (rng.standard_normal(X.shape)
                             + 1j * rng.standard_normal(X.shape))
        prm8['w%d' % l] = np.concatenate([pack_X(X8e), pack_X(X8o)])
    # ---- t: pad ----
    t4 = np.array(prm4['t'])          # [re(2), im(2)]
    prm8['t'] = np.concatenate([t4[:2], [0.0, 0.0], t4[2:], [0.0, 0.0]])
    return prm8

# ---------------- full numpy energy (independent path) ----------------

def numpy_energy(prm, chi):
    us, ws, t = tensors_from_params({k: jnp.array(v) for k, v in prm.items()}, N, chi)
    m = MERA(N, chi, tensors=([np.asarray(x) for x in us],
                              [np.asarray(x) for x in ws], np.asarray(t)))
    e = 0.0
    for a in range(N):
        e += np.real(np.trace(m.rdm([a, (a + 1) % N]) @ _H2))
        e += np.real(np.trace(m.rdm([a]) @ _H1))
    return float(e / N)

def main():
    t0 = time.time()
    d = np.load('/home/z/my-project/scripts/mera_params_N128_chi4_seed0.npz')
    prm4 = {k: d[k] for k in d.files}
    us4, ws4, _t4 = tensors_from_params({k: jnp.array(v) for k, v in prm4.items()}, N, 4)
    us4 = [np.asarray(x) for x in us4]
    ws4 = [np.asarray(x) for x in ws4]
    print('loaded chi=4 Phase-A params:', {k: v.shape for k, v in prm4.items()}, flush=True)

    # exact embedding (coupling = 0)
    prm8_exact = embed_params(prm4, us4, ws4, coupling=0.0)
    # warm start (small normalized coupling; complement sector activates under
    # gradients — at zero coupling the complement directions are exactly
    # stationary, by block preservation)
    prm8_warm = embed_params(prm4, us4, ws4, coupling=0.03, seed=7)

    print('computing chi=4 full energy (numpy, 256 terms)...', flush=True)
    e4 = numpy_energy(prm4, 4)
    print(f'  E4/N = {e4:.9f}', flush=True)
    print('computing chi=8 embedded full energy (numpy, 256 terms)...', flush=True)
    e8 = numpy_energy(prm8_exact, 8)
    print(f'  E8(exact embedding)/N = {e8:.9f}   |E8 - E4| = {abs(e8 - e4):.2e}', flush=True)
    ok = abs(e8 - e4) < 1e-6
    print('EXACT EMBEDDING VALIDATED' if ok else '!!! EMBEDDING MISMATCH', flush=True)

    print('computing chi=8 warm-start energy (coupling 0.02)...', flush=True)
    e8w = numpy_energy(prm8_warm, 8)
    print(f'  E8(warm)/N = {e8w:.9f}   (shift vs exact: {e8w - e8:+.2e})', flush=True)

    np.savez('/home/z/my-project/scripts/mera_params_N128_chi8_embed_exact.npz',
             **prm8_exact)
    np.savez('/home/z/my-project/scripts/mera_params_N128_chi8_warmstart.npz',
             **prm8_warm)
    n4 = sum(int(np.prod(v.shape)) for v in prm4.values())
    n8 = sum(int(np.prod(v.shape)) for v in prm8_warm.values())
    out = dict(E4=e4, E8_exact=e8, E8_exact_diff=float(abs(e8 - e4)),
               E8_warm=e8w, validated=bool(ok),
               params_chi4=n4, params_chi8=n8,
               coupling=0.03, seconds=time.time() - t0)
    with open('/home/z/my-project/scripts/mera_embed8_result.json', 'w') as f:
        json.dump(out, f, indent=1)
    print(f'saved exact + warmstart params ({n4} -> {n8} parameters) '
          f'({time.time()-t0:.0f}s)')


if __name__ == '__main__':
    main()
