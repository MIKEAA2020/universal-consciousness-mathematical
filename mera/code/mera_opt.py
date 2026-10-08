"""mera_opt.py — variational optimization of the binary MERA for the critical
transverse-field Ising model (h=1), plus exact reference energies.

Model: H = -sum_i (sx_i sx_{i+1} + sz_i), periodic. Exact E0/N -> -4/pi.

Z2 symmetry: the ansatz is restricted to states invariant under P = prod_i sz_i.
Each leg's Hilbert space splits into even/odd charge sectors; disentanglers are
block-diagonal in total pair parity, isometries map pair parity to the out-leg
charge, and the top vector lives in the even sector. Without this, the
optimizer falls into symmetry-broken, low-entanglement local minima whose
energy is nearly degenerate with the true ground state at criticality.
"""
import sys
sys.path.insert(0, '/home/z/my-project/scripts')
import numpy as np
import jax
jax.config.update('jax_platform_name', 'cpu')
import jax.numpy as jnp
from jax.scipy.linalg import expm
import mera_real
from mera_real import MERA

def dims_of(N, chi):
    L = int(np.log2(N))
    return [2] + [min(chi, 4)] + [chi] * (L - 1)

# ------------------------------------------------------------------ blocks

def herm(params, m):
    M = (params[:m * m] + 1j * params[m * m:2 * m * m]).reshape(m, m)
    return (M + M.conj().T) / 2.0

def unitary_of(params, m):
    return expm(1j * herm(params, m))

def isometry_of(params, big, small):
    X = (params[:big * small] + 1j * params[big * small:2 * big * small]).reshape(big, small)
    w, V = jnp.linalg.eigh(X.conj().T @ X)
    winv = 1.0 / jnp.sqrt(jnp.maximum(w, 1e-12))
    Sinv = V @ jnp.diag(winv) @ V.conj().T
    return X @ Sinv

def sector_split(d):
    """Z2 charge split of a d-dim leg: even sector first (de), odd last (do)."""
    return (d + 1) // 2, d // 2

def pair_parity_order(d, de):
    """Permutation of pair indices i1*d+i2 grouping total-parity-even first."""
    ev = [i1 * d + i2 for i1 in range(d) for i2 in range(d)
          if ((i1 >= de) + (i2 >= de)) % 2 == 0]
    od = [i1 * d + i2 for i1 in range(d) for i2 in range(d)
          if ((i1 >= de) + (i2 >= de)) % 2 == 1]
    return ev + od

def perm_matrix(d, order):
    P = jnp.zeros((d, d)).at[jnp.arange(d), jnp.array(order)].set(1.0)
    return P

def assemble_u(params, d):
    """Z2-invariant disentangler: block-diagonal (even, odd) in the
    pair-parity basis, permuted back to natural order."""
    de, do = sector_split(d)
    ne, no = de * de + do * do, 2 * de * do
    Ue = unitary_of(params[:2 * ne * ne], ne)
    Uo = unitary_of(params[2 * ne * ne:2 * ne * ne + 2 * no * no], no)
    order = pair_parity_order(d, de)
    P = perm_matrix(d * d, order)
    Ub = jnp.block([[Ue, jnp.zeros((ne, no), dtype=Ue.dtype)],
                    [jnp.zeros((no, ne), dtype=Ue.dtype), Uo]])
    return P.T @ Ub @ P          # U[order[i], order[j]] = Ub[i, j]

def assemble_w(params, d_in, d_out):
    """Z2-invariant isometry: in-pair total parity = out-leg charge.
    Even-parity rows map into the out-leg's even sector (cols 0..de_o),
    odd-parity rows into the odd sector (cols de_o..d_out)."""
    de_i, do_i = sector_split(d_in)
    de_o, do_o = sector_split(d_out)
    ne, no = de_i * de_i + do_i * do_i, 2 * de_i * do_i
    We = isometry_of(params[:2 * ne * de_o], ne, de_o)
    Wo = isometry_of(params[2 * ne * de_o:2 * ne * de_o + 2 * no * do_o], no, do_o)
    order = pair_parity_order(d_in, de_i)
    P = perm_matrix(d_in * d_in, order)
    Wb = jnp.zeros((ne + no, d_out), dtype=We.dtype)
    Wb = Wb.at[:ne, :de_o].set(We)
    Wb = Wb.at[ne:, de_o:].set(Wo)
    return P.T @ Wb               # W[order[i], :] = Wb[i, :]

def u_param_count(d):
    de, do = sector_split(d)
    ne, no = de * de + do * do, 2 * de * do
    return 2 * ne * ne + 2 * no * no

def w_param_count(d_in, d_out):
    de_i, do_i = sector_split(d_in)
    de_o, do_o = sector_split(d_out)
    ne, no = de_i * de_i + do_i * do_i, 2 * de_i * do_i
    return 2 * ne * de_o + 2 * no * do_o

# ------------------------------------------------------------------ params

def params_init(N, chi, seed, scale=0.35):
    rng = np.random.default_rng(seed)
    dims = dims_of(N, chi)
    L = int(np.log2(N))
    prm = {}
    for l in range(L):
        prm['u%d' % l] = jnp.array(rng.standard_normal(u_param_count(dims[l])) * scale)
        prm['w%d' % l] = jnp.array(rng.standard_normal(w_param_count(dims[l], dims[l + 1])) * scale)
    de_top, _ = sector_split(dims[L])
    prm['t'] = jnp.array(rng.standard_normal(2 * de_top) * 0.5)
    return prm

def tensors_from_params(prm, N, chi):
    dims = dims_of(N, chi)
    L = int(np.log2(N))
    us, ws = [], []
    for l in range(L):
        us.append(assemble_u(prm['u%d' % l], dims[l]))
        ws.append(assemble_w(prm['w%d' % l], dims[l], dims[l + 1]))
    de_top, do_top = sector_split(dims[L])
    t_even = prm['t'][:de_top] + 1j * prm['t'][de_top:2 * de_top]
    t_even = t_even / jnp.linalg.norm(t_even)
    t = jnp.concatenate([t_even, jnp.zeros((do_top,), dtype=t_even.dtype)])
    return us, ws, t

# ------------------------------------------------------------------ energy

_SX = np.array([[0, 1], [1, 0]], dtype=complex)
_SZ = np.array([[1, 0], [0, -1]], dtype=complex)
_H2 = -np.kron(_SX, _SX).astype(complex)
_H1 = -_SZ

def energy_fn(prm, N, chi, bonds=None, sites=None):
    """E/N = mean over bonds <h2> + mean over sites <h1>. The finite MERA tree
    has a root which breaks exact translation invariance, so we average over
    ALL N bonds and ALL N sites (or provided subsets for stochastic steps)."""
    us, ws, t = tensors_from_params(prm, N, chi)
    m = MERA(N, chi, tensors=(us, ws, t))
    h2 = jnp.asarray(_H2)
    h1 = jnp.asarray(_H1)
    if bonds is None:
        bonds = [(i, (i + 1) % N) for i in range(N)]
    if sites is None:
        sites = list(range(N))
    e_b = 0.0
    for (a, b) in bonds:
        e_b = e_b + jnp.trace(m.rdm([a, b]) @ h2)
    e_b = jnp.real(e_b) / len(bonds)
    e_s = 0.0
    for s in sites:
        e_s = e_s + jnp.trace(m.rdm([s]) @ h1)
    e_s = jnp.real(e_s) / len(sites)
    return e_b + e_s

# ------------------------------------------------------------------ references

def tfim_exact_energy(N, h=1.0):
    """GS energy of the periodic TFIM via free fermions, both JW sectors."""
    def esum(ks):
        ks = np.asarray(ks)
        eps = 2.0 * np.sqrt(1.0 + h * h - 2.0 * h * np.cos(ks))
        return -0.5 * np.sum(eps)
    E_ap = esum([(2 * m + 1) * np.pi / N for m in range(N)])
    E_p = esum([2 * m * np.pi / N for m in range(N)])
    return min(E_ap, E_p)

def tfim_ed(N):
    """Exact diagonalization (sparse) of the periodic TFIM. For N <= 16."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as sla
    dim = 2 ** N
    sx = sp.csr_matrix(np.array([[0, 1], [1, 0]], dtype=float))
    sz = sp.csr_matrix(np.array([[1, 0], [0, -1]], dtype=float))
    I = sp.identity(2, format='csr')
    def kron_op(op, site):
        out = sp.identity(1, format='csr')
        for s in range(N):
            out = sp.kron(out, op if s == site else I, format='csr')
        return out
    H = sp.csr_matrix((dim, dim))
    for i in range(N):
        H = H - kron_op(sx, i) @ kron_op(sx, (i + 1) % N)
        H = H - kron_op(sz, i)
    vals, vecs = sla.eigsh(H, k=1, which='SA')
    return float(vals[0]), np.asarray(vecs[:, 0])

# ------------------------------------------------------------------ optimizer

def optimize(N, chi, seed=0, steps=1500, lr=0.03, n_sample=None):
    mera_real.set_backend('jax')
    rng = np.random.default_rng(seed + 999)
    use_all = n_sample is None or N <= 32
    if use_all:
        efun = jax.jit(lambda p: energy_fn(p, N, chi))
        gfun = jax.jit(jax.grad(lambda p: energy_fn(p, N, chi)))
    else:
        efun = jax.jit(lambda p, bidx, sidx: energy_fn(
            p, N, chi,
            bonds=[(i, (i + 1) % N) for i in bidx], sites=list(sidx)))
        gfun = jax.jit(jax.grad(lambda p, bidx, sidx: energy_fn(
            p, N, chi,
            bonds=[(i, (i + 1) % N) for i in bidx], sites=list(sidx))))
    prm = params_init(N, chi, seed)
    m = {k: jnp.zeros_like(v) for k, v in prm.items()}
    v = {k: jnp.zeros_like(v) for k, v in prm.items()}
    b1, b2, eps = 0.9, 0.999, 1e-8
    hist = []
    for t in range(1, steps + 1):
        lr_t = lr * (0.5 * (1 + np.cos(np.pi * min(t / steps, 1.0)))) ** 2
        if use_all:
            g = gfun(prm)
        else:
            bidx = jnp.array(rng.choice(N, size=n_sample, replace=False))
            sidx = jnp.array(rng.choice(N, size=n_sample, replace=False))
            g = gfun(prm, bidx, sidx)
        for k in prm:
            m[k] = b1 * m[k] + (1 - b1) * g[k]
            v[k] = b2 * v[k] + (1 - b2) * g[k] ** 2
            mh = m[k] / (1 - b1 ** t)
            vh = v[k] / (1 - b2 ** t)
            prm[k] = prm[k] - lr_t * mh / (jnp.sqrt(vh) + eps)
        if t % max(1, steps // 10) == 0 or t == 1:
            if use_all:
                e = float(efun(prm))
            else:
                bidx = jnp.array(list(range(N)))
                e = float(efun(prm, bidx, bidx))
            hist.append((t, e))
            print(f'  [N={N} chi={chi} seed={seed}] step {t}: E/N = {e:.6f}', flush=True)
    if use_all:
        e_final = float(efun(prm))
    else:
        bidx = jnp.array(list(range(N)))
        e_final = float(efun(prm, bidx, bidx))
    mera_real.set_backend('np')
    return prm, e_final, hist

if __name__ == '__main__':
    import time
    N, chi, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 0
    steps = int(sys.argv[4]) if len(sys.argv) > 4 else 1500
    t0 = time.time()
    prm, e, hist = optimize(N, chi, seed, steps)
    ex = tfim_exact_energy(N) / N
    print(f'FINAL: N={N} chi={chi} seed={seed}  E/N = {e:.6f}  exact = {ex:.6f}  '
          f'err = {abs(e - ex):.2e}  rel = {abs(e - ex) / abs(ex):.2e}  ({time.time() - t0:.0f}s)')
    np.savez(f'/home/z/my-project/scripts/mera_params_N{N}_chi{chi}_seed{seed}.npz',
             **{k: np.asarray(v) for k, v in prm.items()})
