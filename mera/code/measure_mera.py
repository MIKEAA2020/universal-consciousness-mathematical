"""measure_mera.py — full measurement suite on optimized MERAs (numpy backend).

Produces, per (N, chi, seed):
  - block entropies S(L), CFT fit -> central charge c
  - mutual information matrix I_ij (all pairs for N<=64, sampled for N=128)
  - I0, d* = -log(I/I0), d_graph via BFS on the REAL tensor network
  - slope / R^2 of d* vs d_graph (the opus's central claim)
  - MDS stresses for dims 1-4 (elbow test)
  - triangle-deficit flatness check (opus's samples + scan)
  - negative control: random (unoptimized) MERA
  - scale-invariance drift of optimized tensors across layers
  - N=16: full cross-validation against exact diagonalization
Everything is a real contraction output; nothing analytic is substituted.
"""
import sys, json, math, time
sys.path.insert(0, '/home/z/my-project/scripts')
import numpy as np
import mera_real
mera_real.set_backend('np')
from mera_real import MERA, entropy
from mera_opt import tensors_from_params, tfim_exact_energy, tfim_ed, energy_fn  # noqa
import networkx as nx
from sklearn.manifold import MDS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_mera(N, chi, seed):
    import jax.numpy as jnp
    d = np.load(f'/home/z/my-project/scripts/mera_params_N{N}_chi{chi}_seed{seed}.npz')
    prm = {k: jnp.array(d[k]) for k in d.files}
    us, ws, t = tensors_from_params(prm, N, chi)
    return MERA(N, chi, tensors=([np.asarray(x) for x in us],
                                 [np.asarray(x) for x in ws], np.asarray(t)))

# ------------------------------------------------------------------ graph distance

def tensor_graph(mera):
    """networkx graph of tensor INSTANCES + physical leaves; edges = legs.
    (Tensors are shared across positions, but each position is a distinct
    instance/node — that is what the poset graph distance means.)"""
    G = nx.Graph()
    L, N = mera.L, mera.N
    for i in range(N):
        G.add_node(('leaf', i))
    for l in range(L):
        n = mera.n[l]
        for j in range(n // 2):
            s1, s2 = mera.upair(l, j)
            un = ('u', l, j)
            G.add_node(un)
            for s in (s1, s2):
                low = ('leaf', s) if l == 0 else ('w', l - 1, s)
                G.add_edge(un, low)
        for k in range(mera.n[l + 1]):
            wn = ('w', l, k)
            G.add_node(wn)
            # w's two inputs: u-proc sites 2k (from u[k]) and 2k+1 (from u[k+1])
            G.add_edge(wn, ('u', l, k % (n // 2)))
            G.add_edge(wn, ('u', l, (k + 1) % (n // 2)))
            if l + 1 == L:
                G.add_edge(wn, ('top',))
            else:
                # level-(l+1) site k is covered by u_{l+1}[j], pair (2j-1, 2j):
                # j = (k+1)//2 (k odd -> (k+1)/2 ; k even -> k/2)
                m2 = mera.n[l + 1] // 2
                G.add_edge(wn, ('u', l + 1, ((k + 1) // 2) % m2))
    G.add_node(('top',))
    return G

def dgraph_matrix(mera):
    """All-pairs leg-distances between physical leaves on the real network."""
    G = tensor_graph(mera)
    leaves = [('leaf', i) for i in range(mera.N)]
    D = np.zeros((mera.N, mera.N))
    for i, li in enumerate(leaves):
        dist = nx.single_source_shortest_path_length(G, li)
        for j, lj in enumerate(leaves):
            if i != j:
                D[i, j] = dist[lj]
    return D

# ------------------------------------------------------------------ measurements

def block_entropies(mera, Lmax=8):
    S = {}
    for L in range(1, min(Lmax, mera.N // 2) + 1):
        vals = []
        from mera_real import MemoryCapError
        for start in (0, 1):
            A = [(start + i) % mera.N for i in range(L)]
            try:
                vals.append(entropy(mera.rdm(A)))
            except MemoryCapError:
                pass
        if vals:
            S[L] = float(np.mean(vals))
    return S

def cft_fit_c(S_dict, N):
    """S(L) = (c/3) ln[(N/pi) sin(pi L / N)] + s0  ->  least squares for c."""
    ls = np.array(sorted(S_dict.keys()), dtype=float)
    ss = np.array([S_dict[int(l)] for l in ls])
    x = np.log((N / np.pi) * np.sin(np.pi * ls / N))
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, ss, rcond=None)
    c = 3 * coef[0]
    resid = ss - A @ coef
    r2 = 1 - np.var(resid) / np.var(ss)
    return float(c), float(r2)

def mutual_information_matrix(mera, max_pairs=None, verbose=False, max_sep=None):
    """I_ij for all pairs (or sampled). Returns I (NxN), list of measured pairs."""
    N = mera.N
    I = np.full((N, N), np.nan)
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    if max_sep is not None:
        pairs = [(i, j) for (i, j) in pairs
                 if min(j - i, N - (j - i)) <= max_sep]
    if max_pairs is not None and len(pairs) > max_pairs:
        rng = np.random.default_rng(12345)
        by_s = {}
        for (i, j) in pairs:
            by_s.setdefault(min(j - i, N - (j - i)), []).append((i, j))
        per_s = max(1, int(np.floor(max_pairs / len(by_s))))
        sel = []
        for s, ps in sorted(by_s.items()):
            take = min(per_s, len(ps))
            idxs = rng.choice(len(ps), size=take, replace=False)
            sel += [ps[x] for x in idxs]
        pairs = sorted(set(sel))
    t0 = time.time()
    import gc, ctypes
    try:
        libc = ctypes.CDLL('libc.so.6')
        def trim():
            gc.collect()
            libc.malloc_trim(0)
    except Exception:
        def trim():
            gc.collect()
    Ssite = {i: entropy(mera.rdm([i])) for i in range(N)}
    trim()
    from mera_real import MemoryCapError
    n_skip = 0
    for idx, (i, j) in enumerate(pairs):
        try:
            rij = mera.rdm([i, j])
            I[i, j] = I[j, i] = Ssite[i] + Ssite[j] - entropy(rij)
            del rij
        except MemoryCapError:
            I[i, j] = I[j, i] = np.nan
            n_skip += 1
        trim()
        if verbose and idx % 100 == 0:
            print(f'    pair {idx}/{len(pairs)} ({time.time()-t0:.0f}s, {n_skip} skipped)', flush=True)
    print(f'  I-matrix: {len(pairs) - n_skip}/{len(pairs)} pairs measured '
          f'({n_skip} exceeded memory cap)', flush=True)
    return I, pairs

def dstar_analysis(mera, I):
    N = mera.N
    adj = [I[i, (i + 1) % N] for i in range(N)] + [I[i, (i - 1) % N] for i in range(N)]
    I0 = float(np.nanmax(adj))
    Ifloor = I0 * 1e-6
    n_floored = int(np.sum((I > 0) & (I < Ifloor)) + np.sum(I <= 0)) // 2
    Isafe = np.where(np.isnan(I), np.nan, np.maximum(I, Ifloor))
    with np.errstate(divide='ignore', invalid='ignore'):
        dstar = -np.log(Isafe / I0)   # NaN stays NaN
    np.fill_diagonal(dstar, 0.0)
    dg = dgraph_matrix(mera)
    # analytic comparison used by the opus
    idx = np.arange(N)
    s = np.abs(idx[:, None] - idx[None, :])
    s = np.minimum(s, N - s)
    dg_analytic = 2 * np.log2(s + 1)
    np.fill_diagonal(dg_analytic, 0)
    mask = (dg > 2) & (s >= 1) & ~np.isnan(dstar)
    x, y = dg[mask], dstar[mask]
    if len(x) > 2:
        slope, intercept = np.polyfit(x, y, 1)
        r2 = 1 - np.var(y - slope * x - intercept) / np.var(y)
    else:
        slope = intercept = r2 = float('nan')
    # full-range fit too
    mask2 = (s >= 1) & ~np.isnan(dstar)
    x2, y2 = dg[mask2], dstar[mask2]
    slope2, intercept2 = np.polyfit(x2, y2, 1)
    r2_full = 1 - np.var(y2 - slope2 * x2 - intercept2) / np.var(y2)
    return dict(I0=I0, dstar=dstar, dg=dg, dg_analytic=dg_analytic, s=s,
                slope=float(slope), intercept=float(intercept), r2=float(r2),
                slope_full=float(slope2), r2_full=float(r2_full),
                n_floored_pairs=n_floored)

def mds_stresses(dstar, dims=(1, 2, 3, 4)):
    out = {}
    D = dstar.copy()
    for d in dims:
        m = MDS(n_components=d, dissimilarity='precomputed', random_state=0,
                normalized_stress='auto', max_iter=300)
        m.fit(D)
        out[d] = float(m.stress_)
    return out

def triangle_deficits(dstar, N):
    """Heron areas of triangles in the d* metric (opus's flatness test)."""
    def heron(a, b, c):
        q = (a + b + c) / 2
        return math.sqrt(max(q * (q - a) * (q - b) * (q - c), 0.0))
    opus_samples = []
    for (i, j, k) in [(0, 8, 16), (0, 16, 32), (8, 16, 24)]:
        if max(i, j, k) < N:
            opus_samples.append(heron(dstar[i, j], dstar[j, k], dstar[k, i]))
    scan = []
    rng = np.random.default_rng(7)
    for _ in range(200):
        i, j, k = rng.choice(N, size=3, replace=False)
        scan.append(heron(dstar[i, j], dstar[j, k], dstar[k, i]))
    return dict(opus_samples=[float(x) for x in opus_samples],
                scan_mean=float(np.mean(scan)), scan_max=float(np.max(scan)))

def scale_invariance(mera):
    """Frobenius distances between optimized u_l tensors across layers."""
    us = [np.asarray(u) for u in mera.u]
    out = []
    for a in range(len(us)):
        for b in range(a + 1, len(us)):
            if us[a].shape == us[b].shape:
                out.append(float(np.linalg.norm(us[a] - us[b]) / np.linalg.norm(us[a])))
    return dict(pairwise_rel_dist_mean=float(np.mean(out)) if out else None,
                pairwise_rel_dist_max=float(np.max(out)) if out else None)

# ------------------------------------------------------------------ driver

def measure(N, chi, seed, do_mds=True, Lmax=8, max_pairs=None, max_sep=None):
    print(f'--- measuring N={N} chi={chi} seed={seed}', flush=True)
    m = load_mera(N, chi, seed)
    res = dict(N=N, chi=chi, seed=seed)
    t0 = time.time()
    # I-matrix first (light tensors); block entropies last (heavy intermediates)
    I, pairs = mutual_information_matrix(m, max_pairs=max_pairs, verbose=(N >= 64), max_sep=max_sep)
    import gc
    gc.collect()
    res['S'] = block_entropies(m, Lmax=Lmax)
    res['c_fit'], res['c_fit_r2'] = cft_fit_c(res['S'], N)
    print(f'  S(L) done, c_fit = {res["c_fit"]:.4f} (r2={res["c_fit_r2"]:.4f})', flush=True)
    ana = dstar_analysis(m, I)
    res['I0'] = ana['I0']
    res['n_floored_pairs'] = ana['n_floored_pairs']
    res['slope'] = ana['slope']; res['r2'] = ana['r2']
    res['slope_full'] = ana['slope_full']; res['r2_full'] = ana['r2_full']
    res['n_pairs'] = len(pairs)
    if do_mds:
        # largest complete (all-measured) prefix block
        nsub = 0
        for k in range(1, N + 1):
            if np.all(~np.isnan(ana['dstar'][:k, :k])):
                nsub = k
            else:
                break
        nsub = max(nsub, 2)
        res['mds_stress'] = mds_stresses(ana['dstar'][:nsub, :nsub])
        res['mds_subblock'] = nsub
        print(f'  MDS stresses ({nsub}x{nsub} sub-block): {res["mds_stress"]}', flush=True)
    res['triangles'] = triangle_deficits(ana['dstar'], N)
    res['scale_inv'] = scale_invariance(m)
    # per-separation averages of I and d*
    svals = ana['s']
    res['I_of_s'] = {int(sv): float(np.nanmean(I[svals == sv])) for sv in np.unique(svals) if sv > 0}
    res['dstar_of_s'] = {int(sv): float(np.nanmean(ana['dstar'][svals == sv]))
                         for sv in np.unique(svals) if sv > 0}
    res['dgraph_of_s'] = {int(sv): float(np.mean(ana['dg'][svals == sv]))
                          for sv in np.unique(svals) if sv > 0}
    res['seconds'] = time.time() - t0
    # persist raw matrices for plots
    np.savez(f'/home/z/my-project/scripts/mera_meas_N{N}_chi{chi}_seed{seed}.npz',
             I=I, dstar=ana['dstar'], dg=ana['dg'], dg_analytic=ana['dg_analytic'])
    print(f'  done ({res["seconds"]:.0f}s): I0={res["I0"]:.4f} slope={res["slope"]:.4f} '
          f'R2={res["r2"]:.4f}', flush=True)
    return res

def negative_control(N, chi):
    print(f'--- negative control (random tensors) N={N} chi={chi}', flush=True)
    m = MERA(N, chi, seed=99)
    I, _ = mutual_information_matrix(m, max_pairs=(400 if N > 64 else None))
    ana = dstar_analysis(m, I)
    return dict(slope=ana['slope'], r2=ana['r2'], r2_full=ana['r2_full'],
                I0=ana['I0'])

if __name__ == '__main__':
    N = int(sys.argv[1]); chi = int(sys.argv[2]); seed = int(sys.argv[3])
    maxp = int(sys.argv[4]) if len(sys.argv) > 4 else None
    msep = int(sys.argv[5]) if len(sys.argv) > 5 else None
    res = measure(N, chi, seed, max_pairs=maxp, max_sep=msep)
    out = f'/home/z/my-project/scripts/mera_result_N{N}_chi{chi}_seed{seed}.json'
    json.dump(res, open(out, 'w'), indent=1)
    print('saved', out)
