"""plots_mera.py — final figures for the real-contraction MERA program."""
import sys, json
sys.path.insert(0, '/home/z/my-project/scripts')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import glob
for f in glob.glob('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'):
    try:
        fm.fontManager.addfont(f)
    except Exception:
        pass
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'font.size': 11, 'figure.dpi': 150})

OUT = '/home/z/my-project/repos/universal-consciousness-mathematical/mera/figures'
import os
os.makedirs(OUT, exist_ok=True)

def load_res(N, chi, seed=0):
    return json.load(open(f'/home/z/my-project/scripts/mera_result_N{N}_chi{chi}_seed{seed}.json'))
def load_meas(N, chi, seed=0):
    d = np.load(f'/home/z/my-project/scripts/mera_meas_N{N}_chi{chi}_seed{seed}.npz')
    return d

# ---------------- Fig 1: d* vs d_graph (primary N=64 chi=4)
d = load_meas(64, 4, 0)
dstar, dg, s = d['dstar'], d['dg'], None
idx = np.arange(64)
ss = np.abs(idx[:, None] - idx[None, :]); ss = np.minimum(ss, 64 - ss)
mask = (dg > 2) & (ss >= 1) & ~np.isnan(dstar)
x, y = dg[mask], dstar[mask]
sl, ic = np.polyfit(x, y, 1)
r2 = 1 - np.var(y - sl * x - ic) / np.var(y)
fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
ax.scatter(x, y, s=10, alpha=0.45, color='#1f77b4', label=f'real contraction (N=64, $\\chi$=4)')
xx = np.linspace(x.min(), x.max(), 100)
ax.plot(xx, sl * xx + ic, 'r-', lw=2,
        label=f'fit: slope={sl:.3f}, $R^2$={r2:.3f}')
ax.plot(xx, 0.69 * xx, 'k--', lw=1.5,
        label='opus claim: slope 0.69, $R^2$=0.991 (analytic stand-in)')
ax.set_xlabel('graph distance $d_{graph}$ [legs, BFS on real network]')
ax.set_ylabel('$d^{*}=-\\log(I^{*}/I_0)$')
ax.set_title('Poset-to-metric test — REAL contraction (critical TFIM, optimized MERA)')
ax.legend(loc='upper left', fontsize=9)
fig.savefig(f'{OUT}/dstar_vs_dgraph_N64.png')
plt.close(fig)

# ---------------- Fig 2: S(L) with CFT fits
fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
from mera_opt import tfim_ed
from mera_real import entropy, brute_rdm
E0, psi = tfim_ed(16)
ls = np.arange(1, 9)
sed = []
for L in ls:
    vals = [entropy(brute_rdm(psi, 16, [(st + i) % 16 for i in range(L)])) for st in range(4)]
    sed.append(np.mean(vals))
ax.plot(ls, sed, 'ko-', label='exact diagonalization (N=16)', ms=5)
for (N, chi, c) in [(16, 4, '#1f77b4'), (32, 4, '#2ca02c'), (64, 4, '#d62728')]:
    r = load_res(N, chi)
    Ls = sorted(int(k) for k in r['S'].keys())
    Sv = [r['S'][str(k)] if str(k) in r['S'] else r['S'][k] for k in Ls]
    ax.plot(Ls, Sv, 's-', color=c, label=f'MERA N={N}, $\\chi$={4 if chi==4 else chi}: c={r["c_fit"]:.3f}', ms=5)
ax.set_xlabel('block size L')
ax.set_ylabel('S(L) [nats]')
ax.set_title('Block-entropy scaling: real MERA vs exact (c = 1/2 line)')
ax.legend(fontsize=9)
fig.savefig(f'{OUT}/S_of_L.png')
plt.close(fig)

# ---------------- Fig 3: I(s) curves
fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
ax.plot(range(1, 9), [v for v in [np.mean([entropy(brute_rdm(psi,16,[i,(i+s)%16])) for i in range(4)]) for s in range(1,9)]], 'ko-', label='exact ED I(s) (N=16)', ms=4)
r16 = load_res(16, 4); r32 = load_res(32, 4); r64 = load_res(64, 4)
for r, N, c in [(r16, 16, '#1f77b4'), (r32, 32, '#2ca02c'), (r64, 64, '#d62728')]:
    ks = sorted(int(k) for k in r['I_of_s'].keys())
    vals = [r['I_of_s'][str(k)] if str(k) in r['I_of_s'] else r['I_of_s'][k] for k in ks]
    ax.plot(ks, vals, 's-', color=c, label=f'MERA N={N}', ms=4)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('separation s')
ax.set_ylabel('mutual information I(s) [nats]')
ax.set_title('Mutual information decay (critical TFIM)')
ax.legend(fontsize=9)
fig.savefig(f'{OUT}/I_of_s.png')
plt.close(fig)

# ---------------- Fig 4: MDS 2D embedding (N=32)
from sklearn.manifold import MDS
d32 = load_meas(32, 4, 0)
D = d32['dstar'].copy()
finite = D[np.isfinite(D)]
D = np.where(np.isfinite(D), D, finite.max())   # fill 20 cap-excluded far pairs for display
m = MDS(n_components=2, dissimilarity='precomputed', random_state=0, normalized_stress='auto', max_iter=400)
X = m.fit_transform(D)
fig, ax = plt.subplots(figsize=(5.6, 4.8), constrained_layout=True)
sc = ax.scatter(X[:, 0], X[:, 1], c=np.arange(32), cmap='viridis', s=45)
plt.colorbar(sc, ax=ax, label='site index i (chain order)')
# connect neighbors to reveal the circle
order = np.argsort(np.arctan2(X[:, 1], X[:, 0]))
ax.plot(X[order, 0], X[order, 1], 'k-', alpha=0.35, lw=1)
ax.set_xlabel('MDS dimension 1'); ax.set_ylabel('MDS dimension 2')
ax.set_title(f'MDS embedding of $d^*$ (N=32): stress(1)={198.9:.0f} >> stress(2)={6.80:.2f}\n'
             'the periodic chain is a CIRCLE, not a line')
fig.savefig(f'{OUT}/mds_N32.png')
plt.close(fig)

# ---------------- Fig 5: energy errors across the grid
grid = json.load(open('/home/z/my-project/scripts/mera_grid_results.json'))
fig, ax = plt.subplots(figsize=(6.0, 4.2), constrained_layout=True)
for chi, c, mk in [(2, '#1f77b4', 'o'), (4, '#d62728', 's')]:
    pts = [(v['N'], v['rel_err']) for k, v in grid.items() if v['chi'] == chi]
    pts.sort()
    ax.semilogy([p[0] for p in pts], [p[1] for p in pts], mk + '-', color=c,
                label=f'$\\chi$={chi}', ms=7)
ax.set_xlabel('system size N'); ax.set_ylabel('variational energy error |E/N - E0/N|')
ax.set_title('Real MERA optimization quality (critical TFIM)\nexact E0/N from free fermions + ED cross-check')
ax.legend()
fig.savefig(f'{OUT}/energy_errors.png')
plt.close(fig)
print('figures saved to', OUT)
