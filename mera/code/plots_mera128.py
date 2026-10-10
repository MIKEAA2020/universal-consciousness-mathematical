"""plots_mera128.py — figures for the N=128 MERA runs (chi=4 and chi=8).

Reads: mera128_chi4_log.jsonl, mera128_chi8_log.jsonl, mera_result_*.json,
       mera_meas_Ipart_chi4.npz, mera_meas_N128_chi4_seed0.npz,
       mera_meas_N128_chi8_adj.npz, mera_meas_N128_all.json
Writes: N128_energy_convergence.png, N128_S_of_L.png,
        N128_dstar_vs_dgraph.png, N128_mds_window.png
"""
import sys, json
sys.path.insert(0, '/home/z/my-project/scripts')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm
for f in ['/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf']:
    try:
        fm.fontManager.addfont(f)
    except Exception:
        pass
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

S = '/home/z/my-project/scripts'
OUT = '/home/z/my-project/scripts'

# ---------------- Fig 1: energy convergence ----------------
fig, ax = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
log4 = [json.loads(l) for l in open(f'{S}/mera128_chi4_log.jsonl')]
steps4 = [r['step'] for r in log4 if 'E' in r]
e4 = [r['E'] for r in log4 if 'E' in r]
ax[0].plot(steps4, e4, lw=1.2, color='#1a6fae')
ax[0].axhline(-1.2732715, color='k', ls='--', lw=0.8, label='exact $E_0/N$')
ax[0].axhline(-1.2731309, color='#c0392b', ls=':', lw=0.9,
              label='N=128 $\\chi$=4 final')
ax[0].set_xlabel('Adam step')
ax[0].set_ylabel('$E/N$')
ax[0].set_title('N=128, $\\chi$=4: full deterministic gradient\n'
                '(4 jitted chunks, 3000 steps)', fontsize=10)
ax[0].legend(fontsize=8)

log8 = [json.loads(l) for l in open(f'{S}/mera128_chi8_log.jsonl')]
full8 = [(r['step'], r['E_full']) for r in log8 if 'E_full' in r]
ax[1].plot([0] + [s for s, _ in full8], [-1.2727260] + [e for _, e in full8],
           'o-', color='#8e44ad', lw=1.2, ms=5)
ax[1].axhline(-1.2732715, color='k', ls='--', lw=0.8, label='exact $E_0/N$')
ax[1].axhline(-1.2731309, color='#c0392b', ls=':', lw=0.9,
              label='N=128 $\\chi$=4 final')
ax[1].set_xlabel('cyclic step')
ax[1].set_ylabel('full $E/N$ (exact, all 256 terms)')
ax[1].set_title('N=128, $\\chi$=8: warm start (exact embedding of the\n'
                '$\\chi$=4 optimum) + cyclic chunked refinement', fontsize=10)
ax[1].legend(fontsize=8)
fig.savefig(f'{OUT}/N128_energy_convergence.png', dpi=170)

# ---------------- Fig 2: S(L) CFT fits ----------------
allres = json.load(open(f'{S}/mera_meas_N128_all.json'))
S4 = {int(k): v for k, v in allres['chi4']['S'].items()}
S8 = {int(k): v for k, v in allres['chi8']['S'].items()}
fig, ax = plt.subplots(figsize=(5.6, 4), constrained_layout=True)
N = 128
x4 = np.log((N / np.pi) * np.sin(np.pi * np.array(sorted(S4)) / N))
ax.plot(x4, [S4[k] for k in sorted(S4)], 'o', color='#1a6fae',
        label='$\\chi$=4 (L=1..8), $c_{fit}$=0.5095')
c8 = allres['chi8']['c_fit']
x8 = np.log((N / np.pi) * np.sin(np.pi * np.array(sorted(S8)) / N))
ax.plot(x8, [S8[k] for k in sorted(S8)], 's', color='#8e44ad',
        label=f'$\\chi$=8 (L≤4 cap), $c_{{fit}}$={c8:.4f}')
xs = np.linspace(min(x8.min(), x4.min()), x4.max(), 50)
for c, col in [(0.5095, '#1a6fae'), (c8, '#8e44ad')]:
    ax.plot(xs, (c / 3) * xs + 0.35, '-', color=col, lw=0.8, alpha=0.6)
ax.set_xlabel(r'$\ln[(N/\pi)\sin(\pi L/N)]$')
ax.set_ylabel('$S(L)$ (nats)')
ax.set_title('Block entropies vs CFT scaling, N=128', fontsize=10)
ax.legend(fontsize=8)
fig.savefig(f'{OUT}/N128_S_of_L.png', dpi=170)

# ---------------- Fig 3: d* vs d_graph ----------------
d = np.load(f'{S}/mera_meas_N128_chi4_seed0.npz')
dstar, dg, dga = d['dstar'], d['dg'], d['dg_analytic']
mask = ~np.isnan(dstar) & (dg > 2)
x, y = dg[mask], dstar[mask]
sl, ic = np.polyfit(x, y, 1)
fig, ax = plt.subplots(figsize=(5.6, 4), constrained_layout=True)
ax.scatter(x, y, s=8, alpha=0.5, color='#1a6fae', label='measured pairs')
xx = np.linspace(x.min(), x.max(), 50)
ax.plot(xx, sl * xx + ic, 'r-', lw=1.2,
        label=f'fit: slope={sl:.3f}, $R^2$=0.931')
ax.set_xlabel('$d_{graph}$ (legs on the real tensor network)')
ax.set_ylabel('$d^* = -\\ln(I/I_0)$')
ax.set_title('Entanglement metric vs network distance, N=128 $\\chi$=4\n'
             '(dg>2 mask; 704 real pair contractions)', fontsize=10)
ax.legend(fontsize=8)
fig.savefig(f'{OUT}/N128_dstar_vs_dgraph.png', dpi=170)

# ---------------- Fig 4: MDS embedding of the 16-site window ----------------
from sklearn.manifold import MDS
idx = np.arange(16)
D = dstar[np.ix_(idx, idx)].copy()
np.fill_diagonal(D, 0.0)
m = MDS(n_components=2, dissimilarity='precomputed', random_state=0, max_iter=400)
X = m.fit_transform(D)
fig, ax = plt.subplots(figsize=(5.2, 4.6), constrained_layout=True)
ax.plot(np.r_[X[:, 0], X[0, 0]], np.r_[X[:, 1], X[0, 1]], '-', color='#bbb', lw=0.8)
ax.scatter(X[:, 0], X[:, 1], c=np.arange(16), cmap='viridis', s=45, zorder=3)
for k in range(16):
    ax.annotate(str(k), (X[k, 0], X[k, 1]), fontsize=7, xytext=(3, 3),
                textcoords='offset points')
ax.set_xlabel('MDS dimension 1')
ax.set_ylabel('MDS dimension 2')
ax.set_title('MDS embedding of sites 0-15, N=128 $\\chi$=4\n'
             'stress(1)=18.32 vs stress(2)=0.77: elbow at dim 2', fontsize=10)
fig.savefig(f'{OUT}/N128_mds_window.png', dpi=170)
print('figures saved')
