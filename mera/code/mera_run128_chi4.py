"""mera_run128_chi4.py — Phase A: N=128, chi=4, FULL deterministic optimization.

Closes the grid hole "(128,4) not executed: jit compile of the full 256-rdm
energy exceeds the session budget". The fix is chunked jitting: the energy is
split into 4 jitted chunk functions (32 bonds + 32 sites each); the full
deterministic gradient is the exact sum of chunk gradients — identical
mathematics to mera_opt.energy_fn (mean over ALL N bonds and ALL N sites),
same ansatz family (Z2-restricted), same Adam + cosine lr schedule.

Also provides the warm start for the chi=8 run (Phase B, mera_run128_chi8.py).

Usage:  python3 mera_run128_chi4.py [steps] [--resume]
Log:    mera128_chi4_log.jsonl   Checkpoints: mera_params_ckpt128_chi4_*.npz
"""
import sys, os, time, json
sys.path.insert(0, '/home/z/my-project/scripts')
os.environ.setdefault('JAX_COMPILATION_CACHE_DIR',
                      '/home/z/my-project/scripts/.jax_cache')
import numpy as np
import jax
jax.config.update('jax_platform_name', 'cpu')
import jax.numpy as jnp
import mera_real
mera_real.set_backend('jax')
from mera_real import MERA
from mera_opt import params_init, tensors_from_params, tfim_exact_energy

N, CHI = 128, 4
STEPS = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 3000
RESUME = '--resume' in sys.argv
WALL = 540
for i, a in enumerate(sys.argv):
    if a == '--wall' and i + 1 < len(sys.argv):
        WALL = int(sys.argv[i + 1])
LOGF = '/home/z/my-project/scripts/mera128_chi4_log.jsonl'
CKPT = '/home/z/my-project/scripts/mera_params_ckpt128_chi4.npz'
FINAL = '/home/z/my-project/scripts/mera_params_N128_chi4_seed0.npz'
RESULT = '/home/z/my-project/scripts/mera_result_N128_chi4_seed0.json'

_H2 = jnp.asarray(-np.kron(np.array([[0, 1], [1, 0]]), np.array([[0, 1], [1, 0]])).astype(complex))
_H1 = jnp.asarray(-np.array([[1, 0], [0, -1]], dtype=complex))

NCHUNK = 4
SIZE = N // NCHUNK  # 32 bonds + 32 sites per chunk


def chunk_terms(c):
    bonds = [(a, (a + 1) % N) for a in range(c * SIZE, (c + 1) * SIZE)]
    sites = list(range(c * SIZE, (c + 1) * SIZE))
    return bonds, sites


def chunk_energy(prm, bonds, sites):
    us, ws, t = tensors_from_params(prm, N, CHI)
    m = MERA(N, CHI, tensors=(us, ws, t))
    e = 0.0
    for (a, b) in bonds:
        e = e + jnp.trace(m.rdm([a, b]) @ _H2)
    for s in sites:
        e = e + jnp.trace(m.rdm([s]) @ _H1)
    # exact contribution of this chunk to E = (1/N)(sum bonds + sum sites)
    return jnp.real(e) / N


print(f'Phase A: N={N} chi={CHI} steps={STEPS} chunks={NCHUNK}x{2*SIZE}terms', flush=True)
t0 = time.time()
CHUNK_VG = []
for c in range(NCHUNK):
    bonds, sites = chunk_terms(c)
    CHUNK_VG.append(jax.jit(jax.value_and_grad(
        lambda p, bs=bonds, ss=sites: chunk_energy(p, bs, ss))))
print(f'jitted {NCHUNK} chunk vg functions in {time.time()-t0:.0f}s', flush=True)


def full_energy_and_grad(prm):
    e, g = 0.0, None
    for vg in CHUNK_VG:
        v, gi = vg(prm)
        e = e + v
        if g is None:
            g = {k: gi[k] for k in gi}
        else:
            for k in g:
                g[k] = g[k] + gi[k]
    return e, g


prm = params_init(N, CHI, 0)
m = {k: jnp.zeros_like(v) for k, v in prm.items()}
v = {k: jnp.zeros_like(v) for k, v in prm.items()}
step0 = 0
if RESUME and os.path.exists(CKPT):
    d = np.load(CKPT)
    prm = {k: jnp.array(d[k]) for k in d.files if not k.startswith('_')}
    m = {k: jnp.array(d['_m_' + k]) for k in prm}
    v = {k: jnp.array(d['_v_' + k]) for k in prm}
    step0 = int(d['_step'].item())
    print(f'resumed from {CKPT} at step {step0}', flush=True)
    STEPS = max(STEPS, step0 + 250)
b1, b2, eps = 0.9, 0.999, 1e-8
hist = []
t_start = time.time()
DONE = True

for t in range(step0 + 1, STEPS + 1):
    if time.time() - t_start > WALL:
        DONE = False
        print(f'WALL LIMIT {WALL}s reached at step {t-1}; checkpoint saved — '
              f'rerun with --resume to continue', flush=True)
        np.savez(CKPT, _step=t - 1,
                 **{**{k: np.asarray(x) for k, x in prm.items()},
                     **{'_m_' + k: np.asarray(x) for k, x in m.items()},
                     **{'_v_' + k: np.asarray(x) for k, x in v.items()}})
        sys.exit(0)
    lr_t = 0.03 * (0.5 * (1 + np.cos(np.pi * min(t / STEPS, 1.0)))) ** 2
    e, g = full_energy_and_grad(prm)
    for k in prm:
        m[k] = b1 * m[k] + (1 - b1) * g[k]
        v[k] = b2 * v[k] + (1 - b2) * g[k] ** 2
        mh = m[k] / (1 - b1 ** t)
        vh = v[k] / (1 - b2 ** t)
        prm[k] = prm[k] - lr_t * mh / (jnp.sqrt(vh) + eps)
    if t == 1 or t % 100 == 0:
        rec = dict(step=t, E=float(e), lr=float(lr_t),
                   elapsed=time.time() - t_start)
        hist.append(rec)
        print(f"  step {t}: E/N = {float(e):.6f} ({rec['elapsed']:.0f}s)", flush=True)
        with open(LOGF, 'a') as f:
            f.write(json.dumps(rec) + '\n')
    if t % 500 == 0 or t == STEPS:
        np.savez(CKPT, _step=t,
                 **{**{k: np.asarray(x) for k, x in prm.items()},
                     **{'_m_' + k: np.asarray(x) for k, x in m.items()},
                     **{'_v_' + k: np.asarray(x) for k, x in v.items()}})

e_final = float(full_energy_and_grad(prm)[0])
ex = tfim_exact_energy(N) / N
dt = time.time() - t_start
print(f'FINAL jax: E/N = {e_final:.6f}  exact = {ex:.6f}  '
      f'err = {abs(e_final-ex):.2e}  rel = {abs(e_final-ex)/abs(ex):.2e}  ({dt:.0f}s)', flush=True)

np.savez(FINAL, **{k: np.asarray(x) for k, x in prm.items()})

# ---- independent numpy cross-check of the final energy (complex128) ----
mera_real.set_backend('np')
us, ws, tt = tensors_from_params(prm, N, CHI)
mm = MERA(N, CHI, tensors=([np.asarray(x) for x in us],
                           [np.asarray(x) for x in ws], np.asarray(tt)))
_H2n = -np.kron(np.array([[0, 1], [1, 0]]), np.array([[0, 1], [1, 0]])).astype(complex)
_H1n = -np.array([[1, 0], [0, -1]], dtype=complex)
e_np = 0.0
for a in range(N):
    e_np += np.real(np.trace(mm.rdm([a, (a + 1) % N]) @ _H2n))
    e_np += np.real(np.trace(mm.rdm([a]) @ _H1n))
e_np /= N
print(f'numpy cross-check (c128, all {2*N} terms): E/N = {e_np:.6f}  '
      f'|jax - numpy| = {abs(e_np - e_final):.2e}', flush=True)

res = dict(N=N, chi=CHI, seed=0, steps=STEPS, E_over_N=e_final,
           E_over_N_numpy=float(e_np), exact=ex,
           rel_err=abs(e_final - ex) / abs(ex), seconds=dt,
           method='chunked-jit full deterministic gradient (4 chunks x 64 terms), '
                  'same objective/ansatz/optimizer as mera_opt.optimize',
           history=hist)
with open(RESULT, 'w') as f:
    json.dump(res, f, indent=1)
print('saved', RESULT, 'and', FINAL)
