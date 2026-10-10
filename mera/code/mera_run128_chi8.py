"""mera_run128_chi8.py — Phase B v3: N=128, chi=8, per-step CYCLIC chunks.

v1/v2 lessons (recorded for the report):
  - 16 simultaneously-live jitted chunk programs need ~3 GB -> OOM. Fix:
    build ONE chunk program per step, then del + gc.collect() (RSS ~0.9 GB).
  - block-coordinate phases (K consecutive steps on one chunk) OVERFIT the
    16 terms of the phase and degrade the rest (measured: subset positions
    improved ~0.1 while the full energy degraded +6.1e-3). Fix: per-step
    cyclic schedule — every step uses the next chunk; Adam's momentum
    smooths across consecutive chunks, approximating the full gradient.

Warm start = exact embedding of the optimized N=128 chi=4 state (+0.03
normalized coupling, +1e-4 degeneracy guard) — see mera_embed8.py.

Design:
  - per step t: chunk c = (t-1) mod 16 (8 bonds + 8 sites), exact gradient;
  - Adam moments persist across steps (checkpointed);
  - full exact numpy (complex128) energy over all 256 terms at milestones,
    with best-checkpointing keyed on the FULL energy (best initialized to the
    warm-start energy, so the delivered state is never worse than the start);
  - wall-clock limited per process; --resume continues from checkpoint.

Usage: python3 mera_run128_chi8.py [steps] [--resume] [--finish] [--wall S]
"""
import sys, os, time, json, gc
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
from mera_opt import tensors_from_params, tfim_exact_energy

N, CHI = 128, 8
STEPS = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 280
RESUME = '--resume' in sys.argv
FINISH = '--finish' in sys.argv or True  # finish logic at STEPS
WALL = 535
for i, a in enumerate(sys.argv):
    if a == '--wall' and i + 1 < len(sys.argv):
        WALL = int(sys.argv[i + 1])
LOGF = '/home/z/my-project/scripts/mera128_chi8_log.jsonl'
CKPT = '/home/z/my-project/scripts/mera_params_ckpt128_chi8.npz'
BEST = '/home/z/my-project/scripts/mera_params_ckpt128_chi8_best.npz'
FINAL = '/home/z/my-project/scripts/mera_params_N128_chi8_seed0.npz'
RESULT = '/home/z/my-project/scripts/mera_result_N128_chi8_seed0.json'
WARM = '/home/z/my-project/scripts/mera_params_N128_chi8_warmstart.npz'
MILESTONE = 70
LR0 = 0.010

_H2 = jnp.asarray(-np.kron(np.array([[0, 1], [1, 0]]), np.array([[0, 1], [1, 0]])).astype(complex))
_H1 = jnp.asarray(-np.array([[1, 0], [0, -1]], dtype=complex))
_H2n = -np.kron(np.array([[0, 1], [1, 0]]), np.array([[0, 1], [1, 0]])).astype(complex)
_H1n = -np.array([[1, 0], [0, -1]], dtype=complex)
WARM_E = -1.272726008   # measured in mera_embed8.py


def chunk_energy(prm, bonds, sites):
    us, ws, t = tensors_from_params(prm, N, CHI)
    m = MERA(N, CHI, tensors=(us, ws, t))
    e = 0.0
    for (a, b) in bonds:
        e = e + jnp.trace(m.rdm([a, b]) @ _H2)
    for s in sites:
        e = e + jnp.trace(m.rdm([s]) @ _H1)
    return jnp.real(e) / N


def full_numpy_energy(prm):
    mera_real.set_backend('np')
    us, ws, tt = tensors_from_params(prm, N, CHI)
    mm = MERA(N, CHI, tensors=([np.asarray(x) for x in us],
                               [np.asarray(x) for x in ws], np.asarray(tt)))
    e = 0.0
    for a in range(N):
        e += np.real(np.trace(mm.rdm([a, (a + 1) % N]) @ _H2n))
        e += np.real(np.trace(mm.rdm([a]) @ _H1n))
    mera_real.set_backend('jax')
    return float(e / N)


def save_ckpt(path, step, best_e, prm, m, v):
    np.savez(path, _step=step, _best=best_e,
             **{**{k: np.asarray(x) for k, x in prm.items()},
                 **{'_m_' + k: np.asarray(x) for k, x in m.items()},
                 **{'_v_' + k: np.asarray(x) for k, x in v.items()}})


print(f'Phase B v3 cyclic: N={N} chi={CHI} steps={STEPS} lr0={LR0} '
      f'milestones every {MILESTONE}', flush=True)
if os.path.exists(CKPT) and RESUME:
    d = np.load(CKPT)
    prm = {k: jnp.array(d[k]) for k in d.files if not k.startswith('_')}
    m = {k: jnp.array(d['_m_' + k]) for k in prm}
    v = {k: jnp.array(d['_v_' + k]) for k in prm}
    step0 = int(d['_step'].item())
    best_e = float(d['_best'].item())
    print(f'resumed at step {step0} (best full E {best_e:.6f})', flush=True)
else:
    d = np.load(WARM)
    prm = {k: jnp.array(d[k]) for k in d.files}
    m = {k: jnp.zeros_like(x) for k, x in prm.items()}
    v = {k: jnp.zeros_like(x) for k, x in prm.items()}
    step0, best_e = 0, WARM_E
    np.savez(BEST, **{k: np.asarray(x) for k, x in prm.items()})
    print(f'fresh start from warm (E = {WARM_E:.6f}); best initialized to warm',
          flush=True)

b1, b2, eps = 0.9, 0.999, 1e-8
hist = []
t_start = time.time()
DONE = False
t = step0
while t < STEPS:
    if time.time() - t_start > WALL:
        print(f'WALL LIMIT at step {t}; checkpoint saved — rerun with --resume',
              flush=True)
        save_ckpt(CKPT, t, best_e, prm, m, v)
        sys.exit(0)
    t += 1
    c = (t - 1) % 16
    bonds = [(a, (a + 1) % N) for a in range(c * 8, (c + 1) * 8)]
    sites = list(range(c * 8, (c + 1) * 8))
    vg = jax.jit(jax.value_and_grad(
        lambda p, bs=bonds, ss=sites: chunk_energy(p, bs, ss)))
    val, g = vg(prm)
    g['u3'].block_until_ready()
    lr_t = LR0 * (0.5 * (1 + np.cos(np.pi * min(t / STEPS, 1.0)))) ** 2
    for k in prm:
        m[k] = b1 * m[k] + (1 - b1) * g[k]
        v[k] = b2 * v[k] + (1 - b2) * g[k] ** 2
        mh = m[k] / (1 - b1 ** t)
        vh = v[k] / (1 - b2 ** t)
        prm[k] = prm[k] - lr_t * mh / (jnp.sqrt(vh) + eps)
    del vg, g
    gc.collect()
    if t <= 2 or t % 25 == 0:
        rec = dict(step=t, chunk=c, E_chunk=float(val), lr=round(float(lr_t), 5),
                   elapsed=round(time.time() - t_start, 1))
        print(f"  step {t} chunk {c}: E_c = {float(val):.6f}", flush=True)
        with open(LOGF, 'a') as f:
            f.write(json.dumps(rec) + '\n')
    if t % MILESTONE == 0 or t == STEPS:
        efull = full_numpy_energy(prm)
        hist.append(dict(step=t, E_full=efull))
        print(f'  >>> step {t}: FULL E/N = {efull:.6f} (best {best_e:.6f})',
              flush=True)
        with open(LOGF, 'a') as f:
            f.write(json.dumps(dict(step=t, E_full=efull)) + '\n')
        if efull < best_e:
            best_e = efull
            np.savez(BEST, **{k: np.asarray(x) for k, x in prm.items()})
            print(f'      new best saved', flush=True)
    if t % 35 == 0:
        save_ckpt(CKPT, t, best_e, prm, m, v)

save_ckpt(CKPT, t, best_e, prm, m, v)

if t >= STEPS:
    e_final = full_numpy_energy(prm)
    ex = tfim_exact_energy(N) / N
    # deliver the best-measured params (never worse than warm start)
    db = np.load(BEST)
    best_params = {k: db[k] for k in db.files}
    e_best = full_numpy_energy({k: jnp.array(x) for k, x in best_params.items()})
    np.savez(FINAL, **best_params)
    res = dict(N=N, chi=CHI, seed=0, steps=STEPS,
               E_over_N_final=e_final, E_over_N_best=float(e_best),
               E_warmstart=WARM_E, E_exact=ex,
               rel_err_best=abs(e_best - ex) / abs(ex),
               rel_err_final=abs(e_final - ex) / abs(ex),
               full_energy_history=hist,
               method='warm start = exact embedding of the optimized N=128 '
                      'chi=4 state (+0.03 normalized coupling, 1e-4 degeneracy '
                      'guard); per-step CYCLIC chunked gradients (16 chunks x '
                      '16 terms; one live program per step); Adam lr 0.010 '
                      'cosine-sq; exact full numpy c128 energy at milestones; '
                      'best-checkpoint keyed on full energy (init = warm)')
    with open(RESULT, 'w') as f:
        json.dump(res, f, indent=1)
    print(f'FINAL full E/N = {e_final:.6f}; BEST = {e_best:.6f}; exact = {ex:.6f}; '
          f'rel(best) = {abs(e_best-ex)/abs(ex):.2e}', flush=True)
    print('saved', FINAL, 'and', RESULT, flush=True)
