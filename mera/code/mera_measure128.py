"""mera_measure128.py — measurement suite for the N=128 MERA runs.

Hardware reality (disclosed): a two-site RDM at separation s needs a causal
cone of ~6 sites per un-merged level; at chi=8 that is 12 legs of dim 8 =
8^12 entries > memory cap -> only s <= 2 pairs are contractible at chi=8 on
this machine. The full suite (all separations) therefore runs on the SOURCE
N=128 chi=4 optimized state (the lineage of the delivered chi=8 state: its
exact embedding + refinement); chi=8-specific adjacent observables (I0, I(s),
S(L) up to the cap, scale invariance) are measured on the DELIVERED chi=8
state directly.

Sampling for the chi=4 full suite (8128 pairs unaffordable; all measured
entries are real contractions):
  - ALL 128 adjacent pairs (s=1): I0 exact;
  - all 120 pairs within the contiguous window [0..15]: full 16x16 block for
    MDS;
  - 8 stratified random pairs per separation s = 2..64 (504 pairs).

Negative control: random (unoptimized) chi=4 MERA at N=128, same analysis.
"""
import sys, os, json, time, gc, ctypes
sys.path.insert(0, '/home/z/my-project/scripts')
import numpy as np
import jax.numpy as jnp
import mera_real
mera_real.set_backend('np')
from mera_real import MERA, entropy, MemoryCapError
from mera_opt import tensors_from_params
from measure_mera import (dstar_analysis, mds_stresses, triangle_deficits,
                          block_entropies, cft_fit_c, scale_invariance)

N = 128


def trim():
    gc.collect()
    try:
        ctypes.CDLL('libc.so.6').malloc_trim(0)
    except Exception:
        pass


def load_mera(path, chi):
    d = np.load(path)
    prm = {k: jnp.array(d[k]) for k in d.files}
    us, ws, t = tensors_from_params(prm, N, chi)
    return MERA(N, chi, tensors=([np.asarray(x) for x in us],
                                 [np.asarray(x) for x in ws], np.asarray(t)))


def pair_plan(rng):
    pairs = set()
    for i in range(N):
        pairs.add((i, (i + 1) % N))
    for i in range(16):
        for j in range(i + 1, 16):
            pairs.add((i, j))
    for s in range(2, N // 2 + 1):
        cand = [(i, (i + s) % N) for i in range(N)]
        idx = rng.choice(len(cand), size=8, replace=False)
        for x in idx:
            i, j = cand[x]
            pairs.add((min(i, j), max(i, j)))
    return sorted(pairs)


def imatrix(m, pairs, verbose_every=200, part=None, I_init=None, start=0):
    """Resumable: saves partial I to mera_meas_Ipart_{part}.npz every 150 pairs."""
    I = np.full((N, N), np.nan) if I_init is None else I_init
    ck = f'/home/z/my-project/scripts/mera_meas_Ipart_{part}.npz'
    if start == 0:
        Ssite = {i: entropy(m.rdm([i])) for i in range(N)}
        np.savez(ck, I=I, site_S=np.array([Ssite[i] for i in range(N)]), upto=0)
    else:
        d = np.load(ck)
        ss = d['site_S']
        Ssite = {i: float(ss[i]) for i in range(N)}
    trim()
    t0, n_skip = time.time(), 0
    for idx in range(start, len(pairs)):
        i, j = pairs[idx]
        try:
            rij = m.rdm([i, j])
            I[i, j] = I[j, i] = Ssite[i] + Ssite[j] - entropy(rij)
            del rij
        except MemoryCapError:
            I[i, j] = I[j, i] = np.nan
            n_skip += 1
        trim()
        if (idx + 1) % 150 == 0:
            np.savez(ck, I=I, site_S=np.array([Ssite[i] for i in range(N)]),
                     upto=idx + 1)
            print(f'    [saved partial to {idx+1}/{len(pairs)}]', flush=True)
        if idx % verbose_every == 0:
            print(f'    pair {idx}/{len(pairs)} ({time.time()-t0:.0f}s, '
                  f'{n_skip} skipped)', flush=True)
    np.savez(ck, I=I, site_S=np.array([Ssite[i] for i in range(N)]), upto=len(pairs))
    print(f'  I-matrix: complete, {len(pairs)} pairs processed', flush=True)
    return I


def analyze(m, I, tag, do_mds=True):
    res = dict(tag=tag)
    S = block_entropies(m, Lmax=8)
    res['S'] = S
    res['c_fit'], res['c_fit_r2'] = cft_fit_c(S, N)
    print(f'  S(L) -> c_fit = {res["c_fit"]:.4f} (r2={res["c_fit_r2"]:.4f})',
          flush=True)
    ana = dstar_analysis(m, I)
    res['I0'] = ana['I0']
    res['n_floored_pairs'] = ana['n_floored_pairs']
    res['slope'] = ana['slope']; res['r2'] = ana['r2']
    res['slope_full'] = ana['slope_full']; res['r2_full'] = ana['r2_full']
    print(f'  I0 = {res["I0"]:.4f}  slope = {res["slope"]:.4f} (r2={res["r2"]:.4f}) '
          f'slope_full = {res["slope_full"]:.4f} (r2={res["r2_full"]:.4f})', flush=True)
    if do_mds:
        nsub = 0
        for k in range(1, N + 1):
            if np.all(~np.isnan(ana['dstar'][:k, :k])):
                nsub = k
            else:
                break
        nsub = max(nsub, 2)
        res['mds_stress'] = mds_stresses(ana['dstar'][:nsub, :nsub])
        res['mds_subblock'] = nsub
        print(f'  MDS stresses ({nsub}x{nsub}): {res["mds_stress"]}', flush=True)
    res['triangles'] = triangle_deficits(ana['dstar'], N)
    res['scale_inv'] = scale_invariance(m)
    svals = ana['s']
    res['I_of_s'] = {int(sv): float(np.nanmean(I[svals == sv]))
                     for sv in np.unique(svals) if sv > 0}
    res['dstar_of_s'] = {int(sv): float(np.nanmean(ana['dstar'][svals == sv]))
                         for sv in np.unique(svals) if sv > 0}
    res['dgraph_of_s'] = {int(sv): float(np.mean(ana['dg'][svals == sv]))
                          for sv in np.unique(svals) if sv > 0}
    np.savez(f'/home/z/my-project/scripts/mera_meas_N128_{tag}.npz',
             I=I, dstar=ana['dstar'], dg=ana['dg'], dg_analytic=ana['dg_analytic'])
    return res


def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else 'all'
    out = {}

    if stage in ('all', 'chi8'):
        # ---------- chi=8 delivered state: adjacent + s=2 observables ----------
        print('=== delivered N=128 chi=8 state: s<=2 pairs + S(L) + scale inv ===',
              flush=True)
        m8 = load_mera('/home/z/my-project/scripts/mera_params_N128_chi8_seed0.npz', 8)
        t0 = time.time()
        pairs8 = sorted(set([(i, (i + 1) % N) for i in range(N)] +
                            [(i, (i + 2) % N) for i in range(N)]))
        I8 = imatrix(m8, pairs8, verbose_every=128, part='chi8')
        S8 = block_entropies(m8, Lmax=8)
        c8, c8r2 = cft_fit_c(S8, N)
        adj8 = [I8[i, (i + 1) % N] for i in range(N)]
        s2_8 = [I8[i, (i + 2) % N] for i in range(N)]
        out['chi8'] = dict(
            tag='chi8_seed0', S=S8, c_fit=c8, c_fit_r2=c8r2,
            I0=float(np.nanmax(adj8)), I0_mean=float(np.nanmean(adj8)),
            I_s1_mean=float(np.nanmean(adj8)), I_s2_mean=float(np.nanmean(s2_8)),
            n_pairs=len(pairs8), scale_inv=scale_invariance(m8),
            note='chi=8 cones: only s<=2 pair RDMs fit the memory cap; S(L) '
                 'limited to L<=4 by the same cap (code skips overflowing L)',
            seconds=time.time() - t0)
        print(f'  chi8: I0 = {out["chi8"]["I0"]:.4f}  I(s=2) = '
              f'{out["chi8"]["I_s2_mean"]:.4f}  c_fit = {c8:.4f} (r2={c8r2:.4f})',
              flush=True)
        np.savez('/home/z/my-project/scripts/mera_meas_N128_chi8_adj.npz', I=I8)
        json.dump(out['chi8'], open('/home/z/my-project/scripts/'
                'mera_meas_N128_chi8_part.json', 'w'), indent=1)

    if stage in ('all', 'chi4pairs'):
        # ---------- chi=4 source state: pair contractions (resumable) ----------
        print('=== source N=128 chi=4 state: pair contractions ===', flush=True)
        m4 = load_mera('/home/z/my-project/scripts/mera_params_N128_chi4_seed0.npz', 4)
        rng = np.random.default_rng(12345)
        pairs = pair_plan(rng)
        print(f'pair plan: {len(pairs)} pairs', flush=True)
        ck = '/home/z/my-project/scripts/mera_meas_Ipart_chi4.npz'
        start = int(np.load(ck)['upto']) if os.path.exists(ck) else 0
        I_init = np.load(ck)['I'] if start > 0 else None
        if start >= len(pairs):
            print('  already complete', flush=True)
        else:
            imatrix(m4, pairs, part='chi4', I_init=I_init, start=start)

    if stage in ('all', 'analyze'):
        # ---------- analysis + negative control ----------
        out = {}
        try:
            out['chi8'] = json.load(open('/home/z/my-project/scripts/'
                                         'mera_meas_N128_chi8_part.json'))
        except Exception:
            pass
        print('=== chi=4 analysis ===', flush=True)
        m4 = load_mera('/home/z/my-project/scripts/mera_params_N128_chi4_seed0.npz', 4)
        rng = np.random.default_rng(12345)
        pairs = pair_plan(rng)
        I4 = np.load('/home/z/my-project/scripts/mera_meas_Ipart_chi4.npz')['I']
        t0 = time.time()
        out['chi4'] = analyze(m4, I4, 'chi4_seed0')
        out['chi4']['n_pairs'] = len(pairs)
        out['chi4']['seconds'] = time.time() - t0

        print('=== negative control: random chi=4 tensors at N=128 ===', flush=True)
        mr = MERA(N, 4, seed=99)
        ctrl_pairs = set([(i, (i + 1) % N) for i in range(N)])
        rng2 = np.random.default_rng(777)
        for s in range(2, N // 2 + 1):
            cand = [(i, (i + s) % N) for i in range(N)]
            for x in rng2.choice(len(cand), size=4, replace=False):
                i, j = cand[x]
                ctrl_pairs.add((min(i, j), max(i, j)))
        ck2 = '/home/z/my-project/scripts/mera_meas_Ipart_neg.npz'
        start = int(np.load(ck2)['upto']) if os.path.exists(ck2) else 0
        I_init = np.load(ck2)['I'] if start > 0 else None
        if start >= len(ctrl_pairs):
            Ic = np.load(ck2)['I']
        else:
            Ic = imatrix(mr, sorted(ctrl_pairs), part='neg', I_init=I_init, start=start)
        anac = dstar_analysis(mr, Ic)
        out['negctl'] = dict(slope=anac['slope'], r2=anac['r2'],
                             r2_full=anac['r2_full'], I0=anac['I0'],
                             n_pairs=len(ctrl_pairs))
        print(f'  negctl: I0 = {out["negctl"]["I0"]:.4f}  slope = '
              f'{out["negctl"]["slope"]:.4f}  r2 = {out["negctl"]["r2"]:.4f}',
              flush=True)

        with open('/home/z/my-project/scripts/mera_meas_N128_all.json', 'w') as f:
            json.dump(out, f, indent=1)
        print('saved mera_meas_N128_all.json')
        print('DONE')


if __name__ == '__main__':
    main()
