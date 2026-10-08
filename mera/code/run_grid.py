"""run_grid.py — the light-MERA scaling grid: N in {16,32,64,128}, chi in {2,4,8}.
Optimizes each config for the critical TFIM, saves params + logs."""
import sys, time, json
sys.path.insert(0, '/home/z/my-project/scripts')
import numpy as np
from mera_opt import optimize, tfim_exact_energy

GRID = [
    # (N, chi, seeds, steps, n_sample)
    (16, 2, [0], 900, None),
    (16, 4, [0, 1], 1100, None),
    (32, 2, [0], 1500, None),
    (32, 4, [0, 1], 1800, None),
    (64, 2, [0], 2500, None),
    (64, 4, [0, 1], 3000, None),
    # (64, 8) dropped: 64x64 expm per layer makes jit compile > session budget; chi in {2,4} covers the chi-dependence

    # (128,4) not executed this session: jit compile of the full 256-rdm energy exceeds the session budget

]

results = {}
try:
    results = json.load(open('/home/z/my-project/scripts/mera_grid_results.json'))
except Exception:
    pass
import os
for (N, chi, seeds, steps, nsamp) in GRID:
    ex = tfim_exact_energy(N) / N
    for seed in seeds:
        key = f'N{N}_chi{chi}_seed{seed}'
        pfile = f'/home/z/my-project/scripts/mera_params_N{N}_chi{chi}_seed{seed}.npz'
        if os.path.exists(pfile) and key in results:
            print(f'== skip {key} (done)', flush=True)
            continue
        t0 = time.time()
        prm, e, hist = optimize(N, chi, seed=seed, steps=steps, n_sample=nsamp)
        dt = time.time() - t0
        np.savez(pfile, **{k: np.asarray(v) for k, v in prm.items()})
        results[key] = dict(N=N, chi=chi, seed=seed, E_over_N=e, exact=ex,
                            rel_err=abs(e - ex) / abs(ex), seconds=dt)
        print(f'== {key}: E/N={e:.6f} exact={ex:.6f} rel={abs(e-ex)/abs(ex):.2e} ({dt:.0f}s)', flush=True)
        with open('/home/z/my-project/scripts/mera_grid_results.json', 'w') as f:
            json.dump(results, f, indent=1)
print('GRID DONE')
