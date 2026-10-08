"""test_mera.py — validate exact cone contraction against brute force."""
import sys, time
sys.path.insert(0, '/home/z/my-project/scripts')
import numpy as np
from mera_real import MERA, full_state, brute_rdm, entropy, set_backend

def check(N, chi, seed, As):
    m = MERA(N, chi, seed=seed)
    # isometry / unitarity checks
    for l in range(m.L):
        U = np.asarray(m.u[l])
        err_u = np.abs(U @ U.conj().T - np.eye(U.shape[0])).max()
        W = np.asarray(m.w[l])
        err_w = np.abs(W.conj().T @ W - np.eye(W.shape[1])).max()
        assert err_u < 1e-12 and err_w < 1e-12, (l, err_u, err_w)
    psi = full_state(m)
    nrm = np.linalg.norm(psi)
    maxerr = 0.0
    for A in As:
        r1 = m.rdm(A)
        r2 = brute_rdm(psi, N, A)
        e = np.abs(r1 - r2).max()
        maxerr = max(maxerr, e)
        tr = np.trace(r1).real
        herm = np.abs(r1 - r1.conj().T).max()
        assert abs(tr - 1) < 1e-10, (A, tr)
        assert herm < 1e-10, (A, herm)
    print(f"N={N} chi={chi} seed={seed}: ||psi||={nrm:.6f}  max|rdm_cone - rdm_brute| = {maxerr:.2e}")
    return maxerr

t0 = time.time()
As8 = [[0], [1], [0, 1], [1, 2], [0, 3], [2, 5], [0, 1, 2], [1, 2, 3],
       [0, 1, 2, 3], [0, 5], [3, 7], [0, 7], [2, 4], [0, 2, 4, 6], [5, 6, 7]]
e1 = check(8, 2, 1, As8)
e2 = check(8, 2, 7, As8)
As16 = [[0], [1], [0, 1], [1, 2], [0, 5], [3, 11], [0, 1, 2], [0, 1, 2, 3],
        [0, 1, 2, 3, 4], [7, 8], [2, 13], [0, 15], [4, 5, 6, 7, 8], [0, 3], [1, 14]]
e3 = check(16, 2, 3, As16)
e4 = check(16, 4, 5, As16[:10])
assert max(e1, e2, e3, e4) < 1e-9, 'CONE CONTRACTION MISMATCH'
print(f'ALL PASSED in {time.time()-t0:.1f}s — cone contraction is exact')
