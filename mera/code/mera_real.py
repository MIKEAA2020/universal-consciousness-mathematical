"""
mera_real.py — Binary MERA with REAL tensor contraction (repairs fabricated results).

Every number produced here is the output of an exact contraction of real tensors.

Layout (brick-layered binary MERA on a periodic chain, N = 2^L sites):
  levels 0..L; level l has n_l = N>>l sites. Site dims: level 0 -> 2 (spins);
  level 1 -> d1 = min(chi,4); levels >= 2 -> chi.
  Layer l (isometry V_l: level l+1 -> level l):
    u_l[j], j = 0..n_l/2-1 : unitary on raw site pair (2j-1 mod n_l, 2j mod n_l)
      (stored as a d^2 x d^2 matrix, d = site dim at level l)
    w_l[k], k = 0..n_{l+1}-1 : isometry consuming the u-processed sites
      (2k, 2k+1); w_l[k] eats u_l[k].legB (u-proc site 2k) and
      u_l[k+1].legA (u-proc site 2k+1). Stored as a (d*d) x d_{l+1} matrix
      with orthonormal columns (w^dag w = I).
  Top: |t> on the single level-L site. |Psi> = V_0 V_1 ... V_{L-1} |t>.

RDMs by exact causal-cone contraction (descending):
  cones: C^0 = A (target); C^{l+1} = {k: {2k-1,2k,2k+1,2k+2} mod n_l meets C^l}.
  A^L = |t><t| ; A^l = (cone-net_l) A^{l+1} (cone-net_l)^dag, boundary traced.
  Closure facts used (provable from the layout): every u with pair touching C^l
  has both its consumers w[j-1], w[j] in the cone; every k in C^{l+1} has w[k]
  in the cone; fully-traced u (w) reduces to delta on its outputs by unitarity
  (isometry). Validated against brute-force full-state contraction.
"""
import math
import numpy as np

_BACKEND = 'np'

def set_backend(b):
    global _BACKEND
    assert b in ('np', 'jax')
    _BACKEND = b

def _einsum(subs, *arrs):
    import opt_einsum
    if _BACKEND == 'jax':
        return opt_einsum.contract(subs, *arrs, optimize='optimal', backend='jax')
    return opt_einsum.contract(subs, *arrs, optimize='optimal')

_LETTERS = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
_MEM_CAP = 2 ** 25   # ~536MB per intermediate tensor (machine has 4GB)

class MemoryCapError(RuntimeError):
    pass

def named_einsum(ops, out_names):
    """ops: list of (array, [legname,...]). A leg appearing exactly twice is
    contracted; output keeps out_names in order.

    Scheduled by an exact size-greedy pairwise rule: at each step contract
    the pair of operands whose result has the fewest entries. This avoids the
    pathological contraction paths that opt_einsum greedy/optimal produce on
    these multi-leg cone networks."""
    dims = {}
    tensors = []
    for arr, lg in ops:
        shp = list(arr.shape)
        assert len(shp) == len(lg), (shp, lg)
        for name, sz in zip(lg, shp):
            if name in dims:
                assert dims[name] == sz, (name, dims[name], sz)
            else:
                dims[name] = sz
        tensors.append([arr, list(lg)])
    want = list(out_names)
    while len(tensors) > 1:
        best = None
        for i in range(len(tensors)):
            for j in range(i + 1, len(tensors)):
                li, lj = tensors[i][1], tensors[j][1]
                shared = set(li) & set(lj)
                if not shared:
                    continue
                outl = [x for x in li if x not in shared] + \
                       [x for x in lj if x not in shared]
                size = 1
                for x in outl:
                    size *= dims[x]
                if best is None or size < best[0]:
                    best = (size, i, j, outl)
        if best is None:
            raise RuntimeError('disconnected contraction graph')
        if best[0] > _MEM_CAP:
            raise MemoryCapError(f'contraction requires {best[0]} entries > cap {_MEM_CAP}')
        _, i, j, outl = best
        a, la = tensors[i]
        b, lb = tensors[j]
        r = _einsum2(a, la, b, lb, outl)
        for k in sorted((i, j), reverse=True):
            tensors.pop(k)
        tensors.append([r, outl])
    arr, lg = tensors[0]
    assert set(lg) == set(want), (set(lg), set(want))
    if list(lg) != want:
        perm = [lg.index(x) for x in want]
        arr = arr.transpose(perm) if _BACKEND == 'jax' else np.transpose(arr, perm)
    return arr

def _einsum2(a, la, b, lb, outl):
    letters = {}
    def let(n):
        if n not in letters:
            letters[n] = _LETTERS[len(letters)]
        return letters[n]
    sa = ''.join(let(x) for x in la)
    sb = ''.join(let(x) for x in lb)
    so = ''.join(let(x) for x in outl)
    return _einsum(sa + ',' + sb + '->' + so, a, b)

# ------------------------------------------------------------------ random tensors

def rand_unitary(m, rng):
    z = rng.standard_normal((m, m)) + 1j * rng.standard_normal((m, m))
    q, _ = np.linalg.qr(z)
    return q

def rand_isometry(big, small, rng):
    z = rng.standard_normal((big, small)) + 1j * rng.standard_normal((big, small))
    q, _ = np.linalg.qr(z)
    return q[:, :small]

# ------------------------------------------------------------------ MERA

class MERA:
    def __init__(self, N, chi, tensors=None, seed=0):
        assert N >= 8 and (N & (N - 1)) == 0 and chi >= 2
        self.N, self.chi = N, chi
        self.L = int(math.log2(N))
        self.n = [N >> l for l in range(self.L + 1)]
        self.dims = [2] + [min(chi, 4)] + [chi] * (self.L - 1)
        if tensors is not None:
            self.u, self.w, self.t = tensors
            return
        rng = np.random.default_rng(seed)
        self.u, self.w = [], []
        for l in range(self.L):
            d = self.dims[l]
            self.u.append(rand_unitary(d * d, rng))
            self.w.append(rand_isometry(d * d, self.dims[l + 1], rng))
        t = rng.standard_normal(self.dims[self.L]) + 1j * rng.standard_normal(self.dims[self.L])
        self.t = t / np.linalg.norm(t)

    def upair(self, l, j):
        n = self.n[l]
        return ((2 * j - 1) % n, (2 * j) % n)

    def cones_for(self, A):
        C = [set(int(x) for x in A)]
        for l in range(self.L):
            n, cup = self.n[l], C[-1]
            nxt = set()
            for k in range(self.n[l + 1]):
                if ((2 * k - 1) % n in cup or (2 * k) % n in cup
                        or (2 * k + 1) % n in cup or (2 * k + 2) % n in cup):
                    nxt.add(k)
            if not nxt:
                nxt = {0}
            C.append(nxt)
        C[self.L] = {0}
        return C

    # ---------------- exact cone RDM ----------------
    def rdm(self, A):
        """Exact rho_A; legs (ket over sorted A) x (bra over sorted A), dim 2 each."""
        A = sorted(set(int(x) for x in A))
        C = self.cones_for(A)
        tp = self.t
        env = (tp[:, None] * tp.conj()[None, :], ['Tk'], ['Tb'])  # |t><t|
        for l in range(self.L - 1, -1, -1):
            env = self._descend(l, C[l + 1], C[l], env)
        tensor = env[0]
        sq = 2 ** len(A)
        return np.asarray(tensor).reshape(sq, sq) if _BACKEND == 'np' else tensor.reshape(sq, sq)

    def _descend(self, l, C_up, C_low, env):
        n, n2 = self.n[l], self.n[l + 1]
        d = self.dims[l]
        cup_sorted = sorted(C_up)
        ops = []
        # environment tensor over C_up (legs renamed e/E)
        env_legs = ['e%d' % k for k in cup_sorted] + ['E%d' % k for k in cup_sorted]
        ops.append((env[0], env_legs))
        # involved u's: j in {k, k+1 : k in C_up}; include only those whose
        # in-pair meets C_low (fully traced u's reduce to delta on outputs)
        J = set()
        for k in C_up:
            J.add(k % n2)
            J.add((k + 1) % n2)
        u_included = {}
        for j in sorted(J):
            s1, s2 = self.upair(l, j)
            inc = (s1 in C_low) or (s2 in C_low)
            u_included[j] = inc
            if not inc:
                continue
            U = self.u[l].reshape(d, d, d, d)      # shared unitary, 4 legs
            i1k = 'a%d' % s1 if s1 in C_low else 'x%d' % s1
            i1b = 'A%d' % s1 if s1 in C_low else 'x%d' % s1
            i2k = 'a%d' % s2 if s2 in C_low else 'x%d' % s2
            i2b = 'A%d' % s2 if s2 in C_low else 'x%d' % s2
            o1k, o1b = 'p%d' % s1, 'P%d' % s1     # u-proc site s1=2j-1 -> w[j-1]
            o2k, o2b = 'q%d' % s2, 'Q%d' % s2     # u-proc site s2=2j   -> w[j]
            ops.append((U, [i1k, i2k, o1k, o2k]))
            ops.append((U.conj(), [i1b, i2b, o1b, o2b]))
        # w's for k in C_up
        for k in sorted(C_up):
            W = self.w[l]      # shared isometry
            Wt = W.reshape(d, d, self.dims[l + 1])
            # input site 2k produced by u[k] (legB); site 2k+1 by u[k+1] (legA)
            s_even, s_odd = (2 * k) % n, (2 * k + 1) % n
            inc1 = u_included.get(k % n2, False)          # u[k] included?
            inc2 = u_included.get((k + 1) % n2, False)    # u[k+1] included?
            if inc1:
                i1k, i1b = 'q%d' % s_even, 'Q%d' % s_even
            else:
                i1k = i1b = 'zq%d' % s_even
            if inc2:
                i2k, i2b = 'p%d' % s_odd, 'P%d' % s_odd
            else:
                i2k = i2b = 'zp%d' % s_odd
            ops.append((Wt, [i1k, i2k, 'e%d' % k]))
            ops.append((Wt.conj(), [i1b, i2b, 'E%d' % k]))
        low_sorted = sorted(C_low)
        outk = ['a%d' % s for s in low_sorted]
        outb = ['A%d' % s for s in low_sorted]
        T = named_einsum(ops, outk + outb)
        return (T, outk, outb)

# ------------------------------------------------------------------ brute force

def full_state(mera):
    """Dense state vector (small N only), built top-down. Returns numpy."""
    psi = np.asarray(mera.t).reshape([mera.dims[mera.L]])
    for l in range(mera.L - 1, -1, -1):
        n, n2 = mera.n[l], mera.n[l + 1]
        d = mera.dims[l]
        # step 1: w-expansion: top legs t_k -> u-proc legs c_{2k}, c_{2k+1}
        ops = [(psi, ['t%d' % k for k in range(n2)])]
        outlegs = []
        for k in range(n2):
            Wt = np.asarray(mera.w[l]).reshape(d, d, mera.dims[l + 1])
            ops.append((Wt, ['c%d' % (2 * k), 'c%d' % (2 * k + 1), 't%d' % k]))
            outlegs += ['c%d' % (2 * k), 'c%d' % (2 * k + 1)]
        T = named_einsum(ops, outlegs)
        # step 2: u's on staggered raw pairs (2j-1, 2j): c-legs -> r-legs
        ops = [(T, ['c%d' % s for s in range(n)])]
        for j in range(n // 2):
            s1, s2 = mera.upair(l, j)
            U = np.asarray(mera.u[l]).reshape(d, d, d, d)
            ops.append((U, ['r%d' % s1, 'r%d' % s2, 'c%d' % s1, 'c%d' % s2]))
        psi = named_einsum(ops, ['r%d' % s for s in range(n)])
    return psi.reshape(-1)

def brute_rdm(psi, N, A):
    """rho_A from the dense state vector by partial trace. Legs (ket, bra)."""
    A = sorted(set(int(x) for x in A))
    psi = np.asarray(psi).reshape([2] * N)
    ket_legs = ['k%d' % s if s in A else 't%d' % s for s in range(N)]
    bra_legs = ['b%d' % s if s in A else 't%d' % s for s in range(N)]
    out = ['k%d' % s for s in A] + ['b%d' % s for s in A]
    rho = named_einsum([(psi, ket_legs), (psi.conj(), bra_legs)], out)
    sq = 2 ** len(A)
    return rho.reshape(sq, sq)

# ------------------------------------------------------------------ observables

def entropy(rho, tol=1e-14):
    v = np.linalg.eigvalsh(rho)
    v = v[v > tol]
    return float(-np.sum(v * np.log(v)))

def mutual_info(rho_i, rho_j, rho_ij):
    return entropy(rho_i) + entropy(rho_j) - entropy(rho_ij)

def vne(rho):
    return entropy(rho)
