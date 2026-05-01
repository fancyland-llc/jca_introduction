"""
sigma_tau_braid.py — Does V4 = (sigma-parity) x (tau-parity) qubit?

Two involutions on (Z/m)^x that both commute with each other:
  P_sigma:  r -> -r mod m   (additive reflection)
  P_tau:    r -> r^-1 mod m (multiplicative inversion)

C = [D_sym, P_tau] satisfies:
  [C, P_sigma] = 0     => V4 invariant under P_sigma
  {C, P_tau}   = 0     => V4 invariant under P_tau (since C^T C commutes)

[P_sigma, P_tau] = 0 algebraically. Therefore P_sigma|V4 and P_tau|V4 are
two commuting Hermitian involutions on a 4D space, simultaneously
diagonalizable, splitting V4 into four 1D sectors labeled by
(sigma-parity, tau-parity) in {+/-1}^2.

Test: does V4 actually decompose this way? If yes, the rank-4 carrier
IS a 2-bit arithmetic qubit with sigma and tau as the two register bits.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator, svds
import time


def factor_squarefree(m):
    primes = []
    d, n = 2, int(m)
    while d * d <= n:
        if n % d == 0:
            primes.append(d); n //= d
        else:
            d += 1
    if n > 1: primes.append(n)
    return primes


def _pow_vec_mod(base, exp, mod):
    result = np.ones_like(base)
    cur = base % mod
    e = exp
    while e > 0:
        if e & 1: result = (result * cur) % mod
        cur = (cur * cur) % mod
        e >>= 1
    return result


def _inv_R_crt(R, m, primes):
    inv = np.zeros_like(R)
    for p in primes:
        m_over_p = m // p
        m_over_p_inv = pow(int(m_over_p % p), -1, int(p))
        Mi = (m_over_p * m_over_p_inv) % m
        r_mod_p = (R % p).astype(np.int64)
        inv_mod_p = _pow_vec_mod(r_mod_p, p - 2, p)
        inv = (inv + inv_mod_p * Mi) % m
    return inv


def build_substrate(m):
    arr = np.arange(m, dtype=np.int64)
    R = arr[np.gcd(arr, m) == 1]
    n = len(R)
    R_index = np.full(m, -1, dtype=np.int64)
    R_index[R] = np.arange(n)

    primes = factor_squarefree(m)
    inv_R = _inv_R_crt(R, m, primes)
    tau_perm = R_index[inv_R]                      # P_tau: r -> r^-1
    sigma_perm = R_index[(m - R) % m]              # P_sigma: r -> -r

    f_kernel = np.minimum(arr, m - arr).astype(np.float64)
    rfft_kernel = np.fft.rfft(f_kernel)
    return R, tau_perm, sigma_perm, rfft_kernel, n


def make_C(R, tau_perm, rfft_kernel, n, m):
    def D_sym(v):
        v_full = np.zeros(m, dtype=np.float64)
        v_full[R] = v
        return np.fft.irfft(np.fft.rfft(v_full) * rfft_kernel, n=m)[R]
    def C(v):
        v = np.asarray(v).ravel()
        return D_sym(v[tau_perm]) - D_sym(v)[tau_perm]
    return C


def run(m, label):
    print(f"\n{'='*72}")
    print(f"  m = {m:,}  ({label})")
    print(f"{'='*72}")
    t_total = time.time()

    R, tau_perm, sigma_perm, rfft_kernel, n = build_substrate(m)
    print(f"  n = phi(m) = {n:,}")

    # Sanity: P_sigma and P_tau commute as permutations
    sigma_then_tau = tau_perm[sigma_perm]
    tau_then_sigma = sigma_perm[tau_perm]
    print(f"  [P_sigma, P_tau] = 0 check (perm equality): "
          f"{np.array_equal(sigma_then_tau, tau_then_sigma)}")

    # Top-4 SVD
    Cmv = make_C(R, tau_perm, rfft_kernel, n, m)
    C_op = LinearOperator(shape=(n, n), matvec=Cmv,
                          rmatvec=lambda v: -Cmv(v), dtype=np.float64)
    t0 = time.time()
    U, s, Vt = svds(C_op, k=4, which='LM', tol=1e-10)
    idx = np.argsort(-s); s = s[idx]; Vt = Vt[idx, :]
    V4 = Vt.T
    print(f"  top-4 SVD: {time.time()-t0:.2f}s,  "
          f"sigma_0..3 = {[f'{x:.3e}' for x in s]}")

    # Project P_sigma, P_tau onto V4
    Psigma_4D = V4.T @ V4[sigma_perm, :]
    Ptau_4D   = V4.T @ V4[tau_perm, :]

    print(f"\n  --- (sigma, tau) braid on V4 ---")
    print(f"    Psigma_4D - Psigma_4D^T : {np.linalg.norm(Psigma_4D - Psigma_4D.T):.2e}  (should be 0; symmetric)")
    print(f"    Ptau_4D   - Ptau_4D^T   : {np.linalg.norm(Ptau_4D - Ptau_4D.T):.2e}")
    print(f"    Psigma_4D^2 - I         : {np.linalg.norm(Psigma_4D @ Psigma_4D - np.eye(4)):.2e}  (should be 0; involution)")
    print(f"    Ptau_4D^2 - I           : {np.linalg.norm(Ptau_4D @ Ptau_4D - np.eye(4)):.2e}")
    print(f"    [Psigma_4D, Ptau_4D]    : {np.linalg.norm(Psigma_4D @ Ptau_4D - Ptau_4D @ Psigma_4D):.2e}  (should be 0; commute)")

    es = np.linalg.eigvalsh(Psigma_4D)
    et = np.linalg.eigvalsh(Ptau_4D)
    print(f"    Psigma_4D eigenvalues   : {[f'{x:+.6f}' for x in sorted(es)]}  (should be +/-1)")
    print(f"    Ptau_4D   eigenvalues   : {[f'{x:+.6f}' for x in sorted(et)]}")

    # Simultaneous diagonalization
    es_vals, es_vecs = np.linalg.eigh(Psigma_4D)
    Ptau_in_sigma = es_vecs.T @ Ptau_4D @ es_vecs
    print(f"\n    P_tau in sigma-eigenbasis (block-diagonal expected):")
    for row in np.round(Ptau_in_sigma, 6):
        print(f"      [{', '.join(f'{x:+.4f}' for x in row)}]")

    # Build qubit basis: simultaneously diagonalize
    qubit_basis = np.zeros((4, 4))
    sigma_assign, tau_assign = [], []
    eps = 1e-6
    col = 0
    # Group sigma-eigenvectors by eigenvalue, then diagonalize tau within each block
    for sigma_val_target in [-1.0, +1.0]:
        mask = np.abs(es_vals - sigma_val_target) < eps
        if not mask.any(): continue
        block = es_vecs[:, mask]
        Ptau_block = block.T @ Ptau_4D @ block
        et_v, et_w = np.linalg.eigh(Ptau_block)
        for k in range(et_v.shape[0]):
            qubit_basis[:, col] = block @ et_w[:, k]
            sigma_assign.append(int(np.round(sigma_val_target)))
            tau_assign.append(int(np.round(et_v[k])))
            col += 1

    print(f"\n    Qubit basis (sigma, tau) labels:")
    for k in range(4):
        print(f"      |q{k}>:  sigma={sigma_assign[k]:+d}, tau={tau_assign[k]:+d}")

    # Verify diagonalization
    sigma_diag = np.diag(sigma_assign).astype(np.float64)
    tau_diag = np.diag(tau_assign).astype(np.float64)
    err_s = np.linalg.norm(Psigma_4D @ qubit_basis - qubit_basis @ sigma_diag)
    err_t = np.linalg.norm(Ptau_4D @ qubit_basis - qubit_basis @ tau_diag)
    print(f"    ||Psigma U - U diag(s)|| = {err_s:.2e}")
    print(f"    ||Ptau U - U diag(t)||   = {err_t:.2e}")

    # C in the qubit basis
    CV4 = np.column_stack([Cmv(V4[:, k]) for k in range(4)])
    C_4D = V4.T @ CV4
    C_qubit = qubit_basis.T @ C_4D @ qubit_basis
    print(f"\n    C in (sigma, tau) qubit basis (rows/cols ordered as listed above):")
    for row in np.round(C_qubit, 4):
        print(f"      [{', '.join(f'{x:+10.4f}' for x in row)}]")

    # Predicted block structure: C anti-commutes with P_tau, commutes with P_sigma
    # => Block-diagonal in sigma-sectors, off-diagonal in tau within each sector
    print(f"\n    Predicted: block-diag in sigma-sectors, off-diag in tau within sector.")
    sigma_signs = np.array(sigma_assign)
    tau_signs = np.array(tau_assign)
    # Off-sigma elements should be zero
    off_sigma_mask = sigma_signs[:, None] != sigma_signs[None, :]
    same_tau_mask = tau_signs[:, None] == tau_signs[None, :]
    off_sigma_norm = np.linalg.norm(C_qubit[off_sigma_mask])
    diag_tau_norm = np.linalg.norm(C_qubit[same_tau_mask & ~np.eye(4, dtype=bool).astype(bool)])
    diag_norm = np.linalg.norm(np.diag(C_qubit))
    print(f"      ||C entries crossing sigma-sectors||: {off_sigma_norm:.4e}  (should be ~0)")
    print(f"      ||C entries with same tau (off-diag)||: {diag_tau_norm:.4e}  (should be ~0)")
    print(f"      ||C diagonal||:                       {diag_norm:.4e}  (should be ~0; C is skew)")

    print(f"\n  TOTAL: {time.time()-t_total:.2f}s")


PRIMORIAL_LADDER = [
    (30,      "p = 2*3*5  (sanity)"),
    (210,     "p = 2*3*5*7"),
    (2310,    "5th primorial"),
    (30030,   "6th primorial"),
    (510510,  "7th primorial"),
]

for m, label in PRIMORIAL_LADDER:
    try:
        run(m, label)
    except Exception as e:
        print(f"\n!!! FAILED at m={m}: {type(e).__name__}: {e}")
        import traceback; traceback.print_exc()
        break