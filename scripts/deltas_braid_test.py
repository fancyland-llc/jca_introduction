"""
deltas_braid_sweep.py
Sweep the deltas-only interferometric braid up the primorial tower using
matrix-free FFT-accelerated operators (Holographic Eigensolver pattern).

For each m:
  1. Build C = [D_sym, P_tau] as a scipy LinearOperator. D_sym matvec is
     O(m log m) via rfft on the circulant tent kernel restricted to (Z/m)^x.
     P_tau matvec is a pure permutation.
  2. Top-8 singular values + right singular vectors via ARPACK svds.
  3. Spectral signatures: paired-doubled, overtone ratio, fundamental amp.
  4. V4 invariance under H = iC measured matrix-free.
  5. 4x4 deltas-only braid evolution.
  6. Report wall time per stage.

Skip dense full-register reconstruction: validated at m=2310 in the prior
test (matched to 1.12e-11). At larger m, the Hilbert space cannot be
reconstructed densely on any reasonable machine; the matrix-free path is
the only path.
"""
import numpy as np
from math import gcd as _gcd
from scipy.sparse.linalg import LinearOperator, svds
from scipy.linalg import expm
import time
import sys


PRIMORIAL_LADDER = [
    (2310,     "5th primorial 2*3*5*7*11"),
    (30030,    "6th primorial 2*3*5*7*11*13"),
    (510510,   "7th primorial 2*3*5*7*11*13*17"),
    (9699690,  "8th primorial 2*3*5*7*11*13*17*19"),
]


def factor_squarefree(m):
    """Return the prime tuple of a squarefree positive integer m."""
    primes = []
    d, n = 2, int(m)
    while d * d <= n:
        if n % d == 0:
            primes.append(d)
            n //= d
        else:
            d += 1
    if n > 1:
        primes.append(n)
    return primes


def _pow_vec_mod(base_arr, exp, mod):
    """Vectorized base_arr^exp mod mod via binary exponentiation.

    base_arr is a numpy int64 array, exp/mod are Python ints. All ops stay
    inside int64 because mod * mod fits comfortably for the small primes
    p <= 19 of the primorial ladder.
    """
    result = np.ones_like(base_arr)
    cur = base_arr % mod
    e = exp
    while e > 0:
        if e & 1:
            result = (result * cur) % mod
        cur = (cur * cur) % mod
        e >>= 1
    return result


def _inv_R_via_crt(R, m, primes):
    """Vectorized r^-1 mod m for all r in R via CRT over the prime tuple.

    For squarefree m = prod_i p_i, Garner-style CRT gives
        r^-1 mod m  ==  sum_i [(r mod p_i)^(p_i - 2) mod p_i] * M_i  (mod m)
    where M_i = (m/p_i) * ((m/p_i)^-1 mod p_i) mod m.

    Replaces the per-element pow(r, -1, m) Python loop. At the 8th primorial
    (n = 1,658,880) this drops the inv_perm build from ~10 s of Python-loop
    to a fraction of a second of vectorized numpy.
    """
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
    """Vectorized build of (R, inv_perm, rfft_kernel) for (Z/m)^x."""
    arr = np.arange(m, dtype=np.int64)
    R = arr[np.gcd(arr, m) == 1]
    n = len(R)

    R_index = np.full(m, -1, dtype=np.int64)
    R_index[R] = np.arange(n)

    # P_tau permutation: r -> r^-1 mod m, mapped through R indexing.
    # Vectorized CRT replaces the Python-loop pow(r, -1, m) — required for
    # the 8th primorial, where n = 1,658,880 makes the per-element loop
    # dominate the substrate-build wall time.
    primes = factor_squarefree(m)
    inv_R = _inv_R_via_crt(R, m, primes)
    inv_perm = R_index[inv_R]

    # Circulant tent kernel f(k) = min(k, m-k); store its rfft once.
    f_kernel = np.minimum(arr, m - arr).astype(np.float64)
    rfft_kernel = np.fft.rfft(f_kernel)
    return R, inv_perm, rfft_kernel, n


def make_C_operator(R, inv_perm, rfft_kernel, n, m):
    """C = [D_sym, P_tau], real skew, FFT-accelerated matvec."""

    def D_sym_matvec(v):
        v_full = np.zeros(m, dtype=np.float64)
        v_full[R] = v
        u_full = np.fft.irfft(np.fft.rfft(v_full) * rfft_kernel, n=m)
        return u_full[R]

    def C_matvec(v):
        # C v = D_sym (P_tau v) - P_tau (D_sym v)
        return D_sym_matvec(v[inv_perm]) - D_sym_matvec(v)[inv_perm]

    def C_rmatvec(v):
        return -C_matvec(v)  # real skew: C^T = -C

    def C_matmat(M):
        # Apply column-by-column; svds may pass blocks.
        out = np.empty_like(M)
        for j in range(M.shape[1]):
            out[:, j] = C_matvec(M[:, j])
        return out

    return LinearOperator(
        shape=(n, n),
        matvec=C_matvec,
        rmatvec=C_rmatvec,
        matmat=C_matmat,
        dtype=np.float64,
    )


def run_sweep(m, label, do_dense_check=False):
    print(f"\n{'='*74}")
    print(f"  m = {m:,}  ({label})")
    print(f"{'='*74}")

    t0 = time.time()
    R, inv_perm, rfft_kernel, n = build_substrate(m)
    t_build = time.time() - t0
    print(f"  n = phi(m) = {n:,}")
    print(f"  substrate build:       {t_build:7.3f}s")

    C_op = make_C_operator(R, inv_perm, rfft_kernel, n, m)

    # Optional dense cross-check at small m.
    if do_dense_check and n <= 1000:
        t0 = time.time()
        I = np.eye(n)
        C_dense = np.array([C_op.matvec(I[:, k]) for k in range(n)]).T
        skew = np.linalg.norm(C_dense + C_dense.T) / np.linalg.norm(C_dense)
        s_dense = np.linalg.svd(C_dense, compute_uv=False)
        t_dense = time.time() - t0
        print(f"  [dense xcheck]:        {t_dense:7.3f}s")
        print(f"    skew check:                {skew:.2e}")
        print(f"    sigma_0 (dense): {s_dense[0]:.6e}")

    # ---- top-8 SVD via ARPACK Lanczos ----
    t0 = time.time()
    U, s_arpack, Vt = svds(C_op, k=8, which='LM', maxiter=2000, tol=1e-10)
    # ARPACK returns ascending; flip to descending.
    idx = np.argsort(-s_arpack)
    s = s_arpack[idx]
    Vt = Vt[idx, :]
    t_svd = time.time() - t0
    print(f"  matrix-free top-8 SVD: {t_svd:7.3f}s")

    # ---- spectral signatures ----
    print(f"\n  Spectral signatures:")
    print(f"    sigma_0..7 = {[f'{x:.3e}' for x in s]}")
    print(f"    pair (0,1)        ratio: {s[1]/s[0]:.10f}  (target 1.0)")
    print(f"    pair (2,3)        ratio: {s[3]/s[2]:.10f}  (target 1.0)")
    print(f"    pair (4,5)        ratio: {s[5]/s[4]:.10f}  (target 1.0)")
    print(f"    pair (6,7)        ratio: {s[7]/s[6]:.10f}  (target 1.0)")

    # Triangle-wave overtone: sigma_4/sigma_0 -> 1/9
    r_4_0 = s[4] / s[0]
    print(f"    sigma_4/sigma_0:         {r_4_0:.6f}  "
          f"(target 1/9 = {1/9:.6f}, "
          f"err {abs(r_4_0 - 1/9):.2e})")

    # Fundamental amplitude: sigma_0 ~ m * phi(m) / pi^2
    fund_target = m * n / np.pi**2
    print(f"    sigma_0 / (m phi/pi^2):  {s[0] / fund_target:.6f}  "
          f"(target 1.0)")

    # ---- V4 invariance under H = iC, matrix-free ----
    V4 = Vt[:4, :].T  # (n, 4), real
    t0 = time.time()
    CV4 = np.column_stack([C_op.matvec(V4[:, k]) for k in range(4)])
    # H_iC_4D = V4^T (iC) V4 = i (V4^T C V4); 4x4 complex Hermitian
    C_4D = V4.T @ CV4  # 4x4 real (skew, since V4 invariant under skew C)
    H_iC_4D = 1j * C_4D
    # Leakage: ||H V4 - V4 H_4D||_F where H = iC, in real terms ||CV4 - V4 (V4^T CV4)||
    leak = np.linalg.norm(CV4 - V4 @ C_4D)
    leak_rel = leak / np.linalg.norm(CV4)
    t_inv = time.time() - t0
    print(f"  V4 invariance check:   {t_inv:7.3f}s")
    print(f"    ||C V4 - V4 (V4^T C V4)||_F: {leak:.4e}")
    print(f"    relative leakage:            {leak_rel:.4e}  "
          f"(carrier sharpens with m)")

    # ---- 4x4 deltas-only braid ----
    t0 = time.time()
    rng = np.random.default_rng(42)
    q_4D = rng.standard_normal(4) + 1j * rng.standard_normal(4)
    q_4D /= np.linalg.norm(q_4D)
    Pi_0 = np.diag([1.0, 1.0, 0.0, 0.0])
    Pi_1 = np.diag([0.0, 0.0, 1.0, 1.0])
    psi_0_4D = Pi_0 @ q_4D
    psi_1_4D = Pi_1 @ q_4D
    # Normalize evolution time so phase advance is O(1) regardless of scale.
    t_evolve = 1.0 / s[0]
    U_4D = expm(-1j * t_evolve * H_iC_4D)
    delta_4D = U_4D @ psi_0_4D - U_4D @ psi_1_4D
    t_braid = time.time() - t0
    print(f"  4x4 braid (Delta(t)):  {t_braid*1e6:7.1f} us")
    print(f"    |delta_4D| = {np.linalg.norm(delta_4D):.6f}")

    total = t_build + t_svd + t_inv + t_braid
    print(f"  TOTAL:                 {total:7.3f}s")

    return {
        "m": m, "n": n,
        "t_build": t_build, "t_svd": t_svd, "t_inv": t_inv, "t_braid": t_braid,
        "sigma": s.tolist(), "leak_rel": leak_rel,
        "ratio_4_0": r_4_0, "fund_ratio": s[0] / fund_target,
    }


def main():
    results = []
    # Run with dense cross-check at m=2310 to confirm matrix-free correctness.
    first = True
    for m, label in PRIMORIAL_LADDER:
        try:
            r = run_sweep(m, label, do_dense_check=first)
            results.append(r)
            first = False
        except MemoryError as e:
            print(f"\n!!! MemoryError at m={m}: {e}")
            break
        except Exception as e:
            print(f"\n!!! FAILED at m={m}: {type(e).__name__}: {e}")
            import traceback; traceback.print_exc()
            break

    print(f"\n{'='*74}")
    print(f"  Summary")
    print(f"{'='*74}")
    print(f"  {'m':>10} {'n':>10} {'build':>8} {'SVD':>9} "
          f"{'leak_rel':>11} {'sig_4/sig_0':>12}")
    for r in results:
        print(f"  {r['m']:>10,} {r['n']:>10,} "
              f"{r['t_build']:>7.2f}s {r['t_svd']:>8.2f}s "
              f"{r['leak_rel']:>11.3e} {r['ratio_4_0']:>12.6f}")
    print(f"\n  triangle-wave target sigma_4/sigma_0 = {1/9:.6f}")
    print(f"  carrier ratio (1 - leak_rel) sharpens as ~ phi(m)^-1 per Theorem 6.1(iv).")


if __name__ == "__main__":
    main()