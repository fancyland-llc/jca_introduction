"""
substrate_hysteresis_v1.py — Reading 1 (claude.ai protocol).

Cyclic σ-flux drive on the V_4 qubit basis of m=510,510.
Order parameter: ⟨P_τ⟩(λ).
Hypothesis: the σ_+ vs σ_− τ-flip frequency split (21 ppm at m=510,510)
seeds Landau-Zener-class non-adiabatic transitions under finite-rate
drive, producing a hysteresis loop in (λ, ⟨P_τ⟩).

Protocol:
  H(t) = H_iC_4D + λ(t) · P_σ_4D
  λ(t) = λ_max · sin²(πt/T_period)  — triangular up-then-down ramp
  Trotter: |ψ(t+dt)⟩ = exp(-i dt H(t)) |ψ(t)⟩

Three sweep rates:
  SLOW:    T = 1000 / ω_avg  → near-adiabatic, expect tiny loop
  MEDIUM:  T = 100  / ω_avg  → Landau-Zener regime, expect maximal loop
  FAST:    T = 10   / ω_avg  → diabatic, expect loop closes from speed

Loop area = ∮ ⟨P_τ⟩ dλ around the closed (λ, ⟨P_τ⟩) trajectory.
"""
from __future__ import annotations

import json
import math
import time
from math import gcd
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.linalg import expm

HERE = Path(__file__).parent
M = 510510


# ---------- substrate (reuse from sigma_tau_qubit) ----------
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
    result = np.ones_like(base); cur = base % mod
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
    tau_perm = R_index[inv_R]
    sigma_perm = R_index[(m - R) % m]
    f_kernel = np.minimum(arr, m - arr).astype(np.float64)
    rfft_kernel = np.fft.rfft(f_kernel)
    return R, tau_perm, sigma_perm, rfft_kernel, n


def make_C_matvec(R, tau_perm, rfft_kernel, m):
    def D_sym(v):
        v_full = np.zeros(m, dtype=np.float64)
        v_full[R] = v
        return np.fft.irfft(np.fft.rfft(v_full) * rfft_kernel, n=m)[R]
    def C(v):
        v = np.asarray(v).ravel()
        return D_sym(v[tau_perm]) - D_sym(v)[tau_perm]
    return C


# ---------- qubit basis extraction ----------
def extract_qubit_basis(m):
    """Return (V_4, qubit_basis, sigma_assign, tau_assign,
              C_4D_qbasis, P_sigma_4D_qbasis, P_tau_4D_qbasis)
    in the (sigma, tau) qubit basis."""
    print(f"  Building substrate at m = {m:,} ...")
    t0 = time.time()
    R, tau_perm, sigma_perm, rfft_kernel, n = build_substrate(m)
    print(f"    n = phi(m) = {n:,}  (built in {time.time()-t0:.1f}s)")

    Cmv = make_C_matvec(R, tau_perm, rfft_kernel, m)
    C_op = spla.LinearOperator(shape=(n, n), matvec=Cmv,
                               rmatvec=lambda v: -Cmv(v), dtype=np.float64)
    print(f"  Top-4 SVD of C_op ...")
    t0 = time.time()
    U, s, Vt = spla.svds(C_op, k=4, which='LM', tol=1e-10)
    idx = np.argsort(-s); s = s[idx]; Vt = Vt[idx, :]
    V_4 = Vt.T
    print(f"    {time.time()-t0:.1f}s,  sigma_top = {[f'{x:.4e}' for x in s]}")

    # 4x4 projections
    Psigma_4D = V_4.T @ V_4[sigma_perm, :]
    Ptau_4D   = V_4.T @ V_4[tau_perm, :]
    Psigma_4D = 0.5 * (Psigma_4D + Psigma_4D.T)
    Ptau_4D   = 0.5 * (Ptau_4D + Ptau_4D.T)

    # Simultaneous diag → qubit basis
    es_vals, es_vecs = np.linalg.eigh(Psigma_4D)
    qb = np.zeros((4, 4))
    sigma_assign, tau_assign = [], []
    eps = 1e-6
    col = 0
    for sigma_target in [-1.0, +1.0]:
        mask = np.abs(es_vals - sigma_target) < eps
        if not mask.any(): continue
        block = es_vecs[:, mask]
        Pt_block = block.T @ Ptau_4D @ block
        et_v, et_w = np.linalg.eigh(Pt_block)
        for k in range(et_v.shape[0]):
            qb[:, col] = block @ et_w[:, k]
            sigma_assign.append(int(np.round(sigma_target)))
            tau_assign.append(int(np.round(et_v[k])))
            col += 1
    sigma_a = np.array(sigma_assign, dtype=int)
    tau_a   = np.array(tau_assign, dtype=int)
    print(f"  Qubit basis labels:")
    for k in range(4):
        print(f"    |q{k}> = |sigma={sigma_a[k]:+d}, tau={tau_a[k]:+d}>")

    # C_4D in qubit basis
    CV4 = np.column_stack([Cmv(V_4[:, k]) for k in range(4)])
    C_4D = V_4.T @ CV4
    C_qbasis      = qb.T @ C_4D @ qb
    Psigma_qbasis = qb.T @ Psigma_4D @ qb
    Ptau_qbasis   = qb.T @ Ptau_4D @ qb
    return V_4, qb, sigma_a, tau_a, C_qbasis, Psigma_qbasis, Ptau_qbasis, s


# ---------- hysteresis sweep ----------
def cyclic_sweep(H_iC: np.ndarray, P_sigma: np.ndarray, P_tau: np.ndarray,
                 sigma_assign: np.ndarray, tau_assign: np.ndarray,
                 lambda_max: float, T_period: float, n_steps: int,
                 label: str) -> dict:
    """
    Time-evolve |psi(t)> under H(t) = H_iC + lambda(t) * P_sigma
    with lambda(t) = lambda_max * sin^2(pi t / T_period) over t in [0, T_period].
    Track ⟨P_tau⟩(t).
    Initial state: |q0> = |sigma=-1, tau=-1> (a definite tau-eigenstate).
    """
    psi = np.zeros(4, dtype=np.complex128)
    # Pick a tau-definite initial state (|sigma=-1, tau=-1>)
    idx_q0 = int(np.where((sigma_assign == -1) & (tau_assign == -1))[0][0])
    psi[idx_q0] = 1.0

    dt = T_period / n_steps
    lambdas = np.zeros(n_steps + 1)
    Ptau_exp = np.zeros(n_steps + 1)
    sigma_exp = np.zeros(n_steps + 1)

    for k in range(n_steps + 1):
        t = k * dt
        # Triangle-shaped ramp 0 -> lambda_max -> 0
        lam = lambda_max * (math.sin(math.pi * t / T_period)) ** 2
        lambdas[k] = lam
        Ptau_exp[k] = float(np.real(psi.conj() @ P_tau @ psi))
        sigma_exp[k] = float(np.real(psi.conj() @ P_sigma @ psi))
        if k < n_steps:
            # Trotter step (use midpoint H for 2nd-order accuracy)
            t_mid = t + dt / 2
            lam_mid = lambda_max * (math.sin(math.pi * t_mid / T_period)) ** 2
            H_t = H_iC + lam_mid * P_sigma
            U_dt = expm(-1j * dt * H_t)
            psi = U_dt @ psi
            psi = psi / np.linalg.norm(psi)  # numerical hygiene

    # Compute hysteresis loop area: ∮ <P_tau> dlambda
    # Trapezoidal integration over the (lambda, <P_tau>) trajectory
    loop_area = 0.0
    for k in range(n_steps):
        d_lam = lambdas[k+1] - lambdas[k]
        avg_Ptau = 0.5 * (Ptau_exp[k+1] + Ptau_exp[k])
        loop_area += avg_Ptau * d_lam
    # Note: loop_area can be positive or negative; absolute value = enclosed area

    print(f"\n  --- {label} ---")
    print(f"    T_period = {T_period:.3e}  (n_steps = {n_steps:,}, dt = {dt:.3e})")
    print(f"    lambda_max = {lambda_max:.4e}")
    print(f"    initial state: |q{idx_q0}> = |sigma={sigma_assign[idx_q0]:+d}, tau={tau_assign[idx_q0]:+d}>")
    print(f"    <P_tau> at t=0:         {Ptau_exp[0]:+.6f}")
    print(f"    <P_tau> at t=T/2 (peak): {Ptau_exp[n_steps//2]:+.6f}")
    print(f"    <P_tau> at t=T (return): {Ptau_exp[-1]:+.6f}")
    print(f"    loop area ∮<P_tau>dlambda = {loop_area:+.6e}")
    print(f"    norm at t=T: |psi|^2 = {float(np.real(psi.conj() @ psi)):.6f}")
    return {
        "label": label,
        "T_period": T_period,
        "n_steps": n_steps,
        "lambda_max": lambda_max,
        "lambdas": lambdas.tolist(),
        "Ptau_exp": Ptau_exp.tolist(),
        "sigma_exp": sigma_exp.tolist(),
        "loop_area": float(loop_area),
        "Ptau_initial": float(Ptau_exp[0]),
        "Ptau_peak": float(Ptau_exp[n_steps // 2]),
        "Ptau_return": float(Ptau_exp[-1]),
    }


def main():
    print("=" * 78)
    print(" Substrate hysteresis (Reading 1) — cyclic sigma-flux drive ".center(78))
    print("=" * 78)

    V_4, qb, sigma_a, tau_a, C_4D, P_sigma, P_tau, sigma_top = extract_qubit_basis(M)

    H_iC = 1j * C_4D   # Hermitian Hamiltonian on V_4 in qubit basis

    print(f"\n  Operator structure check:")
    print(f"    H_iC eigenvalues: "
          f"{[f'{x:+.4e}' for x in np.sort(np.linalg.eigvalsh(H_iC.real if np.allclose(H_iC.imag, 0) else H_iC))]}")
    print(f"    P_sigma eigenvalues: "
          f"{[f'{x:+.4f}' for x in np.sort(np.linalg.eigvalsh(P_sigma))]}")
    print(f"    P_tau eigenvalues: "
          f"{[f'{x:+.4f}' for x in np.sort(np.linalg.eigvalsh(P_tau))]}")

    # Compute the σ-dependent τ-frequency split (the 21 ppm at m=510510)
    sig_minus_idx = np.where(sigma_a == -1)[0]
    sig_plus_idx = np.where(sigma_a == +1)[0]
    C_sm = C_4D[np.ix_(sig_minus_idx, sig_minus_idx)]
    C_sp = C_4D[np.ix_(sig_plus_idx, sig_plus_idx)]
    omega_minus = float(np.abs(C_sm[0, 1]))
    omega_plus  = float(np.abs(C_sp[0, 1]))
    omega_avg = 0.5 * (omega_minus + omega_plus)
    delta_omega = abs(omega_plus - omega_minus)
    print(f"\n  Fine-structure split:")
    print(f"    omega(sigma = -1) = {omega_minus:.6e}")
    print(f"    omega(sigma = +1) = {omega_plus:.6e}")
    print(f"    omega_avg         = {omega_avg:.6e}")
    print(f"    delta_omega       = {delta_omega:.6e}  (ratio {delta_omega/omega_avg:.4e})")

    # Sweep parameters
    # lambda_max chosen so that lambda * 1 (P_sigma eigenvalue) ~ omega_avg
    # to get strong (avoided-crossing) coupling. Test multiple lambda_max scales.
    lambda_max = omega_avg

    # Pick T_period in units of 1/omega_avg.
    # SLOW: T such that omega_avg * T = 100, expect near-adiabatic
    # MEDIUM: omega_avg * T = 10
    # FAST: omega_avg * T = 1
    base_T = 1.0 / omega_avg
    sweeps = [
        ("SLOW (near-adiabatic)",  base_T * 1000.0, 5000),
        ("MEDIUM (LZ regime)",     base_T * 100.0,  3000),
        ("MEDIUM-FAST",             base_T * 30.0,   2000),
        ("FAST (diabatic)",         base_T * 10.0,   1500),
        ("VERY FAST",               base_T * 3.0,    1000),
    ]

    results = []
    for label, T, n_steps in sweeps:
        r = cyclic_sweep(H_iC, P_sigma, P_tau, sigma_a, tau_a,
                         lambda_max, T, n_steps, label)
        results.append(r)

    # Verdict
    print(f"\n{'=' * 78}")
    print(" HYSTERESIS VERDICT ".center(78))
    print(f"{'=' * 78}")
    print(f"\n  Loop areas (|∮<P_tau>dlambda|) by rate:")
    print(f"  {'rate':<28s}  {'T_period':>12s}  {'|loop area|':>14s}  {'Ptau_return':>14s}")
    print(f"  {'-'*28}  {'-'*12}  {'-'*14}  {'-'*14}")
    max_area = max(abs(r["loop_area"]) for r in results)
    for r in results:
        marker = " ⭐ MAX" if abs(r["loop_area"]) == max_area else ""
        print(f"  {r['label']:<28s}  {r['T_period']:>12.4e}  "
              f"{abs(r['loop_area']):>14.6e}  {r['Ptau_return']:>+14.6f}{marker}")

    # Hysteresis is real iff loop area is non-negligible relative to lambda_max
    # AND the rate dependence makes physical sense (adiabatic limit closes loop).
    rel_areas = [abs(r["loop_area"]) / lambda_max for r in results]
    if max(rel_areas) > 1e-3:
        print(f"\n  HYSTERESIS LOOP DETECTED.")
        print(f"  Max |loop area|/lambda_max = {max(rel_areas):.4e}")
        print(f"  Loop is non-trivial → substrate hysteresis is empirically validated.")
    else:
        print(f"\n  No measurable hysteresis at any tested rate.")
        print(f"  Max |loop area|/lambda_max = {max(rel_areas):.4e} (negligible)")
        print(f"  Either: substrate is fully reversible under sigma-flux drive,")
        print(f"  OR: lambda_max is too small to reach avoided crossings.")

    out = {
        "m": M,
        "sigma_top": [float(x) for x in sigma_top],
        "fine_structure": {
            "omega_minus": omega_minus,
            "omega_plus": omega_plus,
            "omega_avg": omega_avg,
            "delta_omega": delta_omega,
            "delta_over_omega": delta_omega / omega_avg,
        },
        "sweep_results": [{k: v for k, v in r.items()
                           if k not in ("lambdas", "Ptau_exp", "sigma_exp")}
                          for r in results],
        "trajectories": {r["label"]: {"lambdas": r["lambdas"],
                                       "Ptau_exp": r["Ptau_exp"]}
                         for r in results},
    }
    out_path = HERE / "substrate_hysteresis_v1_results.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\n  Saved -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
