"""
substrate_hysteresis_v3_full_register.py — discriminating test for §3.8.

v1/v2 measured cyclic σ-flux drive on the V_4 4×4 qubit-basis projection
of the rank-4 carrier. They saw open curves at every rate.

Three readings of "doesn't close":
  (1) Quasi-periodic torus dynamics on V_4 (substrate-intrinsic, recurrent
      only on Poincaré timescales, but technically reversible).
  (2) Truncation artifact: V_4 is invariant under iC alone (proved bit-
      precision in §3.5); but iC + λ P_σ does NOT preserve V_4. Drive
      couples V_4 to bulk overtones; the 4×4 projection drops the
      leakage. Open curve in 4×4 might be artifact of the projection,
      not a substrate property.
  (3) Genuine open-system dynamics. Requires dissipation we haven't
      added (no Lindblad term in v1/v2).

This script discriminates (1) vs (2) by evolving the FULL n=92,160-dim
substrate register and projecting to V_4 only at measurement time.

Hamiltonian: H(t) = iC + λ(t) P_σ on R^n (Hermitian; iC since C real-skew
gives iC = Hermitian; P_σ real-symmetric permutation; λ real).

Trotter step (Strang split, second order):
  |ψ(t+dt)⟩ = exp(-i λ(t+dt/2) P_σ dt/2) · exp(-i (iC) dt) · exp(-i λ(t+dt/2) P_σ dt/2) |ψ(t)⟩
            = exp(-i λ_mid P_σ dt/2)     · exp(C dt)       · exp(-i λ_mid P_σ dt/2) |ψ⟩

  exp(-iα P_σ) v = cos(α) v - i sin(α) v[sigma_perm]   (P_σ involution)
  exp(C dt)    v = scipy.sparse.linalg.expm_multiply(dt * C_op, v)

Initial state: |q_0⟩ = first column of qubit basis = |σ=-1, τ=-1⟩, lifted
into the full n-dim register by V_4.

Observables (computed at every recorded time):
  ⟨P_τ⟩_full = ⟨ψ_full | P_τ | ψ_full⟩       (full register)
  ⟨P_τ⟩_V4   = ⟨ψ_V4 | P_τ_4D | ψ_V4⟩         (V_4 projection)
  V_4 occupancy = |V_4^T ψ_full|^2            (how much state lives in V_4)

If ⟨P_τ⟩_V4 traces the same open curve as v1/v2 AND ⟨P_τ⟩_full traces a
closed (or near-closed) loop, the v1/v2 open curve was a truncation
artifact (Reading 2). If both are open, the open-curve is a substrate
property (Reading 1, quasi-periodic torus).
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

HERE = Path(__file__).parent
M = 510510


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


def make_C_LinearOperator(R, tau_perm, rfft_kernel, n, m):
    """C = D_sym P_τ - P_τ D_sym, real skew, FFT-accelerated."""
    def D_sym(v):
        v_full = np.zeros(m, dtype=np.float64)
        v_full[R] = v
        return np.fft.irfft(np.fft.rfft(v_full) * rfft_kernel, n=m)[R]
    def C(v):
        v = np.asarray(v).ravel()
        return D_sym(v[tau_perm]) - D_sym(v)[tau_perm]
    return spla.LinearOperator(shape=(n, n), matvec=C,
                               rmatvec=lambda v: -C(v),
                               dtype=np.float64)


def extract_V4_and_qubit_basis(C_op, n, sigma_perm, tau_perm):
    """Return V_4 (n,4 columns), qubit_basis (4,4), sigma_assign,
    tau_assign, plus 4×4 P_sigma and P_tau."""
    print(f"  Top-4 SVD of C ...", flush=True)
    t0 = time.time()
    U, s, Vt = spla.svds(C_op, k=4, which='LM', tol=1e-10)
    idx = np.argsort(-s); s = s[idx]; Vt = Vt[idx, :]
    V_4 = Vt.T
    print(f"    {time.time()-t0:.1f}s, sigma_top = {[f'{x:.3e}' for x in s]}",
          flush=True)

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

    # V_4 columns rotated to qubit basis: V_4_qb has columns = qubit basis vectors in n-dim space
    V_4_qb = V_4 @ qb     # (n, 4) where column k = |q_k⟩ in n-dim register

    Psigma_qb = qb.T @ Psigma_4D @ qb
    Ptau_qb   = qb.T @ Ptau_4D @ qb
    return V_4_qb, sigma_a, tau_a, Psigma_qb, Ptau_qb, s


def main():
    print("=" * 78, flush=True)
    print(" Substrate hysteresis v3 — full-register discriminator ".center(78), flush=True)
    print("=" * 78, flush=True)

    print(f"\n  Building substrate at m = {M:,} ...", flush=True)
    t0 = time.time()
    R, tau_perm, sigma_perm, rfft_kernel, n = build_substrate(M)
    print(f"    n = phi(m) = {n:,}  ({time.time()-t0:.1f}s)", flush=True)

    C_op = make_C_LinearOperator(R, tau_perm, rfft_kernel, n, M)

    V_4_qb, sigma_a, tau_a, Psigma_4D_qb, Ptau_4D_qb, sigma_top = \
        extract_V4_and_qubit_basis(C_op, n, sigma_perm, tau_perm)

    omega_avg = 0.5 * (sigma_top[0] + sigma_top[2])
    print(f"\n  omega_avg = {omega_avg:.4e}", flush=True)

    # Initial state: |q_0> = |sigma=-1, tau=-1>, in n-dim register
    idx_q0 = int(np.where((sigma_a == -1) & (tau_a == -1))[0][0])
    psi = np.zeros(n, dtype=np.complex128)
    psi[:] = V_4_qb[:, idx_q0]   # column = |q_0⟩ in full n-dim register
    psi = psi / np.linalg.norm(psi)
    print(f"  Initial state: |q{idx_q0}⟩ = |sigma={sigma_a[idx_q0]:+d}, "
          f"tau={tau_a[idx_q0]:+d}⟩ lifted to full register",
          flush=True)
    init_V4_occ = float(np.abs(V_4_qb.conj().T @ psi).sum() ** 2)
    print(f"  V_4 occupancy at t=0: |V_4^T psi|^2 = {init_V4_occ:.6f}",
          flush=True)

    # Sweep: MEDIUM rate (T*omega = 100), 1 cycle, fewer steps for speed
    T_omega = 100.0
    T_period = T_omega / omega_avg
    n_steps = 800   # Strang Trotter steps per cycle
    lambda_max = omega_avg
    dt = T_period / n_steps

    print(f"\n  Drive: lambda(t) = lambda_max * sin^2(pi t / T)", flush=True)
    print(f"  T_period = {T_period:.4e},  n_steps = {n_steps},  dt = {dt:.4e}",
          flush=True)
    print(f"  lambda_max = {lambda_max:.4e}", flush=True)

    # Recording cadence: every ~step_record steps
    step_record = max(1, n_steps // 100)
    record_t, record_lambda, record_Ptau_full, record_Ptau_V4, record_V4_occ = \
        [], [], [], [], []

    def measure(psi_full, t, lam):
        """Measure ⟨P_τ⟩_full and ⟨P_τ⟩_V4 (projection)."""
        # Full register: ⟨ψ| P_τ |ψ⟩ where (P_τ ψ)[i] = ψ[tau_perm[i]]
        Ptau_full = float(np.real(psi_full.conj() @ psi_full[tau_perm]))
        # V_4 projection: ⟨ψ| (V_4 V_4^T P_τ V_4 V_4^T) |ψ⟩
        # = ⟨V_4^T ψ| P_τ_4D |V_4^T ψ⟩
        psi_V4 = V_4_qb.conj().T @ psi_full     # (4,)
        norm_V4 = np.linalg.norm(psi_V4)
        if norm_V4 > 0:
            psi_V4_normalized = psi_V4 / norm_V4
            Ptau_V4 = float(np.real(psi_V4_normalized.conj() @ Ptau_4D_qb @ psi_V4_normalized))
        else:
            Ptau_V4 = 0.0
        V4_occ = float(np.real(psi_V4.conj() @ psi_V4))
        return Ptau_full, Ptau_V4, V4_occ

    print(f"\n  Time-evolving (Strang Trotter, expm_multiply per step) ...",
          flush=True)
    t_evo_start = time.time()

    for k in range(n_steps + 1):
        t = k * dt
        lam = lambda_max * (math.sin(math.pi * t / T_period)) ** 2

        if k % step_record == 0 or k == n_steps:
            P_full, P_V4, V4_occ = measure(psi, t, lam)
            record_t.append(t)
            record_lambda.append(lam)
            record_Ptau_full.append(P_full)
            record_Ptau_V4.append(P_V4)
            record_V4_occ.append(V4_occ)
            if k % (step_record * 10) == 0:
                elapsed = time.time() - t_evo_start
                eta = elapsed * (n_steps - k) / max(k, 1)
                print(f"    step {k:>4d}/{n_steps}  t={t:.3e}  lam={lam:.3e}  "
                      f"P_τ_full={P_full:+.4f}  P_τ_V4={P_V4:+.4f}  "
                      f"V4_occ={V4_occ:.4f}  (wall {elapsed:.0f}s, ETA {eta:.0f}s)",
                      flush=True)

        if k < n_steps:
            t_mid = t + dt / 2
            lam_mid = lambda_max * (math.sin(math.pi * t_mid / T_period)) ** 2
            # Strang split:
            # 1. exp(-i lam_mid P_sigma dt/2): cos(lam dt/2) ψ - i sin(lam dt/2) ψ[sigma_perm]
            half_angle = lam_mid * dt / 2.0
            cos_h = math.cos(half_angle)
            sin_h = math.sin(half_angle)
            psi = cos_h * psi - 1j * sin_h * psi[sigma_perm]
            # 2. exp(-i (iC) dt) = exp(C dt) — real-orthogonal
            psi_real = psi.real.copy()
            psi_imag = psi.imag.copy()
            psi_real_new = spla.expm_multiply(dt * C_op, psi_real)
            psi_imag_new = spla.expm_multiply(dt * C_op, psi_imag)
            psi = psi_real_new + 1j * psi_imag_new
            # 3. exp(-i lam_mid P_sigma dt/2) again
            psi = cos_h * psi - 1j * sin_h * psi[sigma_perm]
            # Numerical hygiene
            psi = psi / np.linalg.norm(psi)

    t_evo = time.time() - t_evo_start
    print(f"\n  Time evolution complete: {t_evo:.0f}s", flush=True)

    # ---- Loop area ----
    # ∮ ⟨P_τ⟩ dλ on the (lambda, P_tau) trajectory for both full and V_4
    rec_lam = np.array(record_lambda)
    rec_Pf = np.array(record_Ptau_full)
    rec_Pv = np.array(record_Ptau_V4)
    rec_occ = np.array(record_V4_occ)
    loop_full = float(np.trapezoid(rec_Pf, rec_lam))
    loop_V4   = float(np.trapezoid(rec_Pv, rec_lam))

    print(f"\n{'=' * 78}", flush=True)
    print(" RESULT ".center(78), flush=True)
    print(f"{'=' * 78}", flush=True)
    print(f"\n  ⟨P_τ⟩_full(t=0) = {rec_Pf[0]:+.4f},  ⟨P_τ⟩_full(t=T) = {rec_Pf[-1]:+.4f}",
          flush=True)
    print(f"  ⟨P_τ⟩_V4(t=0)   = {rec_Pv[0]:+.4f},  ⟨P_τ⟩_V4(t=T)   = {rec_Pv[-1]:+.4f}",
          flush=True)
    print(f"  V_4 occupancy(t=0) = {rec_occ[0]:.4f}, (t=T) = {rec_occ[-1]:.4f},  "
          f"min over cycle = {rec_occ.min():.4f}", flush=True)
    print(f"\n  Loop area ∮⟨P_τ⟩_full dλ  = {loop_full:+.4e}", flush=True)
    print(f"  Loop area ∮⟨P_τ⟩_V4 dλ    = {loop_V4:+.4e}", flush=True)
    print(f"\n  Compare to v1/v2 MEDIUM (T*omega=100) result on 4x4-only: "
          f"~ +1.21e+06", flush=True)

    # ---- Discriminate ----
    # If V_4 occupancy stays close to 1 throughout, no leakage → 4×4 was OK
    # → Reading 1 (quasi-periodic torus) wins.
    # If V_4 occupancy drops significantly, leakage is real → Reading 2.
    occ_min = float(rec_occ.min())
    occ_drop = 1.0 - occ_min
    print(f"\n  Max V_4 leakage during cycle: {occ_drop:.4f} "
          f"(1 - min V_4 occupancy)", flush=True)

    if occ_drop < 0.01:
        verdict = "Reading 1 (quasi-periodic torus): V_4 occupancy stays > 99%, no significant leakage. The open curve is substrate-intrinsic."
    elif occ_drop > 0.10:
        verdict = "Reading 2 (truncation artifact): V_4 occupancy drops by >10% during cycle. The drive couples V_4 to bulk; 4x4 projection misses it."
    else:
        verdict = "Mixed / inconclusive: V_4 occupancy drops between 1% and 10%. Some leakage, but small. Need finer test."
    print(f"\n  VERDICT: {verdict}", flush=True)

    # ---- Save ----
    out = {
        "m": M, "n": int(n),
        "T_period": T_period,
        "n_steps": n_steps,
        "lambda_max": lambda_max,
        "sigma_top": [float(x) for x in sigma_top],
        "omega_avg": omega_avg,
        "trajectory": {
            "t": [float(x) for x in record_t],
            "lambda": [float(x) for x in record_lambda],
            "Ptau_full": [float(x) for x in record_Ptau_full],
            "Ptau_V4":   [float(x) for x in record_Ptau_V4],
            "V4_occupancy": [float(x) for x in record_V4_occ],
        },
        "loop_area_full": loop_full,
        "loop_area_V4": loop_V4,
        "V4_occupancy_min": occ_min,
        "V4_leakage_max": occ_drop,
        "verdict": verdict,
        "wall_time_s": t_evo,
    }
    out_path = HERE / "substrate_hysteresis_v3_full_register_results.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\n  Saved -> {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
