"""
substrate_hysteresis_v2.py — bulletproofing v1.

v1 found ⟨P_tau⟩(t=T) != ⟨P_tau⟩(t=0) under cyclic sigma-flux drive at
multiple rates, with single-cycle "loop areas" scaling 4.5 orders of
magnitude with drive rate. v1 cannot distinguish between:

  (A) TRUE CLOSED-LOOP HYSTERESIS: every cycle traces the same
      (lambda, ⟨P_tau⟩) curve, area is per-cycle constant
      after settling, like a magnetic hysteresis loop. The
      mechanism would be Landau-Zener-mediated state-locking on
      the 21-ppm sigma-dependent tau-frequency split.

  (B) COHERENT OSCILLATION (Lissajous): every cycle traces a
      DIFFERENT curve as phase accumulates. Single-cycle "area"
      is non-zero but the trajectory is quasi-periodic on a
      4-torus, no steady-state loop.

v2 protocol:
  Drive 5 consecutive cycles at each rate.
  For each cycle k = 1..5, record:
    - per-cycle loop area  A_k = oint_k <P_tau> d_lambda
    - state at end of cycle  psi_k
    - <P_tau>(t = k*T)
    - <sigma>(t = k*T)

Verdict criterion:
  if std(A_1..A_5) / mean(|A|) < 0.05   -> TRUE HYSTERESIS  (converged loop)
  if std(A_1..A_5) / mean(|A|) > 0.30   -> COHERENT OSCILLATION
  in between                            -> partially-locked, report numbers

If TRUE HYSTERESIS at any rate, paper-2 has the substrate-level loop.
If COHERENT OSCILLATION at all rates, the substrate is unitarily
reversible-on-average and the Kessler hysteresis must come from
the embedding/dissipation, not the substrate alone.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as spla
from scipy.linalg import expm

# Reuse v1 substrate construction (paste of the helpers - keeps script
# standalone)

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


def make_C_matvec(R, tau_perm, rfft_kernel, m):
    def D_sym(v):
        v_full = np.zeros(m, dtype=np.float64)
        v_full[R] = v
        return np.fft.irfft(np.fft.rfft(v_full) * rfft_kernel, n=m)[R]
    def C(v):
        v = np.asarray(v).ravel()
        return D_sym(v[tau_perm]) - D_sym(v)[tau_perm]
    return C


def extract_qubit_basis(m):
    print(f"  Building substrate at m = {m:,} ...", flush=True)
    t0 = time.time()
    R, tau_perm, sigma_perm, rfft_kernel, n = build_substrate(m)
    print(f"    n = phi(m) = {n:,}  ({time.time()-t0:.1f}s)", flush=True)

    Cmv = make_C_matvec(R, tau_perm, rfft_kernel, m)
    C_op = spla.LinearOperator(shape=(n, n), matvec=Cmv,
                               rmatvec=lambda v: -Cmv(v), dtype=np.float64)
    print(f"  Top-4 SVD ...", flush=True)
    t0 = time.time()
    U, s, Vt = spla.svds(C_op, k=4, which='LM', tol=1e-10)
    idx = np.argsort(-s); s = s[idx]; Vt = Vt[idx, :]
    V_4 = Vt.T
    print(f"    {time.time()-t0:.1f}s, sigma_top = {[f'{x:.3e}' for x in s]}",
          flush=True)

    Psigma_4D = V_4.T @ V_4[sigma_perm, :]
    Ptau_4D   = V_4.T @ V_4[tau_perm, :]
    Psigma_4D = 0.5 * (Psigma_4D + Psigma_4D.T)
    Ptau_4D   = 0.5 * (Ptau_4D + Ptau_4D.T)

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

    CV4 = np.column_stack([Cmv(V_4[:, k]) for k in range(4)])
    C_4D = V_4.T @ CV4
    C_qbasis      = qb.T @ C_4D @ qb
    Psigma_qbasis = qb.T @ Psigma_4D @ qb
    Ptau_qbasis   = qb.T @ Ptau_4D @ qb
    return sigma_a, tau_a, C_qbasis, Psigma_qbasis, Ptau_qbasis, s


def multi_cycle_sweep(H_iC, P_sigma, P_tau, sigma_a, tau_a,
                      lambda_max, T_period, n_steps_per_cycle, n_cycles, label):
    """Drive n_cycles consecutive cycles. Track per-cycle loop area + state."""
    psi = np.zeros(4, dtype=np.complex128)
    idx_q0 = int(np.where((sigma_a == -1) & (tau_a == -1))[0][0])
    psi[idx_q0] = 1.0

    dt = T_period / n_steps_per_cycle

    # Storage per cycle
    cycle_areas = np.zeros(n_cycles)
    cycle_end_Ptau = np.zeros(n_cycles)
    cycle_end_sigma = np.zeros(n_cycles)
    cycle_end_psi_overlap_q0 = np.zeros(n_cycles)

    # Also store full trajectory for cycle 1 and cycle N (for shape compare)
    traj_lams_c1 = np.zeros(n_steps_per_cycle + 1)
    traj_Ptau_c1 = np.zeros(n_steps_per_cycle + 1)
    traj_lams_cN = np.zeros(n_steps_per_cycle + 1)
    traj_Ptau_cN = np.zeros(n_steps_per_cycle + 1)

    for c in range(n_cycles):
        lams = np.zeros(n_steps_per_cycle + 1)
        Pts = np.zeros(n_steps_per_cycle + 1)

        for k in range(n_steps_per_cycle + 1):
            t_in_cycle = k * dt
            lam = lambda_max * (math.sin(math.pi * t_in_cycle / T_period)) ** 2
            lams[k] = lam
            Pts[k] = float(np.real(psi.conj() @ P_tau @ psi))

            if k < n_steps_per_cycle:
                t_mid = t_in_cycle + dt / 2
                lam_mid = lambda_max * (math.sin(math.pi * t_mid / T_period)) ** 2
                H_t = H_iC + lam_mid * P_sigma
                U_dt = expm(-1j * dt * H_t)
                psi = U_dt @ psi
                psi = psi / np.linalg.norm(psi)

        # Per-cycle loop area
        area = 0.0
        for k in range(n_steps_per_cycle):
            d_lam = lams[k+1] - lams[k]
            avg = 0.5 * (Pts[k+1] + Pts[k])
            area += avg * d_lam
        cycle_areas[c] = area
        cycle_end_Ptau[c] = Pts[-1]
        cycle_end_sigma[c] = float(np.real(psi.conj() @ P_sigma @ psi))
        # store overlap of psi with initial state q0
        e_q0 = np.zeros(4, dtype=np.complex128); e_q0[idx_q0] = 1.0
        cycle_end_psi_overlap_q0[c] = float(np.abs(np.vdot(e_q0, psi))**2)

        if c == 0:
            traj_lams_c1[:] = lams
            traj_Ptau_c1[:] = Pts
        if c == n_cycles - 1:
            traj_lams_cN[:] = lams
            traj_Ptau_cN[:] = Pts

    # Convergence statistics
    abs_areas = np.abs(cycle_areas)
    mean_area = float(np.mean(abs_areas))
    std_area = float(np.std(cycle_areas))
    rel_std = std_area / max(mean_area, 1e-30)
    drift = float(cycle_end_Ptau[-1] - cycle_end_Ptau[0])

    # Trajectory-shape distance: how far apart cycle 1 and cycle N traverse
    # the same (lambda, P_tau) plane?  Compare at matched lambda.
    # Both cycles have lams[k] following the same sin^2 schedule, so we
    # can compare Ptau at the same k.
    traj_diff = np.max(np.abs(traj_Ptau_c1 - traj_Ptau_cN))

    print(f"\n  --- {label} ---", flush=True)
    print(f"    T_period = {T_period:.3e}  (n_steps_per_cycle = {n_steps_per_cycle})", flush=True)
    print(f"    lambda_max = {lambda_max:.3e}", flush=True)
    print(f"    Per-cycle loop areas:", flush=True)
    for c in range(n_cycles):
        print(f"      cycle {c+1}: area = {cycle_areas[c]:+.4e}, "
              f"P_tau(end) = {cycle_end_Ptau[c]:+.6f}, "
              f"|<q0|psi>|^2 = {cycle_end_psi_overlap_q0[c]:.4f}",
              flush=True)
    print(f"    mean |area|        = {mean_area:.4e}", flush=True)
    print(f"    std (area)         = {std_area:.4e}", flush=True)
    print(f"    rel std (area)     = {rel_std:.4f}", flush=True)
    print(f"    P_tau drift cyc1->cycN = {drift:+.4f}", flush=True)
    print(f"    max |traj_c1 - traj_cN| (in P_tau at matched lambda) = {traj_diff:.4f}",
          flush=True)

    return {
        "label": label,
        "T_period": T_period,
        "lambda_max": lambda_max,
        "cycle_areas": cycle_areas.tolist(),
        "cycle_end_Ptau": cycle_end_Ptau.tolist(),
        "cycle_end_sigma": cycle_end_sigma.tolist(),
        "cycle_end_psi_overlap_q0": cycle_end_psi_overlap_q0.tolist(),
        "mean_abs_area": mean_area,
        "std_area": std_area,
        "rel_std_area": rel_std,
        "Ptau_drift": drift,
        "traj_max_diff": float(traj_diff),
        "traj_c1_lambdas": traj_lams_c1.tolist(),
        "traj_c1_Ptau": traj_Ptau_c1.tolist(),
        "traj_cN_lambdas": traj_lams_cN.tolist(),
        "traj_cN_Ptau": traj_Ptau_cN.tolist(),
    }


def main():
    print("=" * 78, flush=True)
    print(" Substrate hysteresis v2 — multi-cycle bulletproof test ".center(78), flush=True)
    print("=" * 78, flush=True)

    sigma_a, tau_a, C_4D, P_sigma, P_tau, sigma_top = extract_qubit_basis(M)
    H_iC = 1j * C_4D

    sig_minus_idx = np.where(sigma_a == -1)[0]
    sig_plus_idx = np.where(sigma_a == +1)[0]
    omega_minus = float(np.abs(C_4D[np.ix_(sig_minus_idx, sig_minus_idx)][0, 1]))
    omega_plus  = float(np.abs(C_4D[np.ix_(sig_plus_idx,  sig_plus_idx)][0, 1]))
    omega_avg = 0.5 * (omega_minus + omega_plus)
    delta_omega = abs(omega_plus - omega_minus)
    print(f"\n  omega_avg     = {omega_avg:.4e}", flush=True)
    print(f"  delta_omega   = {delta_omega:.4e}  (rel split {delta_omega/omega_avg:.4e})",
          flush=True)

    lambda_max = omega_avg
    base_T = 1.0 / omega_avg
    n_cycles = 5

    sweeps = [
        ("SLOW (T*omega = 1000)",    base_T * 1000.0, 2000),
        ("MEDIUM (T*omega = 100)",   base_T * 100.0,  1500),
        ("MEDIUM-FAST (T*omega = 30)", base_T * 30.0, 1500),
        ("FAST (T*omega = 10)",      base_T * 10.0,   1500),
        ("VERY FAST (T*omega = 3)",  base_T * 3.0,    1000),
    ]

    results = []
    for label, T, n_steps in sweeps:
        r = multi_cycle_sweep(H_iC, P_sigma, P_tau, sigma_a, tau_a,
                              lambda_max, T, n_steps, n_cycles, label)
        results.append(r)

    # Verdict
    print(f"\n{'=' * 78}", flush=True)
    print(" VERDICT: hysteresis vs coherent oscillation ".center(78), flush=True)
    print(f"{'=' * 78}", flush=True)
    print(f"\n  {'rate':<28s}  {'mean|A|':>10s}  {'rel_std':>9s}  "
          f"{'drift':>9s}  {'traj_diff':>10s}  {'verdict':>20s}",
          flush=True)
    print(f"  {'-'*28}  {'-'*10}  {'-'*9}  {'-'*9}  {'-'*10}  {'-'*20}",
          flush=True)

    for r in results:
        if r["rel_std_area"] < 0.05 and r["traj_max_diff"] < 0.05:
            verdict = "TRUE HYSTERESIS"
        elif r["rel_std_area"] > 0.30 or r["traj_max_diff"] > 0.50:
            verdict = "COHERENT OSC."
        else:
            verdict = "partially-locked"
        print(f"  {r['label']:<28s}  {r['mean_abs_area']:>10.3e}  "
              f"{r['rel_std_area']:>9.4f}  {r['Ptau_drift']:>+9.4f}  "
              f"{r['traj_max_diff']:>10.4f}  {verdict:>20s}",
              flush=True)

    print(f"\n  Diagnostic legend:", flush=True)
    print(f"    rel_std    = std(A_1..A_N) / mean|A|   "
          f"-> 0 if loop closes per cycle", flush=True)
    print(f"    drift      = P_tau(N*T) - P_tau(T)     "
          f"-> 0 if state cycles back", flush=True)
    print(f"    traj_diff  = max |Ptau_cycle1(lam) - Ptau_cycleN(lam)|  "
          f"-> 0 if same path each cycle", flush=True)

    out = {
        "m": M,
        "n_cycles": n_cycles,
        "omega_avg": omega_avg,
        "delta_omega": delta_omega,
        "lambda_max": lambda_max,
        "sweeps": [{k: v for k, v in r.items()
                    if k not in ("traj_c1_lambdas", "traj_c1_Ptau",
                                 "traj_cN_lambdas", "traj_cN_Ptau")}
                   for r in results],
        "trajectories_c1_cN": {r["label"]: {
            "lambdas_c1": r["traj_c1_lambdas"], "Ptau_c1": r["traj_c1_Ptau"],
            "lambdas_cN": r["traj_cN_lambdas"], "Ptau_cN": r["traj_cN_Ptau"],
        } for r in results},
    }
    out_path = Path(__file__).parent / "substrate_hysteresis_v2_results.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\n  Saved -> {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())