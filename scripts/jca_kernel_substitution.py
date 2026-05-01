"""
JCA kernel substitution test — axiom-vs-emergent for the rank-4 collapse.

Question: Is the rank-4 boundary block of C(m) = [D_sym, P_tau] a property of:
  (A) the tent kernel D_sym specifically (rank-4 is an AXIOM of the JCA), or
  (B) the bi-involutive symmetry alone (rank-4 is EMERGENT from any K that
      commutes with P_sigma; the geometry of K is irrelevant)?

Test: substitute D_sym with alternative kernels K that all strictly commute
with P_sigma (additive reflection), build C_alt = [K, P_tau], compute the
top-4 singular-value Frobenius-energy fraction.

Outcome A (rank-4 holds only for tent): rank-4 is AXIOMATIC.
Outcome B (rank-4 holds for all kernels): rank-4 is EMERGENT — universal
theorem under the JCA symmetry on primorial lattices.

Run: python test_kernel_substitution.py
"""
from __future__ import annotations
import numpy as np
from math import gcd
import time
import json
from pathlib import Path


# ---------- substrate construction ----------

def coprime_residues(m: int) -> np.ndarray:
    """Sorted residues r in {1,...,m-1} with gcd(r,m)=1."""
    return np.array([r for r in range(1, m) if gcd(r, m) == 1], dtype=np.int64)


def build_perms(residues: np.ndarray, m: int) -> tuple[np.ndarray, np.ndarray]:
    """Build P_sigma (r -> -r mod m) and P_tau (r -> r^-1 mod m) as index permutations."""
    res_to_idx = {int(r): i for i, r in enumerate(residues)}
    sigma_perm = np.empty(len(residues), dtype=np.int64)
    tau_perm = np.empty(len(residues), dtype=np.int64)
    for i, r in enumerate(residues):
        r_int = int(r)
        sigma_perm[i] = res_to_idx[(m - r_int) % m]
        tau_perm[i] = res_to_idx[pow(r_int, -1, m)]
    return sigma_perm, tau_perm


# ---------- kernel constructors (all must commute with P_sigma) ----------

def tent_kernel(residues: np.ndarray, m: int) -> np.ndarray:
    """D_sym(r,s) = min(|r-s|, m-|r-s|). The CREM baseline."""
    diff = np.abs(residues[:, None] - residues[None, :]).astype(np.float64)
    return np.minimum(diff, m - diff)


def squared_tent_kernel(residues: np.ndarray, m: int) -> np.ndarray:
    """Tent squared. Still commutes with P_sigma; smoother profile."""
    return tent_kernel(residues, m) ** 2


def cosine_kernel(residues: np.ndarray, m: int) -> np.ndarray:
    """K_cos(r,s) = cos(2*pi*(r-s)/m). Harmonic; rank-1 on the full Z/m lattice."""
    diff = (residues[:, None] - residues[None, :]).astype(np.float64)
    return np.cos(2 * np.pi * diff / m)


def gaussian_kernel(residues: np.ndarray, m: int, sigma_frac: float = 0.1) -> np.ndarray:
    """K_gauss(r,s) = exp(-D_sym^2 / (2 sigma^2)) with sigma = sigma_frac * m."""
    sigma = sigma_frac * m
    d = tent_kernel(residues, m)
    return np.exp(-(d ** 2) / (2 * sigma ** 2))


def inverse_kernel(residues: np.ndarray, m: int) -> np.ndarray:
    """K_inv(r,s) = 1/(1+D_sym(r,s)). Long-range, slow decay."""
    d = tent_kernel(residues, m)
    return 1.0 / (1.0 + d)


def random_p_sigma_symmetric_kernel(residues: np.ndarray, m: int, seed: int = 0) -> np.ndarray:
    """A random kernel that respects P_sigma symmetry (control: no special geometry).

    Build a random symmetric matrix R, then symmetrize under P_sigma:
        K = (R + R[sigma_perm,:][:,sigma_perm]) / 2
    This is the MOST GENERIC P_sigma-commuting kernel, with no harmonic / metric structure.
    """
    rng = np.random.default_rng(seed)
    n = len(residues)
    R = rng.standard_normal((n, n))
    R = (R + R.T) / 2  # symmetric
    sigma_perm, _ = build_perms(residues, m)
    R_sigma = R[sigma_perm, :][:, sigma_perm]
    K = (R + R_sigma) / 2
    return K


# ---------- commutator + JCA-symmetry verification ----------

def build_commutator(K: np.ndarray, tau_perm: np.ndarray) -> np.ndarray:
    """C = K @ P_tau - P_tau @ K. P_tau acts by index permutation on the columns/rows."""
    K_Ptau = K[:, tau_perm]   # K @ P_tau
    Ptau_K = K[tau_perm, :]   # P_tau @ K
    return K_Ptau - Ptau_K


def verify_jca_symmetries(C: np.ndarray, sigma_perm: np.ndarray,
                           tau_perm: np.ndarray) -> tuple[float, float]:
    """Check [C, P_sigma] = 0 (commute) and {C, P_tau} = 0 (anti-commute).

    Returns relative Frobenius residuals.
    """
    fro_C = np.linalg.norm(C, ord='fro') + 1e-30

    C_Psigma = C[:, sigma_perm]
    Psigma_C = C[sigma_perm, :]
    comm_sigma = np.linalg.norm(C_Psigma - Psigma_C, ord='fro') / fro_C

    C_Ptau = C[:, tau_perm]
    Ptau_C = C[tau_perm, :]
    anticomm_tau = np.linalg.norm(C_Ptau + Ptau_C, ord='fro') / fro_C

    return comm_sigma, anticomm_tau


def rank_k_fraction(C: np.ndarray, k: int = 4) -> tuple[float, np.ndarray]:
    """Fraction of ||C||_F^2 captured by top-k singular values."""
    s = np.linalg.svd(C, compute_uv=False)
    total = np.sum(s ** 2)
    topk = np.sum(s[:k] ** 2)
    return float(topk / total), s[:max(8, k)]


# ---------- experiment driver ----------

def run_one(m: int, kernel_name: str, kernel_fn) -> dict:
    residues = coprime_residues(m)
    n = len(residues)
    sigma_perm, tau_perm = build_perms(residues, m)

    t0 = time.time()
    K = kernel_fn(residues, m)
    t_kernel = time.time() - t0

    t0 = time.time()
    C = build_commutator(K, tau_perm)
    t_comm = time.time() - t0

    comm_sigma, anticomm_tau = verify_jca_symmetries(C, sigma_perm, tau_perm)

    t0 = time.time()
    frac4, top_s = rank_k_fraction(C, k=4)
    t_svd = time.time() - t0

    fro_C = float(np.linalg.norm(C, ord='fro'))

    return {
        'm': m,
        'n': n,
        'kernel': kernel_name,
        'fro_C': fro_C,
        'jca_check_commute_sigma': float(comm_sigma),
        'jca_check_anticomm_tau': float(anticomm_tau),
        'rank4_fraction': frac4,
        'top8_singular_values': [float(x) for x in top_s[:8]],
        'wall_kernel_s': t_kernel,
        'wall_commutator_s': t_comm,
        'wall_svd_s': t_svd,
    }


def main():
    primorials = [
        (30, '3rd'),     # phi=8, sanity check
        (210, '4th'),    # phi=48
        (2310, '5th'),   # phi=480
        (30030, '6th'),  # phi=5760
    ]

    kernels = [
        ('tent (CREM baseline)', tent_kernel),
        ('squared tent', squared_tent_kernel),
        ('cosine harmonic', cosine_kernel),
        ('Gaussian (sigma=0.1*m)', lambda r, m: gaussian_kernel(r, m, 0.1)),
        ('inverse 1/(1+d)', inverse_kernel),
        ('random P_sigma-symmetric', lambda r, m: random_p_sigma_symmetric_kernel(r, m, seed=42)),
    ]

    print("=" * 130)
    print("JCA Kernel Substitution Test — Axiom vs Emergent")
    print("=" * 130)
    print(f"{'m':>8} {'n=phi(m)':>10} {'kernel':>30} "
          f"{'[C,P_sigma]/||C||':>20} {'{C,P_tau}/||C||':>20} {'rank4_frac':>14}")
    print("-" * 130)

    all_results = []
    for m, label in primorials:
        for kname, kfn in kernels:
            try:
                r = run_one(m, kname, kfn)
                all_results.append(r)
                # Highlight rank-4 fraction
                frac_str = f"{r['rank4_fraction']:.6f}"
                if r['rank4_fraction'] >= 0.98:
                    frac_str = "*" + frac_str + "*"  # collapse holds
                print(f"{r['m']:>8} {r['n']:>10} {r['kernel']:>30} "
                      f"{r['jca_check_commute_sigma']:>20.3e} "
                      f"{r['jca_check_anticomm_tau']:>20.3e} "
                      f"{frac_str:>14}")
            except Exception as e:
                print(f"{m:>8} {'-':>10} {kname:>30} FAILED: {e}")
        print("-" * 130)

    # Save raw results
    out = Path(__file__).parent / "jca_kernel_substitution_results.json"
    with open(out, 'w', encoding='utf-8') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to: {out}")

    # Verdict
    print("\n" + "=" * 130)
    print("VERDICT")
    print("=" * 130)
    largest_m = max(p[0] for p in primorials)
    final_results = [r for r in all_results if r['m'] == largest_m]
    tent_frac = next((r['rank4_fraction'] for r in final_results
                      if 'tent (CREM' in r['kernel']), None)
    nontent_fracs = [(r['kernel'], r['rank4_fraction']) for r in final_results
                     if 'tent (CREM' not in r['kernel']]

    print(f"\nAt m = {largest_m} (largest tested):")
    print(f"  CREM baseline (tent):  rank4_frac = {tent_frac:.6f}")
    for name, frac in nontent_fracs:
        verdict = "COLLAPSES" if frac >= 0.98 else "DOES NOT COLLAPSE"
        print(f"  {name:>30}: rank4_frac = {frac:.6f}  -> {verdict}")

    all_collapse = all(frac >= 0.98 for _, frac in nontent_fracs)
    if all_collapse:
        print("\n  --> ALL non-tent kernels exhibit rank-4 collapse.")
        print("      Rank-4 is EMERGENT from the JCA bi-involutive symmetry alone.")
        print("      This is a universal theorem on primorial lattices.")
    else:
        non_collapsing = [name for name, frac in nontent_fracs if frac < 0.98]
        print(f"\n  --> {len(non_collapsing)} non-tent kernel(s) FAIL to collapse to rank-4:")
        for name in non_collapsing:
            print(f"        - {name}")
        print("      Rank-4 is AXIOMATIC of the JCA — it is a property of the tent")
        print("      kernel's geometry, not of the bi-involutive symmetry alone.")
        print("      The JCA must be defined as: bi-involutive + rank-4 invariant subspace.")


if __name__ == '__main__':
    main()
