"""
d_leo_spectrum.py — what's the natural spectral structure of D_LEO?

The Schur-compliance test (schur_compliance_leo.py) showed that
sigma_4/sigma_0 of [D_LEO, P_tau] is 0.559 ± 0.0004 for ANY involution
P_tau (inclination, random) — a graph property, not a CREM signature.

The right next probe: top-8 eigenvalues of D_LEO directly. D_LEO is
real symmetric, so eigenvalues are real. CREM-class graphs have a
clean rank-4 boundary block with sigma_0/sigma_4 ≈ 9 (Triangle-Wave
Theorem: the tent kernel's 1/k^2 Fourier decay puts k=1 at 9x more
energy than k=3). If D_LEO has this structure, some involution might
yet pull out CREM signatures.

This probe also computes:
  - Top-8 eigenvalues of D_LEO (symmetric)
  - Top-8 eigenvalues of D_LEO normalized by sqrt(degree) (correlation
    matrix style)
  - Comparison to randomized D_LEO (preserve sparsity, shuffle weights)
"""
from __future__ import annotations

import csv
import math
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).parent

# Reuse the conjunction-graph builder
from schur_compliance_leo import (
    load_leo_catalog, build_conjunction_graph,
)


def top_k_eigh(A: sp.csr_matrix, k: int = 8) -> np.ndarray:
    """Top-k eigenvalues of symmetric A by absolute magnitude."""
    n = A.shape[0]
    eigs, _ = spla.eigsh(A, k=k, which="LM", tol=1e-10, maxiter=5000)
    # eigsh returns ascending magnitude; sort by absolute value descending
    idx = np.argsort(-np.abs(eigs))
    return eigs[idx]


def randomize_weights(A: sp.csr_matrix, seed: int) -> sp.csr_matrix:
    """Randomize edge weights while preserving sparsity pattern.
    Permutes the upper-triangular weights, then symmetrizes.
    """
    rng = np.random.default_rng(seed)
    A_coo = sp.triu(A, k=1).tocoo()
    perm = rng.permutation(len(A_coo.data))
    new_data = A_coo.data[perm]
    new_upper = sp.coo_matrix((new_data, (A_coo.row, A_coo.col)), shape=A.shape)
    return (new_upper + new_upper.T).tocsr()


def main():
    print("=" * 78)
    print(" D_LEO direct spectrum probe ".center(78))
    print("=" * 78)

    csv_path = HERE / "discosweb_leo_summary.csv"
    print(f"\n  Loading active LEO catalog ...")
    rows = load_leo_catalog(csv_path)
    print(f"    physical-LEO rows: {len(rows):,}")

    print(f"\n  Building D_LEO conjunction graph (k-NN prune top-80) ...")
    t0 = time.time()
    D = build_conjunction_graph(rows)
    print(f"    D_LEO: {D.shape[0]:,} x {D.shape[0]:,}, {D.nnz:,} nnz, "
          f"{time.time() - t0:.1f}s")

    # ---- Top-8 eigenvalues of D_LEO ----
    print(f"\n  --- Top-8 eigenvalues of D_LEO ---")
    t0 = time.time()
    e_D = top_k_eigh(D, k=8)
    print(f"    eigsh: {time.time() - t0:.1f}s")
    print(f"    eigenvalues |sorted desc|: {[f'{x:+.4e}' for x in e_D]}")
    print(f"    abs ratios |e_i / e_0|:")
    e_abs = np.abs(e_D)
    for i in range(8):
        print(f"      i={i}: |e_{i}|/|e_0| = {e_abs[i]/e_abs[0]:.4f}  "
              f"(CREM target i=4: 1/9 = {1/9:.4f})")

    # CREM signature on D directly: paired-doubled? — only happens if D
    # itself has a Z/2 symmetry. For symmetric D this isn't algebraic.
    pair_devs_D = []
    for i in range(0, 7, 2):
        ratio = e_abs[i + 1] / e_abs[i] if e_abs[i] > 0 else float("nan")
        pair_devs_D.append(abs(ratio - 1.0))
        print(f"    pair (e_{i}, e_{i+1}): |e_{i+1}|/|e_{i}| = {ratio:.6f}")
    print(f"    pair max-dev-from-1.0 (only meaningful if D has hidden Z/2): "
          f"{max(pair_devs_D):.4f}")

    sigma_4_over_0_D = e_abs[4] / e_abs[0]
    print(f"    e_4/e_0 (CREM target = 1/9 = {1/9:.4f}): {sigma_4_over_0_D:.4f}")

    # ---- Same on a degree-normalized correlation matrix ----
    print(f"\n  --- Top-8 eigenvalues of degree-normalized D_LEO ---")
    deg = np.array(D.sum(axis=1)).ravel()
    deg_inv_sqrt = np.where(deg > 0, 1.0 / np.sqrt(deg), 0.0)
    D_inv_sqrt = sp.diags(deg_inv_sqrt)
    D_norm = D_inv_sqrt @ D @ D_inv_sqrt
    D_norm = (D_norm + D_norm.T) * 0.5  # symmetry hygiene

    t0 = time.time()
    e_Dn = top_k_eigh(D_norm, k=8)
    print(f"    eigsh: {time.time() - t0:.1f}s")
    print(f"    eigenvalues: {[f'{x:+.4e}' for x in e_Dn]}")
    print(f"    abs ratios |e_i / e_0|:")
    e_Dn_abs = np.abs(e_Dn)
    for i in range(8):
        print(f"      i={i}: |e_{i}|/|e_0| = {e_Dn_abs[i]/e_Dn_abs[0]:.4f}")
    sigma_4_over_0_Dn = e_Dn_abs[4] / e_Dn_abs[0]
    print(f"    e_4/e_0 (CREM target = 1/9 = {1/9:.4f}): {sigma_4_over_0_Dn:.4f}")

    # ---- Randomized-weight baseline ----
    print(f"\n  --- Randomized-weight baseline (preserve sparsity, shuffle weights) ---")
    n_rand = 10
    rand_e4_e0 = []
    rand_top_e0 = []
    t0 = time.time()
    for k in range(n_rand):
        D_rand = randomize_weights(D, seed=k * 17 + 3)
        e_rand = top_k_eigh(D_rand, k=8)
        e_rand_abs = np.abs(e_rand)
        rand_e4_e0.append(e_rand_abs[4] / e_rand_abs[0])
        rand_top_e0.append(e_rand_abs[0])
    print(f"    {n_rand} randomizations: {time.time() - t0:.1f}s")
    rand_e4_e0 = np.array(rand_e4_e0)
    rand_top_e0 = np.array(rand_top_e0)
    print(f"    e_4/e_0 distribution: {rand_e4_e0.mean():.4f} +/- {rand_e4_e0.std():.4f}")
    print(f"    top eigenvalue e_0:   {rand_top_e0.mean():.4e} +/- {rand_top_e0.std():.4e}")
    print(f"    real D_LEO e_4/e_0:   {sigma_4_over_0_D:.4f}  "
          f"(z-score vs random: {(rand_e4_e0.mean() - sigma_4_over_0_D)/max(rand_e4_e0.std(), 1e-12):+.2f})")

    # ---- Verdict ----
    print(f"\n{'=' * 78}")
    print(" VERDICT ".center(78))
    print(f"{'=' * 78}")

    is_crem_rank4 = abs(sigma_4_over_0_D - 1/9) < 0.02
    is_distinguishable = abs((rand_e4_e0.mean() - sigma_4_over_0_D)
                              / max(rand_e4_e0.std(), 1e-12)) > 2.0
    print(f"\n  CREM-class rank-4 dominance on D_LEO directly:")
    print(f"    e_4/e_0 = {sigma_4_over_0_D:.4f}, CREM target 1/9 = {1/9:.4f}")
    print(f"    near target (within 0.02): {'YES' if is_crem_rank4 else 'NO'}")
    print(f"    distinguishable from randomized D_LEO: {'YES' if is_distinguishable else 'NO'}")
    print()

    if is_crem_rank4:
        print(f"  D_LEO has CREM-class rank-4 dominance. Some involution should")
        print(f"  pull out the full CREM signatures — re-test with alternative P_tau.")
    elif is_distinguishable:
        print(f"  D_LEO has structural rank-4 spectrum (different from random graph)")
        print(f"  but NOT the CREM 1/9 target. The graph has its own non-trivial")
        print(f"  spectral profile. Substrate-side framing needs adjustment.")
    else:
        print(f"  D_LEO's rank-4 ratio is indistinguishable from a randomized")
        print(f"  graph with same sparsity. The structure we see in the conjunction")
        print(f"  matrix is essentially that of a random sparse graph weighted by")
        print(f"  degree distribution — no special spectral structure beyond what")
        print(f"  comes from the degree heterogeneity (max degree 7,401, median 85).")
        print(f"  ")
        print(f"  Implication: the CREM substrate's transfer to LEO requires either")
        print(f"  (a) a different conjunction-weight definition that imparts a CREM-")
        print(f"      class kernel structure (e.g., RAAN-cyclic tent kernel), or")
        print(f"  (b) acceptance that LEO's conjunction operator is structurally")
        print(f"      distinct from CREM, and Hypothesis 10.3 needs reformulation.")
    print()


if __name__ == "__main__":
    main()
