"""
sigma_bit_search_leo.py — find the σ-bit of the σ⊗τ structure on V_4(LEO).

Per Tony's BVP framing: the τ-bit is set by P_τ_inc (inclination-rank
antipode) via the universal anticommutation {C_LEO, P_τ} = 0. The σ-bit
is whatever physical ℤ/2 on the LEO catalog COMMUTES with C_LEO on V_4
(equivalently, commutes with P_τ on V_4 and has eigenvalues ±1).

The natural Klein-4 group of orbital dynamics is {I, T, P, TP}:
  T = time reversal = inclination flip i → 180-i  (P_τ_inc, our τ-bit)
  P = orbital-plane reflection = nodal flip Ω → Ω+180  (the cleanest σ-bit candidate;
      requires RAAN, which is in raw envelope but not in summary CSV)
  TP = combination

Without RAAN, we test the remaining physically-meaningful candidates from
the catalog data we have:
  σ_alt   = altitude-band parity (perigee above/below median)
  σ_ecc   = eccentricity parity (above/below median)
  σ_mass  = mass parity (above/below median)
  σ_xsect = cross-section parity (above/below median)
  σ_apo   = apogee parity (above/below median, picks up LEO-bounded vs LEO-piercer)

Physical prediction: σ_alt is most likely to commute, because altitude-band
parity is preserved by conjunction physics (same-band conjunctions dominate)
and by drag dynamics on operational timescales.

For each candidate, project to V_4(C_LEO) and test the σ⊗τ algebra:
  (a) (P_σ^{4D})^2 = I  — involution on V_4
  (b) [P_σ^{4D}, P_τ^{4D}] = 0  — commutes with τ on V_4
  (c) eigenvalues of P_σ^{4D} are {-1, -1, +1, +1}  — clean qubit structure
  (d) [C_LEO, P_σ] = 0 on V_4  — equivalent to (b) by construction
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent

# Reuse infrastructure
from schur_compliance_leo import (
    load_leo_catalog, build_conjunction_graph,
    inclination_antipode_perm, make_C_LEO,
)
from scipy.sparse.linalg import svds


# ---------- σ-candidate involutions ----------
def rank_antipode(values: np.ndarray, label: str) -> tuple[np.ndarray, dict]:
    """Sort by `values` ascending, define σ(j) = obj at rank N-1-rank(j).
    Returns (perm, info_dict). True involution.
    """
    n = len(values)
    sort_idx = np.argsort(values)
    inv_rank = np.empty(n, dtype=np.int64)
    inv_rank[sort_idx] = np.arange(n)
    sigma = np.empty(n, dtype=np.int64)
    for j in range(n):
        sigma[j] = sort_idx[n - 1 - inv_rank[j]]
    assert np.array_equal(sigma[sigma], np.arange(n)), \
        f"σ_{label} not an involution"
    n_fixed = int(np.sum(sigma == np.arange(n)))
    return sigma, {"label": label, "n_fixed": n_fixed,
                   "median": float(np.median(values))}


def median_split_involution(values: np.ndarray, label: str,
                            rng_seed: int = 17) -> tuple[np.ndarray, dict]:
    """Pair objects across the median: each below-median obj pairs with the
    closest above-median obj. Cleaner physical interpretation than rank
    antipode (preserves the bin count explicitly).

    Falls back to greedy nearest-neighbor pairing across the median split.
    """
    n = len(values)
    median = float(np.median(values))
    below = np.where(values < median)[0]
    above = np.where(values >= median)[0]
    # Equal-size pairing: take min(|below|, |above|) pairs
    n_pairs = min(len(below), len(above))
    # Sort each half so pairing is monotone (low-low pairs with low-high, etc.)
    below_sorted = below[np.argsort(values[below])]
    above_sorted = above[np.argsort(values[above])]
    sigma = np.arange(n, dtype=np.int64)  # start as identity
    for k in range(n_pairs):
        a = int(below_sorted[k])
        b = int(above_sorted[k])
        sigma[a] = b
        sigma[b] = a
    # Unpaired (if any) stay as fixed points
    assert np.array_equal(sigma[sigma], np.arange(n)), \
        f"median-split for {label} not an involution"
    n_fixed = int(np.sum(sigma == np.arange(n)))
    return sigma, {"label": label, "n_fixed": n_fixed, "median": median,
                   "n_pairs": n_pairs}


# ---------- σ⊗τ algebra check on V_4 ----------
def test_sigma_tau_algebra(V_4: np.ndarray,
                           sigma_perm: np.ndarray,
                           tau_perm: np.ndarray,
                           label: str) -> dict:
    """Project P_σ and P_τ to V_4 (4x4), test the 2-qubit algebra.

    V_4: (n, 4) right singular vectors of C_LEO
    sigma_perm, tau_perm: candidate involutions on the n-object catalog
    """
    Ps = V_4.T @ V_4[sigma_perm, :]   # 4x4 projection
    Pt = V_4.T @ V_4[tau_perm, :]     # 4x4 projection
    I4 = np.eye(4)

    # (a) involution: (P_σ^4D)^2 = I
    inv_err_s = np.linalg.norm(Ps @ Ps - I4)
    inv_err_t = np.linalg.norm(Pt @ Pt - I4)

    # (b) commutation: [P_σ, P_τ] = 0 on V_4
    comm_err = np.linalg.norm(Ps @ Pt - Pt @ Ps)

    # (c) symmetric (skew-Hermitian projections of involutions become symmetric)
    sym_err_s = np.linalg.norm(Ps - Ps.T)
    sym_err_t = np.linalg.norm(Pt - Pt.T)

    # (d) eigenvalues of P_σ^{4D}
    es = sorted(np.linalg.eigvalsh(Ps).tolist())
    et = sorted(np.linalg.eigvalsh(Pt).tolist())
    es_close_pm1 = max(abs(es[0] + 1), abs(es[1] + 1),
                       abs(es[2] - 1), abs(es[3] - 1))
    et_close_pm1 = max(abs(et[0] + 1), abs(et[1] + 1),
                       abs(et[2] - 1), abs(et[3] - 1))

    # PASS criteria
    passes_inv  = inv_err_s < 1e-6
    passes_comm = comm_err < 1e-6
    passes_eig  = es_close_pm1 < 0.05  # within 5% of {±1}
    qubit_ok    = passes_inv and passes_comm and passes_eig

    return {
        "label": label,
        "inv_err_sigma": float(inv_err_s),
        "inv_err_tau": float(inv_err_t),
        "comm_err": float(comm_err),
        "sym_err_sigma": float(sym_err_s),
        "sym_err_tau": float(sym_err_t),
        "eigvals_sigma": es,
        "eigvals_tau": et,
        "eigvals_sigma_dist_to_pm1": float(es_close_pm1),
        "eigvals_tau_dist_to_pm1": float(et_close_pm1),
        "passes_involution": bool(passes_inv),
        "passes_commutation": bool(passes_comm),
        "passes_eigenvalue": bool(passes_eig),
        "qubit_register_ok": bool(qubit_ok),
    }


def main() -> int:
    print("=" * 78)
    print(" σ-bit search on V_4(C_LEO) ".center(78))
    print("=" * 78)

    csv_path = HERE / "discosweb_leo_summary.csv"
    print(f"\n  Loading active LEO catalog ...")
    rows = load_leo_catalog(csv_path)
    n = len(rows)
    print(f"    physical-LEO rows: {n:,}")

    print(f"\n  Building D_LEO + computing V_4 from C_LEO with τ = inc-antipode ...")
    t0 = time.time()
    D = build_conjunction_graph(rows)
    print(f"    D_LEO: {D.shape[0]:,} x {D.shape[0]:,}, {D.nnz:,} nnz, "
          f"{time.time() - t0:.1f}s")

    sigma_inc = inclination_antipode_perm(rows)
    C_op = make_C_LEO(D, sigma_inc)
    t0 = time.time()
    U, s, Vt = svds(C_op, k=4, which="LM", tol=1e-10, maxiter=5000)
    idx = np.argsort(-s); s = s[idx]; Vt = Vt[idx, :]
    V_4 = Vt.T
    print(f"    top-4 SVD: {time.time() - t0:.1f}s,  σ_0..3 = "
          f"{[f'{x:.3e}' for x in s]}")

    # ----- σ-bit candidates from physical orbital quantities -----
    print(f"\n  Building σ-bit candidates from physical orbital quantities ...")

    # Extract values
    perigee = np.array([r["perigee_alt_km"] for r in rows])
    apogee  = np.array([r["apogee_alt_km"] for r in rows])
    ecc     = np.array([r["ecc"] for r in rows])
    inc     = np.array([r["inc_deg"] for r in rows])
    mass    = np.array([r.get("mass_kg") or 0.0 for r in rows])
    xsect   = np.array([r.get("xSectAvg") or 0.0 for r in rows])

    candidates = []

    # Rank-antipode constructions (one per physical quantity)
    for name, vals in [
        ("altitude_perigee_rank", perigee),
        ("altitude_apogee_rank",  apogee),
        ("eccentricity_rank",     ecc),
        ("inclination_value_rank", inc),  # sanity check: inc-rank is the τ
        ("mass_rank",             mass),
        ("xsect_rank",            xsect),
    ]:
        perm, info = rank_antipode(vals, name)
        candidates.append((name, perm, info, "rank_antipode"))

    # Median-split (cleaner physical interpretation: low-half ↔ high-half)
    for name, vals in [
        ("altitude_perigee_median", perigee),
        ("eccentricity_median",     ecc),
    ]:
        perm, info = median_split_involution(vals, name)
        candidates.append((name, perm, info, "median_split"))

    # τ baseline (should give the IDENTITY σ-bit-on-V_4 — pure self-test)
    tau_self_perm = sigma_inc

    print(f"    {len(candidates)} candidate σ-bit involutions built")

    # ----- σ⊗τ algebra test for each -----
    print(f"\n{'=' * 78}")
    print(" σ⊗τ ALGEBRA TEST RESULTS ".center(78))
    print(f"{'=' * 78}\n")
    print(f"  {'candidate':<32s} {'inv':>8s} {'comm':>10s} {'σ-eigs ok':>12s} {'PASS':>6s}")
    print(f"  {'-'*32} {'-'*8} {'-'*10} {'-'*12} {'-'*6}")

    results = []
    for name, perm, info, kind in candidates:
        r = test_sigma_tau_algebra(V_4, perm, sigma_inc, name)
        r["kind"] = kind
        r["info"] = info
        results.append(r)
        print(f"  {name:<32s} {r['inv_err_sigma']:>8.1e} "
              f"{r['comm_err']:>10.1e} "
              f"{r['eigvals_sigma_dist_to_pm1']:>12.4f} "
              f"{('PASS' if r['qubit_register_ok'] else 'fail'):>6s}")

    # ----- Show details of any PASS candidates -----
    passes = [r for r in results if r["qubit_register_ok"]]
    print(f"\n  --- {len(passes)} candidate(s) pass the σ⊗τ algebra test ---\n")

    for r in passes:
        print(f"  σ-candidate: {r['label']}")
        print(f"    eigenvalues of P_σ^4D: {[f'{x:+.4f}' for x in r['eigvals_sigma']]}")
        print(f"    eigenvalues of P_τ^4D: {[f'{x:+.4f}' for x in r['eigvals_tau']]}")
        print(f"    [P_σ, P_τ] error:      {r['comm_err']:.4e}")

    # If none pass, show closest near-misses
    if not passes:
        print(f"  No σ-candidate passes the strict σ⊗τ algebra test.")
        print(f"  Top near-misses by commutation error:\n")
        sorted_res = sorted(results, key=lambda r: r["comm_err"])
        for r in sorted_res[:5]:
            print(f"  {r['label']:<32s}  comm_err={r['comm_err']:.4e}  "
                  f"σ-eig-dist={r['eigvals_sigma_dist_to_pm1']:.4f}  "
                  f"σ-eigs={[f'{x:+.3f}' for x in r['eigvals_sigma']]}")

    # Save raw results
    out_path = HERE / "sigma_bit_search_leo_results.json"
    out_path.write_text(json.dumps({
        "n_objects": n,
        "n_edges": int(D.nnz / 2),
        "tau_bit": "inclination_rank_antipode",
        "V_4_singular_values": s.tolist(),
        "sigma_candidates": [
            {k: v for k, v in r.items() if k != "info"} | {"info": r["info"]}
            for r in results
        ],
        "n_passes": len(passes),
        "pass_labels": [r["label"] for r in passes],
    }, indent=2), encoding="utf-8")
    print(f"\n  Saved -> {out_path}")
    return 0 if passes else 2


if __name__ == "__main__":
    raise SystemExit(main())
