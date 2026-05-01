"""
sigma_bit_search_v32_rank32.py

The momentum-aware conjunction kernel concentrates 83% of energy in
rank-32 (above CREM's 81% threshold). 32 = 2^5 is a clean 5-bit binary
register width.

Hypothesis: the LEO momentum-aware conjunction operator may admit a
rank-32 invariant subspace with a multi-bit binary register structure
analogous to (but distinct from) the rank-4 σ⊗τ register of CREM. This
script tests whether such a structure exists.

Search procedure:
  1. Build D_LEO momentum kernel (from v4)
  2. Build C_LEO_mom = [D_mom, P_tau_inc]
  3. Top-32 SVD → V_32 (the rank-32 invariant subspace, 83% of energy)
  4. For each candidate physical involution P_i on the LEO catalog,
     project to V_32 (32x32 matrix) and test:
       (a) (P_i^32D)^2 = I
       (b) (P_i^32D)^T = P_i^32D (Hermitian)
       (c) eigenvalues = {-1, ..., -1, +1, ..., +1} with multiplicity 16/16
  5. For each commuting pair (and triple, ..., 5-tuple), test
     [P_i^32D, P_j^32D] = 0 on V_32
  6. Find the largest commuting set → that is the bit-count of the LEO
     binary register on V_32.

If 5 commuting binary involutions exist → 5-bit binary register on LEO
(32-state subspace, JCA-class energy concentration on V_32).

**Result: this hypothesis is FALSIFIED.** No 5-bit register exists; the
LEO momentum operator does not natively organize into a multi-bit binary
register structure on V_32 under any tested physical involution. The
negative result motivates the Shor-class τ-equivariant embedding pivot
of §8 in the preprint.

Candidates tested:
  P_tau_inc:  inclination-rank antipode (the established τ-bit)
  P_perigee:  altitude (perigee-rank) antipode
  P_apogee:   apogee-rank antipode
  P_ecc:      eccentricity-rank antipode
  P_mass:     mass-rank antipode
  P_xsect:    cross-section-rank antipode
"""
from __future__ import annotations

import csv
import json
import math
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).parent
EARTH_R_KM = 6378.137
MU_EARTH_KM3_S2 = 398600.4418
K_RANK = 32   # Target rank-K invariant subspace dimension


# ---------- catalog ----------
def load_leo_catalog(csv_path):
    rows = []
    with csv_path.open("r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            for k in ("mass_kg", "diameter_m", "xSectAvg", "sma_m", "ecc",
                     "inc_deg", "perigee_alt_km", "apogee_alt_km"):
                v = r.get(k, "")
                r[k] = float(v) if v not in ("", None) else None
            rows.append(r)
    keep = []
    for r in rows:
        if r.get("regime") != "LEO":
            continue
        peri, apo, inc, sma = (r.get("perigee_alt_km"), r.get("apogee_alt_km"),
                               r.get("inc_deg"), r.get("sma_m"))
        if any(x is None for x in (peri, apo, inc, sma)):
            continue
        if not (100.0 <= peri <= 2000.0):
            continue
        if apo < peri:
            continue
        keep.append(r)
    return keep


def build_momentum_kernel(rows, altitude_bin_km=25.0, inc_tol_deg=15.0,
                          min_overlap_km=5.0, k_top=80):
    """v4 momentum kernel: w_ij = (overlap/max_apo) * v_rel^2/(vi*vj) * sqrt(xsect)."""
    n = len(rows)
    perigees = np.array([r["perigee_alt_km"] for r in rows])
    apogees  = np.array([r["apogee_alt_km"] for r in rows])
    incs_deg = np.array([r["inc_deg"] for r in rows])
    incs_rad = np.deg2rad(incs_deg)
    xsects   = np.array([r.get("xSectAvg") or 1.0 for r in rows])

    # Orbital velocity at perigee
    v_orbital = np.zeros(n)
    for k, r in enumerate(rows):
        r_peri_km = r["perigee_alt_km"] + EARTH_R_KM
        a_km = r["sma_m"] / 1000.0
        if a_km > 0 and r_peri_km > 0:
            v_orbital[k] = math.sqrt(max(MU_EARTH_KM3_S2 * (2.0 / r_peri_km - 1.0 / a_km), 0.0))

    bin_idx = (perigees // altitude_bin_km).astype(int)
    bins = {}
    for i, b in enumerate(bin_idx):
        bins.setdefault(int(b), []).append(i)
    max_la = min(int(np.ceil((apogees - perigees).max() / altitude_bin_km)) + 1, 80)

    rows_idx, cols_idx, weights = [], [], []
    for b in sorted(bins.keys()):
        candidates = list(bins[b])
        for db in range(1, max_la + 1):
            if (b + db) in bins:
                candidates.extend(bins[b + db])
        for i in bins[b]:
            for j in candidates:
                if j <= i:
                    continue
                lo, hi = max(perigees[i], perigees[j]), min(apogees[i], apogees[j])
                overlap = hi - lo
                if overlap < min_overlap_km:
                    continue
                d_inc_deg = abs(incs_deg[i] - incs_deg[j])
                d_inc_retro = abs(180.0 - incs_deg[i] - incs_deg[j])
                if min(d_inc_deg, d_inc_retro) > inc_tol_deg:
                    continue
                base = (overlap / max(apogees[i], apogees[j])) * math.sqrt(xsects[i] * xsects[j])
                cos_th = math.cos(incs_rad[i] - incs_rad[j])
                vi, vj = v_orbital[i], v_orbital[j]
                if vi > 0 and vj > 0:
                    v_rel_sq = vi * vi + vj * vj - 2.0 * vi * vj * cos_th
                    mom = v_rel_sq / (vi * vj)
                else:
                    mom = 1.0
                w = base * mom
                rows_idx.append(i); cols_idx.append(j); weights.append(w)

    rows_arr = np.array(rows_idx + cols_idx, dtype=np.int32)
    cols_arr = np.array(cols_idx + rows_idx, dtype=np.int32)
    w_arr = np.array(weights + weights, dtype=np.float64)
    A_full = sp.csr_matrix((w_arr, (rows_arr, cols_arr)), shape=(n, n))

    A_csr = A_full.tocsr()
    rk, ck, wk = [], [], []
    for i in range(n):
        s, e = A_csr.indptr[i], A_csr.indptr[i + 1]
        if e - s <= k_top:
            rk.extend([i] * (e - s)); ck.extend(A_csr.indices[s:e]); wk.extend(A_csr.data[s:e])
        else:
            local = A_csr.data[s:e]
            top = np.argpartition(-local, k_top)[:k_top]
            rk.extend([i] * k_top); ck.extend(A_csr.indices[s:e][top]); wk.extend(local[top])
    A_pr = sp.csr_matrix((np.array(wk), (np.array(rk, dtype=np.int32),
                                          np.array(ck, dtype=np.int32))),
                         shape=(n, n))
    return A_pr.maximum(A_pr.T).tocsr()


def build_C_sparse(D, sigma):
    return (D[:, sigma] - D[sigma, :]).tocsr()


def rank_antipode_perm(values, n):
    sort_idx = np.argsort(values)
    inv_rank = np.empty(n, dtype=np.int64)
    inv_rank[sort_idx] = np.arange(n)
    sigma = np.empty(n, dtype=np.int64)
    for j in range(n):
        sigma[j] = sort_idx[n - 1 - inv_rank[j]]
    assert np.array_equal(sigma[sigma], np.arange(n))
    return sigma


def main():
    print("=" * 78)
    print(" σ-bit search on V_32 of LEO momentum kernel ".center(78))
    print(" (5-bit binary register hypothesis on V_32) ".center(78))
    print("=" * 78)

    csv_path = HERE / "discosweb_leo_summary.csv"
    rows = load_leo_catalog(csv_path)
    n = len(rows)
    print(f"\n  physical-LEO rows: {n:,}")

    print(f"\n  Building momentum kernel D_LEO_mom ...")
    t0 = time.time()
    D = build_momentum_kernel(rows)
    print(f"    {time.time()-t0:.1f}s, nnz = {D.nnz:,}")

    # Build P_tau_inc and C_mom
    incs = np.array([r["inc_deg"] for r in rows])
    sigma_inc = rank_antipode_perm(incs, n)
    print(f"\n  Building C_LEO_mom = [D_mom, P_tau_inc] ...")
    C_mom = build_C_sparse(D, sigma_inc)

    # Top-32 SVD → V_32
    print(f"\n  Computing top-{K_RANK} SVD of C_mom ...")
    t0 = time.time()
    U, s, Vt = spla.svds(C_mom, k=K_RANK, which="LM", tol=1e-7,
                          maxiter=10000, solver="arpack")
    idx = np.argsort(-s)
    s = s[idx]
    Vt = Vt[idx, :]
    V_K = Vt.T  # (n, K_RANK)
    print(f"    {time.time()-t0:.1f}s")
    print(f"    σ_0..15 = {[f'{x:.3e}' for x in s[:16]]}")
    print(f"    σ_16..31 = {[f'{x:.3e}' for x in s[16:32]]}")
    C_frob_sq = float((C_mom.data * C_mom.data).sum())
    rank_K_share = float(np.sum(s ** 2) / C_frob_sq)
    print(f"    rank-{K_RANK} share = {rank_K_share:.4f} of ||C||^2")

    # ---- Build candidate involutions ----
    print(f"\n  Building candidate physical involutions on the LEO catalog ...")
    candidates = {}
    candidates["tau_inc"]    = sigma_inc
    candidates["alt_perigee"] = rank_antipode_perm(np.array([r["perigee_alt_km"] for r in rows]), n)
    candidates["alt_apogee"]  = rank_antipode_perm(np.array([r["apogee_alt_km"] for r in rows]), n)
    candidates["ecc"]         = rank_antipode_perm(np.array([r["ecc"] for r in rows]), n)
    candidates["mass"]        = rank_antipode_perm(np.array([r.get("mass_kg") or 0.0 for r in rows]), n)
    candidates["xsect"]       = rank_antipode_perm(np.array([r.get("xSectAvg") or 0.0 for r in rows]), n)
    print(f"    {len(candidates)} candidates")

    # ---- Project each candidate to V_K ----
    print(f"\n  Projecting candidates to V_{K_RANK} (32x32) and testing involution algebra ...")
    P_proj = {}
    cand_results = {}
    for name, perm in candidates.items():
        P = V_K.T @ V_K[perm, :]      # (K, K) = V_K^T * P_perm * V_K
        # Symmetrize
        P_sym = 0.5 * (P + P.T)
        sym_err = np.linalg.norm(P - P.T)
        # (P_sym)^2 should be I
        I_K = np.eye(K_RANK)
        inv_err = np.linalg.norm(P_sym @ P_sym - I_K)
        # Eigenvalues
        eigs = np.linalg.eigvalsh(P_sym)
        eigs_sorted = sorted(eigs.tolist())
        # Distance to {±1}
        dist_to_pm1 = max(abs(e + 1) if e < 0 else abs(e - 1) for e in eigs)
        # Multiplicity of +1 vs -1 (with tolerance)
        n_plus = int(sum(1 for e in eigs if e > 0))
        n_minus = int(sum(1 for e in eigs if e < 0))
        is_balanced = (n_plus == n_minus == K_RANK // 2)
        passes_inv = inv_err < 1e-3
        passes_eig = dist_to_pm1 < 0.05
        P_proj[name] = P_sym
        cand_results[name] = {
            "sym_err": float(sym_err),
            "inv_err": float(inv_err),
            "eigs_min_max": (float(eigs.min()), float(eigs.max())),
            "n_plus": n_plus, "n_minus": n_minus,
            "balanced_16_16": is_balanced,
            "dist_to_pm1": float(dist_to_pm1),
            "passes_inv": passes_inv,
            "passes_eig": passes_eig,
            "qubit_ok": passes_inv and passes_eig,
        }

    print(f"\n  {'candidate':<14s} {'sym_err':>10s} {'inv_err':>10s} "
          f"{'eig_min':>9s} {'eig_max':>9s} {'+/-':>8s} {'PASS':>6s}")
    print(f"  {'-'*14} {'-'*10} {'-'*10} {'-'*9} {'-'*9} {'-'*8} {'-'*6}")
    for name, r in cand_results.items():
        eig_min, eig_max = r["eigs_min_max"]
        print(f"  {name:<14s} {r['sym_err']:>10.2e} {r['inv_err']:>10.2e} "
              f"{eig_min:>+9.4f} {eig_max:>+9.4f} "
              f"{r['n_plus']:>3d}/{r['n_minus']:<3d} "
              f"{('PASS' if r['qubit_ok'] else 'fail'):>6s}")

    # ---- Pairwise commutation on V_K among candidates that pass ----
    qubit_candidates = [name for name, r in cand_results.items() if r["qubit_ok"]]
    print(f"\n  Candidates passing involution algebra on V_{K_RANK}: {len(qubit_candidates)}")
    print(f"    {qubit_candidates}")

    if len(qubit_candidates) < 2:
        print(f"\n  Insufficient binary-register-eligible candidates for σ⊗τ algebra test.")
        print(f"  Note: even the τ_inc candidate may fail if K_RANK isn't a clean")
        print(f"  multiple of the inclination-bit dimension.")
        # Save partial
        return 1

    print(f"\n  --- Pairwise commutation on V_{K_RANK} (lower triangle) ---")
    print(f"    {'pair':<28s} {'comm_err':>10s} {'commutes?':>10s}")
    print(f"    {'-'*28} {'-'*10} {'-'*10}")
    pair_comm = {}
    for a, b in combinations(qubit_candidates, 2):
        Pa, Pb = P_proj[a], P_proj[b]
        comm_err = np.linalg.norm(Pa @ Pb - Pb @ Pa)
        commutes = comm_err < 1e-3
        pair_comm[(a, b)] = {"comm_err": float(comm_err), "commutes": bool(commutes)}
        print(f"    {a:>13s} <-> {b:<13s} {comm_err:>10.2e} "
              f"{('YES' if commutes else 'no'):>10s}")

    # ---- Find largest commuting subset ----
    print(f"\n  --- Finding largest mutually-commuting subset ---")
    best_subset = []
    for r in range(len(qubit_candidates), 0, -1):
        for subset in combinations(qubit_candidates, r):
            ok = True
            for a, b in combinations(subset, 2):
                key = (a, b) if (a, b) in pair_comm else (b, a)
                if not pair_comm[key]["commutes"]:
                    ok = False; break
            if ok:
                best_subset = list(subset)
                break
        if best_subset:
            break

    print(f"  Largest mutually-commuting subset on V_{K_RANK}: "
          f"{len(best_subset)} bits")
    print(f"    {best_subset}")
    print(f"  → LEO supports a {len(best_subset)}-bit binary register on V_{K_RANK}")
    if len(best_subset) >= 5:
        print(f"  ✓ FULL 5-bit register reached (32 = 2^5 sectors)")
    elif len(best_subset) >= 2:
        print(f"  Partial: {len(best_subset)}-bit binary register, "
              f"V_{K_RANK} populated with multiplicity 2^({5 - len(best_subset)}) per sector")
    else:
        print(f"  Single-Z/2 only: τ-bit found, no σ-pair commutes on V_{K_RANK}")

    # ---- Save ----
    def jsonable(v):
        if isinstance(v, (bool, np.bool_)): return bool(v)
        if isinstance(v, np.integer): return int(v)
        if isinstance(v, np.floating): return float(v)
        if isinstance(v, np.ndarray): return v.tolist()
        if isinstance(v, dict): return {str(k): jsonable(vv) for k, vv in v.items()}
        if isinstance(v, tuple): return [jsonable(x) for x in v]
        if isinstance(v, list): return [jsonable(x) for x in v]
        return v

    out = {
        "n_objects": n,
        "K_rank": K_RANK,
        "rank_K_share": float(rank_K_share),
        "sigma_top_K": [float(x) for x in s.tolist()],
        "candidate_results": jsonable(cand_results),
        "pairwise_commutation": jsonable({f"{a}__{b}": pc for (a, b), pc in pair_comm.items()}),
        "largest_commuting_subset": best_subset,
        "register_bits": len(best_subset),
    }
    out_path = HERE / "sigma_bit_search_v32_rank32_results.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\n  Saved -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
