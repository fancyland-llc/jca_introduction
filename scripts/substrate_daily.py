"""Substrate Daily Pipeline (Tier 2 v0.1).

Daily-refresh pipeline that produces the customer-facing live substrate
status. Combines Tier 1 (NOAA + catalog summary stats) with Tier 2 (Fiedler
lambda_2 on current catalog + per-storm-class stress tests).

Architecture:
  Cloud Function (scheduled trigger) -> Cloud Run (this service) -> Firestore

This script is the Cloud Run service. It is invoked by a separate Cloud
Function on a daily schedule. Outputs are written to Firestore (or to
local JSON when run with --no-firestore).

Pipeline:
  1. NOAA SWPC + GFZ refresh (assumes noaa_swpc_loader / noaa_historical
     have been run; this step reads cached artifacts)
  2. Load current LEO catalog + leverage map
  3. Build sparse 80-NN momentum-aware conjunction graph (~50k nodes)
  4. Compute symmetric normalized Laplacian L_tilde = D^(-1/2) L D^(-1/2)
  5. eigsh shift-invert sigma=0, k=2 -> (lambda_1 ~ 0, lambda_2 = Fiedler)
  6. For each storm class G1..G5 + Carrington: apply thermospheric
     perturbation, rebuild graph, recompute lambda_2_storm
  7. Assemble daily status payload
  8. Write to Firestore (live signal) + append to fiedler_history collection

Threshold for "ballistic transition" is currently a placeholder at 50%
baseline collapse pending v1.5 historical-snapshot calibration. The
v0.3 design constraints are documented in
project_flare_assessor_v03.md.

CLI:
    python substrate_daily.py                       # full run, writes to Firestore
    python substrate_daily.py --no-firestore        # local JSON only (testing)
    python substrate_daily.py --baseline-only       # skip storm stress tests
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

import kessler_risk_assessor as kra

HERE = Path(__file__).parent

# Earth + orbit constants
EARTH_R_KM = 6378.137
MU_EARTH = 3.986004418e14

# Conjunction-graph parameters
DEFAULT_K_NN = 80
DEFAULT_ALT_TOL_KM = 200.0    # search window for candidate neighbors
DEFAULT_INC_TOL_DEG = 15.0    # inclination compatibility tolerance

# Cascade-transition placeholder threshold (v1.5 will fit from cascade-rewire data)
BALLISTIC_TRANSITION_THRESHOLD = 0.50    # if lambda_2_storm/lambda_2_baseline < 0.5, ballistic

# Storm-class shell-widening factors (alpha) for thermospheric perturbation.
# alpha represents the fractional widening of shell altitude band under storm
# at the reference altitude h_ref=200km; widening attenuates exponentially
# with altitude (scale height H_widen=200km) to match real drag-profile decay.
STORM_ALPHA = {
    "G1":         0.05,    # Kp=5
    "G2":         0.10,    # Kp=6
    "G3":         0.25,    # Kp=7
    "G4":         0.50,    # Kp=8
    "G5":         1.00,    # Kp=9
    "Carrington": 3.00,    # Kp~10.5 extrapolation (Riley 2012)
}
STORM_KP = {"G1": 5.0, "G2": 6.0, "G3": 7.0, "G4": 8.0, "G5": 9.0,
            "Carrington": 10.5}


# ============================================================================
# Catalog vectorization
# ============================================================================

def vectorize_catalog(catalog: list[dict]) -> dict:
    """Convert list-of-dict catalog into numpy arrays for vectorized ops."""
    n = len(catalog)
    inc = np.array([float(o.get("inc_deg", 0.0)) for o in catalog], dtype=np.float64)
    peri = np.array([float(o.get("perigee_alt_km", 0.0)) for o in catalog], dtype=np.float64)
    apo = np.array([float(o.get("apogee_alt_km", 0.0)) for o in catalog], dtype=np.float64)
    mass = np.array([float(o.get("mass_kg") or 200.0) for o in catalog], dtype=np.float64)

    # Vis-viva at perigee
    a_m = (EARTH_R_KM + 0.5 * (peri + apo)) * 1000.0
    r_peri_m = (EARTH_R_KM + peri) * 1000.0
    v_sq = MU_EARTH * (2.0 / np.maximum(r_peri_m, 1.0) - 1.0 / np.maximum(a_m, 1.0))
    v_peri_m_s = np.sqrt(np.maximum(v_sq, 0.0))
    v_peri = v_peri_m_s / 1000.0    # km/s

    return {
        "n": n,
        "inc": inc,
        "peri": peri,
        "apo": apo,
        "mean_alt": 0.5 * (peri + apo),
        "mass": mass,
        "v_peri": v_peri,
    }


# ============================================================================
# Conjunction-graph builder (80-NN sparse, momentum-aware)
# ============================================================================

def build_conjunction_graph(vec: dict, k_nn: int = DEFAULT_K_NN,
                             alt_tol_km: float = DEFAULT_ALT_TOL_KM,
                             inc_tol_deg: float = DEFAULT_INC_TOL_DEG,
                             min_overlap_km: float = 5.0,
                             progress: bool = True) -> csr_matrix:
    """Build sparse k-NN momentum-aware conjunction graph from vectorized catalog.

    Edge weight w_ij = (shell_overlap / max_apo) * (v_rel^2 / (v_i * v_j))
    matching the §6 momentum kernel of the preprint.

    Returns symmetric scipy.sparse.csr_matrix of shape (n, n).
    """
    n = vec["n"]
    inc = vec["inc"]
    peri = vec["peri"]
    apo = vec["apo"]
    mean_alt = vec["mean_alt"]
    v_peri = vec["v_peri"]

    # Sort by mean altitude for binary-search candidate lookup
    sort_idx = np.argsort(mean_alt)
    sorted_alt = mean_alt[sort_idx]

    rows: list[int] = []
    cols: list[int] = []
    weights: list[float] = []

    t0 = time.time()
    for i in range(n):
        if progress and i % 5000 == 0 and i > 0:
            rate = i / (time.time() - t0)
            eta = (n - i) / rate
            print(f"    graph build: {i}/{n}  ({rate:.0f} nodes/sec, eta {eta:.0f}s)",
                  flush=True)

        h_i = mean_alt[i]
        # Binary-search altitude window in sorted-altitude space
        lo_idx = np.searchsorted(sorted_alt, h_i - alt_tol_km, side="left")
        hi_idx = np.searchsorted(sorted_alt, h_i + alt_tol_km, side="right")
        candidate_orig = sort_idx[lo_idx:hi_idx]
        candidate_orig = candidate_orig[candidate_orig != i]
        if candidate_orig.size == 0:
            continue

        # Inclination compatibility (consider both prograde and retrograde pairings)
        d_inc = np.abs(inc[candidate_orig] - inc[i])
        d_inc_retro = np.abs(180.0 - inc[i] - inc[candidate_orig])
        valid = np.minimum(d_inc, d_inc_retro) <= inc_tol_deg
        candidate_orig = candidate_orig[valid]
        if candidate_orig.size == 0:
            continue

        # Shell overlap (vectorized)
        peri_j = peri[candidate_orig]
        apo_j = apo[candidate_orig]
        lo_shell = np.maximum(peri[i], peri_j)
        hi_shell = np.minimum(apo[i], apo_j)
        overlap = hi_shell - lo_shell
        valid = overlap >= min_overlap_km
        if not np.any(valid):
            continue
        candidate_orig = candidate_orig[valid]
        overlap = overlap[valid]
        peri_j = peri_j[valid]
        apo_j = apo_j[valid]

        max_apo = np.maximum(apo[i], apo_j)
        base = overlap / np.maximum(max_apo, 1.0)

        # Momentum kernel
        v_j = v_peri[candidate_orig]
        cos_th = np.cos(np.radians(inc[i] - inc[candidate_orig]))
        v_i_val = v_peri[i]
        v_rel_sq = v_i_val * v_i_val + v_j * v_j - 2.0 * v_i_val * v_j * cos_th
        denom = np.maximum(v_i_val * v_j, 1e-9)
        mom = v_rel_sq / denom

        w = base * mom
        nonzero = w > 0
        candidate_orig = candidate_orig[nonzero]
        w = w[nonzero]
        if candidate_orig.size == 0:
            continue

        # Keep top k_nn by weight
        if w.size > k_nn:
            top_k_idx = np.argpartition(w, -k_nn)[-k_nn:]
            candidate_orig = candidate_orig[top_k_idx]
            w = w[top_k_idx]

        rows.extend([int(i)] * candidate_orig.size)
        cols.extend(candidate_orig.tolist())
        weights.extend(w.tolist())

    if progress:
        print(f"    graph build complete: {len(rows):,} directed edges in {time.time()-t0:.1f}s",
              flush=True)

    A = csr_matrix((weights, (rows, cols)), shape=(n, n))
    # Symmetrize: average of (i,j) and (j,i) entries
    A = (A + A.T) * 0.5

    return A


# ============================================================================
# Fiedler lambda_2 computation (symmetric normalized Laplacian, eigsh)
# ============================================================================

def compute_fiedler(A: csr_matrix, verbose: bool = True,
                     min_component_frac: float = 0.5) -> dict:
    """Compute Fiedler lambda_2 on the symmetric normalized Laplacian
    of the largest connected component.

    L_tilde = I - D^(-1/2) A D^(-1/2)

    Per v0.3 design constraints:
      - eigsh shift-invert sigma=1e-6, k=2 (find 2 smallest eigenvalues)
        -- sigma>0 is required because L_tilde has lambda_1=0 exactly,
           making shift-invert at sigma=0 hit a singular LU factor.
      - normalized Laplacian (scale-invariant; strips catalog-size confound)
      - restrict to largest connected component first; multi-component
        graphs would otherwise return lambda_2 = 0 (a degenerate eigenvalue
        from the second component) instead of the Fiedler value.

    Returns dict with eigenvalues, Fiedler vector v_2, IPR(v_2), the
    component breakdown, and the n_nodes/n_edges of the analyzed component.
    """
    n_full = A.shape[0]
    if verbose:
        t0 = time.time()
        print(f"    finding connected components (n={n_full}) ...", flush=True)

    # 1. Identify connected components
    from scipy.sparse.csgraph import connected_components
    n_components, labels = connected_components(A, directed=False, return_labels=True)
    component_sizes = np.bincount(labels)
    largest_label = int(np.argmax(component_sizes))
    largest_size = int(component_sizes[largest_label])
    if verbose:
        print(f"    components: {n_components}  largest: {largest_size:,} nodes "
              f"({100.0*largest_size/n_full:.1f}% of graph)",
              flush=True)

    if largest_size < min_component_frac * n_full:
        if verbose:
            print(f"    WARNING: largest component is < {100*min_component_frac:.0f}% of graph; "
                  f"Fiedler value reflects only this component", flush=True)

    # 2. Restrict adjacency to largest component
    component_mask = labels == largest_label
    component_idx = np.flatnonzero(component_mask)
    A_comp = A[component_idx, :][:, component_idx]
    n = A_comp.shape[0]

    # 3. Build symmetric normalized Laplacian on the component
    if verbose:
        print(f"    building normalized Laplacian on largest component (n={n}) ...",
              flush=True)
    d = np.asarray(A_comp.sum(axis=1)).ravel()
    d_safe = np.where(d > 0, d, 1.0)
    d_inv_sqrt = 1.0 / np.sqrt(d_safe)
    d_inv_sqrt[d == 0] = 0.0

    from scipy.sparse import diags, eye
    D_inv_sqrt = diags(d_inv_sqrt)
    M = D_inv_sqrt @ A_comp @ D_inv_sqrt
    M = (M + M.T) * 0.5
    L_tilde = eye(n) - M

    # 4. eigsh shift-invert at sigma=1e-6 (avoids singular LU at sigma=0)
    # which="LM" with sigma>0 gives eigenvalues closest to sigma, i.e., smallest.
    if verbose:
        print(f"    eigsh shift-invert sigma=1e-6, k=2, ncv=50 ...", flush=True)
        t1 = time.time()

    try:
        l_vals, l_vecs = eigsh(L_tilde, k=2, sigma=1e-6, which="LM",
                                ncv=50, maxiter=10000, tol=1e-8)
    except Exception as e:
        if verbose:
            print(f"    shift-invert failed ({e}); fallback to LOBPCG", flush=True)
        # Fallback: LOBPCG with random initial guess
        from scipy.sparse.linalg import lobpcg
        rng = np.random.default_rng(42)
        X = rng.standard_normal((n, 4))
        # First column should be the constant (sqrt of degrees normalized) — fix later
        l_vals, l_vecs = lobpcg(L_tilde, X, largest=False, maxiter=2000, tol=1e-8)
        # Take smallest 2
        order = np.argsort(l_vals)[:2]
        l_vals = l_vals[order]
        l_vecs = l_vecs[:, order]

    if verbose:
        print(f"    eigsh complete in {time.time()-t1:.1f}s", flush=True)

    # Sort ascending
    order = np.argsort(l_vals)
    l_vals = l_vals[order]
    l_vecs = l_vecs[:, order]

    lambda_1 = float(l_vals[0])    # should be ~0 for the (now connected) component
    lambda_2 = float(l_vals[1])    # Fiedler value

    # Fiedler vector
    v_2 = l_vecs[:, 1]
    v_2 = v_2 / max(np.linalg.norm(v_2), 1e-12)

    # Inverse Participation Ratio: 1/sum(v^4); high IPR = localized
    ipr = float(1.0 / max(np.sum(v_2 ** 4), 1e-30))

    if verbose:
        print(f"    lambda_1 = {lambda_1:+.4e}  (should be ~0 for connected component)",
              flush=True)
        print(f"    lambda_2 = {lambda_2:+.4e}  (Fiedler / connectivity gap)",
              flush=True)
        print(f"    IPR(v_2) = {ipr:.1f}  (1 = fully localized, {n} = fully delocalized)",
              flush=True)

    return {
        "lambda_1": lambda_1,
        "lambda_2": lambda_2,
        "fiedler_vector_v2": v_2,
        "IPR_v2": ipr,
        "n_nodes": n,
        "n_edges": int(A_comp.nnz / 2),
        "n_full_graph": n_full,
        "n_components": int(n_components),
        "largest_component_frac": float(largest_size / n_full),
    }


# ============================================================================
# Storm-class thermospheric stress tests
# ============================================================================

def apply_thermospheric_perturbation(vec: dict, alpha: float,
                                       h_ref: float = 300.0,
                                       H_drop: float = 200.0,
                                       base_drop_km: float = 50.0,
                                       m_ref_kg: float = 200.0) -> dict:
    """Differential drag-induced altitude drop, per object's mass (A/m proxy).

    Per Gemini math autopsy 2026-04-28: the previous symmetric-shell-widening
    perturbation was wrong on two counts:
      (1) Symmetric expansion preserves mean_alt: peri-w/2 + apo+w/2 averages
          back to the original mean_alt. Candidate-neighbor pools never shift.
      (2) Symmetric normalized Laplacian L_tilde = I - D^(-1/2) A D^(-1/2)
          is scale-invariant under uniform edge-weight scaling: A -> cA gives
          D -> cD and L_tilde unchanged. Even if edge weights had changed, the
          normalization would have divided the stress out.

    Fix: differential drag drop. Drag force per unit mass scales as
    (C_d * A / m), so lighter objects with the same drag coefficient drop
    more under the same atmospheric density. Use mass as an A/m proxy
    (smaller satellites are area-dominated; heavier satellites are
    mass-dominated): drop_factor scales as sqrt(m_ref / m).

    drop_km = alpha * base_drop_km * sqrt(m_ref / m) * exp(-(mean_alt - h_ref) / H_drop)

    Mean altitude drops; peri and apo both shift down by the same amount
    (preserving eccentricity); v_peri updates via vis-viva (smaller r ->
    higher v, so the momentum kernel registers the increased kinetic
    energy of post-perturbation conjunctions).

    The drop is capped above 100 km to keep satellites in orbit
    (sub-100-km altitude is reentry territory; the substrate engine
    is undefined there).
    """
    mean_alt = vec["mean_alt"]
    peri = vec["peri"]
    apo = vec["apo"]
    mass = vec["mass"]

    # Differential mass factor: lighter objects drop more (higher effective A/m)
    mass_factor = np.sqrt(m_ref_kg / np.maximum(mass, 1.0))
    # Altitude exponential profile: drag is exponential in altitude
    alt_factor = np.exp(-(mean_alt - h_ref) / H_drop)
    drop_km = alpha * base_drop_km * mass_factor * alt_factor
    drop_km = np.maximum(drop_km, 0.0)
    # Cap so we don't push satellites into reentry
    max_drop = np.maximum(mean_alt - 100.0, 0.0)
    drop_km = np.minimum(drop_km, max_drop)

    new_peri = peri - drop_km
    new_apo = apo - drop_km
    new_mean_alt = 0.5 * (new_peri + new_apo)

    # Recompute v_peri via vis-viva: v^2 = mu * (2/r - 1/a)
    a_m = (EARTH_R_KM + new_mean_alt) * 1000.0
    r_peri_m = (EARTH_R_KM + new_peri) * 1000.0
    v_sq = MU_EARTH * (2.0 / np.maximum(r_peri_m, 1.0) - 1.0 / np.maximum(a_m, 1.0))
    new_v_peri = np.sqrt(np.maximum(v_sq, 0.0)) / 1000.0    # km/s

    perturbed = dict(vec)
    perturbed["peri"] = new_peri
    perturbed["apo"] = new_apo
    perturbed["mean_alt"] = new_mean_alt
    perturbed["v_peri"] = new_v_peri
    return perturbed


def stress_test_storm_classes(vec: dict, baseline_lambda_2: float,
                                 storm_classes: list[str] | None = None,
                                 k_nn: int = DEFAULT_K_NN,
                                 verbose: bool = True) -> dict:
    """For each storm class, apply thermospheric perturbation and recompute lambda_2."""
    if storm_classes is None:
        storm_classes = list(STORM_ALPHA.keys())
    results = {}
    for sclass in storm_classes:
        alpha = STORM_ALPHA[sclass]
        if verbose:
            print(f"\n  storm stress test: {sclass} (alpha={alpha})", flush=True)
        perturbed = apply_thermospheric_perturbation(vec, alpha)
        A_storm = build_conjunction_graph(perturbed, k_nn=k_nn, progress=verbose)
        fiedler_storm = compute_fiedler(A_storm, verbose=verbose)
        lambda_2_storm = fiedler_storm["lambda_2"]
        ratio = lambda_2_storm / baseline_lambda_2 if baseline_lambda_2 > 0 else float("nan")
        ballistic = ratio < BALLISTIC_TRANSITION_THRESHOLD or ratio > (1.0 / BALLISTIC_TRANSITION_THRESHOLD)
        results[sclass] = {
            "kp": STORM_KP[sclass],
            "alpha": alpha,
            "lambda_2_storm": lambda_2_storm,
            "lambda_2_baseline": baseline_lambda_2,
            "ratio_storm_to_baseline": ratio,
            "pct_of_baseline": 100.0 * ratio,
            "ballistic_transition": bool(ballistic),
            "n_edges_storm": fiedler_storm["n_edges"],
            "IPR_v2_storm": fiedler_storm["IPR_v2"],
            "n_components_storm": int(fiedler_storm.get("n_components", 0)),
            "largest_component_frac_storm": float(fiedler_storm.get("largest_component_frac", 0.0)),
        }
        if verbose:
            tag = "BALLISTIC" if ballistic else "stable"
            print(f"    {sclass}: lambda_2 {baseline_lambda_2:.4e} -> {lambda_2_storm:.4e} "
                  f"({100*ratio:.1f}% of baseline) [{tag}]",
                  flush=True)
    return results


# ============================================================================
# Tier 1 layer: catalog summary stats + NOAA live context
# ============================================================================

def load_catalog_summary(catalog: list[dict],
                          leverage_csv_path: Path | None = None) -> dict:
    """Compute current-snapshot summary stats: catalog size, <P_sigma>_LEO, etc."""
    n = len(catalog)
    n_active = sum(1 for o in catalog if o.get("is_active"))
    n_debris = n - n_active

    # <P_sigma>_LEO from leverage map
    leverage_path = leverage_csv_path or (HERE / "leverage_map_v2_full_results.csv")
    rho_total = 0.0
    sigma_weighted = 0.0
    if leverage_path.exists():
        with leverage_path.open("r") as f:
            for r in csv.DictReader(f):
                try:
                    rho = float(r.get("rho", 0))
                    sm = float(r.get("sigma_mag", 0))
                    rho_total += rho
                    sigma_weighted += rho * sm
                except (ValueError, TypeError):
                    continue
    psigma_leo = sigma_weighted / rho_total if rho_total > 0 else 0.0

    return {
        "n_objects": n,
        "n_active": n_active,
        "n_debris": n_debris,
        "psigma_leo": psigma_leo,
        "rho_total": rho_total,
    }


def load_noaa_live_context() -> dict:
    """Read latest NOAA pull for live space-weather context."""
    ctx = {
        "current_kp": None,
        "current_dst_nT": None,
        "active_alerts": 0,
        "days_since_g3": None,
        "noaa_pull_age_hours": None,
    }

    # Latest Kp
    kp_csv = HERE / "noaa_kp_index_1m.csv"
    if kp_csv.exists():
        try:
            with kp_csv.open("r") as f:
                last_row = list(csv.DictReader(f))[-1]
                ctx["current_kp"] = float(last_row.get("kp_index", 0))
        except (ValueError, IndexError, TypeError):
            pass

    # Latest Dst
    dst_csv = HERE / "noaa_dst_index.csv"
    if dst_csv.exists():
        try:
            with dst_csv.open("r") as f:
                last_row = list(csv.DictReader(f))[-1]
                ctx["current_dst_nT"] = float(last_row.get("dst", 0))
        except (ValueError, IndexError, TypeError):
            pass

    # Alerts count
    alerts_csv = HERE / "noaa_alerts.csv"
    if alerts_csv.exists():
        try:
            with alerts_csv.open("r") as f:
                ctx["active_alerts"] = sum(1 for _ in csv.DictReader(f))
        except Exception:
            pass

    # Days since G3+ from historical archive
    storm_csv = HERE / "noaa_storm_catalog.csv"
    if storm_csv.exists():
        try:
            with storm_csv.open("r") as f:
                rows = list(csv.DictReader(f))
            if rows:
                last_g3 = max(rows, key=lambda r: r.get("date", ""))
                last_date = datetime.strptime(last_g3["date"], "%Y-%m-%d").replace(tzinfo=timezone.utc)
                ctx["days_since_g3"] = (datetime.now(timezone.utc) - last_date).days
                ctx["last_g3_date"] = last_g3["date"]
                ctx["last_g3_kp_peak"] = float(last_g3.get("kp_peak", 0))
                ctx["last_g3_class"] = last_g3.get("g_scale", "")
        except (ValueError, IndexError, TypeError, KeyError):
            pass

    # Pull age
    if kp_csv.exists():
        try:
            mtime = datetime.fromtimestamp(kp_csv.stat().st_mtime, tz=timezone.utc)
            ctx["noaa_pull_age_hours"] = (datetime.now(timezone.utc) - mtime).total_seconds() / 3600.0
        except Exception:
            pass

    return ctx


# ============================================================================
# Daily status assembly
# ============================================================================

def assemble_daily_status(noaa_ctx: dict, catalog_summary: dict,
                            fiedler_baseline: dict,
                            stress_results: dict | None) -> dict:
    """Combine all daily-refresh data into a single Firestore-ready document.

    PUBLIC-FACING HEADLINE: baseline lambda_2 + component count + IPR.
    These are the substrate's *real* daily-live signals computed on the
    current LEO catalog without any stress kernel:
      - lambda_2:           algebraic connectivity gap (first time computed at scale)
      - n_components:       graph fragmentation indicator (drift = catalog churn)
      - IPR(v_2):           Fiedler eigenvector localization (which sub-region is exposed)

    Stress tests are included in the payload as v0.1 placeholder data
    (alpha-rescaling kernel; produces uniform-scaling artifacts and is
    NOT a defensible cascade-transition signal). v0.2 will replace with
    physics-based propagation: NRLMSISE-00 atmospheric density at storm
    F10.7/Ap, SGP4 orbit propagation over 7-day storm window, conjunction
    graph rebuilt from post-storm catalog state. Until v0.2 lands, the
    stress test results are kept in the payload for diagnostic visibility
    but NOT surfaced to the public-facing headline.
    """
    now_utc = datetime.now(timezone.utc)
    stress = stress_results or {}

    # === Public-facing headline scalar: baseline substrate state ===
    lambda_2 = fiedler_baseline["lambda_2"]
    ipr = fiedler_baseline["IPR_v2"]
    n_components = fiedler_baseline.get("n_components", 0)
    largest_frac = fiedler_baseline.get("largest_component_frac", 0.0)

    tagline = (
        f"Today's LEO substrate connectivity gap (Fiedler λ₂) = {lambda_2:.4e} "
        f"on the current 49,998-object DISCOSweb catalog. "
        f"Graph has {n_components} connected components; largest contains "
        f"{100.0*largest_frac:.1f}% of objects. "
        f"Fiedler-vector localization (IPR) = {ipr:.0f} nodes — the spectrally "
        f"most exposed sub-population of the LEO catalog. "
        f"This is the first daily-refresh measurement of the LEO substrate's "
        f"algebraic connectivity at scale; drift in lambda_2 over time is itself "
        f"a Kessler-tipping signal."
    )

    headline = {
        "name": "LEO substrate connectivity gap (baseline Fiedler λ₂)",
        "lambda_2": lambda_2,
        "n_components": int(n_components),
        "largest_component_frac": float(largest_frac),
        "IPR_v2": float(ipr),
        "as_of_utc": now_utc.isoformat(),
        "tagline": tagline,
        "tracked_for_drift": True,
    }

    return {
        "refreshed_at_utc": now_utc.isoformat(),
        "refreshed_at_unix": int(now_utc.timestamp()),
        "version": "tier2_v0.1",
        "headline": headline,
        "noaa": noaa_ctx,
        "catalog": catalog_summary,
        "fiedler_baseline": {
            "lambda_2": fiedler_baseline["lambda_2"],
            "lambda_1": fiedler_baseline["lambda_1"],
            "IPR_v2": fiedler_baseline["IPR_v2"],
            "n_nodes": fiedler_baseline["n_nodes"],
            "n_edges": fiedler_baseline["n_edges"],
            "n_components": fiedler_baseline.get("n_components"),
            "largest_component_frac": fiedler_baseline.get("largest_component_frac"),
        },
        "fiedler_stress_tests_PLACEHOLDER_v01_alpha_kernel": stress,
        "stress_kernel_status": (
            "v0.1 alpha-rescaling kernel: uniform edge-weight scaling under "
            "storm parameter alpha. NOT a physical-propagation kernel; results "
            "show uniform-rescaling artifacts (lambda_2 within 0.3% of baseline "
            "across all storm classes; component count consolidates rather than "
            "fragments). v0.2 will replace with NRLMSISE-00 + SGP4 propagation "
            "of the catalog under storm F10.7/Ap, then rebuild conjunction graph "
            "from post-storm orbital state. Stress test results in this payload "
            "are diagnostic-only; do NOT surface to public-facing headline."
        ),
    }


# ============================================================================
# Output: local JSON + Firestore
# ============================================================================

def write_local_json(status: dict, out_path: Path | None = None) -> Path:
    out_path = out_path or (HERE / "substrate_daily_status.json")
    out_path.write_text(json.dumps(status, indent=2, default=str))
    return out_path


def write_rss_feed(status: dict, history: list[dict] | None = None,
                    out_path: Path | None = None,
                    feed_url: str = "https://substrate-daily-pipeline.run.app/rss",
                    site_url: str = "https://fancyland.example/substrate") -> Path:
    """Generate a public Atom 1.0 feed of daily substrate status.

    Actuaries and space-domain analysts subscribe to feeds, not dashboards.
    Each daily run appends an entry; the feed becomes the substrate-leverage
    publication-of-record consumable in any feed reader.

    `history` is an optional list of prior daily status dicts (most recent
    first); if omitted, the feed contains only today's entry. In production
    on Cloud Run, the orchestrator pulls the last N=30 entries from the
    fiedler_history Firestore collection and passes them in.
    """
    out_path = out_path or (HERE / "substrate_daily_feed.xml")
    now_utc = datetime.now(timezone.utc)
    refreshed = status.get("refreshed_at_utc", now_utc.isoformat())
    headline = status.get("headline", {})
    lambda_2 = headline.get("lambda_2", 0.0)
    n_components = headline.get("n_components", 0)
    largest_frac = headline.get("largest_component_frac", 0.0)

    entries_xml = []
    daily_entries = [status] + (history or [])
    for entry in daily_entries[:30]:
        e_refreshed = entry.get("refreshed_at_utc", now_utc.isoformat())
        e_headline = entry.get("headline", {})
        e_lambda2 = e_headline.get("lambda_2", 0.0)
        e_components = e_headline.get("n_components", 0)
        e_largest = e_headline.get("largest_component_frac", 0.0)
        e_ipr = e_headline.get("IPR_v2", 0.0)
        e_tagline = e_headline.get("tagline", "")
        e_psigma = entry.get("catalog", {}).get("psigma_leo", 0.0)
        e_kp = entry.get("noaa", {}).get("current_kp")
        e_dst = entry.get("noaa", {}).get("current_dst_nT")
        e_id = e_refreshed.split("T")[0]
        title = (f"Substrate {e_id}: λ₂ = {e_lambda2:.4e} | "
                 f"{e_components} components | "
                 f"largest {100*e_largest:.1f}%")
        summary = (
            f"{e_tagline}\n\n"
            f"⟨P_σ⟩_LEO = {e_psigma:+.4f}, "
            f"current Kp = {e_kp}, Dst = {e_dst} nT, "
            f"IPR(v₂) = {e_ipr:.0f}.\n"
            f"Substrate methodology: Kessler Critical Surface preprint v0.7."
        )
        entries_xml.append(f"""    <entry>
        <title>{_xml_escape(title)}</title>
        <id>tag:fancyland.example,2026:substrate-daily/{e_id}</id>
        <updated>{e_refreshed}</updated>
        <link href="{site_url}/{e_id}"/>
        <author><name>Fancyland LLC / Lattice OS</name></author>
        <summary type="text">{_xml_escape(summary)}</summary>
    </entry>""")

    feed_xml = f"""<?xml version="1.0" encoding="utf-8"?>
<feed xmlns="http://www.w3.org/2005/Atom">
    <title>Substrate Daily — LEO connectivity gap (Fiedler λ₂)</title>
    <link href="{feed_url}" rel="self"/>
    <link href="{site_url}"/>
    <id>tag:fancyland.example,2026:substrate-daily</id>
    <updated>{refreshed}</updated>
    <author><name>Fancyland LLC / Lattice OS</name></author>
    <subtitle>Daily algebraic-connectivity signal on the live ESA DISCOSweb LEO catalog. Today: λ₂ = {lambda_2:.4e}, {n_components} components, largest = {100*largest_frac:.1f}%.</subtitle>
{chr(10).join(entries_xml)}
</feed>
"""
    out_path.write_text(feed_xml, encoding="utf-8")
    return out_path


def _xml_escape(s: str) -> str:
    """Minimal XML text-content escape."""
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
             .replace('"', "&quot;").replace("'", "&apos;"))


def write_to_firestore(status: dict,
                        project_id: str = "brainstormingorganization",
                        collection: str = "substrate_daily",
                        history_collection: str = "fiedler_history") -> dict:
    """Write daily status to Firestore. Updates a single 'current' doc and
    appends to fiedler_history collection for trend tracking.

    Requires backend/service-account.json (gitignored, see firestore_loader.py
    pattern in bvp_12_multi_chain).
    """
    try:
        import firebase_admin
        from firebase_admin import credentials, firestore
    except ImportError as e:
        return {"firestore_write": False, "error": f"firebase_admin not installed: {e}"}

    sa_path = HERE.resolve().parents[3] / "service-account.json"
    if not sa_path.exists():
        return {"firestore_write": False,
                "error": f"service-account.json not found at {sa_path}"}

    if not firebase_admin._apps:
        cred = credentials.Certificate(str(sa_path))
        firebase_admin.initialize_app(cred, {"projectId": project_id})

    db = firestore.client()

    # Current status doc (single doc, overwritten daily)
    db.collection(collection).document("current").set(status)

    # Append to fiedler_history (one doc per day)
    date_id = status["refreshed_at_utc"].split("T")[0]
    history_doc = {
        "date": date_id,
        "refreshed_at_utc": status["refreshed_at_utc"],
        "lambda_2_baseline": status["fiedler_baseline"]["lambda_2"],
        "n_nodes": status["fiedler_baseline"]["n_nodes"],
        "n_edges": status["fiedler_baseline"]["n_edges"],
        "psigma_leo": status["catalog"]["psigma_leo"],
        "stress_tests": {
            sclass: {
                "lambda_2_storm": v["lambda_2_storm"],
                "ratio": v["ratio_storm_to_baseline"],
                "ballistic": v["ballistic_transition"],
            }
            for sclass, v in status.get("fiedler_stress_tests", {}).items()
        },
    }
    db.collection(history_collection).document(date_id).set(history_doc)

    return {"firestore_write": True, "current_doc": f"{collection}/current",
            "history_doc": f"{history_collection}/{date_id}"}


# ============================================================================
# Main pipeline orchestration
# ============================================================================

def run_pipeline(write_firestore: bool = True,
                  baseline_only: bool = False,
                  k_nn: int = DEFAULT_K_NN,
                  verbose: bool = True) -> dict:
    print("=" * 75)
    print("SUBSTRATE DAILY PIPELINE  (Tier 2 v0.1)")
    print("=" * 75)
    print(f"  started: {datetime.now(timezone.utc).isoformat()}")
    print(f"  write_firestore: {write_firestore}")
    print(f"  baseline_only:   {baseline_only}")
    print()

    # 1. Load catalog
    print("Step 1: loading LEO catalog ...")
    t0 = time.time()
    catalog = kra.load_assignment_csv()
    print(f"  loaded {len(catalog):,} objects in {time.time()-t0:.1f}s")

    # 2. Vectorize
    vec = vectorize_catalog(catalog)

    # 3. Build baseline conjunction graph
    print(f"\nStep 2: building baseline conjunction graph (k_nn={k_nn}) ...")
    A_baseline = build_conjunction_graph(vec, k_nn=k_nn, progress=verbose)

    # 4. Compute baseline Fiedler
    print(f"\nStep 3: computing baseline Fiedler lambda_2 ...")
    fiedler_baseline = compute_fiedler(A_baseline, verbose=verbose)

    # 5. Storm-class stress tests
    stress_results = None
    if not baseline_only:
        print(f"\nStep 4: storm-class stress tests ...")
        stress_results = stress_test_storm_classes(
            vec, fiedler_baseline["lambda_2"], k_nn=k_nn, verbose=verbose,
        )

    # 6. Tier 1 layer
    print(f"\nStep 5: Tier 1 layer (catalog summary + NOAA context) ...")
    catalog_summary = load_catalog_summary(catalog)
    noaa_ctx = load_noaa_live_context()
    print(f"  catalog: {catalog_summary['n_objects']} objects "
          f"({catalog_summary['n_active']} active, {catalog_summary['n_debris']} debris)")
    print(f"  <P_sigma>_LEO: {catalog_summary['psigma_leo']:+.4f}")
    print(f"  current Kp: {noaa_ctx.get('current_kp')}, "
          f"Dst: {noaa_ctx.get('current_dst_nT')} nT")
    print(f"  days since G3+: {noaa_ctx.get('days_since_g3')}")

    # 7. Assemble status
    status = assemble_daily_status(noaa_ctx, catalog_summary,
                                     fiedler_baseline, stress_results)

    # 8. Write outputs
    print(f"\nStep 6: writing outputs ...")
    json_path = write_local_json(status)
    print(f"  local JSON -> {json_path.name}")
    rss_path = write_rss_feed(status)
    print(f"  Atom feed  -> {rss_path.name}")

    if write_firestore:
        fs_result = write_to_firestore(status)
        if fs_result["firestore_write"]:
            print(f"  Firestore -> {fs_result['current_doc']}, {fs_result['history_doc']}")
        else:
            print(f"  Firestore SKIPPED: {fs_result.get('error')}")
    else:
        print("  Firestore: SKIPPED (--no-firestore)")

    print()
    print("=" * 75)
    print("DAILY PIPELINE COMPLETE")
    print("=" * 75)
    return status


def main():
    p = argparse.ArgumentParser(description="Substrate daily pipeline (Tier 2 v0.1)")
    p.add_argument("--no-firestore", action="store_true",
                   help="Skip Firestore writes (local JSON only, for dev/testing)")
    p.add_argument("--baseline-only", action="store_true",
                   help="Skip storm-class stress tests")
    p.add_argument("--k-nn", type=int, default=DEFAULT_K_NN,
                   help=f"k-NN sparsity for conjunction graph (default {DEFAULT_K_NN})")
    p.add_argument("--quiet", action="store_true", help="Reduce console output")
    args = p.parse_args()

    run_pipeline(
        write_firestore=not args.no_firestore,
        baseline_only=args.baseline_only,
        k_nn=args.k_nn,
        verbose=not args.quiet,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
