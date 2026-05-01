"""
cascade_trigger_proximity.py -- substrate-native proximity, top-K named-target
ranking, and local counterfactual delta for the slot Layer-2 moat.

Per locked decision #19 (BUILD_PLAN 2026-04-28):

  Customer's slot report includes top-K (K <= 5) named cataloged objects
  ranked by substrate-derived L_signed = rho * sigma_mag, filtered by
  conjunction-graph (D_mom) adjacency to the customer's would-be orbit,
  with local counterfactual delta. Substrate IP exposure mitigated by
  rate-limit on named-target queries; full leverage map stays in
  substrate_internal/.

The proximity metric is conjunction-graph adjacency, NOT Euclidean
(inc, alt) distance. Per Claude.AI 2026-04-28: "two satellites can be at
the same altitude and inclination but at completely different RAAN and
never come close to each other in operational reality." D_mom uses the
same momentum-aware k-NN edges already in leverage_map_v2_full.py
build_momentum_graph -- so cascade-trigger proximity in D_mom space is
the metric the substrate algebra natively respects.

Public API:

  compute_query_to_catalog_dmom(...)
      D_mom edge weight from a query orbit to every catalog object.

  rank_triggers_by_proximity(...)
      Top-K named cascade triggers within proximity, ranked by L_signed.

  counterfactual_local_delta(...)
      In-shell SLS recompute with one trigger excluded.

The full leverage map is required as input; the customer-facing surface
returns only top-K names + L_signed + delta. Mosaic-theory enumeration
attacks are mitigated by per-user rate-limit on these queries (enforced
at the backend, not in this module).
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np

EARTH_R = 6378.137
MU_EARTH = 398600.4418


# ---------------------------------------------------------------------------
# Conjunction-graph (D_mom) adjacency: query orbit -> all catalog objects
# ---------------------------------------------------------------------------

def compute_query_to_catalog_dmom(
    query_inc_deg: float,
    query_peri_km: float,
    query_apo_km: float,
    catalog_inc_deg: np.ndarray,
    catalog_peri_km: np.ndarray,
    catalog_apo_km: np.ndarray,
    catalog_xsect: np.ndarray | None = None,
    inc_tol_deg: float = 15.0,
    min_overlap_km: float = 5.0,
) -> np.ndarray:
    """Return D_mom edge weights from a single query orbit to every catalog row.

    The edge construction mirrors leverage_map_v2_full.build_momentum_graph
    so the same proximity definition the substrate engine uses for the
    leverage map applies to customer-query proximity. Edges are zero where:
      - altitude shells do not overlap (no intersection of [peri, apo]),
      - inclination differs by more than inc_tol_deg (with retrograde
        equivalence: |180 - inc_q - inc_i|),
      - the overlap is below min_overlap_km.

    Where edges are nonzero, weight is base * mom where:
      base = (overlap / max(apo)) * sqrt(xsect_q * xsect_i)
      mom  = relative-velocity-squared / (vq * vi)

    Returns an (N,) array of edge weights (zero where no edge exists).
    """
    n = len(catalog_inc_deg)
    if catalog_xsect is None:
        catalog_xsect = np.ones(n)

    # Catalog orbital velocity at perigee
    a_cat = 0.5 * (catalog_peri_km + catalog_apo_km) + EARTH_R
    r_peri_cat = catalog_peri_km + EARTH_R
    v_sq_cat = MU_EARTH * (2.0 / r_peri_cat - 1.0 / a_cat)
    v_cat = np.sqrt(np.maximum(v_sq_cat, 0.0))

    # Query orbital velocity at perigee
    a_q = 0.5 * (query_peri_km + query_apo_km) + EARTH_R
    r_peri_q = query_peri_km + EARTH_R
    v_sq_q = MU_EARTH * (2.0 / r_peri_q - 1.0 / a_q)
    v_q = math.sqrt(max(v_sq_q, 0.0))

    inc_q_rad = math.radians(query_inc_deg)
    inc_cat_rad = np.deg2rad(catalog_inc_deg)

    # Inclination compatibility (with retrograde-equivalent check)
    d_inc = np.abs(catalog_inc_deg - query_inc_deg)
    d_inc_retro = np.abs(180.0 - catalog_inc_deg - query_inc_deg)
    inc_ok = np.minimum(d_inc, d_inc_retro) <= inc_tol_deg

    # Altitude shell overlap
    lo = np.maximum(catalog_peri_km, query_peri_km)
    hi = np.minimum(catalog_apo_km, query_apo_km)
    overlap = hi - lo
    overlap_ok = overlap >= min_overlap_km

    # Base edge weight
    apo_max = np.maximum(catalog_apo_km, query_apo_km)
    apo_safe = np.where(apo_max > 0, apo_max, 1.0)
    xsect_q = 1.0  # no per-query x-section -- treat as 1
    base = (overlap / apo_safe) * np.sqrt(xsect_q * catalog_xsect)

    # Momentum factor
    cos_th = np.cos(inc_q_rad - inc_cat_rad)
    v_prod = v_q * v_cat
    v_prod_safe = np.where(v_prod > 0, v_prod, 1.0)
    v_rel_sq = v_q * v_q + v_cat * v_cat - 2.0 * v_q * v_cat * cos_th
    mom = np.where(v_prod > 0, v_rel_sq / v_prod_safe, 1.0)

    edge = base * mom
    edge = np.where(inc_ok & overlap_ok, edge, 0.0)
    edge = np.where(edge > 0, edge, 0.0)
    return edge


# ---------------------------------------------------------------------------
# Top-K named cascade triggers within proximity
# ---------------------------------------------------------------------------

def rank_triggers_by_proximity(
    query_inc_deg: float,
    query_peri_km: float,
    query_apo_km: float,
    catalog: list[dict],
    leverage_map: list[dict],
    k: int = 5,
    pool_size: int = 50,
    min_dmom: float = 0.0,
) -> list[dict[str, Any]]:
    """Return top-K named cascade triggers within D_mom proximity.

    Args:
      query_*: customer's proposed orbit
      catalog: full per-object list (used to compute D_mom; needs inc_deg,
               perigee_alt_km, apogee_alt_km, name, cosparId, sigma_mag,
               rho, L_signed)
      leverage_map: top-N highest-L_signed objects (the cascade-trigger
                    pool). May equal catalog when scope is unrestricted;
                    in production this is pre-pruned to ~50 per shell at
                    compile time.
      k: number of named triggers to return (default 5).
      pool_size: how many of the leverage_map's top entries to consider
                 before D_mom filtering (default 50). Leftover IP exposure
                 surface; default value chosen to satisfy K <= 5 ceiling
                 with reasonable D_mom filter pass-through.
      min_dmom: minimum D_mom edge weight for a trigger to qualify
                (zero allows any nonzero edge; production may set this to
                a calibrated floor).

    Returns:
      [{
        "name": str, "cosparId": str,
        "inc_deg": float, "peri_km": float, "apo_km": float,
        "sigma_mag": float, "rho": float, "L_signed": float,
        "dmom_to_query": float,
        "rank_in_proximity": int,  # 1..K
      }, ...]
    """
    pool = leverage_map[:pool_size] if pool_size else leverage_map

    incs = np.array([float(r["inc_deg"]) for r in pool])
    peris = np.array([float(r["perigee_alt_km"]) for r in pool])
    apos = np.array([float(r["apogee_alt_km"]) for r in pool])

    edges = compute_query_to_catalog_dmom(
        query_inc_deg, query_peri_km, query_apo_km,
        incs, peris, apos,
    )

    candidates: list[tuple[int, float]] = []
    for i, e in enumerate(edges):
        if e > min_dmom:
            candidates.append((i, float(e)))

    # Within proximity, rank by abs(L_signed) descending
    candidates.sort(
        key=lambda ie: abs(float(pool[ie[0]].get("L_signed", 0.0))),
        reverse=True,
    )

    out: list[dict[str, Any]] = []
    for rank, (idx, e) in enumerate(candidates[:k], start=1):
        r = pool[idx]
        out.append({
            "name": r.get("name", ""),
            "cosparId": r.get("cosparId", ""),
            "inc_deg": float(r["inc_deg"]),
            "peri_km": float(r["perigee_alt_km"]),
            "apo_km": float(r["apogee_alt_km"]),
            "sigma_mag": float(r.get("sigma_mag", 0.0)),
            "rho": float(r.get("rho", 0.0)),
            "L_signed": float(r.get("L_signed", 0.0)),
            "dmom_to_query": e,
            "rank_in_proximity": rank,
        })
    return out


# ---------------------------------------------------------------------------
# Local counterfactual delta: in-shell SLS recompute with one trigger removed
# ---------------------------------------------------------------------------

def counterfactual_local_delta(
    customer_inc_deg: float,
    customer_peri_km: float,
    customer_apo_km: float,
    customer_shell: str,
    trigger_sigma_mag: float,
    trigger_in_same_shell: bool,
    shell_sigma_mag_mean: float,
    shell_n_samples: int,
    rho_at_query: float,
    rho_percentile_at_query: float,
) -> dict[str, float]:
    """Compute the customer's SLS shift if a single trigger is removed.

    The dominant effect on the customer's SLS from removing one trigger
    is the in-shell sigma_mag mean shift (when the trigger is in the
    customer's shell). For sparse shells (N small), one removal is a
    meaningful fraction of the shell mean. For dense shells (Starlink,
    SSO), one removal is ~1/N which is below 1 SLS point in absolute
    terms, but the trigger still contributes structurally to the system
    cascade leverage and is worth surfacing.

    Returns:
      {
        "sls_substrate_old": int,
        "sls_substrate_new": int,
        "sls_total_old":     int,
        "sls_total_new":     int,
        "sls_delta":         int,   # new - old (positive = improvement)
      }
    """
    # SLS_substrate before removal
    sls_substrate_old = round(100.0 * (1.0 - abs(shell_sigma_mag_mean)))
    # SLS_exposure unchanged at this approximation level (trigger removal
    # affects rho only if the trigger is in the customer's exact (inc, alt)
    # cell; modeled as a separate path below if that becomes relevant).
    sls_exposure = round(100.0 * (1.0 - rho_percentile_at_query / 100.0))
    sls_total_old = round(0.5 * (sls_substrate_old + sls_exposure))

    if trigger_in_same_shell and shell_n_samples >= 2:
        # New shell mean with this trigger excluded
        new_shell_mean = (
            (shell_n_samples * shell_sigma_mag_mean - trigger_sigma_mag)
            / (shell_n_samples - 1)
        )
        sls_substrate_new = round(100.0 * (1.0 - abs(new_shell_mean)))
    else:
        # Trigger is in a different shell: customer's in-shell mean is
        # unaffected (this approximation excludes cross-shell coupling).
        sls_substrate_new = sls_substrate_old

    sls_total_new = round(0.5 * (sls_substrate_new + sls_exposure))
    delta = sls_total_new - sls_total_old

    return {
        "sls_substrate_old": int(sls_substrate_old),
        "sls_substrate_new": int(sls_substrate_new),
        "sls_total_old": int(sls_total_old),
        "sls_total_new": int(sls_total_new),
        "sls_delta": int(delta),
    }
