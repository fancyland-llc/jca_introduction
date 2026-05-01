"""Cauchy / SO(3) tumbling-state test rig — v0.4.1 (orbit-averaged baseline).

Per Tony 2026-04-29 morning compute + Claude.AI iterative review:

The geometric mechanism (Cauchy 1841 surface-area theorem) replaces the
weaker dynamical-rectification mechanism. <A>_isotropic = S/4 for any
convex body tumbling ergodically over SO(3). Flat-spin around the major
(long) axis explores a 1-d orbit in SO(3); chaotic tumble explores a
3-d region. Topologically distinct sectors — rotational analog of
Klein-4 sigma-parity in the substrate.

THE RIGOROUS MAJOR-AXIS BASELINE (Claude.AI 2026-04-29 second pass):
    A long-axis-symmetric body spinning around its long axis with the
    long axis fixed in inertial space, while in circular orbit, sees
    velocity vectors sweeping through theta in [0, 2*pi] each orbit.

    Projected cross-section along velocity:
        A(theta) = A_endon * |cos(theta)| + A_broadside * |sin(theta)|

    Time-averaged over one orbit (uniform on theta):
        <A>_majoraxis = (1/2*pi) * int_0^{2*pi} A(theta) d_theta
                      = (1/2*pi) * 4 * (A_endon + A_broadside)
                      = (2/pi) * (A_endon + A_broadside)

    This is the *physically derived* flat-spin baseline. Yesterday's
    "broadside/2" was a heuristic guess; the arithmetic-mid baseline
    (A_e + A_b)/2 is a different shortcut. The rigorous (2/pi) form
    is what falls out of first principles.

    The earlier alpha=1.46-1.81 numbers were inflated by using a
    different baseline. The rigorous numbers are alpha=1.09-1.20 for
    cylinders, alpha=0.99-1.05 for paneled satellites. The methodology
    is more honest and survives peer review at this rigor.

REAL FINDING — paneled satellites can have alpha < 1:
    For Iridium-1 (panel-broadside-dominant geometry), the major-axis-
    spin orbit-averaged baseline EXCEEDS the Cauchy ergodic-tumble S/4.
    Physical reason: spinning around the long axis with broadside-
    dominant panels sweeps every velocity direction through the panel
    face; Cauchy ergodic dilutes this with end-on-like orientations.
    This is the methodology detecting a real physical effect, not a
    bug. The audit-charter row will name it explicitly.

LEVERAGE WEIGHTING:
    Per Claude.AI 2026-04-29: the headline number is the leverage-
    weighted aggregate across the leverage-dominant subset, not the
    unweighted cohort mean. KH-8 reconnaissance cylinders dominate
    TOP100_LEVERAGE_PORTFOLIO (4 of top-5); their +20% T_state
    contributes ~9 percentage points to the +16% weighted total.

L_signed values pinned to:
    kessler/audit/2026-04-28/TOP100_LEVERAGE_PORTFOLIO.json
    kessler/audit/2026-04-28/SUBSTRATE_COUPLING_PATHWAY_starlink.json
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path


# ============================================================================
# Geometry-class baselines
# ============================================================================

def cauchy_isotropic(surface_area_m2: float) -> float:
    """Cauchy 1841: ergodic-tumble cross-section = S/4 for convex body."""
    return surface_area_m2 / 4.0


def major_axis_spin_orbit_avg(A_endon: float, A_broadside: float) -> float:
    """Orbit-averaged projected area for major-axis-spin in circular orbit.

    A(theta) = A_endon |cos(theta)| + A_broadside |sin(theta)|
    <A> = (1/2*pi) * int_0^{2*pi} A(theta) d_theta
        = (2/pi) * (A_endon + A_broadside)
    """
    return (2.0 / math.pi) * (A_endon + A_broadside)


# ============================================================================
# Geometry primitives
# ============================================================================

@dataclass
class CylinderGeometry:
    name: str
    length_m: float
    diameter_m: float

    @property
    def radius_m(self) -> float:
        return self.diameter_m / 2.0

    @property
    def surface_area(self) -> float:
        return (
            2.0 * math.pi * self.radius_m ** 2
            + 2.0 * math.pi * self.radius_m * self.length_m
        )

    @property
    def A_endon(self) -> float:
        return math.pi * self.radius_m ** 2

    @property
    def A_broadside(self) -> float:
        return self.length_m * self.diameter_m

    @property
    def A_isotropic(self) -> float:
        return cauchy_isotropic(self.surface_area)

    @property
    def A_majoraxis(self) -> float:
        return major_axis_spin_orbit_avg(self.A_endon, self.A_broadside)


@dataclass
class PaneledSatelliteGeometry:
    """Bus + solar-panel composite. Cauchy applied to the convex hull
    (additive surface area is conservative; flagged as approximation)."""
    name: str
    bus_face_m2: float           # bus end-on (smallest stable face)
    panel_face_m2: float         # panel-dominant broadside
    total_surface_m2: float      # additive surface across rigid components

    @property
    def A_endon(self) -> float:
        return self.bus_face_m2

    @property
    def A_broadside(self) -> float:
        return self.panel_face_m2

    @property
    def surface_area(self) -> float:
        return self.total_surface_m2

    @property
    def A_isotropic(self) -> float:
        return cauchy_isotropic(self.total_surface_m2)

    @property
    def A_majoraxis(self) -> float:
        return major_axis_spin_orbit_avg(self.A_endon, self.A_broadside)


# ============================================================================
# T_state observable
# ============================================================================

@dataclass
class TumblingStateResult:
    name: str
    cosparId: str
    geometry_class: str
    convexity_assumption: str
    baseline_derivation: str

    surface_area_m2: float
    A_endon_m2: float
    A_broadside_m2: float
    A_isotropic_m2: float
    A_majoraxis_orbit_avg_m2: float

    alpha_isotropic_over_majoraxis: float
    T_state_pct: float

    L_signed_audit: float | None
    notes: str = ""


def compute_t_state(
    geom: CylinderGeometry | PaneledSatelliteGeometry,
    *,
    cosparId: str,
    L_signed_audit: float | None = None,
    notes: str = "",
) -> TumblingStateResult:
    if isinstance(geom, CylinderGeometry):
        cls = "cylinder"
        convex = "exact"
        derivation = (
            f"Cylinder L={geom.length_m:.2f} m, diameter={geom.diameter_m:.2f} m. "
            f"Cauchy S/4 exact (convex). Major-axis-spin orbit-averaged baseline "
            f"<A>_majoraxis = (2/pi) * (A_endon + A_broadside) = "
            f"{geom.A_majoraxis:.2f} m^2."
        )
    elif isinstance(geom, PaneledSatelliteGeometry):
        cls = "satellite_with_panels"
        convex = "convex_hull_approximation"
        derivation = (
            f"Bus + solar-panel composite. Surface area additive across rigid "
            f"components ({geom.surface_area:.1f} m^2). Cauchy S/4 = "
            f"{geom.A_isotropic:.2f} m^2 is an UPPER BOUND for non-convex "
            f"panel/bus geometry. Major-axis-spin orbit-averaged baseline "
            f"<A>_majoraxis = (2/pi) * (A_endon + A_broadside) = "
            f"{geom.A_majoraxis:.2f} m^2."
        )
    else:
        raise TypeError(f"Unsupported geometry: {type(geom).__name__}")

    alpha = geom.A_isotropic / geom.A_majoraxis
    T_pct = 100.0 * (alpha - 1.0)

    return TumblingStateResult(
        name=geom.name,
        cosparId=cosparId,
        geometry_class=cls,
        convexity_assumption=convex,
        baseline_derivation=derivation,
        surface_area_m2=geom.surface_area,
        A_endon_m2=geom.A_endon,
        A_broadside_m2=geom.A_broadside,
        A_isotropic_m2=geom.A_isotropic,
        A_majoraxis_orbit_avg_m2=geom.A_majoraxis,
        alpha_isotropic_over_majoraxis=alpha,
        T_state_pct=T_pct,
        L_signed_audit=L_signed_audit,
        notes=notes,
    )


# ============================================================================
# Canonical cohort
# ============================================================================

# L_signed values pinned to the audit JSONs at kessler/audit/2026-04-28/.
# KH-8 + Firefly: TOP100_LEVERAGE_PORTFOLIO sample_objects_top10.
# Delta II + Shavit: SUBSTRATE_COUPLING_PATHWAY_starlink top_10_foreign_sinks
# (L_signed_foreign — same as the object's global L_signed since the
# leverage scalar is per-object, not per-pair).
L_SIGNED_AUDIT = {
    "1971-033A": 33384.17,    # KH-8 (05171), TOP100 rank 1
    "2023-202B": 20739.36,    # Firefly Alpha 2nd stage, TOP100 rank 5
    "2007-041B": 10133.29,    # Delta II 2nd stage, Pi top-10 sink rank 2
    "2007-025B": 12269.25,    # Shavit 3rd stage, Pi top-10 sink rank 3
}


def canonical_cohort() -> list[TumblingStateResult]:
    return [
        compute_t_state(
            CylinderGeometry(
                name="KH-8 (Gambit-3) reconnaissance bus",
                length_m=15.0, diameter_m=1.5,
            ),
            cosparId="1971-033A",
            L_signed_audit=L_SIGNED_AUDIT["1971-033A"],
            notes="Top-1 |L_signed| in TOP100_LEVERAGE_PORTFOLIO. KH-8 GAMBIT-3 "
                  "reference dimensions ~15 m x 1.5 m diameter (declassified NRO "
                  "Agena-D-derived bus). 4 of top-5 leverage objects are KH-8 "
                  "spysats; this geometry drives the leverage-weighted aggregate.",
        ),
        compute_t_state(
            CylinderGeometry(
                name="Firefly Alpha 2nd stage",
                length_m=6.0, diameter_m=1.8,
            ),
            cosparId="2023-202B",
            L_signed_audit=L_SIGNED_AUDIT["2023-202B"],
            notes="Top-1 foreign sink in Pi-Class coupling pathway "
                  "(SUBSTRATE_COUPLING_PATHWAY_starlink.json). Reference "
                  "geometry from Tony 2026-04-29 compute.",
        ),
        compute_t_state(
            CylinderGeometry(
                name="Delta II 2nd stage (AJ10-118K)",
                length_m=6.0, diameter_m=2.4,
            ),
            cosparId="2007-041B",
            L_signed_audit=L_SIGNED_AUDIT["2007-041B"],
            notes="Aerojet AJ10-118K Delta-II upper stage. "
                  "Rank 2 in Pi-Class coupling pathway top-10 sinks.",
        ),
        compute_t_state(
            CylinderGeometry(
                name="Shavit 3rd stage (RSA-3)",
                length_m=2.0, diameter_m=1.25,
            ),
            cosparId="2007-025B",
            L_signed_audit=L_SIGNED_AUDIT["2007-025B"],
            notes="Israeli Shavit/RSA-3 third stage solid motor. "
                  "5 of top-10 Pi-Class sinks are Shavit-class.",
        ),
        compute_t_state(
            PaneledSatelliteGeometry(
                name="Envisat",
                bus_face_m2=21.16,         # bus end-on
                panel_face_m2=116.0,       # panel-dominant broadside
                total_surface_m2=366.32,   # additive surface
            ),
            cosparId="2002-009A",
            L_signed_audit=None,
            notes="ESA Earth-observation, 8 metric tons, decommissioned 2012. "
                  "Documented tumbling. NOT in TOP100 leverage cohort; included "
                  "as informative outlier for paneled-satellite geometry class.",
        ),
        compute_t_state(
            PaneledSatelliteGeometry(
                name="Iridium-1 (decommissioned, LM-700 bus)",
                bus_face_m2=3.2,
                panel_face_m2=13.4,
                total_surface_m2=41.84,
            ),
            cosparId="IRIDIUM-1-REF",
            L_signed_audit=None,
            notes="Original Iridium constellation member, decommissioned. "
                  "Lockheed Martin LM-700 bus. Reference geometry only. "
                  "NEGATIVE T_state expected: panel-broadside-dominant geometry "
                  "puts orbit-averaged major-axis spin slightly above Cauchy-ergodic.",
        ),
    ]


# ============================================================================
# Leverage-weighted aggregate
# ============================================================================

def leverage_weighted_aggregate(results: list[TumblingStateResult]) -> dict:
    """sum(L_signed * T_state) / sum(L_signed) over the L-cohort subset."""
    weighted = [r for r in results if r.L_signed_audit is not None]
    if not weighted:
        return {"n_in_aggregate": 0}

    total_L = sum(r.L_signed_audit for r in weighted)  # type: ignore
    weighted_T = sum(r.L_signed_audit * r.T_state_pct for r in weighted) / total_L  # type: ignore
    weighted_alpha = 1.0 + weighted_T / 100.0

    contributions = []
    for r in weighted:
        weight = r.L_signed_audit / total_L  # type: ignore
        contrib_pct = r.L_signed_audit * r.T_state_pct / total_L  # type: ignore
        contributions.append({
            "name": r.name,
            "cosparId": r.cosparId,
            "L_signed": r.L_signed_audit,
            "T_state_pct": r.T_state_pct,
            "weight_in_aggregate": weight,
            "contribution_pct": contrib_pct,
        })

    return {
        "n_in_aggregate": len(weighted),
        "L_signed_total": total_L,
        "leverage_weighted_T_state_pct": weighted_T,
        "leverage_weighted_alpha": weighted_alpha,
        "per_object_contributions": contributions,
    }


# ============================================================================
# Main
# ============================================================================

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).parent / "cauchy_so3_results.json",
    )
    args = parser.parse_args()

    results = canonical_cohort()

    print()
    print("=" * 96)
    print(" " * 24 + "Cauchy/SO(3) cohort — orbit-averaged major-axis baseline (v0.4.1)")
    print("=" * 96)
    hdr = (
        f"{'Object':<42} {'class':<24} {'A_iso':>7} {'A_maj':>7} "
        f"{'alpha':>6} {'T%':>6} {'L_sgn':>9}"
    )
    print(hdr)
    print("-" * 96)
    for r in results:
        Lstr = f"{r.L_signed_audit:>9.0f}" if r.L_signed_audit else "       --"
        print(
            f"{r.name[:42]:<42} {r.geometry_class[:23]:<24} "
            f"{r.A_isotropic_m2:>7.2f} {r.A_majoraxis_orbit_avg_m2:>7.2f} "
            f"{r.alpha_isotropic_over_majoraxis:>6.3f} "
            f"{r.T_state_pct:>+5.1f} {Lstr}"
        )
    print("-" * 96)

    alphas = [r.alpha_isotropic_over_majoraxis for r in results]
    T_pcts = [r.T_state_pct for r in results]
    print()
    print(f"alpha (orbit-averaged major-axis baseline) cohort:")
    print(f"  min={min(alphas):.3f}  median={sorted(alphas)[len(alphas)//2]:.3f}  max={max(alphas):.3f}")
    print(f"T_state %% cohort:")
    print(f"  min={min(T_pcts):+.1f}%  median={sorted(T_pcts)[len(T_pcts)//2]:+.1f}%  max={max(T_pcts):+.1f}%")

    aggregate = leverage_weighted_aggregate(results)
    print()
    print(f"Leverage-weighted aggregate (L-cohort subset, n={aggregate['n_in_aggregate']}):")
    print(f"  Sum |L_signed|:                {aggregate['L_signed_total']:>12,.0f}")
    print(f"  Weighted T_state:              {aggregate['leverage_weighted_T_state_pct']:>+12.2f}%")
    print(f"  Weighted alpha:                {aggregate['leverage_weighted_alpha']:>12.3f}")
    print()
    print(f"  Per-object contribution:")
    for c in aggregate['per_object_contributions']:
        print(f"    {c['name'][:42]:<42}  L={c['L_signed']:>10,.0f}  "
              f"T={c['T_state_pct']:>+5.1f}%  w={c['weight_in_aggregate']:.3f}  "
              f"contrib={c['contribution_pct']:>+5.2f}%")

    out = {
        "claim_id": "T_STATE_COHORT_v0.4.1",
        "claim_name": (
            "Cauchy/SO(3) tumbling-state observable, orbit-averaged major-axis "
            "baseline, leverage-weighted aggregate across leverage-dominant cohort."
        ),
        "methodology_reference": (
            "Cauchy 1841 surface-area theorem on convex bodies tumbling "
            "ergodically over SO(3): <A>_isotropic = S/4. Major-axis-spin "
            "orbit-averaged baseline: <A>_majoraxis = (2/pi)*(A_endon + "
            "A_broadside) for circular-orbit line-of-sight precession. "
            "Topologically distinct SO(3) orbit classes: 1-d major-axis-spin "
            "manifold vs 3-d ergodic-tumble manifold. Rotational analog of "
            "Klein-4 sigma-parity in the substrate."
        ),
        "computation_type": "closed_form_geometric",
        "is_monte_carlo": False,
        "n_objects": len(results),
        "n_in_leverage_cohort": aggregate.get("n_in_aggregate", 0),
        "alpha_cohort_stats": {
            "min": min(alphas),
            "max": max(alphas),
            "median": sorted(alphas)[len(alphas) // 2],
            "mean": sum(alphas) / len(alphas),
        },
        "T_state_pct_cohort_stats": {
            "min_pct": min(T_pcts),
            "max_pct": max(T_pcts),
            "median_pct": sorted(T_pcts)[len(T_pcts) // 2],
        },
        "leverage_weighted_aggregate": aggregate,
        "results": [asdict(r) for r in results],
    }
    args.output.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print()
    print(f"Wrote {args.output}  ({args.output.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
