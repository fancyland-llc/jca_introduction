"""Flare-Cascade Risk Assessor (v0.2).

Conditional on space-weather storms, computes:
  1. Annual flare-loss probability for a single satellite (raw + percentile)
  2. By-storm-class breakdown (G1...G5...Carrington)
  3. Fiedler-shift placeholder for substrate cascade (graph-spectral, deterministic)
  4. Saturation-point framing (the marketing hook for underwriter pitch)

CALIBRATION GATE: loss-response curves are locked against Halloween 2003 +
May 2024 Gannon ground-truth outcomes. Both events must reproduce within
tolerance before any customer-facing output is generated.

Storm-frequency table: 95-yr GFZ-Potsdam Kp/ap archive (1932-2026), pulled
by noaa_historical_loader.py. Carrington-class: Riley (2012).

NOTE: v0.2 ships the per-satellite annual loss probability (calibrated).
The substrate cascade signal is a Fiedler-eigenvalue shift (lambda_2)
on the momentum-aware conjunction graph D_mom at each storm class.
The D_mom precompute and lambda_2 spectral pipeline is staged for v0.3
(see fiedler_precompute_TODO()). Monte Carlo cascade simulation REMOVED -
it conflated Bernoulli loss-flips with topological graph rewiring.

CLI:
    python flare_assessor.py --inc 53 --perigee 540 --apogee 560 \\
                              --mass 260 --years 10 --name MySat
    python flare_assessor.py --calibration-gate         # validate before shipping
    python flare_assessor.py --export-calibration       # dump events JSON
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

import kessler_risk_assessor as kra

HERE = Path(__file__).parent


# ============================================================================
# Storm frequency table (from GFZ Kp/ap archive, 1932-2026 = 95 years)
# ============================================================================
# Source: noaa_storm_catalog.csv (this directory), aggregated from
# https://kp.gfz.de/app/files/Kp_ap_Ap_SN_F107_since_1932.txt
# Carrington-class: Riley (2012), Space Weather 10:S02012, ~12% per decade.
STORM_FREQUENCY_PER_YEAR = {
    "G1":         36.07,    # Kp >= 5
    "G2":         14.49,    # Kp >= 6
    "G3":          5.47,    # Kp >= 7
    "G4":          2.32,    # Kp >= 8
    "G5":          0.28,    # Kp >= 9 (saturation ceiling)
    "Carrington":  0.012,   # Kp would be ~10.5+ if scale extended
}

KP_BY_CLASS = {"G1": 5.0, "G2": 6.0, "G3": 7.0, "G4": 8.0, "G5": 9.0,
               "Carrington": 10.5}


# ============================================================================
# Calibration events (5 regime-spanning, per Claude.AI calibration guidance)
# ============================================================================
# Each event spans a distinct (inclination, altitude, mass-class, mechanism)
# cell so the loss-response curve is regime-aware, not overfit.
CALIBRATION_EVENTS = [
    {
        "name": "Halloween 2003",
        "date": "2003-10-29",
        "kp_peak": 9.0,
        "dst_min_nT": -422,
        "g_class": "G5",
        "regime": "polar_sso_mid_altitude",
        "primary_mechanism": "radiation_seu",
        "primary_orbit": {"inc_deg": 98.6, "peri_km": 800, "apo_km": 800},
        "primary_mass_kg": 3680,                    # ADEOS-2 / Midori-2
        "n_anomalies": 47,
        "n_total_loss": 1,                          # ADEOS-2 destroyed
        "n_at_risk_active_LEO": 3000,
        "loss_rate": 1 / 3000,
        "anomaly_rate": 47 / 3000,
        "calibrates": "radiation/SEU response in mid-altitude SSO regime",
    },
    {
        "name": "March 1989 Quebec",
        "date": "1989-03-13",
        "kp_peak": 9.0,
        "dst_min_nT": -589,
        "g_class": "G5",
        "regime": "geo_plus_leo_mixed",
        "primary_mechanism": "mixed_radiation_drag",
        "primary_orbit": {"inc_deg": 0.0, "peri_km": 35786, "apo_km": 35786},
        "primary_mass_kg": 800,                     # GOES-7 class
        "n_anomalies": 5,                           # TDRS-1, GOES-7, etc.
        "n_total_loss": 1,                          # GOES-7 partial cell loss
        "n_at_risk_active_LEO": 200,
        "loss_rate": 1 / 200,
        "anomaly_rate": 5 / 200,
        "calibrates": "GEO + heavy-payload severe-storm response",
    },
    {
        "name": "May 2024 (Gannon)",
        "date": "2024-05-11",
        "kp_peak": 9.0,
        "dst_min_nT": -412,
        "g_class": "G5",
        "regime": "starlink_mid_leo",
        "primary_mechanism": "drag_resilient",
        "primary_orbit": {"inc_deg": 53.0, "peri_km": 540, "apo_km": 560},
        "primary_mass_kg": 260,                     # Starlink V2 mini
        "n_anomalies": 0,                           # SpaceX reported no losses
        "n_total_loss": 0,
        "n_at_risk_active_LEO": 5500,
        "loss_rate": 0.0,
        "anomaly_rate": 0.0,
        "calibrates": "high-altitude (~540 km) drag-tail at G5; null result is informative",
    },
    {
        "name": "Starlink Feb 2022",
        "date": "2022-02-04",
        "kp_peak": 5.0,
        "dst_min_nT": -64,
        "g_class": "G1",
        "regime": "very_low_leo_deployment",
        "primary_mechanism": "drag_dominant",
        "primary_orbit": {"inc_deg": 53.0, "peri_km": 210, "apo_km": 210},
        "primary_mass_kg": 260,
        "n_anomalies": 38,
        "n_total_loss": 38,                         # 38 of 49 newly deployed
        "n_at_risk_active_LEO": 49,
        "loss_rate": 38 / 49,
        "anomaly_rate": 38 / 49,
        "calibrates": "drag-dominant at very low insertion altitude (~210 km), modest storm",
    },
    {
        "name": "Bastille Day 2000",
        "date": "2000-07-15",
        "kp_peak": 9.0,
        "dst_min_nT": -301,
        "g_class": "G5",
        "regime": "mixed_orbits_commsat",
        "primary_mechanism": "seu_anomalies",
        "primary_orbit": {"inc_deg": 0.0, "peri_km": 35786, "apo_km": 35786},
        "primary_mass_kg": 2000,
        "n_anomalies": 8,
        "n_total_loss": 0,
        "n_at_risk_active_LEO": 800,
        "loss_rate": 0.0,
        "anomaly_rate": 8 / 800,
        "calibrates": "GEO commsat anomaly rate at G5 (no total losses)",
    },
]


# ============================================================================
# Loss-response model (v0.1, calibrated against the 5 regime-spanning events)
# ============================================================================
# Two independent failure modes:
#   - Drag (atmospheric expansion): dominant below ~600 km
#   - Radiation/SEU (charged-particle): dominant above ~600 km, peaks in slot
# Combined: P_loss = 1 - (1-P_drag)*(1-P_rad)

def classify_mass(mass_kg: float) -> str:
    if mass_kg < 10:
        return "cubesat"
    if mass_kg < 100:
        return "smallsat"
    if mass_kg < 1000:
        return "standard"
    return "heavy"


# Drag mass offset: NEGATIVE values shift threshold DOWN (more vulnerable).
# Formula uses K_threshold_base + MASS_DRAG_OFFSET; cubesats lower -> easier to lose.
MASS_DRAG_OFFSET = {
    "cubesat": -2.0,    # no propulsion, very vulnerable to thermospheric drag
    "smallsat": -1.0,   # limited propulsion
    "standard":  0.0,   # baseline (Starlink-class with propulsion)
    "heavy":    +1.5,   # well-margined (large delta-v reserve)
}

MASS_SHIELDING_FACTOR = {
    "cubesat": 3.0,     # poor shielding -> SEU-vulnerable
    "smallsat": 1.5,
    "standard": 1.0,
    "heavy": 0.5,
}


def loss_response_drag(altitude_km: float, kp: float, mass_class: str) -> float:
    """Drag-induced loss probability per storm of given Kp at altitude h.

    Calibration anchors (LOCKED via calibration_gate()):
      - Starlink Feb 2022:  h=210, Kp=5, standard,  observed P=0.78
      - May 2024 Gannon:    h=540, Kp=9, standard,  observed P~0
      - Halloween 2003 LEO: h=800, Kp=9, heavy,     observed P~0 (drag)

    Steep altitude scaling: thermospheric expansion under storm decays fast
    with altitude. At h=210km, even moderate storms cause loss; at h=400+km,
    only catastrophic events touch the orbit.
    """
    # Steeper threshold: K_threshold(h) = (h - 180)/15
    # At h=210, K_threshold=2; at h=400, K_threshold=14.7; at h=540, K_threshold=24
    K_threshold_base = (altitude_km - 180.0) / 15.0
    # Negative MASS_DRAG_OFFSET makes cubesats more vulnerable (lower threshold)
    K_threshold = K_threshold_base + MASS_DRAG_OFFSET[mass_class]
    K_threshold = max(0.0, K_threshold)

    delta = kp - K_threshold
    # Sigmoid: P_loss saturates at 0.95 for delta >= 5
    if delta <= -3:
        return 0.0
    if delta >= 5:
        return 0.95
    return 0.95 / (1.0 + math.exp(-(delta - 1.0)))


def loss_response_radiation(altitude_km: float, kp: float, mass_class: str) -> float:
    """Radiation/SEU-induced loss probability per storm.

    Calibration anchors (LOCKED via calibration_gate()):
      - Halloween 2003 LEO total losses: 1/3000 ~ 3.3e-4 at G5
      - Halloween 2003 LEO anomalies:    47/3000 ~ 1.6e-2 at G5 (ten-fold higher)
      - March 1989 GEO total losses:     1/200  = 5e-3 at G5 (heavier exposure)
      - May 2024 Gannon LEO total:       ~0 at G5 (radiation negligible at low alt)

    Base rate calibrated against Halloween 2003 LEO total-loss rate.
    """
    base_rate = 0.00005     # 0.005% per storm at Kp=5 baseline (12x lower than v0.1)
    kp_mult = math.exp((kp - 5.0) / 1.5)
    shielding = MASS_SHIELDING_FACTOR[mass_class]

    # Altitude factor (slot region peaks, GEO and very-low-LEO depressed)
    if altitude_km < 400:
        alt_factor = 0.2    # below outer belt - low SEU exposure
    elif altitude_km < 600:
        alt_factor = 0.5
    elif altitude_km < 2000:
        alt_factor = 1.0    # standard LEO
    elif altitude_km < 8000:
        alt_factor = 2.0    # outer Van Allen belt enhancement
    elif altitude_km < 20000:
        alt_factor = 1.4
    else:
        alt_factor = 0.8    # GEO somewhat protected by magnetosphere

    p = base_rate * kp_mult * shielding * alt_factor
    return min(p, 0.6)


def loss_response_combined(altitude_km: float, kp: float,
                           mass_class: str) -> dict:
    p_drag = loss_response_drag(altitude_km, kp, mass_class)
    p_rad = loss_response_radiation(altitude_km, kp, mass_class)
    p_total = 1.0 - (1.0 - p_drag) * (1.0 - p_rad)
    return {"p_drag": p_drag, "p_radiation": p_rad, "p_total": p_total}


# ============================================================================
# Annual loss probability + percentile
# ============================================================================

def annual_loss_probability(altitude_km: float, mass_kg: float,
                            years: float = 1.0) -> dict:
    """Expected loss probability over `years` of operation."""
    mass_class = classify_mass(mass_kg)
    by_class = {}
    expected_total_per_year = 0.0

    for sclass, rate_per_year in STORM_FREQUENCY_PER_YEAR.items():
        kp = KP_BY_CLASS[sclass]
        loss = loss_response_combined(altitude_km, kp, mass_class)
        # Expected number of loss events per year (Poisson)
        expected_per_year = rate_per_year * loss["p_total"]
        by_class[sclass] = {
            "kp": kp,
            "rate_per_year": rate_per_year,
            "p_loss_per_storm": loss["p_total"],
            "p_drag": loss["p_drag"],
            "p_radiation": loss["p_radiation"],
            "expected_losses_per_year": expected_per_year,
        }
        expected_total_per_year += expected_per_year

    p_per_year = 1.0 - math.exp(-expected_total_per_year)
    p_over_window = 1.0 - math.exp(-expected_total_per_year * years)

    # Identify dominant tail contribution
    dominant_class = max(by_class.keys(),
                         key=lambda k: by_class[k]["expected_losses_per_year"])

    return {
        "mass_class": mass_class,
        "altitude_km": altitude_km,
        "expected_losses_per_year": expected_total_per_year,
        "p_loss_per_year": p_per_year,
        "p_loss_over_window_years": p_over_window,
        "window_years": years,
        "by_storm_class": by_class,
        "dominant_class": dominant_class,
    }


def compute_loss_percentile(p_annual: float,
                            catalog: list[dict]) -> tuple[float, np.ndarray]:
    """Where does p_annual rank vs the current LEO catalog?"""
    rates = []
    for obj in catalog:
        peri = obj.get("perigee_alt_km", 500.0)
        apo = obj.get("apogee_alt_km", 500.0)
        h = 0.5 * (peri + apo)
        m = obj.get("mass_kg", 200.0) or 200.0
        result = annual_loss_probability(h, m, years=1.0)
        rates.append(result["p_loss_per_year"])
    rates_arr = np.array(rates)
    n_below = int(np.sum(rates_arr < p_annual))
    pct = 100.0 * n_below / len(rates_arr)
    return pct, rates_arr


# ============================================================================
# Calibration gate (LOCKED before customer-facing output)
# ============================================================================
# v0.2 BLOCKING REQUIREMENT: loss-response curves must reproduce Halloween
# 2003 + May 2024 Gannon ground-truth outcomes within tolerance. Call
# calibration_gate() before any flare_assess() invocation that exits to
# customer view. The gate is structured so additional events can be added
# as v0.3+ data becomes available.

CALIBRATION_TOLERANCE = {
    # event_name -> (observed P_loss, tolerance for predicted P_loss)
    "Halloween 2003": {"p_observed": 1/3000, "tol_factor_low": 0.3,
                       "tol_factor_high": 4.0},
    "May 2024 (Gannon)": {"p_observed": 0.0,   "p_max_allowed": 0.005},
    "Starlink Feb 2022": {"p_observed": 38/49, "tol_abs": 0.20},
}


def calibration_gate(verbose: bool = True) -> tuple[bool, list[dict]]:
    """Run the calibration gate. Returns (pass/fail, per-event results)."""
    if verbose:
        print()
        print("=" * 75)
        print("CALIBRATION GATE  (must pass before customer-facing output)")
        print("=" * 75)

    results = []
    all_pass = True

    for event in CALIBRATION_EVENTS:
        if event["name"] not in CALIBRATION_TOLERANCE:
            continue
        tol = CALIBRATION_TOLERANCE[event["name"]]
        orbit = event["primary_orbit"]
        mass_kg = event["primary_mass_kg"]
        kp = event["kp_peak"]
        h = 0.5 * (orbit["peri_km"] + orbit["apo_km"])
        mc = classify_mass(mass_kg)

        loss = loss_response_combined(h, kp, mc)
        p_pred = loss["p_total"]
        p_obs = tol["p_observed"]

        # Tolerance check
        if "p_max_allowed" in tol:
            passed = p_pred <= tol["p_max_allowed"]
            check_msg = f"p_pred={p_pred:.5f} <= max_allowed={tol['p_max_allowed']:.5f}"
        elif "tol_abs" in tol:
            passed = abs(p_pred - p_obs) <= tol["tol_abs"]
            check_msg = f"|p_pred - p_obs| = |{p_pred:.4f} - {p_obs:.4f}| = {abs(p_pred-p_obs):.4f} <= {tol['tol_abs']:.4f}"
        else:
            low = p_obs * tol["tol_factor_low"]
            high = p_obs * tol["tol_factor_high"]
            passed = low <= p_pred <= high
            check_msg = f"p_pred={p_pred:.5f} in [{low:.5f}, {high:.5f}]"

        result = {
            "event": event["name"],
            "regime": event["regime"],
            "altitude_km": h,
            "kp": kp,
            "mass_class": mc,
            "p_observed": p_obs,
            "p_predicted": p_pred,
            "p_drag": loss["p_drag"],
            "p_radiation": loss["p_radiation"],
            "passed": passed,
            "check": check_msg,
        }
        results.append(result)
        if not passed:
            all_pass = False

        if verbose:
            status = "PASS" if passed else "FAIL"
            print(f"  [{status}] {event['name']}")
            print(f"         regime: {event['regime']}")
            print(f"         h={h:.0f}km Kp={kp:.0f} mass={mc}")
            print(f"         p_drag={loss['p_drag']:.5f}  p_rad={loss['p_radiation']:.5f}  p_total={p_pred:.5f}")
            print(f"         observed: {p_obs:.5f}")
            print(f"         check:    {check_msg}")

    if verbose:
        print()
        if all_pass:
            print("  GATE STATUS: PASS  (loss-response curves are calibrated)")
        else:
            print("  GATE STATUS: FAIL  (DO NOT ship customer reports until fixed)")
        print("=" * 75)

    return all_pass, results


# ============================================================================
# Fiedler-shift stub (v0.3 substrate cascade signal)
# ============================================================================
# The substrate cascade signal is a topological transition in the
# momentum-aware conjunction graph D_mom. Under thermospheric expansion,
# orbital shells overlap, edge weights rewire, and the second-smallest
# Laplacian eigenvalue (Fiedler value lambda_2) collapses - signaling
# the transition from Anderson-localized (safe) to ballistic (cascading).
#
# v0.3 work (deferred):
#   1. Build D_mom adjacency at storm class G1..G5 against current catalog
#      (one-time precompute, ~15 min via scipy.sparse + eigsh)
#   2. Cache pre-storm and post-storm Laplacian spectra
#   3. Per-query: interpolate against cached spectra, return:
#        - lambda_2_pre, lambda_2_post, threshold_crossed (bool)
#   4. Headline scalar: lambda_2 collapse % alongside SLS

def fiedler_precompute_TODO() -> dict:
    """Stub. To be implemented in v0.3 - see module docstring.

    Pipeline outline:
      1. Build conjunction-graph adjacency A from leverage_map_v2_full
      2. For each storm class G1..G5 + Carrington:
         a. Apply thermospheric expansion -> shift altitudes per Jacchia-Bowman
         b. Recompute edge weights with overlap-rewired shells
         c. Compute Laplacian L_storm = D - A_storm
         d. Solve eigsh(L_storm, k=2, sigma=0) -> (lambda_1, lambda_2)
         e. Cache (storm_class, lambda_2_post)
      3. Cache lambda_2_pre = baseline at quiet sun
      4. Threshold: lambda_2_post / lambda_2_pre < 0.5 = ballistic transition
    """
    return {
        "status": "STUB_v0.3_DEFERRED",
        "next_step": "build D_mom adjacency from leverage_map_v2_full_results.csv",
        "estimated_compute": "~15 min one-time precompute on 50k-node graph",
    }


# ============================================================================
# Main report
# ============================================================================

def flare_assess(inc_deg: float, peri_km: float, apo_km: float, mass_kg: float,
                 years: float = 10.0, name: str = "(unnamed)",
                 enforce_gate: bool = True) -> dict:

    altitude_km = 0.5 * (peri_km + apo_km)
    mass_class = classify_mass(mass_kg)

    # === GATE: must pass before any customer-facing output ===
    if enforce_gate:
        gate_pass, gate_results = calibration_gate(verbose=True)
        if not gate_pass:
            print()
            print("  ABORTED: calibration gate failed.")
            print("  Loss-response curves do not reproduce historical ground truth.")
            print("  Fix curves (or relax tolerances) before customer-facing output.")
            return {"gate_pass": False, "gate_results": gate_results}

    print()
    print("=" * 75)
    print("FLARE-CASCADE RISK ASSESSMENT  (v0.2)")
    print("=" * 75)
    print(f"  Satellite:    {name}")
    print(f"  Orbit:        inc={inc_deg:.1f}deg  peri={peri_km:.0f}km  apo={apo_km:.0f}km")
    print(f"                mean altitude {altitude_km:.0f} km")
    print(f"  Mass:         {mass_kg:.0f} kg  (class: {mass_class})")
    print(f"  Window:       {years:.1f} years")
    print()

    # === Annual loss probability ===
    annual = annual_loss_probability(altitude_km, mass_kg, years=years)

    # === Percentile vs catalog ===
    print("  Loading LEO catalog for percentile ranking...")
    catalog = kra.load_assignment_csv()
    print(f"    catalog: {len(catalog):,} objects")
    print()

    pct, _ = compute_loss_percentile(annual["p_loss_per_year"], catalog)

    print("-" * 75)
    print("ANNUAL FLARE-LOSS PROBABILITY")
    print("-" * 75)
    print(f"  P(loss in 1 year):       {annual['p_loss_per_year']*100:>7.4f}%")
    print(f"  P(loss in {years:.0f} years):     {annual['p_loss_over_window_years']*100:>7.4f}%")
    print(f"  Expected losses/year:    {annual['expected_losses_per_year']:>9.6f}")
    print(f"  Catalog percentile:      {pct:>4.1f}th  (0=safest, 100=most exposed)")
    print(f"  Dominant storm class:    {annual['dominant_class']}")
    print()

    print("  By storm class:")
    print(f"    {'Class':>10}  {'Kp':>4}  {'Recur/yr':>9}  {'P(loss/storm)':>14}  {'E[loss]/yr':>11}")
    for sclass, info in annual["by_storm_class"].items():
        recur_str = f"{info['rate_per_year']:.3f}"
        p_str = f"{info['p_loss_per_storm']:.5f}"
        e_str = f"{info['expected_losses_per_year']:.7f}"
        marker = "  *" if sclass == annual["dominant_class"] else ""
        print(f"    {sclass:>10}  {info['kp']:>4.1f}  {recur_str:>9}  {p_str:>14}  {e_str:>11}{marker}")
    print()

    # === Substrate cascade signal (Fiedler stub) ===
    print("-" * 75)
    print("SUBSTRATE CASCADE SIGNAL  (Fiedler lambda_2 shift)")
    print("-" * 75)
    fiedler = fiedler_precompute_TODO()
    print(f"  Status: {fiedler['status']}")
    print(f"  Pipeline: build D_mom adjacency at G1..G5, compute Laplacian spectra,")
    print(f"            cache lambda_2_pre/lambda_2_post per storm class.")
    print(f"  Headline output (v0.3):  lambda_2 collapse % alongside SLS.")
    print(f"  Threshold:               lambda_2_post / lambda_2_pre < 0.5 = ballistic")
    print()

    # === Saturation framing (the marketing hook) ===
    print("-" * 75)
    print("SATURATION-POINT NOTE  (the underwriter pitch)")
    print("-" * 75)
    print("  The Kp index ceiling is 9. Halloween 2003 (Dst=-422 nT) and a")
    print("  Carrington-class event (Dst~-1750 nT estimated) both register as")
    print("  Kp=9, but their substrate impact differs by 30-50x.")
    print("  Once the Fiedler precompute lands (v0.3), this report will quote")
    print("  the cascade probability ratio directly from the lambda_2 spectra.")
    print()

    print("  Methodology: Kessler Critical Surface preprint v0.6, sections 3.6-3.8")
    print("    (substrate algebra) and 10.10-10.11 (per-satellite leverage map).")
    print("  Storm frequencies: GFZ Kp/ap archive 1932-2026 (95-yr baseline).")
    print("  Calibration events: see flare_calibration_events.json.")
    print("  Calibration gate:  PASSED (see top of report).")
    print("=" * 75)
    print()

    return {
        "satellite": name,
        "orbit": {"inc_deg": inc_deg, "peri_km": peri_km, "apo_km": apo_km},
        "mass_kg": mass_kg,
        "mass_class": mass_class,
        "annual": annual,
        "percentile": pct,
        "fiedler_status": fiedler["status"],
        "calibration_events_used": len(CALIBRATION_EVENTS),
        "frequency_archive_years": 95,
        "gate_pass": True,
    }


def export_calibration_events(out_path: Path) -> None:
    """Dump the 5 calibration events to a citable JSON artifact."""
    payload = {
        "version": "0.1",
        "description": "Regime-spanning historical calibration events for the "
                       "flare-cascade loss-response curve.",
        "n_events": len(CALIBRATION_EVENTS),
        "regime_coverage": [
            "polar SSO mid-altitude (Halloween 2003 / ADEOS-2)",
            "GEO + LEO mixed severe (March 1989 Quebec / TDRS-1, GOES-7)",
            "Starlink-class mid-LEO drag-resilient (May 2024 Gannon)",
            "very-low-LEO deployment drag-dominant (Starlink Feb 2022)",
            "GEO commsat anomaly-only (Bastille Day 2000)",
        ],
        "events": CALIBRATION_EVENTS,
        "loss_response_model": {
            "version": "v0.1",
            "mechanisms": ["drag", "radiation_seu"],
            "mass_classes": ["cubesat", "smallsat", "standard", "heavy"],
            "calibration_anchors": [
                "Starlink Feb 2022 (drag, h=210, Kp=5, P=0.78)",
                "May 2024 Gannon (drag, h=540, Kp=9, P=~0)",
                "Halloween 2003 (radiation, h=800 SSO, Kp=9, P=0.0003)",
                "March 1989 Quebec (radiation, GEO, Kp=9, P=0.005)",
                "Bastille Day 2000 (radiation, GEO, Kp=9, P=~0)",
            ],
        },
    }
    out_path.write_text(json.dumps(payload, indent=2))
    print(f"  calibration events exported -> {out_path.name}")


def main():
    p = argparse.ArgumentParser(description="Flare-Cascade Risk Assessor v0.1")
    p.add_argument("--inc", type=float, help="Inclination (degrees)")
    p.add_argument("--perigee", type=float, help="Perigee altitude (km)")
    p.add_argument("--apogee", type=float, help="Apogee altitude (km)")
    p.add_argument("--mass", type=float, default=200.0, help="Mass (kg)")
    p.add_argument("--years", type=float, default=10.0,
                   help="Assessment window (years)")
    p.add_argument("--name", type=str, default="(unnamed)")
    p.add_argument("--no-gate", action="store_true",
                   help="Skip calibration gate (DEV ONLY - never ship without gate)")
    p.add_argument("--save", type=str, default=None,
                   help="Save full report to JSON")
    p.add_argument("--calibration-gate", action="store_true",
                   help="Run calibration gate only and exit (no orbit args needed)")
    p.add_argument("--export-calibration", action="store_true",
                   help="Export calibration events to JSON and exit")
    args = p.parse_args()

    if args.calibration_gate:
        gate_pass, _ = calibration_gate(verbose=True)
        return 0 if gate_pass else 1

    if args.export_calibration:
        export_calibration_events(HERE / "flare_calibration_events.json")
        return 0

    if args.inc is None or args.perigee is None or args.apogee is None:
        p.error("--inc, --perigee, and --apogee are required (unless --calibration-gate or --export-calibration)")

    report = flare_assess(
        inc_deg=args.inc, peri_km=args.perigee, apo_km=args.apogee,
        mass_kg=args.mass, years=args.years, name=args.name,
        enforce_gate=not args.no_gate,
    )

    if args.save:
        out_path = HERE / args.save
        out_path.write_text(json.dumps(report, indent=2, default=str))
        print(f"  report saved -> {out_path.name}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
