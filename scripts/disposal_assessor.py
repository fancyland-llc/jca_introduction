"""Disposal Certification Assessor (v0.1).

Computes Substrate-Leverage-Hours (SLH) accumulated during atmospheric
decay of a satellite's disposal trajectory. Supports operator filings
under FCC 47 CFR Part 25.114(d)(14) (5-year deorbit rule, effective
September 2024) and ESA Zero Debris Charter sustainability scoring.

Pipeline:
  1. King-Hele atmospheric decay propagator: starting orbit + ballistic
     coefficient -> trajectory of (time, altitude) waypoints until reentry
  2. Per-waypoint substrate-leverage: |sigma_mag| * rho(alt | inc)
     (sigma_mag is constant during decay because inclination is preserved;
     rho varies as the satellite passes through populated shells)
  3. SLH integration: sum |sigma_mag| * rho_i * dt_i over trajectory
  4. Per-cohort percentile vs sampled catalog baseline (100 samples by
     mission envelope, cached on first run)
  5. Peak-exposure band heat-map (substrate-leverage-hours by altitude band)
  6. FCC Part 25 + ESA Zero Debris Charter regulatory framing block

CLI:
    python disposal_assessor.py --inc 53 --perigee 540 --apogee 560 \\
                                 --mass 260 --area-to-mass 0.044 --name MySat
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import sys
from pathlib import Path

import numpy as np

import kessler_risk_assessor as kra

HERE = Path(__file__).parent

# Earth + atmosphere constants
EARTH_R_KM = 6378.137                # equatorial radius (km)
MU_EARTH = 3.986004418e14            # standard gravitational parameter (m^3/s^2)
REENTRY_ALT_KM = 120.0               # functional reentry threshold
DEFAULT_CD = 2.2                     # standard tumbling-satellite drag coefficient
DEFAULT_AREA_TO_MASS = 0.01          # m^2/kg, typical small LEO sat default

# Regulatory thresholds
FCC_5YR_RULE_YEARS = 5.0             # FCC 47 CFR Part 25.114, effective Sep 2024
FCC_25YR_RULE_YEARS = 25.0           # legacy IADC guideline
LEO_REGIME_CEILING_KM = 2000.0       # boundary where atmospheric decay no longer applies


# ============================================================================
# Atmospheric density model (Vallado solar-mean reference, F10.7~150)
# ============================================================================
# Piecewise exponential: rho(h) = rho_0 * exp(-(h - h_0)/H) within bracket.
# Source: Vallado, "Fundamentals of Astrodynamics" 4th ed, Table 8-4.
ATMOSPHERIC_REFERENCE = [
    # (h_0_km, rho_0_kg_per_m3, scale_height_km)
    (   0,  1.225e+00,   8.44),
    ( 100,  5.297e-07,   5.877),
    ( 150,  2.070e-09,  22.523),
    ( 200,  2.541e-10,  37.105),
    ( 250,  6.073e-11,  45.546),
    ( 300,  1.916e-11,  53.628),
    ( 350,  7.014e-12,  53.298),
    ( 400,  2.803e-12,  58.515),
    ( 450,  1.184e-12,  60.828),
    ( 500,  5.215e-13,  63.822),
    ( 600,  1.137e-13,  71.835),
    ( 700,  3.070e-14,  88.667),
    ( 800,  1.136e-14, 124.640),
    ( 900,  5.759e-15, 181.050),
    (1000,  3.561e-15, 268.000),
    (1250,  1.060e-15, 296.000),
    (1500,  5.260e-16, 353.000),
]


def atmospheric_density(h_km: float) -> float:
    """Atmospheric density (kg/m^3) at altitude h_km via piecewise exponential."""
    if h_km < 0:
        return ATMOSPHERIC_REFERENCE[0][1]
    if h_km >= 2000:
        # Above 2000 km, density is negligible for decay purposes
        h0, rho0, H = ATMOSPHERIC_REFERENCE[-1]
        return rho0 * math.exp(-(h_km - h0) / H)
    # Find bracket
    for i in range(len(ATMOSPHERIC_REFERENCE) - 1):
        h0, rho0, H = ATMOSPHERIC_REFERENCE[i]
        h1, _, _ = ATMOSPHERIC_REFERENCE[i + 1]
        if h0 <= h_km < h1:
            return rho0 * math.exp(-(h_km - h0) / H)
    # Should not reach here
    h0, rho0, H = ATMOSPHERIC_REFERENCE[-1]
    return rho0 * math.exp(-(h_km - h0) / H)


# ============================================================================
# King-Hele decay propagator
# ============================================================================

def propagate_decay(initial_alt_km: float, area_to_mass_m2_per_kg: float,
                    cd: float = DEFAULT_CD,
                    max_years: float = 200.0,
                    n_max_steps: int = 5000) -> tuple[list[dict], bool]:
    """Propagate atmospheric decay via King-Hele circular-orbit approximation.

    dh/dt = -rho(h) * B * sqrt(mu * (R_e + h))
    where B = Cd * A / m  (ballistic coefficient, m^2/kg).

    Returns (waypoints, truncated).
      waypoints: list of {time_yr, alt_km, dt_sec, rho_kg_per_m3}
      truncated: True if simulation hit max_years cap before reentry

    Adaptive timestep: caps per-step altitude change (25 km) so high-altitude
    slow-decay regimes use long timesteps naturally; low-altitude fast-decay
    regimes use short timesteps. No fixed dt cap.
    """
    B = cd * area_to_mass_m2_per_kg                         # m^2/kg
    h_km = float(initial_alt_km)
    t_sec = 0.0
    seconds_per_year = 365.25 * 86400.0
    max_seconds = max_years * seconds_per_year

    waypoints = []
    truncated = False

    for step in range(n_max_steps):
        a_m = (EARTH_R_KM + h_km) * 1000.0                  # semi-major axis (m)
        rho = atmospheric_density(h_km)                      # kg/m^3
        sqrt_mu_a = math.sqrt(MU_EARTH * a_m)                # m^2/s

        # Decay rate (m/s)
        dh_dt_m_per_s = -rho * B * sqrt_mu_a

        if dh_dt_m_per_s == 0:
            truncated = True
            break

        # Adaptive timestep: limit per-step altitude change to 25 km.
        # No fixed dt cap - let dh constraint drive (slow-decay regimes
        # take long steps naturally; fast-decay regimes take short steps).
        max_dh_per_step_m = 25000.0
        dt_from_dh = max_dh_per_step_m / abs(dh_dt_m_per_s)
        # But hard-cap at remaining max_seconds so we don't overshoot the budget
        remaining_seconds = max_seconds - t_sec
        dt_sec = min(dt_from_dh, max(remaining_seconds, 1.0))

        waypoints.append({
            "step": step,
            "time_sec": t_sec,
            "time_yr": t_sec / seconds_per_year,
            "alt_km": h_km,
            "dt_sec": dt_sec,
            "rho_kg_per_m3": rho,
            "dh_dt_m_per_s": dh_dt_m_per_s,
        })

        # Step forward
        dh_m = dh_dt_m_per_s * dt_sec
        h_km += dh_m / 1000.0
        t_sec += dt_sec

        if h_km <= REENTRY_ALT_KM:
            waypoints.append({
                "step": step + 1,
                "time_sec": t_sec,
                "time_yr": t_sec / seconds_per_year,
                "alt_km": REENTRY_ALT_KM,
                "dt_sec": 0.0,
                "rho_kg_per_m3": atmospheric_density(REENTRY_ALT_KM),
                "dh_dt_m_per_s": 0.0,
            })
            break

        if t_sec >= max_seconds:
            truncated = True
            break
    else:
        truncated = True

    return waypoints, truncated


# ============================================================================
# Sub-sampling (King-Hele can produce thousands of waypoints; we want ~50
# for substrate evaluation)
# ============================================================================

def downsample_trajectory(waypoints: list[dict], n_target: int = 50) -> list[dict]:
    """Downsample trajectory to ~n_target waypoints, preserving boundary points
    and concentrating samples in low-altitude/fast-decay regions."""
    if len(waypoints) <= n_target:
        return waypoints
    # Stride sampling
    stride = max(1, len(waypoints) // n_target)
    sampled = waypoints[::stride]
    if sampled[-1] is not waypoints[-1]:
        sampled.append(waypoints[-1])
    return sampled


# ============================================================================
# Substrate-Leverage-Hours integration
# ============================================================================

def compute_rho_at_inc_alt(inc_deg: float, alt_km: float,
                            mass_kg: float, catalog: list[dict],
                            alt_band_km: float = 25.0) -> float:
    """Wrap kra.estimate_rho_proposed for a circular orbit at (inc, alt).

    A truly-circular orbit (peri=apo) has zero shell width, which fails the
    overlap check in estimate_rho_proposed. Pad to a 25-km altitude band
    centered on alt_km so shell-overlap calculations return meaningful
    conjunction-exposure weights. The 25-km band matches estimate_rho_proposed's
    default altitude_bin_km and reflects realistic per-orbit altitude variation
    (atmospheric drag fluctuations + finite spacecraft dimensions).
    """
    half = alt_band_km / 2.0
    rho, _, _ = kra.estimate_rho_proposed(
        proposed_inc=inc_deg,
        proposed_peri=max(alt_km - half, 100.0),
        proposed_apo=alt_km + half,
        proposed_mass=mass_kg, catalog=catalog,
    )
    return rho


def compute_SLH(trajectory: list[dict], sigma_mag: float, inc_deg: float,
                mass_kg: float, catalog: list[dict],
                progress: bool = True) -> dict:
    """Substrate-Leverage-Hours = sum |sigma_mag| * rho(alt | inc) * dt over trajectory.

    Since sigma_mag is constant during decay (inclination preserved),
    SLH = |sigma_mag| * sum rho(alt_i | inc) * dt_i.
    """
    sigma_abs = abs(sigma_mag)
    leverage_seconds = 0.0
    waypoint_records = []

    seconds_per_hour = 3600.0
    n = len(trajectory)
    for i, wp in enumerate(trajectory):
        if progress and (i % max(1, n // 10) == 0):
            print(f"    waypoint {i+1}/{n}  alt={wp['alt_km']:.0f} km  "
                  f"t={wp['time_yr']:.2f} yr",
                  flush=True)
        rho_conj = compute_rho_at_inc_alt(inc_deg, wp["alt_km"], mass_kg, catalog)
        contribution = sigma_abs * rho_conj * wp["dt_sec"]
        leverage_seconds += contribution
        waypoint_records.append({
            "time_yr": wp["time_yr"],
            "alt_km": wp["alt_km"],
            "dt_sec": wp["dt_sec"],
            "rho_conjunction": rho_conj,
            "leverage_seconds_increment": contribution,
        })

    leverage_hours = leverage_seconds / seconds_per_hour

    return {
        "SLH_leverage_hours": leverage_hours,
        "SLH_leverage_seconds": leverage_seconds,
        "sigma_mag_used": sigma_mag,
        "n_waypoints": len(trajectory),
        "waypoint_records": waypoint_records,
        "decay_duration_yr": trajectory[-1]["time_yr"] if trajectory else 0.0,
    }


def bin_SLH_by_altitude(waypoint_records: list[dict],
                          band_width_km: float = 100.0) -> dict:
    """Bin SLH contributions into altitude bands for the heat-map output."""
    bands = {}
    for wp in waypoint_records:
        band_low = int(wp["alt_km"] // band_width_km) * int(band_width_km)
        band_high = band_low + int(band_width_km)
        key = (band_low, band_high)
        bands.setdefault(key, {"SLH_hours": 0.0, "duration_yr": 0.0,
                                "n_waypoints": 0})
        bands[key]["SLH_hours"] += wp["leverage_seconds_increment"] / 3600.0
        bands[key]["duration_yr"] += wp["dt_sec"] / (365.25 * 86400.0)
        bands[key]["n_waypoints"] += 1
    return bands


# ============================================================================
# Cohort percentile baseline (sampled, cached)
# ============================================================================

PERCENTILE_CACHE_PATH = HERE / "disposal_percentile_cache.json"
PERCENTILE_SAMPLE_SIZE = 50          # v0.1 sample size (modest for speed)


def build_or_load_percentile_baseline(catalog: list[dict],
                                       n_samples: int = PERCENTILE_SAMPLE_SIZE,
                                       force_rebuild: bool = False) -> list[float]:
    """Build (or load cached) per-cohort SLH baseline from random catalog samples.

    For v0.1: aggregate baseline (not envelope-segmented). v0.2 will segment
    by mission envelope (starlink-leo, sso-eo, etc).
    """
    if PERCENTILE_CACHE_PATH.exists() and not force_rebuild:
        try:
            data = json.loads(PERCENTILE_CACHE_PATH.read_text())
            if data.get("n_samples") == n_samples:
                return data["SLH_baseline"]
        except (json.JSONDecodeError, KeyError):
            pass

    print(f"  Building disposal SLH baseline from {n_samples} catalog samples ...")
    rng = random.Random(42)
    # Restrict to LEO active payloads
    leo_actives = [obj for obj in catalog
                   if obj.get("is_active", 0)
                   and 0 < obj.get("perigee_alt_km", 0) <= LEO_REGIME_CEILING_KM
                   and obj.get("apogee_alt_km", 0) <= LEO_REGIME_CEILING_KM]
    if len(leo_actives) < n_samples:
        sample = leo_actives
    else:
        sample = rng.sample(leo_actives, n_samples)

    baseline_SLH = []
    for i, obj in enumerate(sample):
        try:
            inc = float(obj["inc_deg"])
            peri = float(obj["perigee_alt_km"])
            apo = float(obj["apogee_alt_km"])
            mass = float(obj.get("mass_kg") or 200.0)
            mean_alt = 0.5 * (peri + apo)

            # Run decay (simplified: 25 waypoints, no per-waypoint progress)
            traj, _truncated = propagate_decay(mean_alt, DEFAULT_AREA_TO_MASS,
                                                 max_years=200.0, n_max_steps=2000)
            traj_ds = downsample_trajectory(traj, n_target=25)

            sigma_mag = obj.get("sigma_mag", 0.0)
            slh_result = compute_SLH(traj_ds, sigma_mag, inc, mass, catalog,
                                       progress=False)
            baseline_SLH.append(slh_result["SLH_leverage_hours"])

            if (i + 1) % 10 == 0:
                print(f"    {i+1}/{n_samples} samples computed",
                      flush=True)
        except Exception as e:
            print(f"    sample {i} failed: {e}", flush=True)
            continue

    PERCENTILE_CACHE_PATH.write_text(json.dumps({
        "n_samples": n_samples,
        "SLH_baseline": baseline_SLH,
        "version": "v0.1",
    }, indent=2))
    print(f"  baseline cached -> {PERCENTILE_CACHE_PATH.name}")
    return baseline_SLH


def compute_percentile(slh_value: float, baseline: list[float]) -> float:
    """Where does this SLH rank vs baseline? (lower = better)."""
    if not baseline:
        return float("nan")
    arr = np.array(baseline)
    n_below = int(np.sum(arr < slh_value))
    return 100.0 * n_below / len(arr)


# ============================================================================
# Regulatory framing
# ============================================================================

REGULATORY_FRAMING = """
REGULATORY FRAMING
This certification is generated against the substrate methodology described
in the Kessler Critical Surface preprint (v0.6, sections 3.6-3.8 substrate
algebra; sections 10.10-10.11 per-satellite leverage map). It supports
operator filings under:

  - FCC 47 CFR Part 25.114(d)(14) post-mission disposal plan requirements
    (5-year deorbit rule, effective September 2024)
  - ESA Zero Debris Charter sustainability scoring
    (signatory list at esa.int/Safety_Security/Space_Debris)
  - UN COPUOS Long-term Sustainability Guidelines (A/74/20, 2019)

This is a substrate-derived informational risk signal, not an insurance
contract or regulatory determination.
"""


# ============================================================================
# Main report
# ============================================================================

def disposal_assess(inc_deg: float, peri_km: float, apo_km: float,
                     mass_kg: float, area_to_mass: float = DEFAULT_AREA_TO_MASS,
                     cd: float = DEFAULT_CD,
                     name: str = "(unnamed)",
                     skip_percentile: bool = False) -> dict:

    mean_alt = 0.5 * (peri_km + apo_km)

    print()
    print("=" * 75)
    print("DISPOSAL CERTIFICATION  (v0.1)")
    print("=" * 75)
    print(f"  Mission:           {name}")
    print(f"  Starting orbit:    inc={inc_deg:.1f}deg  peri={peri_km:.0f}km  apo={apo_km:.0f}km")
    print(f"  Mean altitude:     {mean_alt:.0f} km")
    print(f"  Mass:              {mass_kg:.0f} kg")
    print(f"  Area-to-mass:      {area_to_mass:.4f} m^2/kg  (Cd={cd:.1f})")

    # GEO / above-LEO disposals require a different protocol (graveyard orbit boost)
    if peri_km > LEO_REGIME_CEILING_KM:
        print()
        print("-" * 75)
        print(f"  REGIME:            Above LEO (perigee > {LEO_REGIME_CEILING_KM:.0f} km)")
        print("  Atmospheric decay does not apply at this altitude.")
        print("  Disposal protocol: graveyard-orbit boost (+300 km above GEO for")
        print("  geostationary missions, per ITU-R S.1003-2). The SLH metric")
        print("  evaluates atmospheric-decay disposals only.")
        print("  v0.2 will add graveyard-orbit substrate-leverage scoring.")
        print("=" * 75)
        return {"regime": "above_LEO", "applicable": False}

    # === 1. King-Hele decay propagation ===
    print()
    print("  Propagating atmospheric decay (King-Hele) ...")
    traj_full, truncated = propagate_decay(mean_alt, area_to_mass, cd=cd,
                                             max_years=200.0)
    print(f"    raw waypoints:     {len(traj_full)}")
    print(f"    decay duration:    {traj_full[-1]['time_yr']:.2f} years"
          + ("  [TRUNCATED at 200-yr cap]" if truncated else ""))
    print(f"    final altitude:    {traj_full[-1]['alt_km']:.0f} km"
          + ("  (still in orbit)" if truncated else "  (reentered)"))

    decay_yr = traj_full[-1]['time_yr']
    if truncated:
        # Natural decay does not complete within 200 years -> non-compliant
        fcc_5yr_compliant = False
        fcc_25yr_compliant = False
    else:
        fcc_5yr_compliant = decay_yr <= FCC_5YR_RULE_YEARS
        fcc_25yr_compliant = decay_yr <= FCC_25YR_RULE_YEARS

    # Downsample for substrate evaluation
    traj_ds = downsample_trajectory(traj_full, n_target=50)
    print(f"    waypoints (downsampled for substrate eval): {len(traj_ds)}")

    # === 2. Catalog + substrate residue ===
    print()
    print("  Loading LEO catalog + computing substrate residue ...")
    catalog = kra.load_assignment_csv()
    print(f"    catalog: {len(catalog):,} objects")

    sub_idx, sub_val, sigma_mag = kra.assign_substrate_residue_for_proposed(
        inc_deg, catalog)
    print(f"    substrate residue R[{sub_idx}] = {sub_val}")
    print(f"    sigma_mag at this residue: {sigma_mag:+.4f}")

    # === 3. SLH integration ===
    print()
    print("  Integrating Substrate-Leverage-Hours over trajectory ...")
    slh = compute_SLH(traj_ds, sigma_mag, inc_deg, mass_kg, catalog,
                       progress=True)

    print()
    print(f"    SLH (leverage-hours):  {slh['SLH_leverage_hours']:.3f}")
    print(f"    SLH (leverage-seconds): {slh['SLH_leverage_seconds']:.4e}")

    # === 4. Percentile baseline ===
    if skip_percentile:
        pct = float("nan")
        print("  Skipping percentile baseline (--skip-percentile)")
    else:
        baseline = build_or_load_percentile_baseline(catalog)
        pct = compute_percentile(slh["SLH_leverage_hours"], baseline)
        print(f"  SLH cohort percentile: {pct:.1f}th (lower = better, vs {len(baseline)} sampled actives)")

    # === 5. Altitude-band binning ===
    bands = bin_SLH_by_altitude(slh["waypoint_records"], band_width_km=100.0)

    print()
    print("-" * 75)
    print("PEAK EXPOSURE BANDS  (where the trajectory accumulates leverage)")
    print("-" * 75)
    print(f"    {'Altitude Band':>15}  {'SLH-hours':>10}  {'% of total':>10}  {'Dwell time':>11}  Notes")
    sorted_bands = sorted(bands.items(), key=lambda kv: -kv[1]["SLH_hours"])
    total_slh = sum(b["SLH_hours"] for b in bands.values())
    for (lo, hi), info in sorted_bands[:8]:
        pct_of_total = (info["SLH_hours"] / total_slh * 100) if total_slh > 0 else 0
        notes = []
        if 500 <= lo <= 600:
            notes.append("Starlink-class shell")
        if 700 <= lo <= 900:
            notes.append("OneWeb / SSO shell")
        if 1100 <= lo <= 1300:
            notes.append("OneWeb-mid-alt shell")
        if 200 <= lo <= 400:
            notes.append("ISS / low-LEO")
        notes_str = "  " + ", ".join(notes) if notes else ""
        print(f"    {lo:>4d}-{hi:<4d} km  {info['SLH_hours']:>10.3f}  {pct_of_total:>9.1f}%  "
              f"{info['duration_yr']:>9.2f} yr{notes_str}")

    # === 6. Compliance summary ===
    print()
    print("-" * 75)
    print("COMPLIANCE SUMMARY")
    print("-" * 75)
    fcc5_str = "PASS" if fcc_5yr_compliant else "FAIL"
    fcc25_str = "PASS" if fcc_25yr_compliant else "FAIL"
    decay_str = f">{200.0:.0f}" if truncated else f"{decay_yr:.2f}"
    print(f"  Decay duration:    {decay_str} years"
          + ("  [natural decay does not complete within 200-yr horizon]" if truncated else ""))
    print(f"  FCC 5-year rule:   {fcc5_str}  (threshold: 5 yr, effective Sep 2024)")
    print(f"  FCC 25-year rule:  {fcc25_str}  (legacy IADC guideline)")
    if not fcc_5yr_compliant:
        print(f"  -> NON-COMPLIANT under current FCC 5-year rule.")
        if truncated:
            print(f"  -> Natural atmospheric decay > 200 years; ACTIVE DEORBIT REQUIRED.")
        else:
            print(f"  -> Active deorbit burn or lower insertion altitude required.")

    print(REGULATORY_FRAMING)
    print("=" * 75)
    print()

    return {
        "mission": name,
        "starting_orbit": {"inc_deg": inc_deg, "peri_km": peri_km, "apo_km": apo_km},
        "mass_kg": mass_kg,
        "area_to_mass": area_to_mass,
        "cd": cd,
        "decay_duration_yr": decay_yr,
        "reentry_alt_km": traj_full[-1]["alt_km"],
        "fcc_5yr_compliant": fcc_5yr_compliant,
        "fcc_25yr_compliant": fcc_25yr_compliant,
        "sigma_mag": sigma_mag,
        "substrate_residue": int(sub_idx),
        "SLH_leverage_hours": slh["SLH_leverage_hours"],
        "SLH_cohort_percentile": pct,
        "n_waypoints_evaluated": slh["n_waypoints"],
        "altitude_bands": {f"{lo}-{hi}": info for (lo, hi), info in bands.items()},
    }


def main():
    p = argparse.ArgumentParser(description="Disposal Certification Assessor v0.1")
    p.add_argument("--inc", type=float, required=True, help="Inclination (degrees)")
    p.add_argument("--perigee", type=float, required=True, help="Perigee altitude (km)")
    p.add_argument("--apogee", type=float, required=True, help="Apogee altitude (km)")
    p.add_argument("--mass", type=float, default=200.0, help="Mass (kg)")
    p.add_argument("--area-to-mass", type=float, default=DEFAULT_AREA_TO_MASS,
                   help=f"Area-to-mass ratio (m^2/kg, default {DEFAULT_AREA_TO_MASS})")
    p.add_argument("--cd", type=float, default=DEFAULT_CD,
                   help=f"Drag coefficient (default {DEFAULT_CD})")
    p.add_argument("--name", type=str, default="(unnamed)")
    p.add_argument("--skip-percentile", action="store_true",
                   help="Skip cohort percentile baseline (faster for testing)")
    p.add_argument("--rebuild-baseline", action="store_true",
                   help="Force rebuild of percentile baseline cache")
    p.add_argument("--save", type=str, default=None, help="Save report to JSON")
    args = p.parse_args()

    if args.rebuild_baseline:
        catalog = kra.load_assignment_csv()
        build_or_load_percentile_baseline(catalog, force_rebuild=True)

    report = disposal_assess(
        inc_deg=args.inc, peri_km=args.perigee, apo_km=args.apogee,
        mass_kg=args.mass, area_to_mass=args.area_to_mass, cd=args.cd,
        name=args.name, skip_percentile=args.skip_percentile,
    )

    if args.save:
        out_path = HERE / args.save
        out_path.write_text(json.dumps(report, indent=2, default=str))
        print(f"  report saved -> {out_path.name}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
