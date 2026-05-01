"""
discosweb_loader.py — pull LEO catalog data from ESA DISCOSweb

DISCOSweb (Database and Information System Characterising Objects in
Space) is ESA's public catalog of tracked space objects, including
launch records, decay status, RCS / dimensions / masses, and orbital
classification. It does NOT carry real-time conjunction data (that's
Space-Track / 18 SDS / LeoLabs); but it gives the OBJECT CATALOG side
of the Kessler simulator — what's up there, with what mass and
cross-section.

Reads the API token from the DISCOSWEB_TOKEN environment variable
(typically set in a gitignored .env file at the repo root or in
this directory). Token is generated at:
    https://discosweb.esoc.esa.int/tokens

Usage:
    # one-time:  pip install python-dotenv requests
    # in .env:   DISCOSWEB_TOKEN=<your-token>

    python discosweb_loader.py --query basic
    python discosweb_loader.py --query leo --max-pages 3

Outputs:
    discosweb_<query>_<timestamp>.json — raw API response cache
    discosweb_<query>_summary.csv      — flattened table for downstream

API docs: https://discosweb.esoc.esa.int/apidocs/
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).parent

API_BASE = "https://discosweb.esoc.esa.int/api"


def load_token() -> str:
    """Load DISCOSWEB_TOKEN from environment or .env file in this dir."""
    token = os.environ.get("DISCOSWEB_TOKEN")
    if token:
        return token
    env_path = HERE / ".env"
    if env_path.exists():
        for line in env_path.read_text().splitlines():
            line = line.strip()
            if line.startswith("DISCOSWEB_TOKEN="):
                return line.split("=", 1)[1].strip().strip('"').strip("'")
    raise RuntimeError(
        "DISCOSWEB_TOKEN not set. Either export it in your shell, or copy "
        f"{HERE}/.env.example to {HERE}/.env and fill in the token."
    )


def fetch_page(endpoint: str, params: dict, token: str) -> dict:
    import requests
    headers = {
        "Authorization": f"Bearer {token}",
        "DiscosWeb-Api-Version": "2",
    }
    url = f"{API_BASE}/{endpoint.lstrip('/')}"
    resp = requests.get(url, headers=headers, params=params, timeout=30)
    if resp.status_code == 429:
        # rate limited
        retry = int(resp.headers.get("Retry-After", "60"))
        print(f"    rate-limited; sleeping {retry}s")
        time.sleep(retry)
        return fetch_page(endpoint, params, token)
    if not resp.ok:
        raise RuntimeError(f"DISCOSweb API {resp.status_code}: {resp.text[:300]}")
    return resp.json()


def query_basic(token: str) -> list[dict]:
    """Smallest possible query — one page of objects, default ordering.

    Useful as a 'does the token work + does the API respond' smoke test.
    """
    print("  Querying DISCOSweb /api/objects (1 page, default sort)...")
    data = fetch_page("objects", {"page[size]": 30}, token)
    return data.get("data", [])


def fetch_page_with_meta(endpoint: str, params: dict, token: str) -> dict:
    """Like fetch_page but returns the full envelope (data + meta + links)."""
    import requests
    headers = {
        "Authorization": f"Bearer {token}",
        "DiscosWeb-Api-Version": "2",
    }
    url = f"{API_BASE}/{endpoint.lstrip('/')}"
    resp = requests.get(url, headers=headers, params=params, timeout=30)
    if resp.status_code == 429:
        retry = int(resp.headers.get("Retry-After", "60"))
        print(f"    rate-limited; sleeping {retry}s")
        time.sleep(retry)
        return fetch_page_with_meta(endpoint, params, token)
    if not resp.ok:
        raise RuntimeError(f"DISCOSweb API {resp.status_code}: {resp.text[:300]}")
    return resp.json()


EARTH_RADIUS_KM = 6378.137


def _orbit_to_perigee_apogee_km(orbit_attrs: dict) -> tuple[float | None, float | None]:
    """Compute perigee/apogee altitude (km) from sma (meters) + eccentricity."""
    sma_m = orbit_attrs.get("sma")
    ecc = orbit_attrs.get("ecc")
    if sma_m is None or ecc is None:
        return (None, None)
    sma_km = sma_m / 1000.0
    perigee_alt = sma_km * (1.0 - ecc) - EARTH_RADIUS_KM
    apogee_alt = sma_km * (1.0 + ecc) - EARTH_RADIUS_KM
    return (perigee_alt, apogee_alt)


def query_paginated(token: str, filter_expr: str | None, max_pages: int,
                    label: str = "objects") -> tuple[list[dict], dict[str, dict]]:
    """Paginate /objects with optional filter, embedding initialOrbits inline.

    filter_expr=None pulls EVERYTHING on-orbit (active + decayed-not-purged
    + debris + rocket bodies). filter_expr="eq(active,true)" restricts to
    operational payloads.

    Returns (object_records, orbit_by_id_map).
    """
    all_records: list[dict] = []
    all_orbits: dict[str, dict] = {}
    total_pages_known: int | None = None
    total_records_known: int | None = None

    for page_num in range(1, max_pages + 1):
        params = {
            "page[size]": 100,
            "page[number]": page_num,
            "include": "initialOrbits",
        }
        if filter_expr:
            params["filter"] = filter_expr
        envelope = fetch_page_with_meta("objects", params, token)
        records = envelope.get("data", [])
        included = envelope.get("included", [])
        meta = envelope.get("meta", {})
        pag = meta.get("pagination", {})

        # Index orbits by id from the included array
        for inc in included:
            if inc.get("type") == "initialOrbit":
                all_orbits[inc.get("id")] = inc.get("attributes", {})

        # First-page-only metadata report
        if page_num == 1:
            total_pages_known = pag.get("totalPages")
            total_records_known = pag.get("totalRecords")
            print(f"    pagination meta:  totalPages = {total_pages_known}, "
                  f"totalRecords = {total_records_known}")
            if total_pages_known and total_pages_known < max_pages:
                print(f"    (will stop at page {total_pages_known}; "
                      f"fewer pages exist than requested {max_pages})")

        if not records:
            print(f"    page {page_num} empty - stopping")
            break

        print(f"    page {page_num:>3d}/{max_pages}  obj={len(records):>3d}  "
              f"orbits+={len(included):>3d}  "
              f"(running obj={len(all_records) + len(records):>5d}  "
              f"orbits={len(all_orbits):>5d})")
        all_records.extend(records)
        time.sleep(0.5)  # polite spacing under rate limit

        if total_pages_known and page_num >= total_pages_known:
            print(f"    reached last page ({total_pages_known}) - done")
            break

    return all_records, all_orbits


def query_active(token: str, max_pages: int = 1) -> tuple[list[dict], dict[str, dict]]:
    """Active payloads only (eq(active,true) filter)."""
    return query_paginated(token, "eq(active,true)", max_pages, label="active")


def query_full(token: str, max_pages: int = 1) -> tuple[list[dict], dict[str, dict]]:
    """Full catalog: payloads + rocket bodies + debris + everything DISCOSweb
    has on-orbit (no active filter). Roughly 2x the size of the active subset
    once debris is included.
    """
    return query_paginated(token, None, max_pages, label="full")


# Keep query_leo as alias for backwards compat (it now does the right thing)
query_leo = query_active


def flatten_records(records: list[dict], orbits: dict[str, dict] | None = None) -> list[dict]:
    """Pick the operationally useful fields for the Kessler simulator.

    If `orbits` is provided (orbit_id -> attrs map), attach each object's
    earliest/most-recent initial-orbit + computed perigee/apogee altitude
    and a 'regime' classification (LEO / MEO / GEO / HEO / unknown).
    """
    orbits = orbits or {}
    rows = []
    for rec in records:
        attrs = rec.get("attributes", {})
        rels = rec.get("relationships", {})

        # Pick the first associated initial-orbit ID (earliest by epoch heuristic
        # is fine for catalog-size proxy; precision orbit determination not
        # needed for the percolation-density estimate at Layer 4).
        orbit_links = (rels.get("initialOrbits") or {}).get("data", []) or []
        orbit_attrs: dict = {}
        if orbit_links:
            first_orbit_id = orbit_links[0].get("id")
            orbit_attrs = orbits.get(first_orbit_id, {})

        peri_km, apo_km = _orbit_to_perigee_apogee_km(orbit_attrs)

        regime = "unknown"
        if peri_km is not None and apo_km is not None:
            # Standard regime cuts (perigee-altitude based):
            #   LEO : peri < 2000 km
            #   MEO : 2000 <= peri < 35586 (geosynchronous altitude)
            #   GEO : 35586 km +/- 200 km, ecc ~ 0
            #   HEO : eccentric, apogee >> perigee
            if peri_km < 2000.0:
                regime = "LEO"
            elif 35586.0 - 200.0 <= peri_km <= 35586.0 + 200.0 and apo_km - peri_km < 200.0:
                regime = "GEO"
            elif apo_km - peri_km > 5000.0:
                regime = "HEO"
            else:
                regime = "MEO"

        rows.append({
            "id": rec.get("id"),
            "cosparId": attrs.get("cosparId"),
            "satno": attrs.get("satno"),
            "name": attrs.get("name"),
            "objectClass": attrs.get("objectClass"),
            "mass_kg": attrs.get("mass"),
            "diameter_m": attrs.get("diameter"),
            "length_m": attrs.get("length"),
            "span_m": attrs.get("span"),
            "shape": attrs.get("shape"),
            "xSectMin": attrs.get("xSectMin"),
            "xSectMax": attrs.get("xSectMax"),
            "xSectAvg": attrs.get("xSectAvg"),
            "depth_m": attrs.get("depth"),
            "active": attrs.get("active"),
            # orbit-derived
            "sma_m": orbit_attrs.get("sma"),
            "ecc": orbit_attrs.get("ecc"),
            "inc_deg": orbit_attrs.get("inc"),
            "epoch": orbit_attrs.get("epoch"),
            "perigee_alt_km": peri_km,
            "apogee_alt_km": apo_km,
            "regime": regime,
        })
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    import csv
    if not rows:
        print("    no records to write")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", choices=["basic", "leo", "full"], default="basic",
                    help="basic = smoke test; leo = active payloads only; "
                         "full = ALL on-orbit objects (active + debris + RBs)")
    ap.add_argument("--max-pages", type=int, default=1,
                    help="pagination limit for leo / full queries")
    args = ap.parse_args()

    print("=" * 80)
    print(" DISCOSweb Loader - ESA space-object catalog ".center(80))
    print("=" * 80)

    try:
        token = load_token()
        token_summary = (
            token[:8] + "..." + token[-4:] if len(token) > 16 else "(short)"
        )
        print(f"\n  Token loaded: {token_summary}  (length {len(token)})")
    except RuntimeError as e:
        print(f"\n  ERROR: {e}")
        return 1

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

    orbits: dict[str, dict] = {}
    if args.query == "basic":
        records = query_basic(token)
    elif args.query == "leo":
        records, orbits = query_leo(token, max_pages=args.max_pages)
    elif args.query == "full":
        records, orbits = query_full(token, max_pages=args.max_pages)
    else:
        print(f"  unknown query: {args.query}")
        return 1

    print(f"\n  Records returned: {len(records)}  (orbits indexed: {len(orbits)})")
    if not records:
        print("  (empty result - token may lack scope, or filter may have")
        print("   matched zero records. Check the API docs or rerun with --query basic.)")
        return 1

    # Save raw envelope (objects + orbits side by side) + flattened CSV
    raw_path = HERE / f"discosweb_{args.query}_{timestamp}.json"
    raw_path.write_text(
        json.dumps({"objects": records, "orbits": orbits}, indent=2),
        encoding="utf-8",
    )
    print(f"  Raw  -> {raw_path}")

    rows = flatten_records(records, orbits)
    csv_path = HERE / f"discosweb_{args.query}_summary.csv"
    write_csv(rows, csv_path)
    print(f"  CSV  -> {csv_path}")

    # Quick stats
    masses = [r["mass_kg"] for r in rows if r["mass_kg"] is not None]
    diams = [r["diameter_m"] for r in rows if r["diameter_m"] is not None]
    classes: dict[str, int] = {}
    regimes: dict[str, int] = {}
    leo_masses = []
    leo_xsect = []
    leo_perigee = []
    for r in rows:
        c = r.get("objectClass") or "(unknown)"
        classes[c] = classes.get(c, 0) + 1
        rg = r.get("regime") or "unknown"
        regimes[rg] = regimes.get(rg, 0) + 1
        if rg == "LEO":
            if r.get("mass_kg") is not None:
                leo_masses.append(r["mass_kg"])
            if r.get("xSectAvg") is not None:
                leo_xsect.append(r["xSectAvg"])
            if r.get("perigee_alt_km") is not None:
                leo_perigee.append(r["perigee_alt_km"])

    print()
    print(f"  Quick stats:")
    print(f"    n records:           {len(rows)}")
    print(f"    n with mass:         {len(masses)}  "
          f"({'min='+str(min(masses))+' max='+str(max(masses)) if masses else ''})")
    print(f"    n with diameter:     {len(diams)}  "
          f"({'min='+str(min(diams))+' max='+str(max(diams)) if diams else ''})")
    print(f"    object class breakdown:")
    for c, n in sorted(classes.items(), key=lambda kv: -kv[1]):
        print(f"      {c:>30s}: {n}")
    print(f"    orbital regime breakdown:")
    for c, n in sorted(regimes.items(), key=lambda kv: -kv[1]):
        print(f"      {c:>30s}: {n}")
    if leo_masses:
        leo_masses_sorted = sorted(leo_masses)
        n = len(leo_masses_sorted)
        print(f"    LEO mass distribution (kg, n={n}):")
        print(f"      min={leo_masses_sorted[0]:.1f}  "
              f"p25={leo_masses_sorted[n//4]:.1f}  "
              f"median={leo_masses_sorted[n//2]:.1f}  "
              f"p75={leo_masses_sorted[3*n//4]:.1f}  "
              f"max={leo_masses_sorted[-1]:.1f}")
    if leo_xsect:
        leo_xs_sorted = sorted(leo_xsect)
        n = len(leo_xs_sorted)
        print(f"    LEO xSectAvg (m^2, n={n}):")
        print(f"      min={leo_xs_sorted[0]:.3f}  "
              f"median={leo_xs_sorted[n//2]:.3f}  "
              f"max={leo_xs_sorted[-1]:.3f}")
    if leo_perigee:
        leo_p_sorted = sorted(leo_perigee)
        n = len(leo_p_sorted)
        print(f"    LEO perigee altitude (km, n={n}):")
        print(f"      min={leo_p_sorted[0]:.0f}  "
              f"median={leo_p_sorted[n//2]:.0f}  "
              f"max={leo_p_sorted[-1]:.0f}")
    print()
    print(f"  Token works. Catalog reachable. Layer 3 of Kessler build plan: unblocked.")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
