"""Probe DISCOSweb to discover the schema we need for LEO classification.

Three questions:
  1. What attributes does an /initial-orbits record carry?
  2. Does filter syntax work on those attributes (e.g., lt(perigee,X))?
  3. Does /objects?include=initialOrbits embed orbits inline?
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from discosweb_loader import load_token, fetch_page_with_meta


def probe_initial_orbits_schema(token: str) -> None:
    print("\n--- 1. /initial-orbits — sample attributes ---")
    env = fetch_page_with_meta("initial-orbits", {"page[size]": 3}, token)
    pag = env.get("meta", {}).get("pagination", {})
    print(f"  totalPages={pag.get('totalPages')}  totalRecords={pag.get('totalRecords')}")
    for rec in env.get("data", [])[:2]:
        print(f"\n  id={rec.get('id')}  type={rec.get('type')}")
        print(f"  attributes:")
        for k, v in (rec.get("attributes") or {}).items():
            print(f"    {k:>20s} = {v!r}")
        print(f"  relationships keys: {list((rec.get('relationships') or {}).keys())}")


def probe_filter_on_orbits(token: str) -> None:
    print("\n--- 2. filter syntax on /initial-orbits ---")
    # Probe LEO-shaped filter expressions; report which ones the API accepts.
    candidates = [
        ("lt(sma,8378)", "sma in km, LEO cutoff ~2000km altitude"),
        ("lt(perigee,8378)", "perigee in km, LEO cutoff"),
        ("lt(aPer,8378)", "aPer = altitude of perigee?"),
        ("lt(altPerigee,2000)", "altitude of perigee in km"),
    ]
    import requests
    headers = {"Authorization": f"Bearer {token}", "DiscosWeb-Api-Version": "2"}
    for filt, desc in candidates:
        r = requests.get(
            "https://discosweb.esoc.esa.int/api/initial-orbits",
            headers=headers,
            params={"filter": filt, "page[size]": 1},
            timeout=30,
        )
        if r.ok:
            n = len(r.json().get("data", []))
            tot = r.json().get("meta", {}).get("pagination", {}).get("totalRecords")
            print(f"  OK  filter={filt:<28s}  -> records returned (page1)={n}  total={tot}  ({desc})")
        else:
            err = r.json().get("errors", [{}])[0]
            print(f"  ERR filter={filt:<28s}  -> {err.get('code')} {err.get('meta', {})}  ({desc})")


def probe_include(token: str) -> None:
    print("\n--- 3. /objects?include=initialOrbits — does it embed inline? ---")
    env = fetch_page_with_meta(
        "objects",
        {"page[size]": 3, "filter": "eq(active,true)", "include": "initialOrbits"},
        token,
    )
    print(f"  data records:  {len(env.get('data', []))}")
    print(f"  included recs: {len(env.get('included', []))}")
    if env.get("included"):
        first_inc = env["included"][0]
        print(f"  sample included resource type: {first_inc.get('type')}")
        if first_inc.get("type") == "InitialOrbit" or "orbit" in str(first_inc.get("type", "")).lower():
            print(f"  attributes:")
            for k, v in (first_inc.get("attributes") or {}).items():
                print(f"    {k:>20s} = {v!r}")


def main() -> int:
    token = load_token()
    print(f"Token loaded (len={len(token)})")
    probe_initial_orbits_schema(token)
    probe_filter_on_orbits(token)
    probe_include(token)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
