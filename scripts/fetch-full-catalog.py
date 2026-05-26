#!/usr/bin/env python3
"""
fetch-full-catalog.py
=====================

Downloads the Baumgardt & Vasiliev Milky Way globular cluster catalog and
converts it to the CSV format used by The Furniture of the Galaxy tool.

Source pages:
    https://people.smp.uq.edu.au/HolgerBaumgardt/globular/

The catalog is split across two files:
    orbits_table.txt     — 6D phase space (RA, Dec, distance, PMs, RV)
    combined_table.txt   — structural parameters including mass and tidal radius

The structural-parameter file was previously named parameter.txt; as of 2026
Baumgardt's site uses combined_table.txt, which folds positional and structural
data into one file. We continue to read the kinematic data from orbits_table.txt
and only the mass + tidal radius from combined_table.txt, joined on cluster name.

Cluster names in both files use underscores (e.g. "NGC_104"). The script
normalizes them to space-separated form ("NGC 104") for consistency with the
bundled data/clusters.csv and with the per-cluster description lookup in
index.html.

This script downloads both files, joins them on cluster name, and emits a CSV
matching the format read by index.html and data/clusters.csv.

Requires:
    Python 3.7+ (stdlib only — no pip install)

Usage:
    python3 scripts/fetch-full-catalog.py > data/clusters.csv
    # or, with explicit output path:
    python3 scripts/fetch-full-catalog.py --out data/clusters_full.csv

If the Baumgardt site is unreachable or the table format has shifted, the
script prints diagnostics to stderr and exits non-zero. You can also save
the raw downloaded files to inspect them manually:

    python3 scripts/fetch-full-catalog.py --save-raw raw/

WARNING
-------
Table column layouts on the Baumgardt site occasionally change between
catalog versions. If this script produces nonsense values or skips most
clusters, inspect the first few lines of the downloaded files (--save-raw)
and adjust the COLUMN_GUESS dictionaries below.
"""

import argparse
import sys
import os
import re
import urllib.request
import urllib.error
from typing import Dict, List, Optional

ORBITS_URL = "https://people.smp.uq.edu.au/HolgerBaumgardt/globular/orbits_table.txt"
PARAMS_URL = "https://people.smp.uq.edu.au/HolgerBaumgardt/globular/combined_table.txt"

# -----------------------------------------------------------------------------
# Column guesses.
#
# Baumgardt's tables are whitespace-delimited with a header line that starts
# with '#'. We identify columns by header keywords rather than fixed positions,
# which is robust to layout shifts. The FALLBACK positions below are used
# only when header detection fails.
#
# As of writing, combined_table.txt has 37 columns starting:
#   Cluster RA DEC R_Sun DRSun R_GC DRGC N_RV N_PM Mass DM V Delta_V
#   M/L_V DM/L rc rh,l rh,m rt rho_c rho_h,m ...
# We only extract: Cluster (name), Mass, rt (tidal radius).
# -----------------------------------------------------------------------------

ORBITS_HEADER_KEYS = {
    "name":  ["Cluster", "Name", "ID"],
    "ra":    ["RA"],
    "dec":   ["DEC", "Dec"],
    "dist":  ["R_Sun", "RSun", "Rsun", "Dist"],
    "pmra":  ["mualpha", "muRA", "mu_alpha", "PM_RA", "pmRA"],
    "pmdec": ["mudelta", "muDec", "mu_delta", "PM_Dec", "pmDec"],
    "rv":    ["<RV>", "RV", "Vlos", "Vr"],
}

ORBITS_FALLBACK_COLS = {
    "name": 0, "ra": 1, "dec": 2, "dist": 3, "pmra": 5, "pmdec": 7, "rv": 9,
}

PARAMS_HEADER_KEYS = {
    "name":   ["Cluster", "Name", "ID"],
    "mass":   ["Mass", "M_total", "Mtot"],
    "rtidal": ["rt", "r_t", "rtidal", "r_tidal"],
}

# Positions in combined_table.txt data rows (no leading '#'):
#   0: Cluster, 9: Mass, 18: rt
PARAMS_FALLBACK_COLS = {
    "name": 0, "mass": 9, "rtidal": 18,
}


# -----------------------------------------------------------------------------
# Networking
# -----------------------------------------------------------------------------

def fetch(url: str, timeout: int = 60) -> str:
    """Download a URL with a polite User-Agent."""
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "furniture-of-the-galaxy/1.1 (+https://github.com/catpea/cluster-encounters)"},
    )
    sys.stderr.write(f"Fetching {url} ...\n")
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        data = resp.read()
    sys.stderr.write(f"  got {len(data)} bytes\n")
    return data.decode("utf-8", errors="replace")


# -----------------------------------------------------------------------------
# Parsing — header-driven, with fallback to fixed positions
# -----------------------------------------------------------------------------

_HEADER_HASH_RE = re.compile(r"^\s*#+\s*")


def find_header_line(lines: List[str]) -> Optional[int]:
    """Locate the header row (the line listing column names)."""
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped:
            continue
        # Heuristic: header lines contain column-name keywords and few digits
        digits = sum(c.isdigit() for c in stripped)
        letters = sum(c.isalpha() for c in stripped)
        if letters > digits * 2 and any(
            k.lower() in stripped.lower() for k in ("RA", "DEC", "Cluster", "Name")
        ):
            return i
    return None


def resolve_columns(header_tokens: List[str], keyspec: Dict[str, List[str]]) -> Dict[str, Optional[int]]:
    """Map our logical fields to column indices, using header keywords."""
    out: Dict[str, Optional[int]] = {}
    lowered = [t.lower().strip("<>(){}[]") for t in header_tokens]
    for field, candidates in keyspec.items():
        out[field] = None
        for cand in candidates:
            cl = cand.lower()
            for i, tok in enumerate(lowered):
                if tok == cl or tok.startswith(cl):
                    out[field] = i
                    break
            if out[field] is not None:
                break
    return out


def parse_table(text: str, keyspec: Dict[str, List[str]], fallback: Dict[str, int],
                label: str) -> List[Dict[str, str]]:
    """Parse a whitespace-delimited table into a list of dicts keyed by field."""
    lines = text.splitlines()
    header_idx = find_header_line(lines)

    if header_idx is not None:
        # Strip the leading '#' (and any whitespace around it) so that header
        # column indices align with data row indices. Without this, the '#'
        # would occupy index 0 in the header but not in data rows, throwing
        # every column off by one.
        raw_header = lines[header_idx]
        header_line = _HEADER_HASH_RE.sub("", raw_header)
        header = header_line.split()
        cols = resolve_columns(header, keyspec)
        sys.stderr.write(f"[{label}] Detected header at line {header_idx + 1}: " +
                         ", ".join(f"{k}={v}" for k, v in cols.items()) + "\n")
        data_start = header_idx + 1
    else:
        sys.stderr.write(f"[{label}] WARNING: could not find header — using fallback column positions\n")
        cols = dict(fallback)
        data_start = 0

    # If header resolution missed required fields, fall back for those
    for field, idx in cols.items():
        if idx is None:
            cols[field] = fallback.get(field)
            sys.stderr.write(f"[{label}]   field '{field}' not found in header, using fallback col {cols[field]}\n")

    rows: List[Dict[str, str]] = []
    for raw in lines[data_start:]:
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith("-") or line.startswith("="):
            continue
        # Cluster names in Baumgardt tables use underscores (NGC_104) or
        # hyphens (2MASS-GC01); they do not contain spaces. Plain whitespace
        # split is correct.
        tokens = line.split()
        if len(tokens) < 3:
            continue

        row = {}
        try:
            for field, idx in cols.items():
                if idx is None or idx >= len(tokens):
                    row[field] = ""
                else:
                    row[field] = tokens[idx]
            rows.append(row)
        except (IndexError, ValueError):
            continue

    sys.stderr.write(f"[{label}] Parsed {len(rows)} rows\n")
    return rows


# -----------------------------------------------------------------------------
# Sanity check & coercion
# -----------------------------------------------------------------------------

def coerce_float(s: str) -> Optional[float]:
    if not s:
        return None
    try:
        return float(s.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def normalize_name(raw: str) -> str:
    """Convert Baumgardt-style 'NGC_104' to display form 'NGC 104'."""
    return raw.strip().replace("_", " ")


def sanity_check_orbit(row: Dict[str, str]) -> bool:
    ra = coerce_float(row.get("ra", ""))
    dec = coerce_float(row.get("dec", ""))
    d = coerce_float(row.get("dist", ""))
    if ra is None or dec is None or d is None:
        return False
    if not (0 <= ra <= 360):
        return False
    if not (-90 <= dec <= 90):
        return False
    if not (0 < d < 500):  # kpc; further than this is extragalactic
        return False
    return True


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="-", help="Output path (default: stdout)")
    ap.add_argument("--save-raw", metavar="DIR", help="Save raw downloaded tables to DIR for inspection")
    args = ap.parse_args()

    try:
        orbits_text = fetch(ORBITS_URL)
        params_text = fetch(PARAMS_URL)
    except urllib.error.URLError as e:
        sys.stderr.write(f"ERROR: could not reach the Baumgardt catalog server: {e}\n")
        sys.stderr.write("       Try again later, or save the files manually and adapt the script.\n")
        return 2

    if args.save_raw:
        os.makedirs(args.save_raw, exist_ok=True)
        with open(os.path.join(args.save_raw, "orbits_table.txt"), "w") as f:
            f.write(orbits_text)
        with open(os.path.join(args.save_raw, "combined_table.txt"), "w") as f:
            f.write(params_text)
        sys.stderr.write(f"Raw files saved to {args.save_raw}/\n")

    orbits = parse_table(orbits_text, ORBITS_HEADER_KEYS, ORBITS_FALLBACK_COLS, "orbits")
    params = parse_table(params_text, PARAMS_HEADER_KEYS, PARAMS_FALLBACK_COLS, "params")

    # Build a lookup of mass / tidal radius by normalized cluster name
    params_lookup: Dict[str, Dict[str, str]] = {}
    for p in params:
        name = normalize_name(p.get("name", ""))
        if name:
            params_lookup[name] = p

    # Emit CSV
    out_stream = sys.stdout if args.out == "-" else open(args.out, "w")
    try:
        out_stream.write("# Full Milky Way globular cluster catalog from Baumgardt & Vasiliev.\n")
        out_stream.write("# Downloaded by scripts/fetch-full-catalog.py.\n")
        out_stream.write("# name, RA(deg), Dec(deg), distance(kpc), pmRA*(mas/yr), pmDec(mas/yr), RV(km/s), mass(M_sun), tidal_radius(pc)\n")
        emitted = 0
        skipped = 0
        for row in orbits:
            if not sanity_check_orbit(row):
                skipped += 1
                continue
            name = normalize_name(row.get("name", ""))
            p = params_lookup.get(name, {})

            mass = coerce_float(p.get("mass", ""))
            rtid = coerce_float(p.get("rtidal", ""))

            mass_str = f"{mass:.3e}" if mass else ""
            rtid_str = f"{rtid:.1f}" if rtid else ""

            out_stream.write(",".join([
                name,
                row.get("ra", ""),
                row.get("dec", ""),
                row.get("dist", ""),
                row.get("pmra", ""),
                row.get("pmdec", ""),
                row.get("rv", ""),
                mass_str,
                rtid_str,
            ]) + "\n")
            emitted += 1

        sys.stderr.write(f"Wrote {emitted} clusters ({skipped} rejected by sanity check)\n")
        if emitted < 50:
            sys.stderr.write("WARNING: fewer than 50 clusters parsed. The table format may have shifted.\n")
            sys.stderr.write("         Run with --save-raw raw/ and inspect the downloaded files.\n")
            return 1
        return 0
    finally:
        if out_stream is not sys.stdout:
            out_stream.close()


if __name__ == "__main__":
    sys.exit(main())
