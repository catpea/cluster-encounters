#!/usr/bin/env python3
"""
fetch-stellar-encounters.py
===========================

Downloads the Bailer-Jones (2022, ApJL 935 L9; arXiv:2207.06258) stellar
close-encounter catalog from his MPIA personal page, cross-matches against
Gaia DR3 for the proper-motion components that aren't in the combined table,
and emits a CSV in the format used by The Furniture of the Galaxy tool.

Source: https://bailer-jones.www3.mpia.de/stellar_encounters_gdr3/encounters_combinedtable.csv

DATA FLOW
---------
The Bailer-Jones combined table provides, per encountering star:

    source_id, tphmed, tph5, tph95   (perihelion time and 5/95% CI, in kyr)
    dphmed, dph5, dph95              (perihelion distance, in pc)
    vphmed, vph5, vph95              (relative velocity at perihelion, km/s)
    par, epar                        (parallax in mas, error)
    pm, epm                          (total proper motion magnitude only, mas/yr)
    rv, erv                          (radial velocity in km/s, error)
    abs_rv_to_vtrans                 (sanity ratio)
    phot_g_mean_mag, bp_rp, m_g      (apparent G, color, absolute G)
    ruwe                             (astrometric quality, < 1.4 is clean)
    astrometric_params_solved        (Gaia solution code)
    ipd_frac_multi_peak              (image multiplicity flag)
    rv_expected_sig_to_noise         (RV SNR)
    rv_nb_transits                   (number of RV measurements)
    l, b                             (Galactic longitude/latitude in degrees)

To build a row the tool can ingest we need RA, Dec, distance, pmRA*, pmDec,
RV, and mass. From the combined table:

    - distance (kpc)   = 1 / parallax_mas
    - RA, Dec (deg)    = inverse Galactic-to-ICRS rotation of (l, b)
    - RV (km/s)        = rv (directly)
    - mass (M_sun)     = piecewise mass-luminosity relation from m_g

The proper-motion COMPONENTS aren't in the combined table — only the magnitude.
For each source_id we therefore issue a Gaia DR3 TAP query to fetch pmra,
pmdec. The TAP queries are batched (50 sources per query) to stay under URL
length limits.

USAGE
-----
    # Default — fetch from canonical URL, cross-match, emit to stdout:
    python3 scripts/fetch-stellar-encounters.py > data/stellar-encounters.csv

    # Save to file with quality filter (default ruwe < 1.4):
    python3 scripts/fetch-stellar-encounters.py --out data/stellar-encounters.csv

    # Include all 61 stars, not just the high-confidence ~42:
    python3 scripts/fetch-stellar-encounters.py --no-quality-filter

    # Skip the Gaia cross-match (faster but loses proper-motion components):
    python3 scripts/fetch-stellar-encounters.py --no-gaia-crossmatch

    # Convert a previously-downloaded copy:
    python3 scripts/fetch-stellar-encounters.py --input encounters_combinedtable.csv

Requires:
    Python 3.7+ (stdlib only)
"""

import argparse
import sys
import os
import math
import urllib.request
import urllib.error
import urllib.parse
from typing import Dict, List, Optional

# Canonical references
REFERENCE = "Bailer-Jones 2022, ApJL 935, L9 (arXiv:2207.06258)"
DOI_LINK = "https://doi.org/10.3847/2041-8213/ac816a"
SOURCE_URL = "https://bailer-jones.www3.mpia.de/stellar_encounters_gdr3/encounters_combinedtable.csv"
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"

# ICRS-to-Galactic rotation matrix (J2000), same as in index.html.
# We transpose it below to get Galactic-to-ICRS.
R_ICRS_GAL = (
    (-0.054875560, -0.873437090, -0.483835016),
    ( 0.494109428, -0.444829630,  0.746982244),
    (-0.867666149, -0.198076373,  0.455983776),
)
# Transpose: rows become columns
R_GAL_ICRS = tuple(tuple(R_ICRS_GAL[r][c] for r in range(3)) for c in range(3))


# -----------------------------------------------------------------------------
# Networking
# -----------------------------------------------------------------------------

def fetch(url: str, timeout: int = 60) -> Optional[str]:
    try:
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "furniture-of-the-galaxy/1.1 (+https://github.com/catpea/cluster-encounters)"},
        )
        sys.stderr.write(f"Fetching {url}...\n")
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = resp.read()
        sys.stderr.write(f"  got {len(data)} bytes\n")
        return data.decode("utf-8", errors="replace")
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:
        sys.stderr.write(f"  failed: {e}\n")
        return None


def gaia_tap_crossmatch(source_ids: List[str], chunk_size: int = 50) -> Dict[str, Dict[str, str]]:
    """
    Cross-match Gaia DR3 source_ids against the main catalog via the ESA TAP
    sync service. Returns {source_id: {ra, dec, parallax, pmra, pmdec, rv, ...}}.
    """
    out: Dict[str, Dict[str, str]] = {}
    n_total = len(source_ids)
    n_done = 0

    for i in range(0, n_total, chunk_size):
        chunk = source_ids[i:i + chunk_size]
        id_list = ",".join(str(sid) for sid in chunk)
        adql = (
            "SELECT source_id, ra, dec, parallax, pmra, pmdec, "
            "radial_velocity, phot_g_mean_mag, bp_rp "
            "FROM gaiadr3.gaia_source "
            f"WHERE source_id IN ({id_list})"
        )
        params = urllib.parse.urlencode({
            "REQUEST": "doQuery",
            "LANG": "ADQL",
            "FORMAT": "csv",
            "QUERY": adql,
        })
        sys.stderr.write(f"  Gaia DR3 TAP query: {len(chunk)} sources ({n_done + len(chunk)}/{n_total})...\n")
        text = fetch(GAIA_TAP + "?" + params, timeout=90)
        if not text:
            sys.stderr.write("  WARNING: Gaia TAP query failed; continuing without these sources\n")
            n_done += len(chunk)
            continue

        # Parse CSV response from Gaia
        lines = text.strip().split("\n")
        if len(lines) < 2:
            n_done += len(chunk)
            continue
        header = [h.strip() for h in lines[0].split(",")]
        try:
            sid_col = header.index("source_id")
        except ValueError:
            n_done += len(chunk)
            continue
        for row in lines[1:]:
            cells = [c.strip() for c in row.split(",")]
            if len(cells) != len(header):
                continue
            sid = cells[sid_col]
            out[sid] = dict(zip(header, cells))
        n_done += len(chunk)

    sys.stderr.write(f"  Cross-match returned {len(out)} of {n_total} sources\n")
    return out


# -----------------------------------------------------------------------------
# Coordinate conversion
# -----------------------------------------------------------------------------

def galactic_to_icrs(l_deg: float, b_deg: float) -> tuple:
    """Convert Galactic (l, b) in degrees to ICRS (RA, Dec) in degrees."""
    l_rad = math.radians(l_deg)
    b_rad = math.radians(b_deg)
    cb = math.cos(b_rad)
    # Galactic unit vector
    gx = cb * math.cos(l_rad)
    gy = cb * math.sin(l_rad)
    gz = math.sin(b_rad)
    # Rotate Galactic → ICRS
    x = R_GAL_ICRS[0][0]*gx + R_GAL_ICRS[0][1]*gy + R_GAL_ICRS[0][2]*gz
    y = R_GAL_ICRS[1][0]*gx + R_GAL_ICRS[1][1]*gy + R_GAL_ICRS[1][2]*gz
    z = R_GAL_ICRS[2][0]*gx + R_GAL_ICRS[2][1]*gy + R_GAL_ICRS[2][2]*gz
    ra = math.degrees(math.atan2(y, x))
    if ra < 0:
        ra += 360.0
    dec = math.degrees(math.asin(max(-1.0, min(1.0, z))))
    return ra, dec


# -----------------------------------------------------------------------------
# Mass estimation from absolute G magnitude
# -----------------------------------------------------------------------------

def estimate_mass_from_mg(m_g: float) -> float:
    """
    Rough piecewise mass-luminosity relation from absolute Gaia G magnitude.
    Calibrated against MIST/PARSEC isochrones for ~solar-metallicity main
    sequence. Accurate to ~factor of 2; sufficient for ranking by tidal
    impulse, insufficient for precision dynamics.
    """
    if m_g is None or m_g != m_g:  # None or NaN
        return 0.5  # default M dwarf
    if m_g < 1.5:    return 2.5     # early A
    if m_g < 3.0:    return 1.5     # F
    if m_g < 4.5:    return 1.1     # G
    if m_g < 6.0:    return 0.9     # K0-K3
    if m_g < 8.0:    return 0.7     # K4-K7 (Gliese 710 falls here)
    if m_g < 10.0:   return 0.5     # M0-M3
    if m_g < 13.0:   return 0.25    # M4-M6
    if m_g < 17.0:   return 0.10    # late M / L
    return 0.08


# -----------------------------------------------------------------------------
# Parsing the combined table
# -----------------------------------------------------------------------------

EXPECTED_COLUMNS = [
    "source_id", "tphmed", "dphmed", "vphmed",
    "par", "epar", "pm", "epm", "rv", "erv",
    "m_g", "bp_rp", "ruwe", "l", "b",
]


def parse_combined_table(text: str) -> List[Dict[str, str]]:
    """Parse Bailer-Jones encounters_combinedtable.csv. Plain CSV with header row."""
    rows: List[Dict[str, str]] = []
    lines = [ln for ln in text.split("\n") if ln.strip()]
    if not lines:
        return rows

    header = [h.strip() for h in lines[0].split(",")]
    missing = [c for c in EXPECTED_COLUMNS if c not in header]
    if missing:
        sys.stderr.write(f"  WARNING: expected columns missing from header: {missing}\n")
        sys.stderr.write(f"  Header found: {header}\n")

    for line in lines[1:]:
        cells = [c.strip() for c in line.split(",")]
        if len(cells) != len(header):
            continue
        rows.append(dict(zip(header, cells)))

    sys.stderr.write(f"Parsed {len(rows)} encounter rows from combined table\n")
    return rows


def to_float(s: Optional[str]) -> Optional[float]:
    if s is None or not s or s.lower() in ("nan", "null", "--", "none"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


# -----------------------------------------------------------------------------
# Naming
# -----------------------------------------------------------------------------

WELL_KNOWN_NAMES = {
    "4270814637616488064": "Gliese 710 (HIP 89825)",
    "3133367744202841344": "WISE 0720-08 (Scholz's Star)",
    # The first few well-known encounters; everything else gets "Gaia DR3 {id}"
}


def name_for(source_id: str) -> str:
    if source_id in WELL_KNOWN_NAMES:
        return WELL_KNOWN_NAMES[source_id]
    return f"Gaia DR3 {source_id}"


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="-", help="Output path (default: stdout)")
    ap.add_argument("--input", help="Path to a previously-downloaded copy of encounters_combinedtable.csv")
    ap.add_argument("--no-gaia-crossmatch", action="store_true",
                    help="Skip Gaia DR3 cross-match (faster but no pmRA/pmDec; output unusable for orbit integration)")
    ap.add_argument("--no-quality-filter", action="store_true",
                    help="Include all 61 stars (default: filter ruwe < 1.4, ~42 high-confidence)")
    ap.add_argument("--ruwe-cut", type=float, default=1.4,
                    help="Max RUWE (Gaia astrometric quality flag) to include (default: 1.4)")
    args = ap.parse_args()

    sys.stderr.write(f"Source: {REFERENCE}\n")
    sys.stderr.write(f"DOI:    {DOI_LINK}\n")
    sys.stderr.write(f"URL:    {SOURCE_URL}\n\n")

    # Get the combined table
    if args.input:
        sys.stderr.write(f"Reading {args.input}...\n")
        try:
            with open(args.input, encoding="utf-8") as f:
                text = f.read()
        except OSError as e:
            sys.stderr.write(f"ERROR: could not read {args.input}: {e}\n")
            return 2
    else:
        text = fetch(SOURCE_URL)
        if text is None:
            sys.stderr.write(
                "\nERROR: could not reach the Bailer-Jones MPIA page.\n"
                "Download manually with a browser or curl:\n"
                f"    curl -O {SOURCE_URL}\n"
                "then re-run with --input encounters_combinedtable.csv\n"
            )
            return 3

    rows = parse_combined_table(text)
    if not rows:
        sys.stderr.write("ERROR: no rows parsed. Inspect the file format.\n")
        return 4

    # Quality filter (skip stars with bad astrometric solutions)
    if not args.no_quality_filter:
        before = len(rows)
        filtered = []
        for r in rows:
            ruwe = to_float(r.get("ruwe"))
            if ruwe is not None and ruwe < args.ruwe_cut:
                filtered.append(r)
        rows = filtered
        sys.stderr.write(f"Quality filter (ruwe < {args.ruwe_cut}): {len(rows)} of {before} pass\n")

    # Cross-match against Gaia DR3 for pmRA, pmDec components
    gaia_data: Dict[str, Dict[str, str]] = {}
    if not args.no_gaia_crossmatch:
        source_ids = [r["source_id"] for r in rows if r.get("source_id")]
        if source_ids:
            sys.stderr.write(f"\nCross-matching {len(source_ids)} sources against Gaia DR3...\n")
            gaia_data = gaia_tap_crossmatch(source_ids)

    # Emit CSV
    out_stream = sys.stdout if args.out == "-" else open(args.out, "w")
    try:
        out_stream.write("# Stellar close encounters from Bailer-Jones 2022 / Gaia DR3\n")
        out_stream.write(f"# Source: {REFERENCE}\n")
        out_stream.write(f"# Downloaded from: {SOURCE_URL}\n")
        out_stream.write("# Generated by scripts/fetch-stellar-encounters.py\n")
        out_stream.write("#\n")
        out_stream.write("# RECOMMENDED TOOL SETTINGS for stellar-encounter analysis:\n")
        out_stream.write("#   Integration time:  10 Ma     (encounters lie within +/- 6 Myr)\n")
        out_stream.write("#   Timestep:          0.005 Ma  (sub-Myr resolution required)\n")
        out_stream.write("#   Match window:      0.5 Ma\n")
        out_stream.write("#\n")
        out_stream.write("# Position from Galactic (l, b) via inverse J2000 rotation.\n")
        out_stream.write("# Distance from parallax: d_kpc = 1/parallax_mas.\n")
        out_stream.write("# pmRA*, pmDec from Gaia DR3 main catalog (cross-matched by source_id).\n")
        out_stream.write("# Mass from absolute G magnitude via piecewise M-L relation (factor-of-2 accuracy).\n")
        out_stream.write("#\n")
        out_stream.write("# name, RA(deg), Dec(deg), distance(kpc), pmRA*(mas/yr), pmDec(mas/yr), RV(km/s), mass(M_sun)\n")

        emitted = 0
        skipped = 0
        for r in rows:
            sid = r.get("source_id", "")
            par = to_float(r.get("par"))
            l_deg = to_float(r.get("l"))
            b_deg = to_float(r.get("b"))
            rv = to_float(r.get("rv"))
            m_g = to_float(r.get("m_g"))

            if not sid or par is None or par <= 0 or l_deg is None or b_deg is None:
                skipped += 1
                continue

            # Distance from parallax (in mas -> kpc directly)
            dist_kpc = 1.0 / par

            # Galactic (l, b) → ICRS (RA, Dec)
            ra, dec = galactic_to_icrs(l_deg, b_deg)

            # Proper motion components from Gaia cross-match
            gaia = gaia_data.get(sid, {})
            pmra = to_float(gaia.get("pmra"))
            pmdec = to_float(gaia.get("pmdec"))
            if pmra is None or pmdec is None:
                # Without components, orbit integration won't produce meaningful results.
                # Skip the row rather than emit a half-defined source.
                skipped += 1
                continue

            # Prefer Gaia DR3's RV; fall back to combined-table RV
            rv_use = to_float(gaia.get("radial_velocity")) or rv
            if rv_use is None:
                skipped += 1
                continue

            # Mass from absolute G magnitude
            mass = estimate_mass_from_mg(m_g)

            name = name_for(sid)
            out_stream.write(",".join([
                name,
                f"{ra:.4f}", f"{dec:.4f}", f"{dist_kpc:.5f}",
                f"{pmra:.3f}", f"{pmdec:.3f}", f"{rv_use:.2f}",
                f"{mass:.3f}",
            ]) + "\n")
            emitted += 1

        sys.stderr.write(f"\nWrote {emitted} stars (skipped {skipped} for incomplete data)\n")
        if emitted < 10:
            sys.stderr.write("WARNING: very few stars emitted. Check Gaia cross-match.\n")
            return 1
        return 0
    finally:
        if out_stream is not sys.stdout:
            out_stream.close()


if __name__ == "__main__":
    sys.exit(main())
