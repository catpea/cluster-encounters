# The Furniture of the Galaxy

> *"We thought we were searching for life among the stars. We were the search.
> The stars we were looking at had been answering for ten billion years."*

A browser-based tool that integrates the Sun's orbit and the orbits of known Milky Way globular clusters backward through galactic time, finds the closest historical approaches, and tests whether any of those encounters coincide with major events in Earth's biological and geological record.

It is a **hypothesis-testing instrument, not a confirmation of anything.** The list of candidate matches it produces is where the looking starts.

**Live demo:** https://catpea.github.io/cluster-encounters/
**Source:** https://github.com/catpea/cluster-encounters

---

## The Question

There are about 150 globular clusters orbiting the Milky Way. Each is a tight ball of 10⁵–10⁶ ancient stars. They are 12–13 billion years old — *older than the galactic disk itself.* They do not live in the disk; they orbit the galactic center on highly inclined, often nearly perpendicular paths, plunging through the plane on cycles of hundreds of millions of years.

The Sun, meanwhile, oscillates vertically through the disk every ~30 Myr while completing its 225 Myr orbit around the galactic center. Two or three times each galactic year, it passes near enough to a globular cluster to feel its gravity.

**The question this tool asks:** do those close approaches correlate with terrestrial events?

The hypothesised causal chain is gravitational, not anthropic:

1. A cluster's mass perturbs the Sun's Oort cloud (Hills mechanism)
2. The perturbation knocks long-period comets onto inward-bound orbits
3. Cometary infall delivers organic chemistry and impactors to the inner solar system
4. The biosphere responds — either *bloom* (radiation, recovery) or *cull* (extinction)

Modern Gaia-era 6D phase space for globular clusters makes backward integration over 200+ Myr quantitatively meaningful for the first time. This tool runs the integration and presents the result.

> ⚠️ **What this is not.** It is not evidence that panspermia happened. The cluster-encounter → comet shower → biological event causal chain has multiple uncertain steps. Backward orbital integration over 200+ Myr is sensitive to the assumed Galactic potential and to every digit in the input data. Treat any correlation as suggestive, never confirmatory.

---

## Quick Start

### In the browser (the easy way)

Visit https://catpea.github.io/cluster-encounters/ — that's it. The tool loads the bundled cluster catalog and event list automatically, and runs entirely in your browser.

### Run locally

```bash
git clone https://github.com/catpea/cluster-encounters.git
cd cluster-encounters
python3 -m http.server 8000
# then visit http://localhost:8000/
```

A local HTTP server is needed for the auto-loaded CSV files to work. If you just double-click `index.html` and open it via `file://`, the tool still runs — it falls back silently to the same data baked into the HTML as defaults.

A full run with the bundled 13-cluster sample and 8 events, integrated 300 Myr backward at 0.2 Myr resolution, finishes in under 100 ms.

---

## Data

### Bundled (loaded automatically)

The tool fetches two CSV files on page load:

- [`data/clusters.csv`](./data/clusters.csv) — 13 well-characterized clusters with approximate Gaia-era 6D phase space, masses, and tidal radii
- [`data/events.csv`](./data/events.csv) — 8 major Phanerozoic biological and geochemical events

Edit either file (or the textareas in the running tool) to customize the analysis. The CSV format is documented in headers of each file.

### Full catalog (~160 clusters)

For serious work, replace the bundled sample with the complete Baumgardt & Vasiliev catalog:

```bash
python3 scripts/fetch-full-catalog.py --out data/clusters.csv
```

This downloads the canonical orbits and structural-parameter tables from
[Holger Baumgardt's site](https://people.smp.uq.edu.au/HolgerBaumgardt/globular/), joins them on cluster name, and emits a CSV in the format the tool expects. Pure Python 3 stdlib — no pip install.

If you want to inspect the raw downloads (useful if the table format on the Baumgardt site has shifted since this was written):

```bash
python3 scripts/fetch-full-catalog.py --save-raw raw/ --out data/clusters.csv
```

The script writes diagnostics to stderr and exits non-zero if it parses fewer than 50 clusters — a sign the column layout needs adjustment. See `scripts/fetch-full-catalog.py` for the column-detection logic and how to override it.

---

## What's Being Computed

### Galactic Potential

Three-component axisymmetric model:

| Component | Form | Parameters |
|-----------|------|------------|
| Bulge | Hernquist | M = 1.2 × 10¹⁰ M☉, a = 0.5 kpc |
| Disk | Miyamoto–Nagai | M = 8.0 × 10¹⁰ M☉, a = 3.0 kpc, b = 0.3 kpc |
| Halo | Logarithmic | v_h = 200 km/s, r_h = 12 kpc |

Produces v_circ ≈ 228 km/s at the solar radius — within ~2% of the measured value.

### Coordinate Transformation

Each cluster's published `(RA, Dec, distance, μα*·cos δ, μδ, RV)` is converted to Galactocentric Cartesian phase space via:

1. ICRS Cartesian unit vector from `(RA, Dec, d)`
2. Rotation by the standard J2000 ICRS-to-Galactic matrix
3. Translation to the Galactocentric frame (Sun at `(−8.122, 0, 0.0208)` kpc)
4. Addition of the solar velocity `(11.1, 245.24, 7.25)` km/s

### Integration

Classical fourth-order Runge–Kutta on pre-allocated `Float64Array` buffers. Zero per-step allocation. The hot path reads more like a C inner loop than idiomatic JavaScript — within ~2× of ANSI C speed in V8/SpiderMonkey. A WebAssembly port is straightforward but unnecessary for the current workload.

### Per-Encounter Output

Five physical numbers come out of each cluster's integration:

| Quantity | Symbol | Units | Meaning |
|----------|--------|-------|---------|
| Minimum distance | b_min | kpc | Closest historical approach |
| Approach time | t_min | Ma | When the closest approach happened |
| Relative velocity | v_rel | km/s | At the moment of closest approach |
| Encounter timescale | τ_enc | Myr | e-folding time = b_min / v_rel |
| Tidal impulse | ΔV_Oort | m/s | Velocity kick to outer Oort comets |

The tidal impulse uses the standard impulse approximation:

```
ΔV_Oort ≈ 2 · G · M · a_Oort / (b_min² · v_rel)
```

with `a_Oort = 50,000 AU`. This is the right physical observable to rank candidate encounters by — it is what actually shakes comets loose. Outer Oort orbital velocities are roughly 100 m/s, so a tidal kick of even 1 m/s is meaningful.

---

## Input Formats

### `data/clusters.csv`

One cluster per line, comma-separated. Lines starting with `#` are ignored.

```
name, RA(deg), Dec(deg), distance(kpc), pmRA*(mas/yr), pmDec(mas/yr), RV(km/s), mass(M☉), tidal_radius(pc)
```

The last two columns are optional. Without mass, the tidal-impulse output shows `—` and ranking falls back to closest-approach distance. With tidal radius, candidate cards display `b_min` as a multiple of `r_tidal`.

Example:

```
NGC 5139, 201.697, -47.480, 5.4, -3.250, -6.755, 232.7, 3.55e6, 64
```

### `data/events.csv`

```
name, age(Ma), type
```

Type is `cull` (extinction, anoxic event) or `bloom` (radiation, recovery, climatic optimum).

```
K-Pg Mass Extinction, 66.0, cull
PETM (Paleocene-Eocene), 56.0, bloom
```

---

## Reading the Output

The results section has two views:

**Candidate cards.** One card per cluster whose closest approach landed within the match window of a real event. Each card shows the cluster's name, common name, brief description, the five encounter numbers, the matched event with `ΔT`, and a distance-vs-time chart with the matched event window shaded amber and the closest approach marked.

**Full ranked table.** Every cluster you supplied, sorted by tidal impulse (highest first). Matched rows are highlighted.

**What to look for.** A genuinely interesting candidate has all three of:

- `b_min < ~1 kpc`
- `ΔV_Oort` of order 1 m/s or more
- A matched event with small `ΔT`

A candidate with one of those but not the others is weaker. A candidate with `b_min > 2 kpc` is probably too distant for the comet-perturbation mechanism to matter.

---

## Honest Limits

The integrator is excellent. The science is hard. Several layers of uncertainty stack:

- **Backward integration is sensitive to the potential.** Switching from a logarithmic halo to NFW shifts individual cluster trajectories at the 5–15% level over 200 Myr. Re-run with alternative potentials before trusting any candidate.
- **The potential is assumed smooth and time-independent.** No spiral arms, no bar, no time-varying perturbations from past mergers. The Milky Way was not exactly today's Milky Way 250 Myr ago.
- **Input uncertainties are not propagated.** A single integration is one realization. A real analysis requires Monte Carlo sampling of each cluster's `(RA, Dec, d, μα*, μδ, RV)` from its published error ellipsoid, ideally 1000+ realizations per cluster.
- **Clusters are treated as point masses.** They have internal structure and tidal tails; close-approach physics within the cluster's tidal radius would require N-body treatment.
- **The cluster → comet → biology chain has its own large uncertainties.** Even a strong perturbation does not deterministically produce an extinction event 50 Myr later.

The tool's purpose is to surface candidates worth that level of follow-up, not to perform the follow-up.

---

## Data Sources

For accurate inputs, the canonical sources are:

- **Baumgardt & Vasiliev (2021)** — *Accurate distances, proper motions, and orbits of Milky Way globular clusters.* Gaia EDR3 6D phase space for ~160 clusters with covariance matrices.
  [people.smp.uq.edu.au/HolgerBaumgardt/globular/](https://people.smp.uq.edu.au/HolgerBaumgardt/globular/)
- **Baumgardt & Hilker (2018)** — Globular cluster masses from N-body fits to surface brightness profiles.
- **Harris (1996, 2010 edition)** — The classical Milky Way globular cluster catalogue. Less precise on proper motions than Baumgardt, but still standard for finder-chart purposes.
- **Schönrich, Binney & Dehnen (2010)** — Solar peculiar motion relative to the LSR.
- **GRAVITY Collaboration / Reid & Brunthaler** — Galactocentric distance `R₀` and angular velocity of the Sun.

For the geological record:

- **Geologic Time Scale 2020** (Gradstein et al.) — Stage and series boundary ages.
- **Paleobiology Database (PBDB)** — Extinction and radiation events.

---

## References

The tool's components rest on standard literature:

- Bovy (2015), *galpy: a Python library for galactic dynamics* — inspiration for MWPotential2014
- Bailer-Jones (2018), *Close encounters of the stellar kind from Gaia DR2* — backward-integration methodology
- Hills (1981), *Comet showers and the steady-state infall of comets from the Oort cloud* — perturbation mechanism
- Heisler & Tremaine (1986), *Influence of the galactic tidal field on the Oort comet cloud*
- Vasiliev (2019), *agama: action-based galaxy modelling architecture*

The framing — that globular clusters are senior structures of the galaxy whose orbits we periodically intersect — owes more to Carl Sagan and Jill Tarter than to any specific paper. The science is real. The wonder is real. The list of candidates is a list of places to look.

---

## Repository Structure

```
cluster-encounters/
├── index.html                          # The tool (entry point for GitHub Pages)
├── README.md                            # You are here
├── LICENSE                              # MIT
├── the_furniture_of_the_galaxy.md      # Popularizer-style essay
├── .nojekyll                            # Tell GitHub Pages to serve files as-is
├── data/
│   ├── clusters.csv                     # 13 sample clusters (auto-loaded)
│   └── events.csv                       # 8 Earth events (auto-loaded)
└── scripts/
    └── fetch-full-catalog.py            # Download Baumgardt's full ~160-cluster catalog
```

---

## GitHub Pages Deployment

This repo is pre-configured for GitHub Pages. To enable:

1. Push to GitHub: `git push origin main`
2. Repository **Settings → Pages**
3. **Source:** Deploy from a branch
4. **Branch:** `main`, folder `/ (root)`
5. Save. The site will be live at `https://catpea.github.io/cluster-encounters/` within a minute.

The `.nojekyll` file ensures GitHub Pages serves the `data/` and `scripts/` directories as-is rather than running them through Jekyll.

---

## Contributing

Useful directions, roughly ordered by leverage:

1. **Load the full Baumgardt & Vasiliev table.** Run `scripts/fetch-full-catalog.py` and commit the resulting `data/clusters.csv`. If the fetch script needs adjustment for the current table format, fix it and PR.
2. **Monte Carlo over input uncertainties.** Wrap the integration in an outer loop that samples each cluster's 6D phase space from its published covariance matrix 1000+ times. Report medians and percentiles for `t_min`, `b_min`, and `ΔV_Oort`.
3. **Alternative potentials.** Add MWPotential2014 (NFW halo + power-law bulge), McMillan (2017), or a barred potential. Agreement between potentials is much stronger signal than any single integration.
4. **Forward integration.** Show each cluster's continuing orbit forward in time too, so the user can see the encounter geometry from both sides.
5. **3D orbital visualization.** A top-down + edge-on view of the encounter geometry, ideally with WebGL.
6. **Periodicity analysis.** FFT or Lomb–Scargle on the encounter density along the Sun's orbit, to test the ~30 Myr disk-crossing hypothesis quantitatively.
7. **WebAssembly port.** Not necessary for current workloads; useful if full-catalog Monte Carlo gets serious.
8. **Cluster description database.** Expand the per-cluster blurbs in `CLUSTER_INFO` (inside `index.html`) to cover the full catalogue.

Pull requests welcome. Please keep the tool itself self-contained — no build step, no bundler, single HTML file.

---

## See Also

- [`the_furniture_of_the_galaxy.md`](./the_furniture_of_the_galaxy.md) — A popularizer-style essay explaining the underlying hypothesis in long form, for readers without an astronomy background.

---

## Citation

If this tool is useful in your work, please cite it as:

```bibtex
@software{furniture_of_the_galaxy,
  title  = {The Furniture of the Galaxy: A Tool for Solar–Globular Cluster Encounter Analysis},
  author = {catpea},
  year   = {2026},
  url    = {https://github.com/catpea/cluster-encounters},
  note   = {Browser-based backward orbital integration with tidal impulse ranking}
}
```

---

## License

MIT. See [`LICENSE`](./LICENSE).

---

## Acknowledgments

The Gaia mission produced the data that makes Gyr-scale backward integration of cluster orbits meaningful at all. Holger Baumgardt, Eugene Vasiliev, William Harris, and the long tradition of globular-cluster catalogers built the foundation.

The intellectual framing — that the old things are exactly where they have always been, and that we are the latecomers passing through their patient geometry — comes from Carl Sagan, Jill Tarter, and everyone who taught us that the question of life elsewhere deserves both rigour and wonder.

> *And we, who thought we were the ones arriving,
> are the long swing of a pendulum
> through a room where the furniture
> has not moved in ten billion years.*
