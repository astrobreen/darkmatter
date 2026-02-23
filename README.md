# Dark Matter in Our Galaxy

**Author:** Dillon Breen (15719601) — Dublin City University, PHA4
**Supervisor:** Dr. Eamonn Cunningham
**Module:** PS451 — Final Year Project (April 2019)

---

## Overview

This project estimates the dark matter mass distribution of the Milky Way by measuring the discrepancy between the observed (flat) stellar rotation curve and the Keplerian decline predicted by visible matter alone. Data are drawn from the ESA Gaia DR2 catalogue.

The pipeline has three stages:

1. **Query** — cylindrical slices of the Milky Way are queried from Gaia DR2 using ADQL, in two opposing directions along the line of sight to Sgr A*
2. **Compute** — stellar positions and velocities are transformed from ICRS to the Galactocentric frame; enclosed masses are derived from M(R) = RV²/G
3. **Visualise** — three analysis plots are produced: stellar density distribution, rotation curve, and dark matter mass distribution

---

## Project Structure

```
rotation_curve/
    __init__.py
    query.py        Step 1 — Gaia DR2 data acquisition
    compute.py      Step 2 — coordinate transforms and kinematic analysis
    visualise.py    Step 3 — plots (imports compute.py)
cli.py              Command-line entry point
pyproject.toml      Project metadata and dependencies (uv)
README.md
.gitignore
```

---

## Requirements

- Python ≥ 3.9
- [uv](https://docs.astral.sh/uv/) (recommended) **or** pip

---

## Installation

```bash
# Install uv if needed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Install project dependencies
uv sync
```

Or with pip:

```bash
pip install astropy astroquery matplotlib numpy pandas
```

---

## Usage

```bash
uv run dm --help
```

### Step 1 — Acquire data

```bash
uv run dm query --type G --interval 20 --range 8000 --angle 4
uv run dm query --type all --interval 20 --range 8000 --angle 4
```

| Flag | Thesis value | Description |
|------|-------------|-------------|
| `--type` | `G` (or `all`) | Spectral type: O B A F G K M or all |
| `--interval` | `20` | Shell thickness (pc). Thesis used 20 pc |
| `--range` | `8000` | Outer radius (pc). Observations extend to 8 kpc |
| `--angle` | `4` | Cone half-angle (degrees). Gives ~280 pc FoV at 8 kpc, matching the thin disc scale height (~400 pc) |

Queries are run towards Sgr A* (RA=266.417°, Dec=−29.008°) and its antipode. Each spectral type produces two CSV files: `SgrA-Type<X>.csv` and `OppSgrA-Type<X>.csv`.

> **Note:** Stellar types K and M require very large query volumes and were excluded from the thesis results due to computational constraints. Types O–G are sufficient for the rotation curve analysis.

### Step 2 — Analyse and plot

```bash
uv run dm analyse --type G --sample all
uv run dm analyse --type all --sample 5000
```

| Flag | Description |
|------|-------------|
| `--type` | Must match the downloaded data |
| `--sample` | Integer to randomly subsample, or `all` (default) |

The `analyse` command produces three sequential plots and prompts for:
- A galactic mass estimate (solar masses) — thesis used **5.17 × 10¹⁰ M☉** (Licquia & Newman 2014, within 25 kpc)
- A reference rotational velocity (km/s) — thesis used **236 km/s**

---

## Key Physical Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| G (astronomical) | 4.30091 × 10⁻³ pc (km/s)² M☉⁻¹ | Standard |
| R☉ (galactocentric distance) | 8300 pc (8.3 kpc) | Karim & Mamajek (2016) |
| z☉ (height above mid-plane) | 17 pc | Karim & Mamajek (2016) |
| V☉ (solar circular velocity) | 236 km/s | Used in coordinate transform |
| Solar peculiar motion (U, W) | 12.4, 7.7 km/s | Used in coordinate transform |
| Analysis radius | 0 – 16 kpc (20 bins) | Thesis §4.1 |
| Density model: bulge radius r_b | 1.5 kpc | Thesis conclusion |
| Density model: disc length l_c | 2.0 kpc | Thesis conclusion |
| Stellar mass within 25 kpc | 5.17 ± 1.11 × 10¹⁰ M☉ | Licquia & Newman (2014) |

---

## Results Summary

| Analysis | Result |
|----------|--------|
| Dark matter rotational velocity contribution | 100.07 ± 59.89 km/s within 16 kpc |
| Dark matter mass within 25 kpc | 1.08 ± 0.58 × 10¹¹ M☉ |
| Total virial mass within 25 kpc | 1.59 ± 0.58 × 10¹¹ M☉ |

These results are consistent with Sofue (2011) who found 8.87 ± 1.37 × 10¹⁰ M☉ for the dark matter halo within 20 kpc, and Posti & Helmi (2018) who estimated 1.37 ± 0.17 × 10¹¹ M☉ within 20 kpc.

---

## Method

### Query methodology (§3.1.2)

Gaia DR2 is queried in concentric cylindrical shells using ADQL cone queries. Each shell has a parallax interval of 20 pc and a field of view of ~280 pc radius at 8 kpc distance (4° cone half-angle), chosen to encompass the thin galactic disc scale height (~400 pc). The total cylindrical observation region extends 8 kpc in each direction along the Sgr A* line of sight. Volume error from this approach is less than 0.25%.

### Coordinate transform (§3.1.3)

Stellar positions and proper motions are transformed from ICRS (heliocentric) to Galactocentric Cartesian coordinates using Astropy, following the rotation matrix formulation of Eadie, Springford & Harris (2017). The Sun is placed at R☉ = 8300 pc, z☉ = 17 pc with circular velocity 236 km/s.

### Mass estimation (§3.5)

Enclosed mass is derived from the circular velocity approximation:

```
M(R) = R · V²(R) / G
```

where R is in pc, V in km/s, and G = 4.30091 × 10⁻³ pc (km/s)² M☉⁻¹. The dark matter contribution is the difference between the observed mass distribution and the Keplerian prediction using the known stellar mass.

---

## Data Source

ESA Gaia DR2 catalogue (reference epoch J2015.5, ICRS), accessed via [astroquery](https://astroquery.readthedocs.io/):

> Gaia Collaboration, Brown A.G.A., et al. (2018). *Gaia Data Release 2: Summary of the contents and survey properties.* A&A 616, A1.
