# ---------------------------------------------------------------------------
# Title:        Observation Data Query
# Author:       Dillon Breen
# Institute:    Dublin City University
# ---------------------------------------------------------------------------

# Step 1 of 2 — download Gaia DR2 observations before running analyse.
# Output CSV files are consumed by visualise.py.
#
# Usage (via CLI):
#   uv run dm query --type G --interval 100 --range 8000 --angle 5

from astroquery.gaia import Gaia
import astropy.units as u
import numpy as np
import pandas as pd
import warnings

# ---------------------------------------------------------------------------
# Stellar type absolute-magnitude bounds (Mv min, Mv max)
# ---------------------------------------------------------------------------

STELLAR_TYPES = {
    "O": (-4.5, -3.6),
    "B": (-3.3,  1.1),
    "A": ( 1.5,  2.4),
    "F": ( 3.0,  4.4),
    "G": ( 4.7,  5.6),
    "K": ( 6.0,  8.1),
    "M": ( 8.7, 14.4),
}

# Galactic centre (SgrA*) and its antipode in ICRS
RA_SGRA,  DEC_SGRA  =  266.41683, -29.00781
RA_OPP,   DEC_OPP   =   86.41683,  29.00781


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def fov(r, theta, h):
    """Return the angular diameter (degrees) of a shell at distance h.

    Parameters
    ----------
    r      : outer scan radius (pc)
    theta  : half-angle of the desired cone (degrees)
    h      : distance to the shell (pc)
    """
    radius = r * np.tan(np.radians(theta / 2))
    return np.round(2 * np.degrees(np.arctan(radius / h)), 6)


def query(radius_deg, ra, dec, g_min, g_max, par_min, par_max):
    """Execute an asynchronous Gaia DR2 cone query with magnitude/parallax cuts.

    Returns a pandas DataFrame of the results.
    """
    warnings.filterwarnings("ignore")
    adql = (
        f"SELECT * FROM gaiadr2.gaia_source "
        f"WHERE CONTAINS("
        f"  POINT('ICRS', gaiadr2.gaia_source.ra, gaiadr2.gaia_source.dec),"
        f"  CIRCLE('ICRS', {ra}, {dec}, {radius_deg})"
        f")=1 "
        f"AND phot_g_mean_mag BETWEEN {g_min} AND {g_max} "
        f"AND parallax BETWEEN {par_min} AND {par_max}"
    )
    job = Gaia.launch_job_async(adql, background=True)
    return job.get_results().to_pandas()


# ---------------------------------------------------------------------------
# Data generation
# ---------------------------------------------------------------------------

def gen_data(interval, r, degree, vmin, vmax, ra, dec, file, stellar_type):
    """Query Gaia DR2 in concentric distance shells and append results to a CSV.

    Parameters
    ----------
    interval     : shell thickness (pc)
    r            : outer observation radius (pc)
    degree       : cone half-angle (degrees)
    vmin, vmax   : absolute-magnitude bounds for the target stellar type
    ra, dec      : pointing direction (ICRS degrees)
    file         : output CSV filename
    stellar_type : label used in progress messages
    """
    # Towards SgrA* we query out to 8 kpc; in the opposite direction out to r
    limit = 8000 if file == f"SgrA-Type{stellar_type}.csv" else r
    values = np.arange(interval, limit + interval, interval)

    frames = []
    for angular_radius, dist in zip(fov(r, degree, values), values):
        # Parallax limits (near edge = larger parallax, far edge = smaller)
        par_far  = np.round((dist            * u.pc).to(u.mas, u.parallax()).value, 6)
        par_near = np.round(((dist + interval) * u.pc).to(u.mas, u.parallax()).value, 6)

        # Apparent magnitude from distance modulus: m = M + 5 log10(d/10)
        m_lo = np.round(vmin + 5 * np.log10(dist / 10), 6)
        m_hi = np.round(vmax + 5 * np.log10(dist / 10), 6)

        print(
            f"  deg={angular_radius}  par=[{par_far}, {par_near}]  "
            f"mag=[{m_lo}, {m_hi}]"
        )
        frames.append(query(angular_radius, ra, dec, m_lo, m_hi, par_far, par_near))

    stars = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    stars.to_csv(file, mode="a", index=False)
    print(f"  → Type {stellar_type} complete: {len(stars)} sources written to {file}")


def query_type(t, interval, r, degree):
    """Run forward and reverse queries for a single spectral type."""
    vmin, vmax = STELLAR_TYPES[t]
    print(f"\nQuerying stellar type {t}:")
    gen_data(interval, r, degree, vmin, vmax, RA_SGRA,  DEC_SGRA,  f"SgrA-Type{t}.csv",    t)
    gen_data(interval, r, degree, vmin, vmax, RA_OPP,   DEC_OPP,   f"OppSgrA-Type{t}.csv", t)


