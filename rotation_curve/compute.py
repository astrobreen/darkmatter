# ---------------------------------------------------------------------------
# Title:        Kinematic Computation
# Author:       Dillon Breen
# Institute:    Dublin City University
# ---------------------------------------------------------------------------

# Pure-computation module 

import numpy as np
import pandas as pd
import astropy.units as u
import astropy.coordinates as coord

# Gravitational constant: pc (km/s)^2 / M_sun
G = 4.30091e-3

# Sun's Galactocentric position and velocity [Karim & Mamajek 2016]
# R_sun = 8300 pc, z_sun = 17 pc, v_sun = (12.4, 236.0, 7.7) km/s
_V_SUN = coord.CartesianDifferential([12.4, 236.0, 7.7] * u.km / u.s)
_GALACTOCENTRIC = coord.Galactocentric(
    galcen_distance=8.3 * u.kpc,
    galcen_v_sun=_V_SUN,
    z_sun=17 * u.pc,
)


# ---------------------------------------------------------------------------
# Coordinate transform
# ---------------------------------------------------------------------------

def gc(ra, dec, p, pmra, pmdec, rv):
    """Transform ICRS astrometry to Galactocentric Cartesian coordinates.

    Accepts pandas Series or array-like inputs.
    When called with error columns it performs first-order error propagation.

    Parameters
    ----------
    ra, dec   : right ascension / declination (degrees)
    p         : parallax (mas)
    pmra      : proper motion in RA cos(dec) (mas/yr)
    pmdec     : proper motion in Dec (mas/yr)
    rv        : radial velocity (km/s)

    Returns
    -------
    astropy Galactocentric object with .x, .y, .z, .v_x, .v_y, .v_z
    """
    icrs = coord.ICRS(
        ra=np.multiply(ra, u.degree),
        dec=np.multiply(dec, u.degree),
        distance=(list(p.values) * u.mas).to(u.pc, u.parallax()),
        pm_ra_cosdec=np.multiply(pmra, u.mas / u.yr),
        pm_dec=np.multiply(pmdec, u.mas / u.yr),
        radial_velocity=np.multiply(rv, u.km / u.s),
    )
    return icrs.transform_to(_GALACTOCENTRIC)


# ---------------------------------------------------------------------------
# Polynomial fitting with fixed points
# ---------------------------------------------------------------------------

def fix_points(n, x, y, xf, yf):
    """Fit a degree-n polynomial to (x, y) constrained to pass through (xf, yf).

    Uses a Lagrange-multiplier formulation to enforce the fixed points exactly.

    Returns
    -------
    params : array of n+1 polynomial coefficients (constant term first)
    """
    size = n + 1 + len(xf)
    mat = np.empty((size, size))
    vec = np.empty(size)

    x_n = x ** np.arange(2 * n + 1)[:, None]
    yx_n = np.sum(x_n[: n + 1] * y, axis=1)
    x_n = np.sum(x_n, axis=1)

    idx = np.arange(n + 1) + np.arange(n + 1)[:, None]
    mat[: n + 1, : n + 1] = np.take(x_n, idx)

    xf_n = xf ** np.arange(n + 1)[:, None]
    mat[: n + 1, n + 1 :] = xf_n / 2
    mat[n + 1 :, : n + 1] = xf_n.T
    mat[n + 1 :, n + 1 :] = 0

    vec[: n + 1] = yx_n
    vec[n + 1 :] = yf

    return np.linalg.solve(mat, vec)[: n + 1]


# ---------------------------------------------------------------------------
# Data loading and processing
# ---------------------------------------------------------------------------

def _add_positions(df):
    """Add Galactocentric (x, y, parsec) and propagated errors to df in place."""
    gal = gc(df["ra"], df["dec"], df["parallax"], df["pmra"], df["pmdec"], df["radial_velocity"])
    df["x"] = gal.x.value
    df["y"] = gal.y.value
    df["parsec"] = np.sqrt(df["x"] ** 2 + df["y"] ** 2)

    gal_err = gc(
        df["ra_error"], df["dec_error"], df["parallax_error"],
        df["pmra_error"], df["pmdec_error"], df["radial_velocity_error"],
    )
    df["x_error"] = gal_err.x.value
    df["y_error"] = gal_err.y.value
    df["parsec_error"] = np.sqrt(df["x_error"] ** 2 + df["y_error"] ** 2)


def gen_data(stellar_type):
    """Load Gaia CSV files for one spectral type, compute kinematics and masses.

    Reads SgrA-Type<X>.csv and OppSgrA-Type<X>.csv from the working directory.

    Returns
    -------
    DataFrame with columns: kiloparsec, kiloparsec_error, rotational_velocity,
    rotational_velocity_error, mass, mass_error (plus intermediate columns).
    """
    cols = [
        "ra", "ra_error", "dec", "dec_error",
        "parallax", "parallax_error",
        "pmra", "pmra_error", "pmdec", "pmdec_error",
        "radial_velocity", "radial_velocity_error",
    ]
    read_opts = dict(usecols=cols, sep=",", low_memory=False, on_bad_lines="skip")

    df1 = pd.read_csv(f"SgrA-Type{stellar_type}.csv",    **read_opts).astype(float)
    df2 = pd.read_csv(f"OppSgrA-Type{stellar_type}.csv", **read_opts).astype(float)

    _add_positions(df1)
    _add_positions(df2)

    df = pd.concat([df1, df2], sort=True, ignore_index=True)

    df["kiloparsec"]       = (np.array(df["parsec"])       * u.pc).to(u.kpc).value
    df["kiloparsec_error"] = (np.array(df["parsec_error"]) * u.pc).to(u.kpc).value

    # Drop rows without radial velocity (required for full velocity vector)
    df = df[np.isfinite(df["radial_velocity"])].copy()

    # 3-D Galactocentric velocity
    gal = gc(df["ra"], df["dec"], df["parallax"], df["pmra"], df["pmdec"], df["radial_velocity"])
    df["v_x"] = gal.v_x.value
    df["v_y"] = gal.v_y.value
    df["v_z"] = gal.v_z.value
    df["rotational_velocity"] = np.sqrt(df["v_x"] ** 2 + df["v_y"] ** 2 + df["v_z"] ** 2)

    gal_err = gc(
        df["ra_error"], df["dec_error"], df["parallax_error"],
        df["pmra_error"], df["pmdec_error"], df["radial_velocity_error"],
    )
    df["v_x_error"] = gal_err.v_x.value
    df["v_y_error"] = gal_err.v_y.value
    df["v_z_error"] = gal_err.v_z.value
    df["rotational_velocity_error"] = np.sqrt(
        df["v_x_error"] ** 2 + df["v_y_error"] ** 2 + df["v_z_error"] ** 2
    )

    # Enclosed mass: M = R V^2 / G  (R in pc, V in km/s → M in M_sun)
    df["mass"] = (df["parsec"] * df["rotational_velocity"] ** 2) / G

    # Mass error via propagation: dM/M = sqrt((dR/R)^2 + (2 dV/V)^2)
    vel_sq_err = 2 * df["rotational_velocity"] * df["rotational_velocity_error"]
    df["mass_error"] = (
        np.abs(df["parsec"] * df["rotational_velocity"] ** 2)
        * np.sqrt(
            (df["parsec_error"] / df["parsec"]) ** 2
            + (vel_sq_err / df["rotational_velocity"] ** 2) ** 2
        )
    ) / G

    df = df.reset_index(drop=True)
    print(df[["kiloparsec", "kiloparsec_error", "rotational_velocity",
              "rotational_velocity_error", "mass", "mass_error"]])
    return df


def gen_types(stellar_type, sample):
    """Return a (optionally sampled) DataFrame for the selected stellar type(s).

    Parameters
    ----------
    stellar_type : one of 'O','B','A','F','G','K','M', or 'all'
    sample       : int (number of rows to sample) or 'all'
    """
    types = ["O", "B", "A", "F", "G", "K", "M"] if stellar_type == "all" else [stellar_type]
    df = pd.concat([gen_data(t) for t in types], ignore_index=True)

    if isinstance(sample, int) and 0 <= sample < len(df):
        df = df.sample(sample)
    elif sample == "all":
        print(f"Total data points: {len(df)} sources.")
    else:
        print(f"Sample size out of range — {len(df)} data points available.")

    return df


def gen_bins(frame):
    """Bin stars by Galactocentric distance (0–16 kpc) and normalise counts.

    Returns a DataFrame with columns: bins, bin_count, bins_error, count_error.
    """
    bin_range = np.linspace(0, 16, 20)
    frame = frame.copy()
    frame["bins"]       = pd.cut(frame["kiloparsec"],       bin_range)
    frame["bins_error"] = pd.cut(frame["kiloparsec_error"], bin_range)

    count = (frame.groupby("bins",       observed=False)["kiloparsec"]
                  .count().reset_index(name="bin_count"))
    error = (frame.groupby("bins_error", observed=False)["kiloparsec_error"]
                  .count().reset_index(name="count_error"))

    mid = [(a + b) / 2 for a, b in zip(bin_range[:-1], bin_range[1:])]
    count["bins"]       = mid
    error["bins_error"] = mid

    # Scale bin position errors by the mean parallax relative uncertainty
    mean_par_rel_err     = (frame["parallax_error"] / frame["parallax"]).mean()
    error["bins_error"] *= mean_par_rel_err

    # Normalise to arbitrary density units
    count["bin_count"]   /= count["bin_count"].max()
    error["count_error"] /= error["count_error"].max()

    result = pd.merge(count, error, how="left", left_index=True, right_index=True)
    print(result)
    return result
