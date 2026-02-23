# ---------------------------------------------------------------------------
# Title:        Visualisation & Interactive Analysis
# Author:       Dillon Breen
# Institute:    Dublin City University
# ---------------------------------------------------------------------------
#
# Step 2 of 2 — compute kinematics and produce plots from the queried CSV files.
#
# Usage (via CLI):
#   uv run dm analyse --type G --sample 1000

import numpy as np
import matplotlib.pyplot as plt

from .compute import G, fix_points, gen_types, gen_bins


# ---------------------------------------------------------------------------
# Reference model
# ---------------------------------------------------------------------------

def ss_model():
    """Plot Solar System orbital velocities as a Keplerian reference."""
    planets    = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
    distances  = [0.39, 0.723, 1, 1.524, 5.203, 9.539, 19.18, 30.06]   # AU
    velocities = [47.36, 35.02, 29.78, 24.077, 13.07, 9.69, 6.81, 5.43]  # km/s

    # Keplerian curve (1 AU = 4.848e-6 pc; total mass ≈ 1.0014 M_sun)
    xc = np.linspace(0.01, 35, 100)
    yc = np.sqrt((G * 1.0014) / (4.84814e-6 * xc))
    plt.plot(xc, yc, color="blue", label="Keplerian model", zorder=1)

    for name, d, v in zip(planets, distances, velocities):
        plt.scatter(d, v, marker="o", s=10, color="black", zorder=2)
        plt.text(d + 0.5, v + 0.5, name, fontsize=10, zorder=3)

    plt.title("Solar System Rotation Curve")
    plt.xlabel("Distance ($AU$)")
    plt.ylabel("Orbital velocity ($km\\ s^{-1}$)")
    plt.xlim(0, 35)
    plt.ylim(0, 55)
    plt.legend(loc="upper right")
    plt.show()


# ---------------------------------------------------------------------------
# Stellar density distribution
# ---------------------------------------------------------------------------

def stellar_distribution_plot(frame, fit):
    """Plot binned stellar density with a constrained parametric model.

    Parameters
    ----------
    frame : output of compute.gen_bins()
    fit   : polynomial degree for the model fit
    """
    plt.errorbar(
        frame["bins"], frame["bin_count"],
        xerr=frame["bins_error"], yerr=frame["count_error"],
        fmt="o", ms=0.5, color="black", ecolor="grey", elinewidth=0.25,
        label="Gaia DR2", zorder=1,
    )

    # Parametric disc density profile S(R) ∝ l_c / sqrt((R - R_b)^2 + l_c^2)
    # r_b = 1.5 kpc (central bulge radius), l_c = 2 kpc (characteristic disk length) [thesis §5]
    s0, l_c, r_b = frame["bin_count"].max(), 2.0, 1.5
    xb = np.linspace(0, 16, 50)
    yb = (s0 * l_c) / np.sqrt((xb - r_b) ** 2 + l_c ** 2)

    # Constrain polynomial through three representative observed points
    xp = np.array([frame["bins"].iloc[0], frame["bins"].iloc[14], frame["bins"].iloc[-1]])
    yp = np.array([frame["bin_count"].iloc[0], frame["bin_count"].iloc[14], frame["bin_count"].iloc[-1]])
    params = fix_points(fit, xb, yb, xp, yp)
    ppoly  = np.polynomial.Polynomial(params)
    pdata  = np.linspace(0, 16, 50)
    plt.plot(pdata, ppoly(pdata), ms=1, color="blue", label="Best-fit density model", zorder=2)

    plt.title("Stellar Density Distribution")
    plt.xlabel("Galactic radius ($kpc$)")
    plt.ylabel("Density (arbitrary units)")
    plt.xlim(0, 16)
    plt.ylim(0, s0 * 1.1)
    plt.legend(loc="upper right")
    plt.show()


def stellar_distribution(stellar_type, sample, fit):
    """Load data, bin, and plot the stellar density distribution."""
    print(f"\nStellar density distribution — type: {stellar_type}, sample: {sample}")
    stellar_distribution_plot(gen_bins(gen_types(stellar_type, sample)), fit)


# ---------------------------------------------------------------------------
# Rotation curve
# ---------------------------------------------------------------------------

def rotational_curvature_plot(stellar_type, sample, fit):
    """Plot the Milky Way rotation curve, Keplerian prediction, and dark matter excess.

    Prompts for a trial galactic mass (M_sun) used to set the Keplerian baseline.
    """
    try:
        gal_mass = float(input("Enter a galactic mass estimate (solar masses): "))
    except ValueError:
        print("Invalid input — please enter a numeric value.")
        return rotational_curvature_plot(stellar_type, sample, fit)

    df    = gen_types(stellar_type, sample)
    pdata = np.linspace(0, 16, 50)

    # Observational data
    plt.errorbar(
        df["kiloparsec"], df["rotational_velocity"],
        fmt="o", ms=0.5, color="black",
        xerr=df["kiloparsec_error"], yerr=df["rotational_velocity_error"],
        elinewidth=0.25, ecolor="grey", label="Gaia DR2", zorder=1,
    )

    # Best-fit rotation curve (constrained through origin)
    params = fix_points(fit, df["kiloparsec"].values, df["rotational_velocity"].values,
                        np.array([0]), np.array([0]))
    ppoly = np.polynomial.Polynomial(params)
    plt.plot(pdata, ppoly(pdata), ms=2, color="blue", label="Best-fit rotation curve", zorder=4)

    # Keplerian prediction: V(R) = sqrt(G M / R_pc)  where R_pc = 1000 * R_kpc
    xkm   = np.linspace(0.1, 16, 50)
    ykm   = np.sqrt((G * gal_mass) / (1000 * xkm))
    params = fix_points(fit, xkm, ykm, np.array([0]), np.array([0]))
    kpoly  = np.polynomial.Polynomial(params)
    kdata  = np.linspace(0, 16, 50)
    plt.plot(kdata, kpoly(kdata), "--", ms=2, color="blue", label="Keplerian prediction", zorder=2)

    # Dark matter contribution: excess of observed curve over Keplerian
    ynk    = ppoly(pdata) - kpoly(pdata)
    params = fix_points(fit, pdata, ynk, np.array([0]), np.array([0]))
    npoly  = np.polynomial.Polynomial(params)
    plt.plot(pdata, npoly(pdata), "--", ms=2, color="red", label="Dark matter contribution", zorder=3)

    # Sun: R_sun = 8.3 kpc, V_sun = 236 km/s (consistent with coordinate transform parameters)
    plt.plot(8.3, 236, marker="*", ms=8, color="yellow", label="Sun", zorder=5)
    plt.title("Milky Way Rotation Curve")
    plt.xlabel("Galactic radius ($kpc$)")
    plt.ylabel("Rotational velocity ($km\\ s^{-1}$)")
    plt.legend(loc="upper right")
    plt.xlim(0, 16)
    plt.ylim(0, 400)
    plt.show()


def rotational_curvature(stellar_type, sample, fit):
    """Wrapper: print info and produce the rotation-curve plot."""
    print(f"\nRotation curve — type: {stellar_type}, sample: {sample}")
    rotational_curvature_plot(stellar_type, sample, fit)


# ---------------------------------------------------------------------------
# Dark matter mass distribution
# ---------------------------------------------------------------------------

def dark_matter_distribution_plot(stellar_type, sample, fit):
    """Plot the inferred enclosed mass distribution and Keplerian prediction.

    Prompts for a reference rotational velocity (km/s) to set the Keplerian baseline.
    """
    try:
        gal_vel = float(input("Enter a reference galactic rotational velocity (km/s): "))
    except ValueError:
        print("Invalid input — please enter a numeric value.")
        return dark_matter_distribution_plot(stellar_type, sample, fit)

    df   = gen_types(stellar_type, sample)
    data = np.linspace(0, 16, 50)

    # Observational data
    plt.errorbar(
        df["kiloparsec"], df["mass"],
        fmt="o", ms=0.5, color="black",
        xerr=df["kiloparsec_error"], yerr=df["mass_error"],
        ecolor="grey", elinewidth=0.25, label="Gaia DR2", zorder=1,
    )

    # Best-fit mass distribution (constrained through origin)
    params = fix_points(fit, df["kiloparsec"].values, df["mass"].values,
                        np.array([0]), np.array([0]))
    poly = np.polynomial.Polynomial(params)
    plt.plot(data, poly(data), color="blue", label="Mass distribution fit", zorder=4)

    # Keplerian prediction: M(R) = R_pc V^2 / G
    xkm   = np.linspace(0.1, 16, 50)
    ykm   = (1000 * xkm * gal_vel ** 2) / G
    params = fix_points(fit, xkm, ykm, np.array([0]), np.array([0]))
    kpoly  = np.polynomial.Polynomial(params)
    kdata  = np.linspace(0, 16, 50)
    plt.plot(kdata, kpoly(kdata), "--", ms=2, color="blue", label="Keplerian prediction", zorder=2)

    # Dark matter contribution
    ynk    = poly(data) - kpoly(data)
    params = fix_points(fit, data, ynk, np.array([0]), np.array([0]))
    npoly  = np.polynomial.Polynomial(params)
    plt.plot(data, npoly(data), "--", ms=2, color="red", label="Dark matter contribution", zorder=3)

    # Reference: Sun at R = 8.3 kpc = 8300 pc, V = 236 km/s (matches coordinate transform)
    sun_mass = (8300 * 236 ** 2) / G
    plt.plot(8, sun_mass, marker="*", ms=8, color="yellow", label="Sun", zorder=5)

    plt.title("Dark Matter Mass Distribution")
    plt.xlabel("Galactic radius ($kpc$)")
    plt.ylabel("Enclosed mass ($M_\\odot$)")
    plt.legend(loc="upper right")
    plt.xlim(0, 16)
    plt.ylim(0, 1e12)
    plt.show()


def dark_matter_distribution(stellar_type, sample, fit):
    """Wrapper: print info and produce the dark matter distribution plot."""
    print(f"\nMass distribution — type: {stellar_type}, sample: {sample}")
    dark_matter_distribution_plot(stellar_type, sample, fit)


