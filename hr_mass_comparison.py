"""HR diagrams comparing the "ours" models across mass.

Two figures are produced: one for the pre-MS (``*_reduction``) runs and one for
the post-MS (``*_reduction_post``) runs.  Every mass is drawn as its own
colour-coded evolutionary track (``log_Teff`` versus ``log_L``) with a matching
legend entry, so the mass of each line is unambiguous.

The empirical Humphreys-Davidson (HD) upper-luminosity limit is drawn on both
figures for the metallicity of these models (solar, ``Z = 0.02``).  The
classical Galactic limit is temperature-independent at
``log(L/L_sun) = 5.8`` for ``log(Teff/K) <= 4.0`` and slopes up to
``log(L/L_sun) ~ 6.8`` at ``log(Teff/K) = 4.6`` for the hot side
(Humphreys & Davidson 1979, ApJ 232, 409; see also Woosley 2014 and
Sabhahit et al. 2021, MNRAS 506, 1).

Usage
-----
    uv run hr_mass_comparison.py
    uv run hr_mass_comparison.py --masses 30 50 100
    uv run hr_mass_comparison.py --phase pre_ms
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np

import utils

# The masses and model variants to compare.
DEFAULT_MASSES = [30, 40, 50, 60, 70, 80, 90, 100]

# phase label -> run-name suffix
PHASES = {
    "pre_ms": "ours_reduction",
    "post_ms": "ours_reduction_post",
}

# One clearly different colour per mass.  The same mapping is used for both
# the pre- and post-MS figures, so a colour always means the same mass.
MASS_COLOURS = {
    30: "#1F77B4",   # blue
    40: "#FF7F0E",   # orange
    50: "#2CA02C",   # green
    60: "#D62728",   # red
    70: "#9467BD",   # purple
    80: "#8C564B",   # brown
    90: "#E377C2",   # pink
    100: "#17BECF",  # cyan
}


# ---------------------------------------------------------------------------
# Humphreys-Davidson limit for solar metallicity (Z = Z_sun = 0.02).
#
# The empirical boundary has a temperature-independent (horizontal) part at
# log(L/L_sun) = 5.8 for cool stars, and a sloping hot part that rises towards
# the O-star luminosity limit.  We join the two at log(Teff/K) = 4.0 (the cool
# edge quoted for the temperature-independent section, e.g. Sabhahit et al.
# 2021) and reach log(L/L_sun) = 6.8 at log(Teff/K) = 4.6 (~40,000 K, the
# maximum luminosity of O stars; Woosley 2014).
#
# These constants are kept here so the limit can easily be adjusted for a
# different metallicity (e.g. log(L/L_sun) = 5.5 for the downward-revised
# Magellanic/Galactic value of Davies et al. 2018 / Davies & Beasor 2020).
# ---------------------------------------------------------------------------
HD_LIMIT = {
    "log_l": 5.8,          # horizontal (temperature-independent) part
    "log_teff_knee": 4.0,  # where the horizontal part ends
    "log_teff_hot": 4.6,   # hottest point of the sloping part
    "log_l_hot": 6.8,      # upper luminosity reached at log_teff_hot
}
HD_LIMIT_COLOUR = "#222222"


def hd_limit_curve(log_teff, hd=HD_LIMIT):
    """Return ``log(L/L_sun)`` of the two-part HD limit at ``log_teff``.

    ``log_teff`` may be a scalar or an array.  Cool stars sit on the
    temperature-independent part; hotter stars follow the sloping segment and
    the limit is capped at ``log_l_hot`` beyond ``log_teff_hot``.
    """
    log_teff = np.asarray(log_teff, dtype=float)
    log_l = hd["log_l"]
    knee = hd["log_teff_knee"]
    hot = hd["log_teff_hot"]
    log_l_hot = hd["log_l_hot"]

    slope = (log_l_hot - log_l) / (hot - knee)
    sloping = log_l + slope * (log_teff - knee)
    return np.where(
        log_teff <= knee,
        log_l,
        np.minimum(sloping, log_l_hot),
    )


def _mass_colour(mass, masses):
    if mass in MASS_COLOURS:
        return MASS_COLOURS[mass]
    # Fallback for masses outside the standard table.
    sample = np.linspace(0.05, 0.95, len(masses))
    idx = list(masses).index(mass)
    return plt.cm.viridis(float(sample[idx]))


def load_track(mass: int, suffix: str):
    """Return (log_Teff, log_L) for one run, or None if it cannot be read."""
    path = f"{mass}m-{suffix}/LOGS/history.data"
    try:
        hist = mr.MesaData(path)
    except Exception as exc:
        print(f"    HR comparison: could not read {path}: {exc}")
        return None

    teff = np.asarray(hist.log_Teff, dtype=float)
    lum = np.asarray(hist.log_L, dtype=float)
    model = np.asarray(hist.model_number, dtype=float)
    good = np.isfinite(teff) & np.isfinite(lum) & (model > 0)
    teff, lum = teff[good], lum[good]
    if teff.size < 2:
        print(f"    HR comparison: no usable history rows in {path}")
        return None
    return teff, lum


def plot_phase(phase: str, masses, out_dir: str = "") -> str | None:
    """Render one HR comparison figure for a phase; return its path."""
    if phase not in PHASES:
        raise ValueError(f"unknown phase {phase!r}; choose from {sorted(PHASES)}")
    suffix = PHASES[phase]
    phase_label = phase.replace("_", "-").replace("ms", "MS")

    fig, ax = plt.subplots(figsize=(10.0, 7.5))
    plotted = 0
    for mass in masses:
        track = load_track(mass, suffix)
        if track is None:
            continue
        teff, lum = track
        colour = _mass_colour(mass, masses)
        ax.plot(
            teff, lum, color=colour, linewidth=1.8, alpha=0.95,
            label=rf"{mass} $M_\odot$",
        )
        # Mark the first and last model of each track.
        ax.plot(teff[0], lum[0], "o", color=colour, markersize=4.5,
                markeredgecolor="white", markeredgewidth=0.6, zorder=4)
        ax.plot(teff[-1], lum[-1], "s", color=colour, markersize=4.5,
                markeredgecolor="white", markeredgewidth=0.6, zorder=4)
        plotted += 1

    if plotted == 0:
        plt.close(fig)
        print(f"    HR comparison: nothing to plot for phase={phase}")
        return None

    # Draw the empirical Humphreys-Davidson limit for solar metallicity.
    # Span the full plotted temperature range (the x-axis is inverted later,
    # which does not affect the data coordinates used here).
    x_lo, x_hi = ax.get_xlim()
    teff_grid = np.linspace(min(x_lo, x_hi), max(x_lo, x_hi), 512)
    hd_line, = ax.plot(
        teff_grid,
        hd_limit_curve(teff_grid),
        color=HD_LIMIT_COLOUR,
        linewidth=2.2,
        linestyle="--",
        alpha=0.9,
        zorder=6,
        label=r"HD limit",
    )

    ax.set_xlabel(r"$\log\,(T_{\mathrm{eff}}/\mathrm{K})$")
    ax.set_ylabel(r"$\log\,(L/L_\odot)$")
    ax.set_title(
        f"HR diagram of our models, {min(masses)}-{max(masses)} "
        rf"$M_\odot$ ({phase_label})"
    )
    # HR diagrams conventionally have the effective temperature decreasing
    # to the right.
    ax.invert_xaxis()
    ax.grid(True, which="major", alpha=0.25, linewidth=0.6)
    mass_legend = ax.legend(
        title=r"Mass", loc="best", fontsize=9, title_fontsize=9,
        framealpha=1.0, ncol=2,
    )
    ax.add_artist(mass_legend)
    ax.legend(
        handles=[hd_line], loc="lower left", bbox_to_anchor=(0.0, 0.045),
        fontsize=8.5, framealpha=1.0,
    )
    ax.text(
        0.01, 0.01, "circles: first model    squares: last model",
        transform=ax.transAxes, fontsize=7.5, color="#555555",
        ha="left", va="bottom",
    )
    fig.tight_layout()

    out_dir = out_dir or f"{utils.PLOTS_ROOT}/hr_mass_comparison"
    os.makedirs(out_dir, exist_ok=True)
    out_path = (
        f"{out_dir}/HR_ours_masses_{min(masses)}-{max(masses)}m_{phase}.png"
    )
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"    Saved: {out_path}")
    return out_path


def plot_all(masses=None, phases=None, out_dir: str = "") -> list[str]:
    masses = sorted(masses or DEFAULT_MASSES)
    phases = phases or list(PHASES)
    paths = []
    for phase in phases:
        path = plot_phase(phase, masses, out_dir=out_dir)
        if path is not None:
            paths.append(path)
    return paths


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="HR diagrams comparing the our-model runs across mass."
    )
    parser.add_argument(
        "--masses", type=int, nargs="+", default=DEFAULT_MASSES,
        help=f"masses to plot (default: {DEFAULT_MASSES})",
    )
    parser.add_argument(
        "--phase", choices=sorted(PHASES), nargs="+", default=None,
        help="only render these phases (default: both)",
    )
    parser.add_argument(
        "--out-dir", default="",
        help=f"output directory (default: {utils.PLOTS_ROOT}/hr_mass_comparison)",
    )
    args = parser.parse_args(argv)

    paths = plot_all(args.masses, args.phase, out_dir=args.out_dir)
    return 0 if paths else 1


if __name__ == "__main__":
    raise SystemExit(main())
