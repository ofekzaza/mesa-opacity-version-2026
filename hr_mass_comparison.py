"""HR diagrams comparing the "ours" models across mass.

Two figures are produced: one for the pre-MS (``*_reduction``) runs and one for
the post-MS (``*_reduction_post``) runs.  Every mass is drawn as its own
colour-coded evolutionary track (``log_Teff`` versus ``log_L``) with a matching
legend entry, so the mass of each line is unambiguous.

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
    ax.legend(
        title=r"Mass", loc="best", fontsize=9, title_fontsize=9,
        framealpha=1.0, ncol=2,
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
