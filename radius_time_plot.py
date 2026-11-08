"""Radius-versus-time diagrams comparing the "ours" models across mass.

For each phase (pre-MS and post-MS) two versions are produced: one with a
logarithmic time axis and one with a linear time axis.  Both use
``log(R / R_sun)`` on the radius axis.  Each mass is drawn as its own
colour-coded track against ``star_age``, using the same colour per mass as the
HR comparison plots.

Usage
-----
    uv run radius_time_plot.py
    uv run radius_time_plot.py --masses 30 50 100
    uv run radius_time_plot.py --phase pre_ms --scale log
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np

import utils
from hr_mass_comparison import DEFAULT_MASSES, MASS_COLOURS, PHASES


def _mass_colour(mass, masses):
    if mass in MASS_COLOURS:
        return MASS_COLOURS[mass]
    sample = np.linspace(0.05, 0.95, len(masses))
    idx = list(masses).index(mass)
    return plt.cm.viridis(float(sample[idx]))


def load_track(mass: int, suffix: str):
    """Return (star_age, log_R) for one run, or None if it cannot be read."""
    path = f"{mass}m-{suffix}/LOGS/history.data"
    try:
        hist = mr.MesaData(path)
    except Exception as exc:
        print(f"    Radius-time: could not read {path}: {exc}")
        return None

    age = np.asarray(hist.star_age, dtype=float)
    log_r = np.asarray(hist.log_R, dtype=float)
    model = np.asarray(hist.model_number, dtype=float)
    good = (
        np.isfinite(age) & np.isfinite(log_r) & (model > 0) & (age > 0.0)
    )
    age, log_r = age[good], log_r[good]
    if age.size < 2:
        print(f"    Radius-time: no usable history rows in {path}")
        return None
    return age, log_r


def plot_phase(
    phase: str, masses, out_dir: str = "", log_scale: bool = True
) -> str | None:
    """Render one radius-time comparison figure for a phase."""
    if phase not in PHASES:
        raise ValueError(f"unknown phase {phase!r}; choose from {sorted(PHASES)}")
    suffix = PHASES[phase]
    phase_label = phase.replace("_", "-").replace("ms", "MS")
    scale_label = "log time" if log_scale else "linear time"
    scale_suffix = "log" if log_scale else "linear"

    fig, ax = plt.subplots(figsize=(10.0, 7.5))
    plotted = 0
    for mass in masses:
        track = load_track(mass, suffix)
        if track is None:
            continue
        age, y = track
        # Both versions use log(R/R_sun) on the radius axis; they differ only
        # in whether the time axis is logarithmic.
        colour = _mass_colour(mass, masses)
        ax.plot(
            age, y, color=colour, linewidth=1.8, alpha=0.95,
            label=rf"{mass} $M_\odot$",
        )
        ax.plot(age[0], y[0], "o", color=colour, markersize=4.5,
                markeredgecolor="white", markeredgewidth=0.6, zorder=4)
        ax.plot(age[-1], y[-1], "s", color=colour, markersize=4.5,
                markeredgecolor="white", markeredgewidth=0.6, zorder=4)
        plotted += 1

    if plotted == 0:
        plt.close(fig)
        print(f"    Radius-time: nothing to plot for phase={phase}")
        return None

    if log_scale:
        ax.set_xscale("log")
    ax.set_xlabel("Star age [yr]")
    ax.set_ylabel(r"$\log\,(R/R_\odot)$")
    ax.set_title(
        f"Radius evolution of our models, {min(masses)}-{max(masses)} "
        rf"$M_\odot$ ({phase_label}, {scale_label})"
    )
    ax.grid(True, which="both", alpha=0.2, linewidth=0.6)
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

    out_dir = out_dir or f"{utils.PLOTS_ROOT}/radius_time"
    os.makedirs(out_dir, exist_ok=True)
    out_path = (
        f"{out_dir}/RADIUS_ours_masses_{min(masses)}-{max(masses)}m_{phase}"
        f"_{scale_suffix}.png"
    )
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"    Saved: {out_path}")
    return out_path


def plot_all(masses=None, phases=None, out_dir: str = "",
             scales=("log", "linear")) -> list[str]:
    masses = sorted(masses or DEFAULT_MASSES)
    phases = phases or list(PHASES)
    paths = []
    for phase in phases:
        for scale in scales:
            path = plot_phase(
                phase, masses, out_dir=out_dir, log_scale=(scale == "log")
            )
            if path is not None:
                paths.append(path)
    return paths


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Radius-versus-time diagrams for the our-model runs."
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
        "--scale", choices=("log", "linear"), nargs="+",
        default=("log", "linear"),
        help="axis scale(s) to render (default: both)",
    )
    parser.add_argument(
        "--out-dir", default="",
        help=f"output directory (default: {utils.PLOTS_ROOT}/radius_time)",
    )
    args = parser.parse_args(argv)

    paths = plot_all(
        args.masses, args.phase, out_dir=args.out_dir, scales=args.scale
    )
    return 0 if paths else 1


if __name__ == "__main__":
    raise SystemExit(main())
