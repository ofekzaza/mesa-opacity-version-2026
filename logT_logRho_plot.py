"""Standalone log T vs log rho diagram for the last profile of each model."""

import os

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Polygon, Rectangle
import mesa_reader as mr
import numpy as np

import utils
from generic_last_log_plot import _convective_mask, _plot_convective_segments


# Burning-stage reference data (kept here, not in the generic plot module).
BURNING_CONDITIONS = {
    "H": {"log_T": 7.1, "log_Rho": (1.0, 2.0), "color": "#1f77b4"},
    "He": {"log_T": 8.0, "log_Rho": (2.0, 4.0), "color": "#2ca02c"},
    "C": {"log_T": 8.8, "log_Rho": (4.5, 5.5), "color": "#ff7f0e"},
    "O": {"log_T": 9.3, "log_Rho": (6.5, 7.5), "color": "#d62728"},
}


def draw_burning_reference_lines(ax, stages=("H", "He", "C", "O")):
    for stage in stages:
        spec = BURNING_CONDITIONS[stage]
        log_t = spec["log_T"]
        rho_lo, rho_hi = spec["log_Rho"]
        color = spec["color"]
        ax.axvline(
            log_t, color=color, linestyle="--", linewidth=1.1,
            alpha=0.85, zorder=1,
            label=f"{stage} burning (log T ≈ {log_t:g}, log ρ ≈ {rho_lo:g}-{rho_hi:g})",
        )


def draw_instability_zones(ax):
    zone_color = "#e74c3c"
    logT_grid = np.linspace(3, 10, 500)

    # EOS boundaries (unified style, single legend entry)
    boundary_color = "#7f8c8d"
    boundary_style = "--"
    ax.plot(logT_grid, 3.0 * logT_grid - 23.32, color=boundary_color,
            linestyle=boundary_style, linewidth=1.2)
    ax.plot(logT_grid, 1.5 * logT_grid - 8.25, color=boundary_color,
            linestyle=boundary_style, linewidth=1.2)
    ax.axhline(y=6.0, color=boundary_color, linestyle=boundary_style,
               linewidth=1.2)
    ax.plot([], [], color=boundary_color, linestyle=boundary_style,
            linewidth=1.2, label="EOS boundaries")

    # H / He ionization strips
    ax.axvspan(3.9, 4.1, ymin=0.0, ymax=0.32, color=zone_color, alpha=0.22)
    ax.axvspan(4.4, 4.6, ymin=0.0, ymax=0.32, color=zone_color, alpha=0.22)

    # Pair instability (ellipse)
    pair_ellipse = Ellipse(xy=(9.15, 2.5), width=0.44, height=4.0,
                          edgecolor=zone_color, facecolor=zone_color,
                          alpha=0.30, linewidth=1.5)
    ax.add_patch(pair_ellipse)

    # Iron photodisintegration (polygon)
    rho_poly = np.linspace(5.5, 10.0, 100)
    T_poly = 0.11 * rho_poly + 9.0
    poly_vertices = np.vstack([
        [10.0, 5.5],
        np.column_stack([T_poly, rho_poly]),
        [10.0, 10.0],
    ])
    iron_polygon = Polygon(poly_vertices, edgecolor=zone_color,
                          facecolor=zone_color, alpha=0.30, linewidth=1.5)
    ax.add_patch(iron_polygon)

    # Electron capture (rectangle)
    ec_rect = Rectangle((8.9, 9.5), 0.4, 1.5, facecolor=zone_color,
                        edgecolor=zone_color, alpha=0.30, linewidth=1.0)
    ax.add_patch(ec_rect)

    # Single legend entry for all instability zones
    ax.plot([], [], color=zone_color, alpha=0.30, linewidth=8,
            label=r"Instability zones ($\Gamma_1<4/3$)")


def load_last_profiles(mass: int, names: list[str]) -> "dict[str, mr.MesaData | None]":
    loaded: dict[str, mr.MesaData | None] = {}
    for name in names:
        path = f"{mass}m-{name}/LOGS"
        try:
            num = utils.get_last_profile_num(path)
            loaded[name] = mr.MesaData(f"{path}/profile{num}.data")
        except Exception as e:
            print(f"    logT_logRho_plot: could not load last profile for {name}: {e}")
            loaded[name] = None
    return loaded


def plot(mass: int, names: list[str], phase: str = "",
         _preloaded: "dict | None" = None):
    plt.figure()
    phase_suffix = f"_{phase}" if phase else ""
    base_colors: dict[str, str] = {}
    conv_plotted = False

    all_winds = all("_reduction" in n for n in names)
    all_post = all(n.endswith("_post") for n in names)

    profiles = _preloaded if _preloaded is not None else load_last_profiles(mass, names)

    data_x_min = data_x_max = None
    data_y_min = data_y_max = None

    for name in names:
        prof = profiles.get(name) if profiles is not None else None
        if prof is None:
            continue
        if not hasattr(prof, "logT") or not hasattr(prof, "logRho"):
            print(f"logT/logRho not in prof for {name}")
            continue

        x_data = getattr(prof, "logT")
        y_data = getattr(prof, "logRho")

        xa = np.asarray(x_data, dtype=float)
        ya = np.asarray(y_data, dtype=float)
        data_x_min = xa.min() if data_x_min is None else min(data_x_min, float(xa.min()))
        data_x_max = xa.max() if data_x_max is None else max(data_x_max, float(xa.max()))
        data_y_min = ya.min() if data_y_min is None else min(data_y_min, float(ya.min()))
        data_y_max = ya.max() if data_y_max is None else max(data_y_max, float(ya.max()))

        kwargs, is_winds, base_name = utils.get_line_kwargs(
            name, base_colors, drop_winds=all_winds, drop_post=all_post
        )
        (line,) = plt.plot(x_data, y_data, **kwargs)
        color = line.get_color()
        if not is_winds:
            base_colors[base_name] = color
        plt.plot(x_data[-1], y_data[-1], "o", markersize=6, color=color)
        plt.plot(x_data[0], y_data[0], "o", markersize=6, color=color)

        if _plot_convective_segments(x_data, y_data, _convective_mask(prof), color):
            conv_plotted = True

    if conv_plotted:
        plt.plot([], [], color="0.25", linewidth=3.2,
                 label="Convection (thick line, conv_vel > 0)")

    ax = plt.gca()
    draw_instability_zones(ax)
    draw_burning_reference_lines(ax, ("H", "He", "C", "O"))

    if data_x_min is not None:
        pad_x = 0.05 * (data_x_max - data_x_min) or 0.5
        ax.set_xlim(data_x_min - pad_x, data_x_max + pad_x)
    if data_y_min is not None:
        pad_y = 0.05 * (data_y_max - data_y_min) or 0.5
        ax.set_ylim(data_y_min - pad_y, data_y_max + pad_y)

    plt.xlabel(r"$\log(T)$ [K]")
    plt.ylabel(r"$\log(\rho)$ [g/cm$^3$]")

    plt.legend()
    plt.tight_layout()

    os.makedirs(f"{utils.PLOTS_ROOT}/plot{mass}m{phase_suffix}", exist_ok=True)
    plt.savefig(
        f"{utils.PLOTS_ROOT}/plot{mass}m{phase_suffix}/Mass_{mass}_logT_vs_logRho.png",
        dpi=300,
    )
    plt.close()