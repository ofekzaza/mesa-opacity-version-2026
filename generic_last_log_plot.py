import glob
import os

import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np

import utils

def _convective_mask(prof):
    """True for zones where convection is active (conv_vel > 0)."""
    if not hasattr(prof, "conv_vel"):
        return None
    return np.array(prof.conv_vel, dtype=float) > 0.0


def _plot_convective_segments(x, y, mask, color):
    """Overlay contiguous convective blocks along an existing profile curve."""
    if mask is None:
        return False
    diff = np.diff(np.concatenate(([False], mask, [False])).astype(int))
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    for s, e in zip(starts, ends):
        s0 = max(0, s - 1)
        e0 = min(len(x), e + 1)
        plt.plot(x[s0:e0], y[s0:e0], color=color, linewidth=3.2, alpha=0.85, zorder=5)
    return len(starts) > 0

y_axises = [
    "L_div_Ledd_effective",
    "log_L_div_Ledd",
    "log_opacity",
    "extra_opacity_factor",
    "logT",
    "entropy",
]


def load_last_profiles(mass: int, names: list[str]) -> "dict[str, mr.MesaData | None]":
    """Load each name's last profile exactly once.

    Returns a dict mapping name -> MesaData (or None if loading failed).
    Pass the result as `_preloaded` to plot() / plot_generic_last_log_plot()
    so that repeated calls for the same names share the already-loaded objects.
    """
    loaded: dict[str, mr.MesaData | None] = {}
    for name in names:
        path = f"{mass}m-{name}/LOGS"
        try:
            num = utils.get_last_profile_num(path)
            loaded[name] = mr.MesaData(f"{path}/profile{num}.data")
        except Exception as e:
            print(f"    generic_last_log_plot: could not load last profile for {name}: {e}")
            loaded[name] = None
    return loaded


def plot(
    mass: int,
    names: list[str],
    x_axis: str,
    x_units: str,
    y_axis: str,
    y_units: str,
    phase: str = "",
    inverse_x: bool = False,
    _preloaded: "dict | None" = None,
):
    """Plot y_axis vs x_axis from the last profile of each name.

    If _preloaded is given (a dict name -> MesaData), profiles are taken from
    there instead of being read from disk.  Pass the output of
    load_last_profiles() to share one disk read across many plot() calls.
    """
    plt.figure()
    phase_suffix = f"_{phase}" if phase else ""
    base_colors: dict[str, str] = {}
    conv_plotted = False

    # Load profiles if not pre-supplied
    profiles = _preloaded if _preloaded is not None else load_last_profiles(mass, names)

    for name in names:
        prof = profiles.get(name)
        if prof is None:
            continue

        if not hasattr(prof, y_axis):
            print(f"{y_axis} is not in prof for {name}")
            continue
        if not hasattr(prof, x_axis):
            print(f"{x_axis} is not in prof for {name}")
            continue

        x_data = getattr(prof, x_axis)
        y_data = getattr(prof, y_axis)

        kwargs, is_winds, base_name = utils.get_line_kwargs(name, base_colors)
        (line,) = plt.plot(x_data, y_data, **kwargs)
        color = line.get_color()
        if not is_winds:
            base_colors[base_name] = color
        plt.plot(x_data[-1], y_data[-1], "o", markersize=6, color=color)
        plt.plot(x_data[0], y_data[0], "o", markersize=6, color=color)

        if _plot_convective_segments(x_data, y_data, _convective_mask(prof), color):
            conv_plotted = True

    if conv_plotted:
        plt.plot([], [], color="0.25", linewidth=3.2, label="Convection (thick line, conv_vel > 0)")

    x_units_str = fr" [{x_units}]" if x_units else ""
    y_units_str = fr" [{y_units}]" if y_units else ""
    plt.xlabel(f"{utils.pretty_axis_name(x_axis)}{x_units_str}")
    plt.ylabel(f"{utils.pretty_axis_name(y_axis)}{y_units_str}")

    if inverse_x:
        plt.gca().invert_xaxis()

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    os.makedirs(f"{utils.PLOTS_ROOT}/plot{mass}m{phase_suffix}", exist_ok=True)
    plt.savefig(f"{utils.PLOTS_ROOT}/plot{mass}m{phase_suffix}/Mass_{mass}_{y_axis}_vs_{x_axis}.png", dpi=300)
    plt.close()


def plot_generic_last_log_plot(
    mass: int,
    names: list[str],
    x_axis: str = "logRho",
    x_units: str = r"g/cm$^3$",
    phase: str = "",
    y_units: str = "",
    inverse_x: bool = False,
    _preloaded: "dict | None" = None,
):
    """Plot all y_axises vs x_axis, sharing one last-profile read per name.

    Loads last profiles once (if _preloaded is not given), then reuses them
    for every y_axis — so N names × M y-axes = N profile reads, not N×M.
    """
    if _preloaded is None:
        _preloaded = load_last_profiles(mass, names)
    for y_axis in y_axises:
        print(y_axis, names)
        plot(
            mass,
            names,
            x_axis=x_axis,
            x_units=x_units,
            y_axis=y_axis,
            y_units=y_units,
            phase=phase,
            inverse_x=inverse_x,
            _preloaded=_preloaded,
        )
