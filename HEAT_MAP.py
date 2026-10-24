import json
import math
import os

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np

import utils


class HeatMapCollector:
    """Streams through profiles once, collecting logRho + one-or-more value columns.

    Only numpy arrays are held in memory — MesaData objects are discarded
    by the caller (stream.stream_profiles) after each consume() call.
    """

    def __init__(self, value_axes: list[str], y_axis: str = "logRho"):
        self.y_axis = y_axis
        self.value_axes = list(value_axes)
        self.model_numbers: list = []
        self.star_ages: list = []
        self.Y_all: list = []
        self.VAL_all: dict[str, list] = {va: [] for va in value_axes}

    def consume(self, prof) -> None:
        if not hasattr(prof, self.y_axis):
            return
        self.model_numbers.append(getattr(prof, "model_number"))
        self.star_ages.append(getattr(prof, "star_age"))
        self.Y_all.append(np.array(getattr(prof, self.y_axis), dtype=float))
        for va in self.value_axes:
            if hasattr(prof, va):
                self.VAL_all[va].append(np.array(getattr(prof, va), dtype=float))
            else:
                self.VAL_all[va].append(None)

    def render_all(
        self,
        mass: int,
        name: str,
        normalize_map: "dict[str, bool] | None" = None,
        y_units: str = "g/cm^3",
        phase: str = "",
    ) -> None:
        if not self.model_numbers:
            raise FileNotFoundError(f"No profiles collected for {mass}m-{name}")
        if normalize_map is None:
            normalize_map = {}
        phase_suffix = f"_{phase}" if phase else ""
        os.makedirs(f"{utils.PLOTS_ROOT}/plot{mass}m{phase_suffix}", exist_ok=True)
        model_numbers = np.array(self.model_numbers)
        star_ages = np.array(self.star_ages)
        for va in self.value_axes:
            _render_heatmap(
                model_numbers=model_numbers,
                star_ages=star_ages,
                Y_all=self.Y_all,
                VAL_list=self.VAL_all[va],
                mass=mass,
                name=name,
                value_axis=va,
                y_axis=self.y_axis,
                y_units=y_units,
                normalize=normalize_map.get(va, True),
                phase=phase,
            )


def _render_heatmap(
    model_numbers,
    star_ages,
    Y_all,
    VAL_list,
    mass,
    name,
    value_axis,
    y_axis,
    y_units,
    normalize,
    phase,
):
    """Pure rendering — no file I/O. All data already extracted from profiles."""
    phase_suffix = f"_{phase}" if phase else ""

    # Filter profiles where this value_axis was missing
    valid_idx = [j for j, v in enumerate(VAL_list) if v is not None]
    if not valid_idx:
        print(f"    HEATMAP: no data for {value_axis} in {mass}m-{name}, skipping")
        return

    mn = model_numbers[valid_idx]
    sa = star_ages[valid_idx]
    Y = [Y_all[j] for j in valid_idx]
    VAL = [VAL_list[j] for j in valid_idx]
    n_v = len(valid_idx)

    # Build common y-grid
    min_y = min(np.min(arr) for arr in Y)
    max_y = max(np.max(arr) for arr in Y)
    n_points = max(len(arr) for arr in Y)
    common_y = np.linspace(min_y, max_y, n_points)

    # Build heatmap matrix via interpolation
    DATA = np.full((n_points, n_v), np.nan)
    for j in range(n_v):
        y_data = Y[j]
        val_data = VAL[j]
        sort_idx = np.argsort(y_data)
        DATA[:, j] = np.interp(
            common_y, y_data[sort_idx], val_data[sort_idx], left=np.nan, right=np.nan
        )

    # Save raw matrix for debugging
    with open("tmp.json", "w") as f:
        json.dump(np.where(np.isnan(DATA), None, DATA).tolist(), f)

    for x_axis in ("model_number", "star_age"):
        _draw_heatmap(
            mn=mn,
            sa=sa,
            Y=Y,
            DATA=DATA,
            common_y=common_y,
            mass=mass,
            name=name,
            value_axis=value_axis,
            y_axis=y_axis,
            y_units=y_units,
            normalize=normalize,
            phase_suffix=phase_suffix,
            x_axis=x_axis,
        )


def _draw_heatmap(
    mn,
    sa,
    Y,
    DATA,
    common_y,
    mass,
    name,
    value_axis,
    y_axis,
    y_units,
    normalize,
    phase_suffix,
    x_axis,
):
    """Render one heatmap image with either model number or star age on the x axis."""
    n_v = DATA.shape[1]

    # Sort columns by the chosen x axis
    if x_axis == "star_age":
        order = np.argsort(sa)
        x = sa[order]
        secondary = mn[order]
        t_end = x[-1]
        x = t_end - x
        x = np.maximum(x, 0.0)
        x_label = "Star Time To Live [years]\n(Model Number)"
        file_suffix = "_xstar_age"
    else:
        order = np.argsort(mn)
        x = mn[order]
        secondary = sa[order]
        x_label = "Model Number\n(Star Age [years])"
        file_suffix = ""

    DATA = DATA[:, order]
    Y_sorted = [Y[k] for k in order]

    # Plot
    plt.figure(figsize=(10, 6))
    ax = plt.gca()

    if value_axis == "extra_opacity_factor":
        data_for_plot = DATA * 100
        cbar_label = utils.pretty_axis_name(value_axis) + " [%]"
    else:
        data_for_plot = DATA
        cbar_label = utils.pretty_axis_name(value_axis)

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color=ax.get_facecolor())

    if normalize:
        min_val = math.floor(np.nanmin(data_for_plot))
        norm = mcolors.PowerNorm(gamma=4.0, vmin=min_val, vmax=100)
    else:
        norm = None

    if x_axis == "star_age":
        secondary_plot = secondary
        x_plot_axis = x
    else:
        x_plot_axis = x
        secondary_plot = secondary

    if x_axis == "star_age":
        ax.invert_xaxis()

    plt.pcolormesh(x_plot_axis, common_y, data_for_plot, cmap=cmap, norm=norm, shading="nearest")

    min_rhos = np.array([np.min(Y_sorted[k]) for k in range(n_v)])
    max_rhos = np.array([np.max(Y_sorted[k]) for k in range(n_v)])
    plt.plot(x_plot_axis, min_rhos, color="red", linewidth=1.5, linestyle="--", label="Star Boundaries")
    plt.plot(x_plot_axis, max_rhos, color="red", linewidth=1.5, linestyle="--")
    plt.legend(loc="upper right", framealpha=0.7)

    cbar = plt.colorbar(label=cbar_label)
    if normalize:
        norm_vals = np.linspace(0, 1, 6)
        tick_vals = min_val + (100 - min_val) * (norm_vals ** (1 / 4.0))
        tick_vals = np.unique(np.round(tick_vals).astype(int))
        cbar.set_ticks(tick_vals)
        cbar.set_ticklabels([str(t) for t in tick_vals])

    xticks = ax.get_xticks()
    if x_axis == "star_age":
        xticks = [t for t in xticks if x_plot_axis.min() <= t <= x_plot_axis.max()]
        new_labels = []
        for t in xticks:
            mn_interp = np.interp(t, x, secondary_plot)
            new_labels.append(f"{t:.2e}\n{int(round(mn_interp))}")
    else:
        xticks = [t for t in xticks if x_plot_axis.min() <= t <= x_plot_axis.max()]
        new_labels = [f"{int(t)}\n{np.interp(t, x_plot_axis, secondary_plot):.2e}" for t in xticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(new_labels)

    plt.xlabel(x_label)

    plt.ylabel(f"{utils.pretty_axis_name(y_axis)} [{y_units}]")
    plt.tight_layout()
    out_path = (
        f"{utils.PLOTS_ROOT}/plot{mass}m{phase_suffix}/"
        f"HEATMAP_{mass}_{name}_{value_axis}_vs_{y_axis}{file_suffix}.png"
    )
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"    Saved: {out_path}")


def plot_batch(
    mass: int,
    name: str,
    value_axes: list[str],
    normalize_map: "dict[str, bool] | None" = None,
    y_axis: str = "logRho",
    y_units: str = "g/cm^3",
    phase: str = "",
):
    """Stream profiles once, render one heatmap per value_axis.

    This is the preferred entry point when calling from main.py.
    """
    from stream import stream_profiles

    logs_path = f"{mass}m-{name}/LOGS"
    print(f"    HEATMAP: streaming {logs_path} for {value_axes}...")
    collector = HeatMapCollector(value_axes=value_axes, y_axis=y_axis)
    count = stream_profiles(logs_path, [collector])
    if count == 0:
        raise FileNotFoundError(f"No profiles found for {mass}m-{name}")
    collector.render_all(mass, name, normalize_map=normalize_map, y_units=y_units, phase=phase)


def plot(
    mass: int,
    name: str,
    model_var: str,
    y_axis: str,
    y_units: str,
    value_axis: str,
    normalize: bool = True,
    phase: str = "",
):
    """Backward-compatible wrapper — delegates to plot_batch."""
    plot_batch(
        mass,
        name,
        value_axes=[value_axis],
        normalize_map={value_axis: normalize},
        y_axis=y_axis,
        y_units=y_units,
        phase=phase,
    )


def plot_heat_map(
    mass: int,
    name: str,
    value_axis: str = "extra_opacity_factor",
    normalize: bool = True,
    phase: str = "",
):
    """Standalone convenience wrapper (single value_axis)."""
    plot_batch(
        mass,
        name,
        value_axes=[value_axis],
        normalize_map={value_axis: normalize},
        phase=phase,
    )


if __name__ == "__main__":
    for mass in [80]:
        plot(
            mass=mass,
            name="ours-post",
            model_var="model_number",
            y_axis="logRho",
            y_units="g/cm^3",
            value_axis="extra_opacity_factor",
            normalize=False,
        )
