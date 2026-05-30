import json
import math
import os

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np

import utils


def plot(
    mass: int,
    name: str,
    model_var: str,
    y_axis: str,
    y_units: str,
    value_axis: str,
    normalize: bool = True,
):
    """
    model_var: the field holding model number (always 'model_number')
    y_axis: the vertical axis variable (e.g. 'logRho')
    value_axis: the heatmap color variable (e.g. 'extra_opacity_factor')
    """

    path = f"{mass}m-{name}/LOGS"

    print(f"Reading profiles from: {path}")

    profiles = []

    # Read all profiles and collect model numbers
    model_numbers = []
    star_ages = []
    Y_all = []
    VAL_all = []

    for i in range(1, 9999):
        f_path = f"{path}/profile{i}.data"
        if not os.path.exists(f_path):
            break

        prof = mr.MesaData(f_path)

        # Ensure required fields exist
        if not hasattr(prof, y_axis):
            print(f"Profile {i} missing required fields, skipping. {y_axis}")
            continue
        if not hasattr(prof, value_axis):
            print(f"Profile {i} missing required fields, skipping. {value_axis}")
            continue

        profiles.append(f_path)
        model_numbers.append(getattr(prof, "model_number"))
        star_ages.append(getattr(prof, "star_age"))

        Y_all.append(np.array(getattr(prof, y_axis), dtype=float))
        VAL_all.append(np.array(getattr(prof, value_axis), dtype=float))

    if not profiles:
        raise FileNotFoundError(f"No profiles found. for {mass}m-{name}")

    # Build common y-grid (vertical axis)
    min_y = min(np.min(arr) for arr in Y_all)
    max_y = max(np.max(arr) for arr in Y_all)
    n_points = max(len(arr) for arr in Y_all)

    common_y = np.linspace(min_y, max_y, n_points)

    # Heatmap matrix
    n_models = len(model_numbers)
    DATA = np.zeros((n_points, n_models))

    # Fill the heatmap interpolating each profile
    for j in range(n_models):
        y_data = Y_all[j]
        val_data = VAL_all[j]

        # Interpolate to common y grid
        # Sort y_data to ensure it is increasing for np.interp
        sort_idx = np.argsort(y_data)
        interp_val = np.interp(
            common_y, y_data[sort_idx], val_data[sort_idx], left=np.nan, right=np.nan
        )
        DATA[:, j] = interp_val

    # Save raw matrix for debugging
    with open("tmp.json", "w") as f:
        # replace np.nan with None for valid JSON
        json.dump(np.where(np.isnan(DATA), None, DATA).tolist(), f)

    # Sort by model number (important!)
    model_numbers = np.array(model_numbers)
    star_ages = np.array(star_ages)
    sort_idx = np.argsort(model_numbers)

    DATA = DATA[:, sort_idx]
    model_numbers = model_numbers[sort_idx]
    star_ages = star_ages[sort_idx]

    # Plot heatmap
    plt.figure(figsize=(10, 6))
    ax = plt.gca()
    data_for_plot = DATA * 100

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color=ax.get_facecolor())

    # Use PowerNorm to expand the contrast near the maximum values (gamma > 1)
    # This prevents the bulk of data near 100 from looking identical, making tiny drops obvious!
    if normalize:
        min_val = math.floor(np.nanmin(data_for_plot))
        norm = mcolors.PowerNorm(gamma=4.0, vmin=min_val, vmax=100)
    else:
        norm = None

    # Use pcolormesh to map non-linear model_numbers accurately
    plt.pcolormesh(
        model_numbers, common_y, data_for_plot, cmap=cmap, norm=norm, shading="nearest"
    )

    # Add boundary lines for min and max logRho of the star
    min_rhos = np.array([np.min(Y_all[k]) for k in range(n_models)])[sort_idx]
    max_rhos = np.array([np.max(Y_all[k]) for k in range(n_models)])[sort_idx]
    plt.plot(
        model_numbers,
        min_rhos,
        color="red",
        linewidth=1.5,
        linestyle="--",
        label="Star Boundaries",
    )
    plt.plot(model_numbers, max_rhos, color="red", linewidth=1.5, linestyle="--")
    plt.legend(loc="upper right", framealpha=0.7)

    cbar = plt.colorbar(label=utils.parse_name_for_plots(value_axis) + " [%]")
    # Manually configure visually evenly-spaced ticks to prevent overlapping labels at the bottom due to PowerNorm
    if normalize:
        norm_vals = np.linspace(0, 1, 6)
        tick_vals = min_val + (100 - min_val) * (norm_vals ** (1 / 4.0))
        tick_vals = np.unique(np.round(tick_vals).astype(int))
        cbar.set_ticks(tick_vals)
        cbar.set_ticklabels([str(t) for t in tick_vals])

    # Update x-axis labels to include star_age
    ax = plt.gca()
    # Get current ticks (auto-selected by matplotlib based on extent)
    # We force a draw to make sure ticks are populated, though usually get_xticks works if we use the current axes
    xticks = ax.get_xticks()

    # Filter xticks to those within valid range
    xticks = [x for x in xticks if model_numbers[0] <= x <= model_numbers[-1]]

    new_labels = []
    for x in xticks:
        # Interpolate to find age for this model number
        # We can use np.interp because model_numbers are sorted
        age = np.interp(x, model_numbers, star_ages)
        new_labels.append(f"{int(x)}\n{age:.2e}")

    ax.set_xticks(xticks)
    ax.set_xticklabels(new_labels)

    plt.xlabel("Model Number\n(Star Age [years])")
    plt.ylabel(f"{utils.parse_name_for_plots(y_axis)} [{y_units}]")
    # plt.title(f"{value_axis} vs {y_axis} across Model Number")

    plt.tight_layout()
    plt.savefig(
        f"plot{mass}m/HEATMAP_{mass}_{name}_{value_axis}_vs_{y_axis}.png", dpi=300
    )
    plt.close()


def plot_heat_map(
    mass: int,
    name: str,
    value_axis: str = "extra_opacity_factor",
    normalize: bool = True,
):
    plot(
        mass=mass,
        name=name,
        model_var="model_number",
        y_axis="logRho",
        y_units="g/cm^3",
        value_axis=value_axis,
        normalize=normalize,
    )


if __name__ == "__main__":
    # for mass in [60]:
    #     plot(
    #         mass=mass,
    #         name = "ours_reduction",
    #         model_var="model_number",
    #         y_axis="logRho",
    #         y_units="g/cm^3",
    #         value_axis="extra_opacity_factor",
    #     )
    for mass in [80]:
        plot(
            mass=mass,
            name="ours-post",
            model_var="model_number",
            y_axis="logRho",
            y_units="g/cm^3",
            value_axis="extra_opacity_factor",  # "L_div_Ledd_effective",
            normalize=False,
        )
