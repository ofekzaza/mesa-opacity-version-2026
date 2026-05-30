import mesa_reader as mr
import matplotlib.pyplot as plt

import os

# List of history file names (without extensions)
names = [
    # "normal",
    # "ours_interval_2",
    # "ours_interval_3",
    # "ours_thermal_time",
    # "ours_interval_5",
    # "ours_interval_10",
    # "mlt++",
    # "supereduction_a=2",
    # "ours_interval_40",
    # "supereduction_a=5",
    "ours",
    "ours_reduction",
]
#'supereduction_a=5' empty it doesnt run
# Initialize the plot
plt.figure()

mass = 80

import utils

# Iterate over each history file
def plot_history_diagrams(names: list[str], mass: int, x_axis: str, y_axis: str, y_scale: str):
    plt.figure()
    for name in names:
        path = f"{mass}m-{name}/LOGS/history.data"
        plot_name = utils.parse_name_for_plots(name)
        try:
            # Load the history data
            history = mr.MesaData(path)
            # Plot log_L vs. log_Teff
            x = getattr(history, x_axis)
            if y_axis == "1":
                y = getattr(history, "time_step") / getattr(history, "kh_timescale")
            else:
                y = getattr(history, y_axis)
            (line,) = plt.plot(x, y, label=plot_name, linewidth=0.8)
            color = line.get_color()
            plt.plot(x[-1], y[-1], "o", color=color)
        except FileNotFoundError:
            print(f"Warning: {path}.data not found. Skipping.")

    # Customize the plot
    plt.xlabel(utils.parse_name_for_plots(x_axis))
    plt.ylabel(utils.parse_name_for_plots(y_axis))
    if y_scale == "log":
        plt.yscale("log")
    # plt.title(f"{mass} Mass Star {y_axis} by {x_axis}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    # Save the plot
    os.makedirs(f"plot{mass}m", exist_ok=True)
    plt.savefig(f"plot{mass}m/Mass_{mass}_{y_axis}_vs_{x_axis}.png", dpi=300)
    plt.close()


todo = [
    ("model_number", "1", "notlog"),
    ("model_number", "time_step_sec", "log"),
    ("model_number", "kh_timescale", "log"),
]
for x_axis, y_axis, scale in todo:
    print(x_axis, y_axis)
    plot_history_diagrams(names, mass, x_axis, y_axis, scale)
