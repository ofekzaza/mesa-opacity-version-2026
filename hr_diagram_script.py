import os

import matplotlib.pyplot as plt
import mesa_reader as mr

# List of history file names (without extensions)
NAMES = [
    # "normal",
    # "ours_interval_2",
    # "ours_interval_3",
    # "ours_thermal_time",
    # "ours_interval_5",
    # "ours_interval_10",
    # "normal",
    "ours",
    "ours_reduction",
    # "mlt++",
    "supereduction_a=2",
    # "supereduction_a=2_post",
    "supereduction_a=2_reduction",
    # "help",
    # "ours_interval_40",
    # "supereduction_a=5",
]
#'supereduction_a=5' empty it doesnt run
# Initialize the plot
plt.figure()

MASS = 60

import utils


def plot_hr_diagram(mass: int, names: list[str], phase: str = ""):
    base_colors = {}
    phase_suffix = f"_{phase}" if phase else ""
    for name in names:
        path = f"{mass}m-{name}/LOGS/history.data"

        try:
            # Load the history data
            history = mr.MesaData(path)

            kwargs, is_winds, base_name = utils.get_line_kwargs(name, base_colors)

            # Plot log_L vs. log_Teff
            (line,) = plt.plot(history.log_Teff, history.log_L, **kwargs)
            color = line.get_color()

            if not is_winds:
                base_colors[base_name] = color

            plt.plot(history.log_Teff[-1], history.log_L[-1], "o", color=color)
        except FileNotFoundError:
            print(f"Warning: {path}.data not found. Skipping.")

    # Customize the plot
    plt.xlabel(r"$\log (T_{\text{eff}}/\text{K})$")
    plt.ylabel(r"$\log (L/L_{\sun})$")
    # plt.title(f"{mass} Mass Star Hertzsprung-Russell Diagram")
    plt.gca().invert_xaxis()  # HR diagrams have Teff decreasing to the right
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    # Save the plot
    os.makedirs(f"plot{mass}m{phase_suffix}", exist_ok=True)
    plt.savefig(f"plot{mass}m{phase_suffix}/mass_{mass}_hr_diagram_comparison.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    plot_hr_diagram(MASS, NAMES)
