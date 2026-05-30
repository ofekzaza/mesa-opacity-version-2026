import mesa_reader as mr
import matplotlib.pyplot as plt
import os

mass = 50

# List of history file names (without extensions)
names = [
    # "normal",
    # "ours_interval_2",
    # "ours_interval_3",
    # "ours_thermal_time",
    # "ours_interval_5",
    # "ours_interval_10",
    "normal",
    "ours",
    "mlt++",
    "supereduction_a=2",
    # "ours_interval_40",
    # "supereduction_a=5",
]
#'supereduction_a=5' empty it doesnt run
# Initialize the plot
plt.figure()

# Iterate over each history file
for name in names:
    path = f"{mass}m-{name}/LOGS"
    try:
        print(path)
        # Load the history data
        history = mr.MesaLogDir(path)
        # Plot log_L vs. log_Teff
        prof = mr.MesaData(f"{path}/profile{history.profile_numbers[-1]}.data")
        (line,) = plt.plot(prof.logRho, prof.logT, "-", linewidth=0.8, label=name)
        color = line.get_color()
        plt.plot(prof.logRho[-1], prof.logT[-1], "o", markersize=6, color=color)
    except FileNotFoundError:
        print(f"Warning: {path}.data not found. Skipping.")

# Customize the plot
plt.xlabel(r"$\log(\rho)\ \,[\mathrm{g/cm^3}]$")
plt.ylabel(r"$\log(T)\ \,[\mathrm{K}]$")
plt.title(f"mass {mass}: Temperature vs Density at onset of Carbon burning")
# plt.gca().invert_xaxis()  # HR diagrams have Teff decreasing to the right
plt.grid(True)
plt.legend()
plt.tight_layout()
os.makedirs(f"plot{mass}m", exist_ok=True)
# Save the plot
plt.savefig(f"plot{mass}m/Mass_{mass}_Temperature_vs_Density.png", dpi=300)
plt.close()
