import json

import generic_last_log_plot
import HEAT_MAP
import hr_diagram_script

folders_for_graphs = {
    # 30: ["ours", "supereduction_a=2", "mlt++", "normal", "supereduction_a=5"],
    # 40: ["ours", "supereduction_a=2", "mlt++", "normal"],
    # 50: ["ours", "supereduction_a=2", "mlt++", "normal"],
    # 60: ["ours", "supereduction_a=2", "mlt++", "ours_reduction", "supereduction_a=2_reduction"],
    80: [
        # "ours",
        # "ours_post",
        # "supereduction_a=2",
        # "supereduction_a=2_post",
        # "mlt_tdc_only",
        # "mlt_tdc_only_post",
        "ours_reduction",
        # "mlt++",
        # "",
        # "ours_post"
        # "supereduction_a=2_post",
        # "ours_reduction",
        # "supereduction_a=2_reduction",
        # "ours_low_alpha",
    ],  # ours and _post run start at the end of main sequence and symolate the end of the life, cause our model have problem running with no winds for 80 mass stars
}

for mass, names in folders_for_graphs.items():
    print("\n")
    print(f"start graphs for {mass}M☉ folders")
    print(json.dumps(names, ensure_ascii=False))
    generic_last_log_plot.plot_generic_last_log_plot(mass, names)
    generic_last_log_plot.plot_generic_last_log_plot(mass, names, x_axis="logR", x_units="R☉")
    hr_diagram_script.plot_hr_diagram(mass, names)

    for name in names:
        if "ours" in name:
            HEAT_MAP.plot_heat_map(mass, name)
        HEAT_MAP.plot_heat_map(
            mass, name, value_axis="L_div_Ledd_effective", normalize=False
        )

    print("\n")
