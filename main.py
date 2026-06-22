import json

import generic_last_log_plot
import HEAT_MAP
import hr_diagram_script
import material_composition_plot


def _get_phase(name: str) -> str:
    if name.endswith("_post"):
        return "post_ms"
    elif name.endswith("_reduction"):
        return "pre_ms"
    return ""


def _group_names(names: list[str]) -> list[tuple[str, str]]:
    groups: list[tuple[str, str]] = []
    for name in names:
        phase = _get_phase(name)
        groups.append((name, phase))
    return groups


folders_for_graphs = {
    10: [
        "ours_reduction",
        "ours_reduction_post",
        "supereduction_a=2_reduction",
        "supereduction_a=2_reduction_post",
    ],
    20: [
        "ours_reduction",
        "ours_reduction_post",
        "supereduction_a=2_reduction",
        "supereduction_a=2_reduction_post",
    ],
    30: [
        "ours_reduction",
        "ours_reduction_post",
        "supereduction_a=2_reduction",
        "supereduction_a=2_reduction_post",
        "mlt++_reduction",
        "mlt++_reduction_post",
        "normal_reduction",
        "normal_reduction_post",
    ],
    40: [
        "ours_reduction",
        "ours_reduction_post",
        "supereduction_a=2_reduction",
        "supereduction_a=2_reduction_post",
    ],
    50: [
        "ours_reduction",
        "ours_reduction_post",
        "supereduction_a=2_reduction",
        "supereduction_a=2_reduction_post",
    ],
    60: [
        "ours_reduction",
        "ours_reduction_post",
        "supereduction_a=2_reduction",
        "supereduction_a=2_reduction_post",
    ],
    70: [
        "ours_reduction",
        "ours_reduction_post",
        "supereduction_a=2_reduction",
        "supereduction_a=2_reduction_post",
    ],
    80: [
        "ours_reduction",
        "ours_reduction_post",
        "supereduction_a=2_reduction",
        "supereduction_a=2_reduction_post",
    ],
    90: [
        "ours_reduction",
        "ours_reduction_post",
        "supereduction_a=2_reduction",
        "supereduction_a=2_reduction_post",
    ],
    100: [
        "ours_reduction",
        "ours_reduction_post",
        "supereduction_a=2_reduction",
        "supereduction_a=2_reduction_post",
    ],
}

for mass, names in folders_for_graphs.items():
    print("\n")
    print(f"start graphs for {mass}M\u2609 folders")
    print(json.dumps(names, ensure_ascii=False))

    grouped = _group_names(names)

    # Group by phase for comparison plots
    by_phase: dict[str, list[str]] = {}
    for name, phase in grouped:
        by_phase.setdefault(phase, []).append(name)

    # Comparison plots: all names per phase together
    for phase, phase_names in by_phase.items():
        if not phase:
            continue
        print(f"  comparison plots for phase={phase}: {phase_names}")

        try:
            generic_last_log_plot.plot_generic_last_log_plot(
                mass, phase_names, phase=phase
            )
        except Exception as e:
            print(f"    generic_last_log_plot failed: {e}")

        try:
            generic_last_log_plot.plot_generic_last_log_plot(
                mass, phase_names, x_axis="logR", x_units="R\u2609", phase=phase
            )
        except Exception as e:
            print(f"    generic_last_log_plot(logR) failed: {e}")

        try:
            hr_diagram_script.plot_hr_diagram(mass, phase_names, phase=phase)
        except Exception as e:
            print(f"    hr_diagram_script failed: {e}")

    # Individual per-folder plots
    for name, phase in grouped:
        print(f"  individual plots for {name} (phase={phase or 'none'})...")

        material_composition_plot.plot(mass, name, phase=phase)

        if "ours" in name:
            try:
                HEAT_MAP.plot_heat_map(mass, name, phase=phase)
            except Exception as e:
                print(f"    HEAT_MAP({name}) failed: {e}")

        try:
            HEAT_MAP.plot_heat_map(
                mass,
                name,
                value_axis="L_div_Ledd_effective",
                normalize=False,
                phase=phase,
            )
        except Exception as e:
            print(f"    HEAT_MAP({name}, L_div_Ledd) failed: {e}")


    print("\n")
