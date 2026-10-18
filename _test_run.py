"""Temporary test runner — runs the main pipeline for masses 10-50 only."""
import json

import evolution_graph
import generic_last_log_plot
import HEAT_MAP
import hr_diagram_script
import material_composition_plot
from stream import stream_profiles


def _get_phase(name: str) -> str:
    if name.endswith("_post"):
        return "post_ms"
    elif name.endswith("_reduction"):
        return "pre_ms"
    return ""


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
}

for mass, names in folders_for_graphs.items():
    print("\n")
    print(f"=== {mass}M\u2609 ===")
    print(json.dumps(names, ensure_ascii=False))

    grouped = [(name, _get_phase(name)) for name in names]

    by_phase: dict[str, list[str]] = {}
    for name, phase in grouped:
        by_phase.setdefault(phase, []).append(name)

    # Comparison plots
    for phase, phase_names in by_phase.items():
        if not phase:
            continue
        print(f"  comparison plots phase={phase}: {phase_names}")

        preloaded = generic_last_log_plot.load_last_profiles(mass, phase_names)

        try:
            generic_last_log_plot.plot(
                mass, phase_names,
                x_axis="R", x_units="", phase=phase,
                y_axis="logRho", y_units="g/cm³",
                _preloaded=preloaded,
            )
            generic_last_log_plot.plot_generic_last_log_plot(
                mass, phase_names, phase=phase, _preloaded=preloaded
            )
        except Exception as e:
            print(f"    generic_last_log_plot failed: {e}")

        try:
            generic_last_log_plot.plot_generic_last_log_plot(
                mass, phase_names,
                x_axis="logR", x_units="", phase=phase,
                inverse_x=True,
                _preloaded=preloaded,
            )
        except Exception as e:
            print(f"    generic_last_log_plot(logR) failed: {e}")

        try:
            hr_diagram_script.plot_hr_diagram(mass, phase_names, phase=phase)
        except Exception as e:
            print(f"    hr_diagram_script failed: {e}")

        try:
            evolution_graph.plot(mass, phase_names, y_axis="mass", phase=phase)
        except Exception as e:
            print(f"    evolution_graph failed: {e}")

    # Individual plots — single streaming pass per name
    for name, phase in grouped:
        print(f"  individual: {name} (phase={phase or 'none'})")

        logs_path = f"{mass}m-{name}/LOGS"
        phase_suffix = f"_{phase}" if phase else ""

        if "ours" in name:
            hm_value_axes = ["extra_opacity_factor", "logT", "L_div_Ledd_effective"]
            hm_normalize_map = {
                "extra_opacity_factor": True,
                "logT": False,
                "L_div_Ledd_effective": False,
            }
        else:
            hm_value_axes = ["L_div_Ledd_effective"]
            hm_normalize_map = {"L_div_Ledd_effective": False}

        hm_collector = HEAT_MAP.HeatMapCollector(value_axes=hm_value_axes)
        evo_collector = material_composition_plot.CompositionEvolutionCollector()

        print(f"    streaming {logs_path}...")
        stream_profiles(logs_path, [hm_collector, evo_collector])

        last_prof = material_composition_plot.load_last_profile(mass, name)
        if last_prof is not None:
            material_composition_plot.plot_snapshots(last_prof, mass, name, phase=phase)

        try:
            evo_collector.render(mass, name, phase=phase)
        except Exception as e:
            print(f"    composition evolution failed: {e}")

        try:
            hm_collector.render_all(mass, name, normalize_map=hm_normalize_map, phase=phase)
        except Exception as e:
            print(f"    HEAT_MAP({name}) failed: {e}")

print("\nDone.")
