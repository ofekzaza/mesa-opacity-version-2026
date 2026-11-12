import json

import matplotlib

# Use a headless backend before any plotting module imports pyplot.  On macOS
# the default ``macosx`` backend starts AppKit/XPC listeners, which emit
# "NSXPCSharedListener ... Connection interrupted" warnings.  This must happen
# before the imports below, since they import ``matplotlib.pyplot``.
matplotlib.use("Agg")

import evolution_graph
import generic_last_log_plot
import HEAT_MAP
import hr_diagram_script
import hr_mass_comparison
import kippenhahn_plot
import logT_logRho_plot
import material_composition_plot
import radius_time_plot
from stream import stream_profiles


def _get_phase(name: str) -> str:
    if name.endswith("_post"):
        return "post_ms"
    elif name.endswith("_suite"):
        return "suite"
    elif name.endswith("_default"):
        return "default"
    elif name.endswith("_reduction"):
        return "pre_ms"
    return ""


def _group_names(names: list[str]) -> list[tuple[str, str]]:
    return [(name, _get_phase(name)) for name in names]


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
        "ours_reduction_suite",
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
        "ours_reduction_default",
        "ours_reduction_post_default",
        "supereduction_a=2_reduction_default",
        "supereduction_a=2_reduction_post_default",
        "mlt++_reduction_default",
        "mlt++_reduction_post_default",
        "normal_reduction_default",
        "normal_reduction_post_default",
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
        # "supereduction_a=2_reduction",
        # "supereduction_a=2_reduction_post",
    ],
    100: [
        "ours_reduction",
        "ours_reduction_post",
        # "supereduction_a=2_reduction",
        # "supereduction_a=2_reduction_post",
    ],
    300: [
        "ours_reduction",
    ],
}

# ----------------------------------------------------------------------
# Cross-mass HR comparison for the "ours" models (30-100 M☉): one figure
# for the pre-MS runs and one for the post-MS runs.
# ----------------------------------------------------------------------
print("Cross-mass HR comparison for the our models (30-100 M\u2609)")
try:
    hr_mass_comparison.plot_all()
except Exception as e:
    print(f"    hr_mass_comparison failed: {e}")

# ----------------------------------------------------------------------
# Cross-mass radius-versus-time comparison for the "ours" models
# (30-100 M☉): one figure for the pre-MS runs and one for the post-MS runs.
# ----------------------------------------------------------------------
print("Cross-mass radius-time comparison for the our models (30-100 M\u2609)")
try:
    radius_time_plot.plot_all()
except Exception as e:
    print(f"    radius_time_plot failed: {e}")



for mass, names in folders_for_graphs.items():
    print("\n")
    print(f"start graphs for {mass}M\u2609 folders")
    print(json.dumps(names, ensure_ascii=False))

    grouped = _group_names(names)

    # Group by phase for comparison plots
    by_phase: dict[str, list[str]] = {}
    for name, phase in grouped:
        by_phase.setdefault(phase, []).append(name)

    # ------------------------------------------------------------------
    # Comparison plots — one last-profile read per name, shared across
    # all generic_last_log_plot / hr_diagram / evolution_graph calls
    # ------------------------------------------------------------------
    for phase, phase_names in by_phase.items():
        if not phase:
            continue
        print(f"  comparison plots for phase={phase}: {phase_names}")

        # Load last profiles once; reuse across all 13 plot() calls below
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
            logT_logRho_plot.plot(
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

        try:
            evolution_graph.plot_percentages(mass, phase_names, phase=phase)
        except Exception as e:
            print(f"    evolution_graph percentages failed: {e}")

        try:
            evolution_graph.plot_taus(mass, phase_names, phase=phase)
        except Exception as e:
            print(f"    evolution_graph taus failed: {e}")

    # ------------------------------------------------------------------
    # Individual per-folder plots
    # One streaming pass per name feeds the HeatMap, CompositionEvolution and
    # Kippenhahn collectors — profiles are read exactly once.
    # ------------------------------------------------------------------
    for name, phase in grouped:
        print(f"  individual plots for {name} (phase={phase or 'none'})...")

        logs_path = f"{mass}m-{name}/LOGS"
        phase_suffix = f"_{phase}" if phase else ""

        # Determine which heatmap value_axes are needed
        if "ours" in name:
            hm_value_axes = [
                "extra_opacity_factor",
                "logT",
                "L_div_Ledd_effective",
            ]
            hm_normalize_map = {
                "extra_opacity_factor": True,
                "logT": False,
                "L_div_Ledd_effective": False,
            }
        else:
            hm_value_axes = ["L_div_Ledd_effective"]
            hm_normalize_map = {"L_div_Ledd_effective": False}

        # Set up streaming collectors
        hm_collector = HEAT_MAP.HeatMapCollector(value_axes=hm_value_axes)
        evo_collector = material_composition_plot.CompositionEvolutionCollector()
        kipp_collector = kippenhahn_plot.KippenhahnCollector()

        # Single pass through all profiles — all collectors consume each profile
        print(f"    streaming {logs_path}...")
        stream_profiles(logs_path, [hm_collector, evo_collector, kipp_collector])

        # Load last profile once for snapshot plots
        last_prof = material_composition_plot.load_last_profile(mass, name)
        if last_prof is not None:
            material_composition_plot.plot_snapshots(last_prof, mass, name, phase=phase)

        # Render composition evolution from collected scalars
        try:
            evo_collector.render(mass, name, phase=phase)
        except Exception as e:
            print(f"    composition evolution({name}) failed: {e}")

        # Render heatmaps from collected arrays
        try:
            hm_collector.render_all(
                mass, name,
                normalize_map=hm_normalize_map,
                phase=phase,
            )
        except Exception as e:
            print(f"    HEAT_MAP({name}) failed: {e}")

        # Render Kippenhahn diagrams from the same streaming pass
        try:
            kipp_collector.render(mass, name, phase=phase)
        except Exception as e:
            print(f"    KIPPENHAHN({name}) failed: {e}")

    print("\n")

