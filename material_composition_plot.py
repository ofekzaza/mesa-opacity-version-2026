import os

import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np

import utils

KEY_ISOS = ["h1", "he4", "c12", "n14", "o16", "ne20", "mg24", "si28", "s32", "ar36", "ca40", "fe52", "fe54", "fe56"]
ISO_LABELS = {
    "h1": r"H",
    "he3": r"$^3$He",
    "he4": r"$^4$He",
    "c12": r"$^{12}$C",
    "n14": r"$^{14}$N",
    "o16": r"$^{16}$O",
    "ne20": r"$^{20}$Ne",
    "mg24": r"$^{24}$Mg",
    "si28": r"$^{28}$Si",
    "s32": r"$^{32}$S",
    "ar36": r"$^{36}$Ar",
    "ca40": r"$^{40}$Ca",
    "fe52": r"$^{52}$Fe",
    "fe54": r"$^{54}$Fe",
    "fe56": r"$^{56}$Fe",
}
ISO_COLORS = {
    "h1": "#1f77b4",
    "he3": "#aec7e8",
    "he4": "#ff7f0e",
    "c12": "#2ca02c",
    "n14": "#d62728",
    "o16": "#9467bd",
    "ne20": "#8c564b",
    "mg24": "#e377c2",
    "si28": "#7f7f7f",
    "s32": "#bcbd22",
    "ar36": "#17becf",
    "ca40": "#9edae5",
    "fe52": "#dbdb8d",
    "fe54": "#c5b0d5",
    "fe56": "#f15854",
}


def _get_valid_isos(prof) -> list[str]:
    valid = []
    for iso in KEY_ISOS:
        try:
            getattr(prof, iso)
            valid.append(iso)
        except Exception:
            pass
    return valid


def _setup_ax(ax):
    ax.set_ylabel("Mass Fraction")
    ax.set_yscale("symlog", linthresh=0.01)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="best")


def _finalize_ax(ax):
    ax.set_ylim(bottom=0, top=1)


# ---------------------------------------------------------------------------
# Streaming collector
# ---------------------------------------------------------------------------

class CompositionEvolutionCollector:
    """Accumulates mass-weighted isotope averages (scalars) during a streaming pass.

    Only scalars are held in memory after each profile — the full arrays are
    reduced immediately inside consume() before the MesaData object is discarded.
    """

    def __init__(self):
        self.model_numbers: list[float] = []
        self.star_ages: list[float] = []
        self.valid_isos: list[str] = []
        self.iso_avgs: dict[str, list[float]] = {}
        self._initialized = False

    def consume(self, prof) -> None:
        # Determine which isotopes are present from the first profile
        if not self._initialized:
            self.valid_isos = [iso for iso in KEY_ISOS if hasattr(prof, iso)]
            self.iso_avgs = {iso: [] for iso in self.valid_isos}
            self._initialized = True

        try:
            mn = float(getattr(prof, "model_number"))
            sa = float(getattr(prof, "star_age"))
            dm = np.array(getattr(prof, "dm"), dtype=float)
            total_mass = float(np.sum(dm))
        except Exception:
            return

        self.model_numbers.append(mn)
        self.star_ages.append(sa)
        for iso in self.valid_isos:
            try:
                iso_data = np.array(getattr(prof, iso), dtype=float)
                self.iso_avgs[iso].append(float(np.sum(iso_data * dm) / total_mass))
            except Exception:
                self.iso_avgs[iso].append(float("nan"))

    def render(self, mass: int, name: str, phase: str = "") -> None:
        if not self.model_numbers:
            return
        out_dir = f"{utils.PLOTS_ROOT}/plot{mass}m_material_composition{'_' + phase if phase else ''}"
        os.makedirs(out_dir, exist_ok=True)
        _plot_abundance_evolution_from_data(
            self.model_numbers,
            self.star_ages,
            self.iso_avgs,
            self.valid_isos,
            mass,
            name,
            phase,
            out_dir,
        )


# ---------------------------------------------------------------------------
# Snapshot plots (operate on a single pre-loaded last profile)
# ---------------------------------------------------------------------------

def _plot_abundance_vs_logr(prof, mass: int, name: str, phase: str, out_dir: str) -> None:
    valid_isos = _get_valid_isos(prof)
    if not valid_isos:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    x_data = getattr(prof, "logR", None)
    if x_data is not None:
        for iso in valid_isos:
            ax.plot(
                x_data,
                getattr(prof, iso),
                label=ISO_LABELS.get(iso, iso),
                color=ISO_COLORS.get(iso, "#333333"),
                linewidth=1.2,
            )
        ax.set_xlabel("log(R/R\u2609)")
        _setup_ax(ax)
        _finalize_ax(ax)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/Mass_{mass}_{name}_abundance_vs_logR.png", dpi=300)
    plt.close()
    print(f"Saved: {out_dir}/Mass_{mass}_{name}_abundance_vs_logR.png")


def _plot_abundance_vs_logrho(prof, mass: int, name: str, phase: str, out_dir: str) -> None:
    valid_isos = _get_valid_isos(prof)
    if not valid_isos:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    x_data = getattr(prof, "logRho", None)
    if x_data is not None:
        for iso in valid_isos:
            ax.plot(
                x_data,
                getattr(prof, iso),
                label=ISO_LABELS.get(iso, iso),
                color=ISO_COLORS.get(iso, "#333333"),
                linewidth=1.2,
            )
        ax.set_xlabel("log(\u03c1) [g/cm\u00b3]")
        _setup_ax(ax)
        _finalize_ax(ax)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/Mass_{mass}_{name}_abundance_vs_logRho.png", dpi=300)
    plt.close()
    print(f"Saved: {out_dir}/Mass_{mass}_{name}_abundance_vs_logRho.png")


def _plot_abundance_evolution_from_data(
    model_numbers, star_ages, iso_avgs, valid_isos, mass, name, phase, out_dir
) -> None:
    model_numbers = np.array(model_numbers)
    star_ages = np.array(star_ages)

    for iso in valid_isos:
        iso_avgs[iso] = np.array(iso_avgs[iso])

    used_isos = [iso for iso in valid_isos if not all(np.isnan(iso_avgs[iso]))]

    sort_idx = np.argsort(model_numbers)
    model_numbers = model_numbers[sort_idx]
    star_ages = star_ages[sort_idx]
    for iso in valid_isos:
        iso_avgs[iso] = iso_avgs[iso][sort_idx]

    fig, ax = plt.subplots(figsize=(8, 6))
    for iso in used_isos:
        ax.plot(
            model_numbers,
            iso_avgs[iso],
            label=ISO_LABELS.get(iso, iso),
            color=ISO_COLORS.get(iso, "#333333"),
            linewidth=1.2,
        )
    ax.set_xlabel("Model Number")
    _setup_ax(ax)
    _finalize_ax(ax)

    if len(model_numbers) > 1:
        sec_ax = ax.secondary_xaxis("top")
        n_ticks = min(5, len(model_numbers))
        idx = np.linspace(0, len(model_numbers) - 1, n_ticks, dtype=int)
        sec_ax.set_xticks(model_numbers[idx])
        sec_ax.set_xticklabels(
            [f"{star_ages[i]:.2e}" for i in idx], fontsize=7, rotation=45
        )
        sec_ax.set_xlabel("Star Age [yr]", fontsize=8)

    plt.tight_layout()
    plt.savefig(f"{out_dir}/Mass_{mass}_{name}_abundance_evolution.png", dpi=300)
    plt.close()
    print(f"Saved: {out_dir}/Mass_{mass}_{name}_abundance_evolution.png")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def load_last_profile(mass: int, name: str) -> "mr.MesaData | None":
    """Load the last profile for a run exactly once."""
    path = f"{mass}m-{name}/LOGS"
    try:
        num = utils.get_last_profile_num(path)
        return mr.MesaData(f"{path}/profile{num}.data")
    except Exception as e:
        print(f"    material_composition: could not load last profile for {name}: {e}")
        return None


def plot_snapshots(prof, mass: int, name: str, phase: str = "") -> None:
    """Render abundance-vs-logR and abundance-vs-logRho from a pre-loaded profile."""
    out_dir = f"{utils.PLOTS_ROOT}/plot{mass}m_material_composition{'_' + phase if phase else ''}"
    os.makedirs(out_dir, exist_ok=True)
    _plot_abundance_vs_logr(prof, mass, name, phase, out_dir)
    _plot_abundance_vs_logrho(prof, mass, name, phase, out_dir)


def plot(mass: int, name: str, phase: str = "") -> None:
    """Standalone entry point — loads last profile once, streams all profiles once."""
    from stream import stream_profiles

    last_prof = load_last_profile(mass, name)
    if last_prof is not None:
        plot_snapshots(last_prof, mass, name, phase=phase)

    evo_collector = CompositionEvolutionCollector()
    path = f"{mass}m-{name}/LOGS"
    stream_profiles(path, [evo_collector])
    evo_collector.render(mass, name, phase=phase)


def plot_all(mass: int, names: list[str]) -> None:
    for name in names:
        try:
            plot(mass, name)
        except Exception as e:
            print(f"Error plotting {name}: {e}")


if __name__ == "__main__":
    for mass in [80]:
        names = ["ours", "ours_reduction"]
        plot_all(mass, names)
