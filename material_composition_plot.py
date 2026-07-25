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


def _get_valid_isos(prof):
    valid = []
    for iso in KEY_ISOS:
        try:
            getattr(prof, iso)
            valid.append(iso)
        except Exception:
            pass
    return valid


def _get_last_profile_num(path):
    try:
        hist_dir = mr.MesaLogDir(path)
        return hist_dir.profile_numbers[-1]
    except Exception:
        import glob

        prof_files = sorted(
            glob.glob(f"{path}/profile*.data"),
            key=lambda x: int(x.split("profile")[-1].split(".data")[0]),
        )
        return int(prof_files[-1].split("profile")[-1].split(".data")[0])


def _setup_ax(ax):
    ax.set_ylabel("Mass Fraction")
    ax.set_yscale("symlog", linthresh=0.01)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="best")


def _finalize_ax(ax):
    ax.set_ylim(bottom=0, top=1)


def _plot_abundance_vs_logr(mass: int, name: str, phase: str = ""):
    path = f"{mass}m-{name}/LOGS"
    out_dir = f"plot{mass}m_material_composition{'_' + phase if phase else ''}"
    os.makedirs(out_dir, exist_ok=True)

    last_profile_num = _get_last_profile_num(path)
    prof = mr.MesaData(f"{path}/profile{last_profile_num}.data")

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
        # ax.set_title(f'{utils.parse_name_for_plots(name)} \u2014 Mass {mass} M\u2609')

    plt.tight_layout()
    plt.savefig(f"{out_dir}/Mass_{mass}_{name}_abundance_vs_logR.png", dpi=300)
    plt.close()
    print(f"Saved: {out_dir}/Mass_{mass}_{name}_abundance_vs_logR.png")


def _plot_abundance_vs_logrho(mass: int, name: str, phase: str = ""):
    path = f"{mass}m-{name}/LOGS"
    out_dir = f"plot{mass}m_material_composition{'_' + phase if phase else ''}"
    os.makedirs(out_dir, exist_ok=True)

    last_profile_num = _get_last_profile_num(path)
    prof = mr.MesaData(f"{path}/profile{last_profile_num}.data")

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
        # ax.set_title(f'{utils.parse_name_for_plots(name)} \u2014 Mass {mass} M\u2609')

    plt.tight_layout()
    plt.savefig(f"{out_dir}/Mass_{mass}_{name}_abundance_vs_logRho.png", dpi=300)
    plt.close()
    print(f"Saved: {out_dir}/Mass_{mass}_{name}_abundance_vs_logRho.png")


def _plot_abundance_evolution(mass: int, name: str, phase: str = ""):
    path = f"{mass}m-{name}/LOGS"
    out_dir = f"plot{mass}m_material_composition{'_' + phase if phase else ''}"
    os.makedirs(out_dir, exist_ok=True)

    last_profile_num = _get_last_profile_num(path)
    prof = mr.MesaData(f"{path}/profile{last_profile_num}.data")
    valid_isos = _get_valid_isos(prof)
    if not valid_isos:
        return

    model_numbers = []
    star_ages = []
    iso_avg = {iso: [] for iso in valid_isos}

    for i in range(1, 9999):
        f_path = f"{path}/profile{i}.data"
        if not os.path.exists(f_path):
            break
        try:
            p = mr.MesaData(f_path)
        except Exception:
            continue

        try:
            model_numbers.append(float(getattr(p, "model_number")))
            star_ages.append(float(getattr(p, "star_age")))
            dm = np.array(getattr(p, "dm"), dtype=float)
            total_mass = np.sum(dm)
        except Exception:
            continue

        for iso in valid_isos:
            try:
                iso_data = np.array(getattr(p, iso), dtype=float)
                iso_avg[iso].append(np.sum(iso_data * dm) / total_mass)
            except AttributeError:
                iso_avg[iso].append(np.nan)

    if not model_numbers:
        return

    model_numbers = np.array(model_numbers)
    star_ages = np.array(star_ages)
    used_isos = [iso for iso in valid_isos if not all(np.isnan(v) for v in iso_avg[iso])]
    for iso in valid_isos:
        iso_avg[iso] = np.array(iso_avg[iso])

    sort_idx = np.argsort(model_numbers)
    model_numbers = model_numbers[sort_idx]
    star_ages = star_ages[sort_idx]
    for iso in valid_isos:
        iso_avg[iso] = iso_avg[iso][sort_idx]

    fig, ax = plt.subplots(figsize=(8, 6))
    for iso in used_isos:
        ax.plot(
            model_numbers,
            iso_avg[iso],
            label=ISO_LABELS.get(iso, iso),
            color=ISO_COLORS.get(iso, "#333333"),
            linewidth=1.2,
        )
    ax.set_xlabel("Model Number")
    _setup_ax(ax)
    _finalize_ax(ax)
    # ax.set_title(f'{utils.parse_name_for_plots(name)} \u2014 Mass {mass} M\u2609')

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


def plot(mass: int, name: str, phase: str = ""):
    _plot_abundance_vs_logr(mass, name, phase=phase)
    _plot_abundance_vs_logrho(mass, name, phase=phase)
    _plot_abundance_evolution(mass, name, phase=phase)


def plot_all(mass: int, names: list[str]):
    for name in names:
        try:
            plot(mass, name)
        except Exception as e:
            print(f"Error plotting {name}: {e}")


if __name__ == "__main__":
    for mass in [80]:
        names = ["ours", "ours_reduction"]
        plot_all(mass, names)
