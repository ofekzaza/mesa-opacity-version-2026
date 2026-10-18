OURS_NAME = "Porosity Reduction"

PLOTS_ROOT = "plotted_results"
SUPEREDUCTION_NAME = "Super Reduction"


REPLACE_DICT = {
    "supereduction": SUPEREDUCTION_NAME,
    "_reduction": "_winds",
    "_": " ",
    "-": " ",
    "a=": "α=",
    "mlt++": "MLT++",
    "Ledd": r"$L_{edd}$",
    "normal": "default",
}


def parse_name_for_plots(name: str) -> str:
    # plot_name = name.lower()
    plot_name = name.replace("ours", OURS_NAME) if "ours" in name else name

    for old, new in REPLACE_DICT.items():
        plot_name = plot_name.replace(old, new)

    plot_name = plot_name.title()
    return plot_name


pretty_axises = {
    "L_div_Ledd_effective": r"$L/L_{edd}$ effective",
    "log_L_div_Ledd": r"$\log(L/L_{edd})$",
    "log_opacity": r"$\log(\kappa)$",
    "extra_opacity_factor": "extra opacity factor",
    "logT": r"$\log(T)$ [K]",
    "entropy": r"$S$",
    "mass": r"$M/M_{\odot}$",
    "logR": r"$\log(R/R_{\odot})$",
    "logRho": r"$\log(\rho)$",
    "logP": r"$\log(P)$",
    "star_age": "Star Age [yr]",
    "star_mass": r"$M/M_{\odot}$",
    "logL": r"$\log(L/L_{\odot})$",
    "R": r"$R/R_{\odot}$",
}


def pretty_axis_name(axis: str) -> str:
    return pretty_axises.get(axis, axis)


def get_last_profile_num(logs_path: str) -> int:
    """Return the last profile number in a LOGS directory."""
    import glob as _glob
    import mesa_reader as _mr
    try:
        return _mr.MesaLogDir(logs_path).profile_numbers[-1]
    except Exception:
        prof_files = sorted(
            _glob.glob(f"{logs_path}/profile*.data"),
            key=lambda x: int(x.split("profile")[-1].split(".data")[0]),
        )
        if not prof_files:
            raise FileNotFoundError(f"No profiles found in {logs_path}")
        return int(prof_files[-1].split("profile")[-1].split(".data")[0])


def get_line_kwargs(name: str, base_colors: dict) -> tuple[dict, bool, str]:
    plot_name = parse_name_for_plots(name)
    is_winds = "_reduction" in name
    base_name = name.replace("_reduction", "")

    kwargs = {"label": plot_name, "linewidth": 0.8}
    if is_winds and base_name in base_colors:
        kwargs["color"] = base_colors[base_name]
        kwargs["linestyle"] = "--"
    elif is_winds:
        kwargs["linestyle"] = "--"
    else:
        kwargs["linestyle"] = "-"

    return kwargs, is_winds, base_name
