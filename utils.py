OURS_NAME = "Porosity Reduction"
SUPEREDUCTION_NAME = "Super Reduction"


REPLACE_DICT = {
    "supereduction": SUPEREDUCTION_NAME,
    "_reduction": "_winds",
    "_": " ",
    "-": " ",
    "a=": "α=",
    "mlt++": "MLT++",
    " div ": "/",
    "Ledd": r"$L_{edd}$",
} 

def parse_name_for_plots(name: str) -> str:
    # plot_name = name.lower()
    plot_name = name.replace("ours", OURS_NAME) if "ours" in name else name

    for old, new in REPLACE_DICT.items():
        plot_name = plot_name.replace(old, new)
        
    if plot_name.startswith("log"):
        idx = plot_name[2:].find("log") or -1
        inner = plot_name[3 + idx:] if idx != -1 else plot_name[3:]
        recursive = parse_name_for_plots(plot_name[3 + idx:]) if idx != -1 else ""
        plot_name = rf"$\log ({inner})$" + recursive
    else:
        plot_name = plot_name.title()
    return plot_name


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
