OURS_NAME = "ours model"

def parse_name_for_plots(name: str) -> str:
    plot_name = name.replace("ours", OURS_NAME) if "ours" in name else name
    plot_name = plot_name.replace("_reduction", "_winds")
    plot_name = plot_name.replace("_", " ")
    plot_name = plot_name.replace("-", " ")
    plot_name = plot_name.capitalize()
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
