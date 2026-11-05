"""Kippenhahn diagrams for the MESA runs in this repository.

The plot is built from the saved profiles (``LOGS/profile*.data``).  It shows
the enclosed mass on the vertical axis and the star time to live on the
horizontal axis.

Encoded regions
---------------
* Nuclear burning, split by the dominant fuel/reaction group
  (H, He, C, Ne, O, Si) using the per-zone energy-generation columns
  (``pp``, ``cno``, ``tri_alpha``, ...).
* Mixing regions from ``mix_type`` (convection, overshoot, semiconvection,
  thermohaline, rotation, ...).
* Temperature inversions (``gradT < 0``).
* Porosity reduction regions (``extra_opacity_factor < 1``), when the run
  has this column.

The heavy use of hatched, unfilled overlays means the burning colours stay
readable even where several regions overlap.

When all burning is confined to a very thin layer (for example the shell of a
nova on a white dwarf), an extra pair of ``_zoom`` plots is produced around
that layer so it does not vanish against the surface line.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

import utils

# Keep the hatched overlays light so they do not hide the burning colours.
plt.rcParams["hatch.linewidth"] = 0.7

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# MESA mixing_type values (see const_def.f90).
MIX_INFO = {
    -1: ("Crystallised", "#9E9E9E", ".."),
    1: ("Convective", "#D62728", "///"),
    2: ("Overshoot", "#FF7F0E", "\\\\\\"),
    3: ("Semiconvection", "#17BECF", "xx"),
    4: ("Thermohaline", "#E377C2", "oo"),
    5: ("Rotation", "#1F77B4", "++"),
    6: ("Rayleigh-Taylor", "#7F7F7F", "**"),
    7: ("Minimum", "#BCBD22", ".."),
    8: ("Anonymous", "#000000", "--"),
    9: ("Leftover convective", "#9467BD", "//"),
    10: ("Phase separation", "#8C564B", "OO"),
}

# Reaction columns are grouped into the standard burning stages.  Columns not
# present in a given profile are simply skipped.
BURN_GROUPS = [
    ("H", ["pp", "cno"], "#4C72B0"),
    ("He", ["tri_alpha", "n_alpha", "c_alpha"], "#DD8452"),
    ("C", ["c12_c12", "c12_o16"], "#55A868"),
    ("Ne", ["ne_alpha", "o_alpha", "na_alpha", "photo", "pnhe4"], "#8172B3"),
    ("O", ["o16_o16"], "#C44E52"),
    (
        "Si",
        [
            "mg_alpha",
            "si_alpha",
            "s_alpha",
            "ar_alpha",
            "ca_alpha",
            "ti_alpha",
            "cr_alpha",
            "fe_co_ni",
            "ni56_co56",
            "co56_fe56",
            "other",
        ],
        "#937860",
    ),
]

# Minimum nuclear energy generation (erg g^-1 s^-1) for a zone to be shown as
# a burning region.  This is deliberately absolute so weak C/O shells still
# appear; raise it if a run looks noisy.
BURN_THRESHOLD = 1e-4

# A zone is part of a porosity-reduction region when the opacity factor is
# below this value.
REDUCTION_LIMIT = 1.0 - 1e-6

# Number of mass points in the interpolated Kippenhahn grid.
DEFAULT_N_Y = 600


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------


def _arr(prof, name):
    """Return a float array with NaN/inf replaced by 0, or None if missing."""
    if not hasattr(prof, name):
        return None
    a = np.asarray(getattr(prof, name), dtype=float)
    return np.nan_to_num(a, nan=0.0, posinf=0.0, neginf=0.0)


class KippenhahnCollector:
    """Collect the quantities needed for a Kippenhahn diagram.

    It is designed to be passed to ``stream.stream_profiles`` together with
    the other collectors, so the profiles are only read once per run.
    """

    def __init__(self, burn_threshold: float = BURN_THRESHOLD):
        self.burn_threshold = burn_threshold
        self.model_numbers: list[float] = []
        self.star_ages: list[float] = []
        self.mass_all: list[np.ndarray] = []
        self.mix_all: list[np.ndarray] = []
        self.burn_all: list[np.ndarray] = []
        self.inv_all: list[np.ndarray] = []
        # ``None`` means this run has no extra_opacity_factor column.
        self.red_all: list[np.ndarray | None] = []
        # Enclosed-mass extent of all burning zones (used for an optional zoom
        # when burning is confined to a very thin layer, e.g. a nova shell).
        self.burn_min: float | None = None
        self.burn_max: float | None = None

    # ------------------------------------------------------------------
    def consume(self, prof) -> None:
        if not hasattr(prof, "mass"):
            return

        mass = np.asarray(getattr(prof, "mass"), dtype=float)
        if mass.size < 2 or not np.all(np.isfinite(mass)):
            return
        if not np.all(np.diff(mass) >= 0):
            order = np.argsort(mass)
            mass = mass[order]
        else:
            order = None

        nz = mass.size

        # --- mixing ---------------------------------------------------
        if hasattr(prof, "mix_type"):
            mix = np.asarray(getattr(prof, "mix_type"), dtype=int)
            if order is not None:
                mix = mix[order]
        else:
            mix = np.zeros(nz, dtype=int)

        # --- nuclear burning -----------------------------------------
        eps = _arr(prof, "eps_nuc")
        if eps is None:
            eps = np.zeros(nz)
        elif order is not None:
            eps = eps[order]

        rates = []
        for _label, cols, _colour in BURN_GROUPS:
            rate = np.zeros(nz)
            for col in cols:
                a = _arr(prof, col)
                if a is not None:
                    if order is not None:
                        a = a[order]
                    rate += np.clip(a, 0.0, None)
            rates.append(rate)
        rates = np.stack(rates)
        dominant = np.argmax(rates, axis=0)
        burn = np.where(
            (eps > self.burn_threshold) & (np.max(rates, axis=0) > 0.0),
            dominant + 1,  # 1-based; 0 means "no burning"
            0,
        ).astype(np.int8)

        burning = burn > 0
        if burning.any():
            bm = mass[burning]
            self.burn_min = (
                float(bm[0]) if self.burn_min is None else min(self.burn_min, float(bm[0]))
            )
            self.burn_max = (
                float(bm[-1]) if self.burn_max is None else max(self.burn_max, float(bm[-1]))
            )

        # --- inversion ------------------------------------------------
        grad_t = _arr(prof, "gradT")
        if grad_t is None:
            inversion = np.zeros(nz, dtype=bool)
        else:
            if order is not None:
                grad_t = grad_t[order]
            inversion = grad_t < 0.0

        # --- porosity reduction --------------------------------------
        factor = _arr(prof, "extra_opacity_factor")
        if factor is None:
            reduction = None
        else:
            if order is not None:
                factor = factor[order]
            # ``_arr`` turns missing values into 0; a real opacity factor is
            # always positive, so this avoids flagging NaNs as reduction.
            reduction = (factor > 0.0) & (factor < REDUCTION_LIMIT)

        self.model_numbers.append(float(getattr(prof, "model_number")))
        self.star_ages.append(float(getattr(prof, "star_age")))
        self.mass_all.append(mass)
        self.mix_all.append(mix)
        self.burn_all.append(burn)
        self.inv_all.append(inversion)
        self.red_all.append(reduction)

    # ------------------------------------------------------------------
    def render(self, mass: int, name: str, phase: str = "", n_y: int = DEFAULT_N_Y) -> None:
        if not self.model_numbers:
            print(f"    KIPPENHAHN: no profiles collected for {mass}m-{name}")
            return

        phase_suffix = f"_{phase}" if phase else ""
        out_dir = f"{utils.PLOTS_ROOT}/plot{mass}m{phase_suffix}"
        os.makedirs(out_dir, exist_ok=True)

        # Sort by model number (the natural time ordering in the profiles).
        order = np.argsort(self.model_numbers)
        mn = np.asarray(self.model_numbers, dtype=float)[order]
        age = np.asarray(self.star_ages, dtype=float)[order]
        mass_all = [self.mass_all[i] for i in order]
        mix_all = [self.mix_all[i] for i in order]
        burn_all = [self.burn_all[i] for i in order]
        inv_all = [self.inv_all[i] for i in order]
        red_all = [self.red_all[i] for i in order]

        # Common enclosed-mass grid.  Regions above a profile's surface are
        # masked out, which naturally shows the shrinking surface if the run
        # loses mass.
        ymax = max(float(np.max(m)) for m in mass_all)
        if not np.isfinite(ymax) or ymax <= 0:
            print(f"    KIPPENHAHN: invalid maximum mass for {mass}m-{name}")
            return
        y = np.linspace(0.0, ymax, n_y)

        surface = np.empty(len(mn))

        def build_grids(y_grid):
            burn = np.full((len(y_grid), len(mn)), np.nan)
            mix = np.full((len(y_grid), len(mn)), np.nan)
            inv = np.zeros((len(y_grid), len(mn)), dtype=bool)
            red = np.zeros((len(y_grid), len(mn)), dtype=bool)
            for j, m in enumerate(mass_all):
                burn[:, j] = _map_to_grid(m, burn_all[j], y_grid)
                mix[:, j] = _map_to_grid(m, mix_all[j], y_grid)
                inv[:, j] = _map_to_grid(m, inv_all[j].astype(float), y_grid) > 0.5
                if red_all[j] is not None:
                    red[:, j] = _map_to_grid(m, red_all[j].astype(float), y_grid) > 0.5
                surface[j] = m[-1]
            return burn, mix, inv, red

        burn_grid, mix_grid, inv_grid, red_grid = build_grids(y)

        core_lines = _load_core_lines(mass, name, mn, ymax)

        self._render_pair(
            mass=mass,
            name=name,
            out_dir=out_dir,
            mn=mn,
            age=age,
            y=y,
            burn_grid=burn_grid,
            mix_grid=mix_grid,
            inv_grid=inv_grid,
            red_grid=red_grid,
            surface=surface,
            core_lines=core_lines,
            y_lim=(0.0, ymax),
            suffix="",
            title_suffix="",
        )

        # If all burning happens in a very thin layer (e.g. the shell of a
        # nova), the linear mass axis makes it disappear against the surface.
        # Add a zoomed pair of plots around that layer.
        if (
            self.burn_min is not None
            and ymax > 0.0
            and 0.0 < (self.burn_max - self.burn_min) < 0.02 * ymax
        ):
            span = self.burn_max - self.burn_min
            lo = max(0.0, self.burn_min - 0.15 * span)
            hi = min(ymax, self.burn_max + 0.15 * span)
            if 0.0 < (hi - lo) < 0.05 * ymax:
                y_zoom = np.linspace(lo, hi, n_y)
                zburn, zmix, zinv, zred = build_grids(y_zoom)
                self._render_pair(
                    mass=mass,
                    name=name,
                    out_dir=out_dir,
                    mn=mn,
                    age=age,
                    y=y_zoom,
                    burn_grid=zburn,
                    mix_grid=zmix,
                    inv_grid=zinv,
                    red_grid=zred,
                    surface=surface,
                    core_lines=core_lines,
                    y_lim=(lo, hi),
                    suffix="_zoom",
                    title_suffix=" (zoom)",
                )

    # ------------------------------------------------------------------
    def _render_pair(
        self,
        *,
        mass,
        name,
        out_dir,
        mn,
        age,
        y,
        burn_grid,
        mix_grid,
        inv_grid,
        red_grid,
        surface,
        core_lines,
        y_lim,
        suffix,
        title_suffix,
    ) -> None:
        # ---- time-to-live variant -----------------------------------
        age_order = np.argsort(age)
        age_sorted = age[age_order]
        mn_sorted = mn[age_order]
        t_end = float(age_sorted[-1])
        if len(mn) < 2 or t_end <= 0.0:
            return
        x_ttl = t_end - age_sorted
        # Log scale cannot show exactly zero; keep the final column visible.
        floor = max(t_end * 1e-6, 1e-30)
        x_ttl = np.maximum(x_ttl, floor)
        # Reverse to ascending order for contouring, then invert the axis so
        # time still runs left -> right.
        x_asc = x_ttl[::-1]
        reorder = age_order[::-1]
        core_ttl = (
            {k: (v[0], v[1], v[2][reorder]) for k, v in core_lines.items()}
            if core_lines
            else {}
        )

        self._render_one(
            x=x_asc,
            y=y,
            burn_grid=burn_grid[:, reorder],
            mix_grid=mix_grid[:, reorder],
            inv_grid=inv_grid[:, reorder],
            red_grid=red_grid[:, reorder],
            surface=surface[reorder],
            core_lines=core_ttl,
            x_label="Star Time To Live [years]",
            log_x=True,
            invert_x=True,
            y_lim=y_lim,
            out_path=(
                f"{out_dir}/KIPPENHAHN_{mass}_{name}_vs_mass_xstar_age{suffix}.png"
            ),
        )

    # ------------------------------------------------------------------
    def _render_one(
        self,
        *,
        x,
        y,
        burn_grid,
        mix_grid,
        inv_grid,
        red_grid,
        surface,
        core_lines,
        x_label,
        log_x,
        invert_x,
        y_lim,
        out_path,
    ) -> None:
        n_burn = len(BURN_GROUPS)
        burn_cmap = ListedColormap([g[2] for g in BURN_GROUPS])
        burn_cmap.set_bad((0.0, 0.0, 0.0, 0.0))
        burn_norm = BoundaryNorm(np.arange(0.5, n_burn + 1.5), n_burn)
        burn_plot = np.where(burn_grid >= 1.0, burn_grid, np.nan)

        fig, ax = plt.subplots(figsize=(12, 7))
        # contourf/contour need at least two time columns.
        can_contour = len(x) >= 2

        # 1. Burning regions: filled colours.
        ax.pcolormesh(
            x, y, burn_plot, cmap=burn_cmap, norm=burn_norm,
            shading="nearest", zorder=1,
        )

        # 2. Mixing regions: hatched, unfilled overlays.
        mix_handles = []
        present_mix = sorted(
            int(v) for v in np.unique(mix_grid[np.isfinite(mix_grid)]) if int(v) != 0
        )
        for code in present_mix if can_contour else []:
            label, colour, hatch = MIX_INFO.get(
                code, (f"Mix type {code}", "#333333", "..")
            )
            if code == -1:
                hatch = ""
                face = (0.6, 0.6, 0.6, 0.15)
            else:
                face = (1.0, 1.0, 1.0, 0.0)
            mask = (mix_grid == code).astype(float)
            cs = ax.contourf(
                x, y, mask, levels=[0.5, 1.5], hatches=[hatch],
                colors=[face], zorder=2,
            )
            cs.set_edgecolor(colour)
            cs.set_linewidth(0.0)
            mix_handles.append(
                Patch(facecolor=face, edgecolor=colour, hatch=hatch, label=label)
            )

        # 3. Temperature inversion: light blue fill + blue hatching.
        inv_handles = []
        if can_contour and inv_grid.any():
            cs = ax.contourf(
                x, y, inv_grid.astype(float), levels=[0.5, 1.5],
                hatches=["\\\\"], colors=[(0.12, 0.47, 0.71, 0.18)], zorder=3,
            )
            cs.set_edgecolor("#1F77B4")
            cs.set_linewidth(0.0)
            ax.contour(
                x, y, inv_grid.astype(float), levels=[0.5],
                colors="#1F77B4", linewidths=1.6, zorder=4,
            )
            inv_handles.append(
                Patch(
                    facecolor=(0.12, 0.47, 0.71, 0.18), edgecolor="#1F77B4",
                    hatch="\\\\", label=r"Inversion ($\nabla_T<0$)",
                )
            )

        # 4. Porosity reduction: subtle grey fill + dotted outline.
        red_handles = []
        if can_contour and red_grid.any():
            cs = ax.contourf(
                x, y, red_grid.astype(float), levels=[0.5, 1.5],
                colors=[(0.25, 0.25, 0.25, 0.10)], zorder=4,
            )
            cs.set_edgecolor("#555555")
            cs.set_linewidth(0.0)
            ax.contour(
                x, y, red_grid.astype(float), levels=[0.5],
                colors="#444444", linestyles=":", linewidths=1.2, zorder=5,
            )
            red_handles.append(
                Patch(
                    facecolor=(0.25, 0.25, 0.25, 0.10), edgecolor="#444444",
                    linestyle=":", linewidth=1.2,
                    label=r"Porosity reduction ($f_\kappa<1$)",
                )
            )

        # 5. Surface and core boundaries.
        ax.plot(x, surface, color="black", linewidth=0.8, zorder=4,
                label="Stellar surface")
        core_handles = []
        for label, colour, values in core_lines.values():
            ax.plot(x, values, color=colour, linewidth=1.0, linestyle="--",
                    zorder=6)
            core_handles.append(
                plt.Line2D([], [], color=colour, linewidth=1.0,
                           linestyle="--", label=label)
            )

        # 6. Axes, legends, output.
        ax.set_xlim(float(np.min(x)), float(np.max(x)))
        ax.set_ylim(float(y_lim[0]), float(y_lim[1]))
        if float(y_lim[1]) - float(y_lim[0]) < 0.01:
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.5f}"))
        if log_x:
            ax.set_xscale("log")
        if invert_x:
            ax.invert_xaxis()
        ax.set_xlabel(x_label)
        ax.set_ylabel(r"Enclosed mass $M_r$ [$M_\odot$]")
        ax.grid(True, which="major", alpha=0.18, linewidth=0.5)

        burn_handles = [
            Patch(facecolor=colour, label=f"{label} burning")
            for i, (label, _cols, colour) in enumerate(BURN_GROUPS)
            if np.any(burn_grid == i + 1)
        ]
        all_handles = (
            burn_handles + inv_handles + red_handles + core_handles
            + [plt.Line2D([], [], color="black", linewidth=0.8, label="Stellar surface")]
            + mix_handles
        )
        leg = ax.legend(
            handles=all_handles, loc="upper right", fontsize=7,
            framealpha=1.0, ncol=2,
        )
        leg.get_frame().set_facecolor("white")
        leg.get_frame().set_edgecolor("#aaaaaa")
        leg.set_zorder(10)

        fig.tight_layout()
        fig.savefig(out_path, dpi=300)
        plt.close(fig)
        print(f"    Saved: {out_path}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _map_to_grid(m, values, y):
    """Nearest-zone map of a profile quantity onto the common mass grid.

    ``m`` must be ascending.  Grid points above the stellar surface become
    NaN so they are transparent in the plot.
    """
    m = np.asarray(m, dtype=float)
    values = np.asarray(values)
    idx = np.searchsorted(m, y, side="left")
    idx = np.clip(idx, 0, len(m) - 1)
    left = np.clip(idx - 1, 0, len(m) - 1)
    choose = np.where(np.abs(m[idx] - y) <= np.abs(m[left] - y), idx, left)
    out = values[choose].astype(float, copy=True)
    out[y > m[-1]] = np.nan
    return out


def _even_ticks(values, n):
    values = np.asarray(values, dtype=float)
    if values.size <= n:
        return values
    return values[np.linspace(0, values.size - 1, n, dtype=int)]


def _load_core_lines(mass: int, name: str, record_mn, ymax: float):
    """Return core-boundary lines, interpolated onto the profile times."""
    path = f"{mass}m-{name}/LOGS/history.data"
    try:
        import mesa_reader as mr

        hist = mr.MesaData(path)
    except Exception as e:
        print(f"    KIPPENHAHN: could not read history for {mass}m-{name}: {e}")
        return {}

    hist_mn = np.asarray(hist.model_number, dtype=float)
    order = np.argsort(hist_mn)
    hist_mn = hist_mn[order]

    definitions = [
        ("He core", "#222222", "he_core_mass"),
        ("CO core", "#006400", "co_core_mass"),
        ("Fe core", "#8B0000", "fe_core_mass"),
        ("O/Ne core", "#8B4513", "one_core_mass"),
        ("Neutron-rich core", "#4B0082", "neutron_rich_core_mass"),
    ]
    lines = {}
    for label, colour, col in definitions:
        if not hasattr(hist, col):
            continue
        vals = np.asarray(getattr(hist, col), dtype=float)[order]
        vals = np.nan_to_num(vals, nan=0.0)
        if np.nanmax(vals) <= 0.0:
            continue
        # Skip a boundary that is effectively at the surface (e.g. a pure WD).
        if np.nanmax(vals) >= 0.97 * ymax or np.nanmax(vals) <= 0.02 * ymax:
            continue
        lines[col] = (label, colour, np.interp(record_mn, hist_mn, vals))
    return lines


# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------


def plot(mass: int, name: str, phase: str = "") -> None:
    """Standalone entry point: stream one run and render both Kippenhahn variants."""
    from stream import stream_profiles

    logs_path = f"{mass}m-{name}/LOGS"
    collector = KippenhahnCollector()
    count = stream_profiles(logs_path, [collector])
    if count == 0:
        print(f"    KIPPENHAHN: no profiles found in {logs_path}")
        return
    collector.render(mass, name, phase=phase)


def plot_all(mass: int, names: list[str]) -> None:
    for name in names:
        try:
            plot(mass, name)
        except Exception as e:
            print(f"    KIPPENHAHN({name}) failed: {e}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Make Kippenhahn diagrams.")
    parser.add_argument("mass", type=int, help="stellar mass, e.g. 10")
    parser.add_argument("name", help="run folder suffix, e.g. ours_reduction")
    parser.add_argument("--phase", default="", help="phase suffix used in the output folder")
    args = parser.parse_args()
    plot(args.mass, args.name, phase=args.phase)
