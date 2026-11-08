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
* Density inversions ("inflation") in the outer envelope, where
  ``logRho`` increases outwards before decreasing again.
* Porosity reduction regions (``extra_opacity_factor < 1``), when the run
  has this column.
* Zones where inflation and porosity reduction coincide keep both fills; their
  shared border is drawn as a dark solid line with an amber dashed overlay so
  both colours stay visible in thin overlap layers.

The heavy use of hatched, unfilled overlays means the burning colours stay
readable even where several regions overlap.
"""

import os

import matplotlib

matplotlib.use("Agg")  # headless backend: avoids macOS AppKit/XPC listeners

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

# Porosity-reduction shading is split into steps of this size (2.5 %).  Each
# step towards a smaller opacity factor is drawn darker.
REDUCTION_STEP = 0.025

# Density-inversion ("inflation") detection.  A zone is marked when logRho
# increases outwards over a contiguous run by at least this many dex and the
# run lies in the outer part of the star.  The small numeric threshold rejects
# profile-to-profile noise while keeping genuine envelope inversions.
INFLATION_MIN_DELTA = 0.05
INFLATION_OUTER_Q = 0.5

# The red star-boundary line is lifted slightly above the true surface and
# drawn in the background (lowest zorder).  Inflation and porosity-reduction
# outlines also run along the outer envelope, so this keeps the thick red
# dashes from hiding them.  The lift is a fraction of the plotted mass range.
SURFACE_OFFSET_FRAC = 0.0025

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
        # Otherwise each entry is the reduction fraction ``1 - f_kappa`` per
        # zone, with 0 where there is no reduction.
        self.red_all: list[np.ndarray | None] = []
        # ``None`` means this run has no logRho column.  Otherwise each entry
        # is a boolean mask of zones inside a significant density inversion.
        self.infl_all: list[np.ndarray | None] = []

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

        # --- inversion ------------------------------------------------
        grad_t = _arr(prof, "gradT")
        if grad_t is None:
            inversion = np.zeros(nz, dtype=bool)
        else:
            if order is not None:
                grad_t = grad_t[order]
            inversion = grad_t < 0.0

        # --- density inversion ("inflation") --------------------------
        logrho = getattr(prof, "logRho", None)
        if logrho is None:
            inflation = None
        else:
            logrho = np.asarray(logrho, dtype=float)
            if order is not None:
                logrho = logrho[order]
            inflation = _density_inversion_mask(mass, logrho)

        # --- porosity reduction --------------------------------------
        factor = _arr(prof, "extra_opacity_factor")
        if factor is None:
            reduction = None
        else:
            if order is not None:
                factor = factor[order]
            # ``_arr`` turns missing values into 0; a real opacity factor is
            # always positive, so this avoids flagging NaNs as reduction.
            reduction = np.where(
                (factor > 0.0) & (factor < REDUCTION_LIMIT),
                1.0 - factor,
                0.0,
            ).astype(float)

        self.model_numbers.append(float(getattr(prof, "model_number")))
        self.star_ages.append(float(getattr(prof, "star_age")))
        self.mass_all.append(mass)
        self.mix_all.append(mix)
        self.burn_all.append(burn)
        self.inv_all.append(inversion)
        self.red_all.append(reduction)
        self.infl_all.append(inflation)

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
        infl_all = [self.infl_all[i] for i in order]

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
            infl = np.zeros((len(y_grid), len(mn)), dtype=bool)
            red = np.zeros((len(y_grid), len(mn)), dtype=float)
            # Thin regions (e.g. the near-surface layers of the 40 m and 100 m
            # runs) can be many orders of magnitude thinner than one grid cell.
            # Give them a minimum on-screen thickness so they remain visible
            # instead of falling between two grid points.
            dy = float(y_grid[1] - y_grid[0]) if len(y_grid) > 1 else 0.0
            min_thickness = 1.5 * dy
            for j, m in enumerate(mass_all):
                burn[:, j] = _map_to_grid(m, burn_all[j], y_grid)
                mix[:, j] = _map_to_grid(m, mix_all[j], y_grid)
                inv[:, j] = _map_to_grid(m, inv_all[j].astype(float), y_grid) > 0.5
                if infl_all[j] is not None:
                    infl[:, j] = _region_to_grid(
                        m, infl_all[j], y_grid,
                        min_thickness=min_thickness,
                        surface=m[-1],
                    )
                if red_all[j] is not None:
                    red[:, j] = _reduction_to_grid(
                        m, red_all[j], y_grid,
                        min_thickness=min_thickness,
                        surface=m[-1],
                    )
                surface[j] = m[-1]
            return burn, mix, inv, infl, red

        burn_grid, mix_grid, inv_grid, infl_grid, red_grid = build_grids(y)

        core_lines = _load_core_lines(mass, name, mn, ymax)

        self._render_plot(
            mass=mass,
            name=name,
            out_dir=out_dir,
            mn=mn,
            age=age,
            y=y,
            burn_grid=burn_grid,
            mix_grid=mix_grid,
            inv_grid=inv_grid,
            infl_grid=infl_grid,
            red_grid=red_grid,
            surface=surface,
            core_lines=core_lines,
            y_lim=(0.0, ymax),
        )

    # ------------------------------------------------------------------
    def _render_plot(
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
        infl_grid,
        red_grid,
        surface,
        core_lines,
        y_lim,
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
            infl_grid=infl_grid[:, reorder],
            red_grid=red_grid[:, reorder],
            surface=surface[reorder],
            core_lines=core_ttl,
            x_label="Star Time To Live [years]",
            log_x=True,
            invert_x=True,
            y_lim=y_lim,
            out_path=f"{out_dir}/KIPPENHAHN_{mass}_{name}_vs_mass_xstar_age.png",
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
        infl_grid,
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
                    hatch="\\\\", label=r"Inversion ($\nabla T<0$)",
                )
            )

        # Inflation and porosity reduction can cover the same outer-envelope
        # zones.  Draw both fills (porosity underneath, translucent amber
        # inflation on top) so the overlap keeps a warm grey tint.  Their
        # outlines use the same small linewidth; in the overlap the dark
        # porosity line is overlaid with an amber dashed line instead of a
        # wider amber halo, so both colours are visible without making either
        # region look larger than it really is.
        overlap = infl_grid & (red_grid > 0.0)

        # 4. Inflation (density inversion): amber fill + continuous outline.
        # These are outer-envelope zones where logRho increases outwards before
        # decreasing again, so they appear as isolated bands in mass.
        infl_handles = []
        infl_colour = "#E69F00"
        infl_face = (0.90, 0.62, 0.00, 0.18)
        if can_contour and infl_grid.any():
            ax.contourf(
                x, y, infl_grid.astype(float), levels=[0.5, 1.5],
                colors=[infl_face], zorder=3.5,
            )
            ax.contour(
                x, y, infl_grid.astype(float), levels=[0.5],
                colors=infl_colour, linestyles="-", linewidths=1.5, zorder=4.5,
            )
            infl_handles.append(
                Patch(
                    facecolor=infl_face, edgecolor=infl_colour,
                    linestyle="-", linewidth=1.5,
                    label="Density Inversion",
                )
            )

        # 5. Porosity reduction: grey bands (darker = stronger reduction) with
        # a continuous outline.  The bands are kept strong enough that a thin,
        # one-cell-tall layer is still visible.  The fill sits just below the
        # inflation fill so both contribute in overlapping zones.
        red_handles = []
        red_mask = (red_grid > 0.0).astype(float)
        if can_contour and red_mask.any():
            red_plot = np.where(red_grid > 0.0, red_grid, np.nan)
            max_red = float(np.nanmax(red_plot))
            n_bins = max(1, int(np.ceil(max_red / REDUCTION_STEP)))
            levels = np.linspace(0.0, n_bins * REDUCTION_STEP, n_bins + 1)
            shades = np.linspace(0.60, 0.15, n_bins)
            alphas = np.linspace(0.35, 0.75, n_bins)
            red_colors = [
                (float(s), float(s), float(s), float(a))
                for s, a in zip(shades, alphas)
            ]
            ax.contourf(
                x, y, red_plot, levels=levels, colors=red_colors, zorder=3.4,
            )
            # Continuous edge around the whole reduction region.
            ax.contour(
                x, y, red_mask, levels=[0.5],
                colors="#2F2F2F", linestyles="-", linewidths=1.3, zorder=5,
            )
            red_handles.append(
                Patch(
                    facecolor=(0.42, 0.42, 0.42, 0.55),
                    edgecolor="#2F2F2F", linestyle="-", linewidth=1.3,
                    label=r"Porosity reduction",
                )
            )

        # 6. Overlap outline.  Draw the dark porosity border along the shared
        # boundary and an amber dashed line on top of it.  The equal on/off
        # dash pattern alternates the two colours along the border, so both
        # remain visible even for a layer that is only one grid row thick.
        # Both lines have the same width, so the overlap does not exaggerate
        # either region.
        if can_contour and overlap.any():
            both = overlap.astype(float)
            ax.contour(
                x, y, both, levels=[0.2],
                colors="#2F2F2F", linestyles="-", linewidths=1.4, zorder=5.2,
            )
            ax.contour(
                x, y, both, levels=[0.2],
                colors=infl_colour, linestyles=[(0, (2, 2))],
                linewidths=1.4, zorder=5.2,
            )

        # 7. Surface and core boundaries.  Use the same red dashed style as the
        # star boundaries in the heat maps.  Lift the surface line a little
        # above the true surface and draw it in the first (background) layer so
        # it cannot hide the inflation and porosity-reduction outlines.
        surface_offset = SURFACE_OFFSET_FRAC * (
            float(y_lim[1]) - float(y_lim[0])
        )
        ax.plot(
            x, surface + surface_offset, color="red", linewidth=1.5,
            linestyle="--", zorder=0, label="Star Boundaries",
        )
        core_handles = []
        for label, colour, values in core_lines.values():
            ax.plot(x, values, color=colour, linewidth=1.0, linestyle="--",
                    zorder=6)
            core_handles.append(
                plt.Line2D([], [], color=colour, linewidth=1.0,
                           linestyle="--", label=label)
            )

        # 8. Axes, legends, output.
        ax.set_xlim(float(np.min(x)), float(np.max(x)))
        # Leave headroom on top for the lifted star-boundary line.
        ax.set_ylim(float(y_lim[0]), float(y_lim[1]) + surface_offset)
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
            burn_handles + inv_handles + infl_handles + red_handles + core_handles
            + [
                plt.Line2D(
                    [], [], color="red", linewidth=1.5, linestyle="--",
                    label="Star Boundaries",
                )
            ]
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


def _density_inversion_mask(m, logrho, min_delta=INFLATION_MIN_DELTA,
                            outer_q=INFLATION_OUTER_Q):
    """Return zones that belong to a significant outer density inversion.

    A density inversion is a contiguous run where ``logRho`` increases
    outwards (``d logRho / dm > 0``) instead of the usual decrease.  Only runs
    in the outer envelope (``m / m_surface >= outer_q``) with a total rise of
    at least ``min_delta`` dex are kept, so profile noise does not fill the
    diagram.  The surrounding zones must decrease again on both sides, which
    selects the down-up-down pattern of an inflated layer.
    """
    m = np.asarray(m, dtype=float)
    logrho = np.asarray(logrho, dtype=float)
    out = np.zeros(m.size, dtype=bool)
    if m.size < 3 or logrho.size != m.size:
        return out

    valid = np.isfinite(logrho)
    q = m / m[-1] if m[-1] > 0.0 else np.zeros_like(m)
    rising = (np.diff(logrho) > 0.0) & valid[:-1] & valid[1:]
    rising &= q[:-1] >= outer_q
    idx = np.flatnonzero(rising)
    if idx.size == 0:
        return out

    runs = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
    n = logrho.size
    for run in runs:
        i0 = int(run[0])
        i1 = int(run[-1]) + 1
        if (logrho[i1] - logrho[i0]) < min_delta:
            continue
        # Require a decrease immediately inside and outside the rising run.
        if i0 > 0 and valid[i0 - 1] and logrho[i0] >= logrho[i0 - 1]:
            continue
        if i1 + 1 < n and valid[i1 + 1] and logrho[i1] <= logrho[i1 + 1]:
            continue
        out[i0:i1 + 1] = True
    return out


def _region_to_grid(m, mask, y, min_thickness=0.0, surface=None):
    """Map a boolean zone mask onto the mass grid without losing thin layers.

    Each contiguous run of True values is projected onto ``y``.  Runs thinner
    than ``min_thickness`` are widened symmetrically around their centre (and
    clipped to the stellar surface) so they still cover at least one grid row.
    """
    out = np.zeros(len(y), dtype=bool)
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return out

    idx = np.flatnonzero(mask)
    runs = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
    top = float(m[-1]) if surface is None else float(surface)
    for run in runs:
        if run.size == 0:
            continue
        lo = float(m[run[0]])
        hi = float(m[run[-1]])
        if hi < lo:
            lo, hi = hi, lo
        lo = max(0.0, min(lo, top))
        hi = max(0.0, min(hi, top))
        if min_thickness > 0.0 and (hi - lo) < min_thickness:
            centre = 0.5 * (lo + hi)
            lo = centre - 0.5 * min_thickness
            hi = centre + 0.5 * min_thickness
            # Keep the full visual thickness inside the star.  A layer sitting
            # on the surface is pushed downwards rather than clipped in half.
            if hi > top:
                hi = top
                lo = max(0.0, hi - min_thickness)
            if lo < 0.0:
                lo = 0.0
                hi = min(top, min_thickness)

        first = int(np.searchsorted(y, lo, side="left"))
        last = int(np.searchsorted(y, hi, side="right"))
        if last <= first:
            # The whole layer is thinner than one grid cell: light up the
            # cell closest to its centre.
            first = int(np.argmin(np.abs(y - 0.5 * (lo + hi))))
            last = first + 1
        out[first:last] = True

    # Never draw a region above the stellar surface.
    out[y > top + 1e-12] = False
    return out


def _reduction_to_grid(m, reduction, y, min_thickness=0.0, surface=None):
    """Map per-zone porosity reduction onto the mass grid.

    ``reduction`` is ``1 - extra_opacity_factor`` indexed like the ascending
    ``m`` array, with 0 where there is no reduction.  Thin layers are widened
    to ``min_thickness`` and assigned the strongest reduction in the layer, so
    they remain visible instead of falling between two grid points.
    """
    reduction = np.asarray(reduction, dtype=float)
    mapped = _map_to_grid(m, reduction, y)
    out = np.where(np.isnan(mapped), 0.0, np.clip(mapped, 0.0, None))
    mask = reduction > 0.0
    if not mask.any():
        return out

    idx = np.flatnonzero(mask)
    runs = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
    top = float(m[-1]) if surface is None else float(surface)
    for run in runs:
        if run.size == 0:
            continue
        lo = float(m[run[0]])
        hi = float(m[run[-1]])
        if hi < lo:
            lo, hi = hi, lo
        lo = max(0.0, min(lo, top))
        hi = max(0.0, min(hi, top))
        peak = float(np.max(reduction[run]))
        if min_thickness > 0.0 and (hi - lo) < min_thickness:
            centre = 0.5 * (lo + hi)
            lo = centre - 0.5 * min_thickness
            hi = centre + 0.5 * min_thickness
            # Keep the full visual thickness inside the star.  A layer sitting
            # on the surface is pushed downwards rather than clipped in half.
            if hi > top:
                hi = top
                lo = max(0.0, hi - min_thickness)
            if lo < 0.0:
                lo = 0.0
                hi = min(top, min_thickness)

            first = int(np.searchsorted(y, lo, side="left"))
            last = int(np.searchsorted(y, hi, side="right"))
            if last <= first:
                # The whole layer is thinner than one grid cell: light up the
                # cell closest to its centre.
                first = int(np.argmin(np.abs(y - 0.5 * (lo + hi))))
                last = first + 1
            out[first:last] = np.maximum(out[first:last], peak)

    # Never draw reduction above the stellar surface.
    out[y > top + 1e-12] = 0.0
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
