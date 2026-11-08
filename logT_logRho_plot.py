"""Standalone log T vs log rho diagram for the last profile of each model.

Background physics is computed in `stellar_regimes.py`. Conventions:

* EOS boundary lines are drawn at pressure equality (P_gas = P_rad,
  ideal gas <-> NR-degenerate electrons, NR <-> ER-degenerate electrons),
  so the three lines meet at one triple point.
* Envelope-like composition: mu = 0.6, mu_e = 1.17. A C/O core (mu_e ~ 2)
  shifts the degenerate lines up by log10(2/1.17) ~ 0.23 dex.
* The degeneracy line uses ideal electron pressure = NR-degenerate pressure
  (zero point -7.62 + log10(mu_e)), not the psi = 0 criterion.
* Gamma_1 < 4/3 zones: pair creation (full e+- EOS contour), H/He partial
  ionization (Saha), Fe photodisintegration and electron capture (schematic).
* Ignition curves are the eps(rho, T) = EPS_THRESHOLD = 10 erg/g/s loci
  (no neutrino losses; electron screening included).
"""

import os

import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
import mesa_reader as mr
import numpy as np

import stellar_regimes as sr
import utils
from generic_last_log_plot import _convective_mask, _plot_convective_segments


BOUNDARY_COLOR = "#7f8c8d"
ZONE_COLOR = "#e74c3c"          # Gamma_1 < 4/3
MARGINAL_COLOR = "#f1c40f"      # Gamma_1 -> 4/3 (radiation / ER-degenerate)
BURN_STAGES = ("H", "He", "C", "Ne", "O", "Si")

_LOG_T_GRID = np.linspace(3.0, 10.5, 800)


# --- helpers -----------------------------------------------------------------
def _in_view(ax, x, y, pad=0.0):
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    return (xlim[0] + pad <= x <= xlim[1] - pad
            and ylim[0] + pad <= y <= ylim[1] - pad)


def _visible_indices(ax, x, y):
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    x, y = np.asarray(x), np.asarray(y)
    return np.where((x >= xlim[0]) & (x <= xlim[1])
                    & (y >= ylim[0]) & (y <= ylim[1]))[0]


def _clamped_display(ax, x, y, margin=10.0):
    """(px, py, bbox): the data point (x, y) in figure pixels, clamped so it
    lies inside the axes frame with a small `margin` gap."""
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    xn = float(np.clip(x, xlim[0], xlim[1]))
    yn = float(np.clip(y, ylim[0], ylim[1]))
    px, py = ax.transData.transform((xn, yn))
    bb = ax.get_window_extent()
    px = min(max(px, bb.x0 + margin), bb.x1 - margin)
    py = min(max(py, bb.y0 + margin), bb.y1 - margin)
    return px, py, bb


def _pick_align(bb, px, py):
    """Alignment that keeps a label inside the frame: text extends away from
    the nearer edge (ha="right" extends left, va="top" extends downward)."""
    cx = 0.5 * (bb.x0 + bb.x1)
    cy = 0.5 * (bb.y0 + bb.y1)
    ha = "left" if px < cx else ("right" if px > cx else "center")
    va = "top" if py > cy else ("bottom" if py < cy else "center")
    return ha, va


def _fit_inside(ax, ann, px, py, max_iter=40):
    """Move an annotation's anchor (figure pixels) toward the axes center until
    the full rendered text box lies inside the axes frame."""
    bb = ax.get_window_extent()
    try:
        ren = ax.figure.canvas.get_renderer()
    except Exception:
        return
    for _ in range(max_iter):
        tb = ann.get_window_extent(ren)
        dx = dy = 0.0
        if tb.x0 < bb.x0:
            dx = bb.x0 - tb.x0
        elif tb.x1 > bb.x1:
            dx = bb.x1 - tb.x1
        if tb.y0 < bb.y0:
            dy = bb.y0 - tb.y0
        elif tb.y1 > bb.y1:
            dy = bb.y1 - tb.y1
        if abs(dx) < 0.5 and abs(dy) < 0.5:
            return
        px += dx
        py += dy
        ann.xy = (px, py)


def _text_inside(ax, x, y, text, color, fontsize=8, **kw):
    """Text pinned at (x, y), forced fully inside the axes frame."""
    px, py, bb = _clamped_display(ax, x, y)
    ha, va = _pick_align(bb, px, py)
    kw.setdefault("ha", ha)
    kw.setdefault("va", va)
    ann = ax.annotate(text, (px, py), xycoords="figure pixels",
                      xytext=(0, 0), textcoords="offset points",
                      fontsize=fontsize, color=color, alpha=0.9,
                      zorder=6, **kw)
    _fit_inside(ax, ann, px, py)
    return ann


def _label_along(ax, x, y, text, frac, color, offset_pts=(0, 4),
                 at_rho=None, angle_deg=None, **kw):
    """Put `text` on a curve at `frac` (0..1) of its visible length, rotated
    to follow the curve on screen, kept fully inside the axes frame.

    If `at_rho` is given, the label is pinned at the visible curve point with
    that log-rho value instead of at `frac`. `angle_deg` overrides the
    automatic tangent angle."""
    idx = _visible_indices(ax, x, y)
    if len(idx) < 3:
        return
    if at_rho is not None:
        k = idx[int(np.argmin(np.abs(np.asarray(y)[idx] - at_rho)))]
    else:
        k = idx[int(np.clip(frac, 0, 1) * (len(idx) - 1))]
    k = min(max(k, 1), len(x) - 2)
    if angle_deg is None:
        p0 = ax.transData.transform((x[k - 1], y[k - 1]))
        p1 = ax.transData.transform((x[k + 1], y[k + 1]))
        angle_deg = np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0]))
        if angle_deg > 90:              # keep the text upright
            angle_deg -= 180
        elif angle_deg < -90:
            angle_deg += 180
    angle = angle_deg
    px, py, bb = _clamped_display(ax, x[k], y[k])
    ha, va = _pick_align(bb, px, py)
    kw.setdefault("ha", ha)
    kw.setdefault("va", va)
    kw.setdefault("zorder", 6)
    ann = ax.annotate(text, (px, py), xycoords="figure pixels",
                      xytext=offset_pts, textcoords="offset points",
                      rotation=angle, rotation_mode="anchor",
                      fontsize=8, color=color, alpha=0.9,
                      **kw)
    _fit_inside(ax, ann, px, py)
    return ann


def _er_boundary(log_t):
    """Lower edge of the ER-degenerate regime (horizontal, then slope 3)."""
    t_tp = sr.log_t_triple_point()
    return np.where(log_t <= t_tp, sr.log_rho_nr_er(),
                    sr.log_rho_ideal_er(log_t))


# --- background layers -------------------------------------------------------
def draw_eos_regimes(ax):
    """EOS boundaries (pressure-equality convention) + marginal-Gamma_1 shading."""
    t = _LOG_T_GRID
    t_tp = sr.log_t_triple_point()
    style = dict(color=BOUNDARY_COLOR, linestyle="--", linewidth=1.2)

    rad = sr.log_rho_prad_eq_pgas(t)
    ax.plot(t, rad, **style)
    left = t <= t_tp
    ax.plot(t[left], sr.log_rho_degeneracy(t[left]), **style)
    ax.plot(t[left], np.full(left.sum(), sr.log_rho_nr_er()), **style)
    right = t >= t_tp
    ax.plot(t[right], sr.log_rho_ideal_er(t[right]), **style)
    ax.axvline(sr.LOG_T_THERMAL_RELATIVISTIC, color=BOUNDARY_COLOR,
               linestyle=":", linewidth=1.0, alpha=0.7)

    # Regimes where Gamma_1 -> 4/3 asymptotically (marginal stability).
    ax.fill_between(t, -40, rad, color=MARGINAL_COLOR, alpha=0.10,
                    linewidth=0, zorder=0.5)
    ax.fill_between(t, _er_boundary(t), 40, color=MARGINAL_COLOR, alpha=0.10,
                    linewidth=0, zorder=0.5)


def draw_gamma1_zones(ax):
    """Zones with Gamma_1 < 4/3."""
    # Pair creation: contour of Gamma_1 from a full e+- EOS (Ye=0.5, A=16).
    log_t, log_rho, g1 = sr.pair_gamma1_map()
    level = 4.0 / 3.0 - 0.003
    ax.contourf(log_t, log_rho, g1.T, levels=[0.0, level],
                colors=[ZONE_COLOR], alpha=0.30, zorder=1)
    ax.contour(log_t, log_rho, g1.T, levels=[level], colors=[ZONE_COLOR],
               linewidths=1.0, alpha=0.8, zorder=1)

    # H / He partial ionization: bands that follow the Saha curves.
    log_rho_i = np.linspace(-12.0, -2.0, 80)
    for lo, hi in sr.ionization_bands(log_rho_i).values():
        ax.fill_betweenx(log_rho_i, lo, hi, color=ZONE_COLOR, alpha=0.22,
                         linewidth=0, zorder=1)

    # Fe photodisintegration (schematic, parametrized edge).
    r0, r1 = sr.FE_PHOTO_LOG_RHO
    rho = np.linspace(r0, r1, 60)
    edge = sr.fe_photodisintegration_edge(rho)
    verts = np.vstack([
        np.column_stack([edge, rho]),
        [sr.FE_PHOTO_LOG_T_MAX, r1],
        [sr.FE_PHOTO_LOG_T_MAX, r0],
    ])
    ax.add_patch(Polygon(verts, closed=True, facecolor=ZONE_COLOR,
                         edgecolor=ZONE_COLOR, alpha=0.30, linewidth=1.0,
                         zorder=1))

    # Electron capture: set by the Fermi energy, so a band across all T.
    lo, hi = sr.ELECTRON_CAPTURE_LOG_RHO
    ax.axhspan(lo, hi, facecolor=ZONE_COLOR, edgecolor=ZONE_COLOR,
               alpha=0.30, linewidth=1.0, zorder=1)


def draw_ignition_curves(ax):
    """eps(rho, T) = EPS_THRESHOLD for each burning stage."""
    curves = {}
    for stage in BURN_STAGES:
        t, r = sr.ignition_curve(stage)
        if len(t) == 0:
            continue
        color = sr.BURNING_STAGES[stage]["color"]
        ax.plot(t, r, color=color, linestyle="--", linewidth=1.1, alpha=0.85,
                zorder=1)
        curves[stage] = (t, r)
    return curves


def _add_labels(ax, curves):
    """Text labels; call after the axis limits are final."""
    t = _LOG_T_GRID
    t_tp = sr.log_t_triple_point()

    for stage, (ct, cr) in curves.items():
        angle_deg = None
        if stage == "H":
            idxv = _visible_indices(ax, ct, cr)
            if len(idxv) > 2:
                k = idxv[int(0.9 * (len(idxv) - 1))]
                if abs(float(cr[k]) + 3.0) < 0.05:   # never pin at rho=-3 end
                    k = idxv[max(1, int(0.9 * (len(idxv) - 1)) - 1)]
                at_rho, frac = float(cr[k]), 1.0
            else:
                at_rho, frac = None, 0.9
        elif stage == "He":
            at_rho, frac = 0.0, 1.0
        elif stage == "Ne":
            at_rho, frac = 2.25, 1.0
        elif stage in ("C", "O"):
            at_rho, frac = None, 0.0
        else:
            at_rho, frac = None, 1.0
        _label_along(ax, ct, cr, f"{stage} burn", frac=frac,
                     color=sr.BURNING_STAGES[stage]["color"],
                     at_rho=at_rho, angle_deg=angle_deg,
                     offset_pts=(-3, -3), zorder=1.5)

    _label_along(ax, t, sr.log_rho_prad_eq_pgas(t), r"$P_{\rm gas}=P_{\rm rad}$",
                 0.15, BOUNDARY_COLOR, offset_pts=(0, -18))
    left = t <= t_tp
    _label_along(ax, t[left], sr.log_rho_degeneracy(t[left]), "degeneracy",
                 0.08, BOUNDARY_COLOR, offset_pts=(0, 8))

    # Regime names, placed in each region (only when inside the view).
    regimes = [
        (7.2, sr.log_rho_nr_er() + 0.7, "ER degenerate"),
    ]
    for x, y, name in regimes:
        if _in_view(ax, x, y, pad=0.3):
            _text_inside(ax, x, y, name, BOUNDARY_COLOR,
                         fontsize=8, style="italic")

    if _in_view(ax, 9.3, -4.0):
        _text_inside(ax, 9.3, -4.0, r"$\Gamma_1 = 4/3$", ZONE_COLOR)
    if _in_view(ax, 10.1, 8.0):
        _text_inside(ax, 10.1, 8.0, "Fe photodis.", ZONE_COLOR)
    ec_mid = 0.5 * sum(sr.ELECTRON_CAPTURE_LOG_RHO)
    if ax.get_ylim()[0] <= ec_mid <= ax.get_ylim()[1]:
        ax.text(0.02, ec_mid, r"$e^-$ capture", fontsize=8, color=ZONE_COLOR,
                transform=ax.get_yaxis_transform(), va="center", zorder=6)
    if _in_view(ax, 3.8, -3.5):
        _text_inside(ax, 3.8, -3.5, "H/He\nionization", ZONE_COLOR)

    # Marginal-stability zones (Gamma_1 -> 4/3): labelled inside the zones.
    if _in_view(ax, 7.5, -5.0):
        _text_inside(ax, 7.5, -5.0,
                     r"$\Gamma_1\to4/3$", MARGINAL_COLOR)
    if _in_view(ax, 9.9, 5.0):
        _text_inside(ax, 9.9, 5.0,
                     r"$\Gamma_1\to4/3$", MARGINAL_COLOR)


# --- main entry points -------------------------------------------------------
def load_last_profiles(mass: int, names: list[str]) -> "dict[str, mr.MesaData | None]":
    loaded: dict[str, mr.MesaData | None] = {}
    for name in names:
        path = f"{mass}m-{name}/LOGS"
        try:
            loaded[name] = mr.MesaData(utils.get_last_profile_path(path))
        except Exception as e:
            print(f"    logT_logRho_plot: could not load last profile for {name}: {e}")
            loaded[name] = None
    return loaded


def plot(mass: int, names: list[str], phase: str = "",
         _preloaded: "dict | None" = None):
    plt.figure()
    phase_suffix = f"_{phase}" if phase else ""
    base_colors: dict[str, str] = {}

    all_winds = all("_reduction" in n for n in names)
    all_post = all(n.endswith("_post") for n in names)

    profiles = _preloaded if _preloaded is not None else load_last_profiles(mass, names)

    data_x_min = data_x_max = None
    data_y_min = data_y_max = None

    for name in names:
        prof = profiles.get(name) if profiles is not None else None
        if prof is None:
            continue
        if not hasattr(prof, "logT") or not hasattr(prof, "logRho"):
            print(f"logT/logRho not in prof for {name}")
            continue

        x_data = getattr(prof, "logT")
        y_data = getattr(prof, "logRho")

        xa = np.asarray(x_data, dtype=float)
        ya = np.asarray(y_data, dtype=float)
        data_x_min = xa.min() if data_x_min is None else min(data_x_min, float(xa.min()))
        data_x_max = xa.max() if data_x_max is None else max(data_x_max, float(xa.max()))
        data_y_min = ya.min() if data_y_min is None else min(data_y_min, float(ya.min()))
        data_y_max = ya.max() if data_y_max is None else max(data_y_max, float(ya.max()))

        kwargs, is_winds, base_name = utils.get_line_kwargs(
            name, base_colors, drop_winds=all_winds, drop_post=all_post
        )
        (line,) = plt.plot(x_data, y_data, **kwargs)
        color = line.get_color()
        if not is_winds:
            base_colors[base_name] = color
        plt.plot(x_data[-1], y_data[-1], "o", markersize=6, color=color)
        plt.plot(x_data[0], y_data[0], "o", markersize=6, color=color)

        _plot_convective_segments(x_data, y_data, _convective_mask(prof), color)

    ax = plt.gca()
    draw_eos_regimes(ax)
    draw_gamma1_zones(ax)
    curves = draw_ignition_curves(ax)

    if data_x_min is not None:
        pad_x = 0.05 * (data_x_max - data_x_min) or 0.5
        ax.set_xlim(data_x_min - pad_x, data_x_max + pad_x)
    if data_y_min is not None:
        pad_y = 0.05 * (data_y_max - data_y_min) or 0.5
        ax.set_ylim(data_y_min - pad_y, data_y_max + pad_y)

    plt.xlabel(r"$\log[T(K)]$")
    plt.ylabel(r"$\log[\rho(\text{g}/\text{cm}^3)]$")

    plt.legend(fontsize=8)
    plt.tight_layout()

    # Lock the render geometry the saved image will use (300 dpi) before
    # placing labels, so the inside-frame clamping matches the final raster.
    ax.figure.set_dpi(300)
    ax.figure.canvas.draw()
    _add_labels(ax, curves)

    os.makedirs(f"{utils.PLOTS_ROOT}/plot{mass}m{phase_suffix}", exist_ok=True)
    plt.savefig(
        f"{utils.PLOTS_ROOT}/plot{mass}m{phase_suffix}/Mass_{mass}_logT_vs_logRho.png",
        dpi=300,
    )
    plt.close()
