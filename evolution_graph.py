import glob
import os

import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np

import utils

options = """
model_number
num_zones
star_age
star_age_sec
star_age_min
star_age_hr
star_age_day
time_step
time_step_sec
time_step_days
log_dt
log_dt_sec
log_dt_days
star_mass
log_xmstar
log_abs_mdot
mass_conv_core
conv_mx1_top
conv_mx1_bot
conv_mx2_top
conv_mx2_bot
mx1_top
mx1_bot
mx2_top
mx2_bot
log_LH
log_LHe
log_LZ
log_Lnuc
pp
cno
tri_alpha
epsnuc_M_1
epsnuc_M_2
epsnuc_M_3
epsnuc_M_4
epsnuc_M_5
epsnuc_M_6
epsnuc_M_7
epsnuc_M_8
he_core_mass
co_core_mass
one_core_mass
fe_core_mass
neutron_rich_core_mass
kh_timescale
effective_T
log_Teff
luminosity
log_L
log_R
log_g
log_L_div_Ledd
lum_div_Ledd
v_surf_km_s
v_div_csound_surf
surf_avg_j_rot
surf_avg_omega
surf_avg_omega_crit
surf_avg_omega_div_omega_crit
surf_avg_v_rot
surf_avg_v_crit
surf_avg_v_div_v_crit
surf_avg_Lrad_div_Ledd
log_rotational_mdot_boost
log_center_T
log_center_Rho
log_cntr_P
log_cntr_Rho
log_cntr_T
center_T
center_Rho
center_mu
center_ye
center_abar
center_entropy
center_omega
center_omega_div_omega_crit
center_h1
center_he4
center_c12
center_o16
surface_c12
surface_o16
total_mass_h1
total_mass_he4
log_max_T
total_energy
log_total_energy
total_energy_foe
rel_E_err
log_rel_E_err
log_rel_run_E_err
log_total_angular_momentum
gamma1_min
num_retries
num_iters
num_solver_iterations
burn_type_1
burn_qtop_1
burn_type_2
burn_qtop_2
burn_type_3
burn_qtop_3
burn_type_4
burn_qtop_4
burn_type_5
burn_qtop_5
burn_type_6
burn_qtop_6
burn_type_7
burn_qtop_7
burn_type_8
burn_qtop_8
burn_type_9
burn_qtop_9
burn_type_10
burn_qtop_10
burn_type_11
burn_qtop_11
burn_type_12
burn_qtop_12
burn_type_13
burn_qtop_13
burn_type_14
burn_qtop_14
burn_type_15
burn_qtop_15
burn_type_16
burn_qtop_16
burn_type_17
burn_qtop_17
burn_type_18
burn_qtop_18
burn_type_19
burn_qtop_19
burn_type_20
burn_qtop_20
burn_type_21
burn_qtop_21
burn_type_22
burn_qtop_22
burn_type_23
burn_qtop_23
burn_type_24
burn_qtop_24
burn_type_25
burn_qtop_25
burn_type_26
burn_qtop_26
burn_type_27
burn_qtop_27
burn_type_28
burn_qtop_28
burn_type_29
burn_qtop_29
burn_type_30
burn_qtop_30
burn_type_31
burn_qtop_31
burn_type_32
burn_qtop_32
burn_type_33
burn_qtop_33
burn_type_34
burn_qtop_34
burn_type_35
burn_qtop_35
burn_type_36
burn_qtop_36
burn_type_37
burn_qtop_37
burn_type_38
burn_qtop_38
burn_type_39
burn_qtop_39
burn_type_40
burn_qtop_40
mix_type_1
mix_qtop_1
mix_type_2
mix_qtop_2
mix_type_3
mix_qtop_3
mix_type_4
mix_qtop_4
mix_type_5
mix_qtop_5
mix_type_6
mix_qtop_6
mix_type_7
mix_qtop_7
mix_type_8
mix_qtop_8
mix_type_9
mix_qtop_9
mix_type_10
mix_qtop_10
mix_type_11
mix_qtop_11
mix_type_12
mix_qtop_12
mix_type_13
mix_qtop_13
mix_type_14
mix_qtop_14
mix_type_15
mix_qtop_15
mix_type_16
mix_qtop_16
mix_type_17
mix_qtop_17
mix_type_18
mix_qtop_18
mix_type_19
mix_qtop_19
mix_type_20
mix_qtop_20
mix_type_21
mix_qtop_21
mix_type_22
mix_qtop_22
mix_type_23
mix_qtop_23
mix_type_24
mix_qtop_24
mix_type_25
mix_qtop_25
mix_type_26
mix_qtop_26
mix_type_27
mix_qtop_27
mix_type_28
mix_qtop_28
mix_type_29
mix_qtop_29
mix_type_30
mix_qtop_30
mix_type_31
mix_qtop_31
mix_type_32
mix_qtop_32
mix_type_33
mix_qtop_33
mix_type_34
mix_qtop_34
mix_type_35
mix_qtop_35
mix_type_36
mix_qtop_36
mix_type_37
mix_qtop_37
mix_type_38
mix_qtop_38
mix_type_39
mix_qtop_39
mix_type_40
mix_qtop_40
"""


def plot(mass: int, names: list[str], y_axis: str, phase: str = ""):
    phase_suffix = f"_{phase}" if phase else ""

    plt.figure()

    base_colors: dict[str, str] = {}
    any_plotted = False

    for name in names:
        path = f"{mass}m-{name}/LOGS"

        history = None
        try:
            history = mr.MesaData(f"{path}/history.data")
        except Exception:
            pass

        if history is not None and hasattr(history, y_axis) and hasattr(history, "star_age"):
            x = getattr(history, "star_age")
            y = getattr(history, y_axis)

            kwargs, is_winds, base_name = utils.get_line_kwargs(name, base_colors)
            (line,) = plt.plot(x, y, **kwargs)
            color = line.get_color()
            if not is_winds:
                base_colors[base_name] = color
            plt.plot(x[-1], y[-1], "o", markersize=6, color=color)
            plt.plot(x[0], y[0], "o", markersize=6, color=color)
            any_plotted = True
        else:
            prof_files = sorted(
                glob.glob(f"{path}/profile*.data"),
                key=lambda x: int(x.split("profile")[-1].split(".data")[0]),
            )
            if not prof_files:
                print(f"    evolution_graph({name}): no profiles found")
                continue
            raw = []
            missing = False
            for f_path in prof_files:
                try:
                    prof = mr.MesaData(f_path)
                except Exception:
                    continue
                if not hasattr(prof, y_axis):
                    if not missing:
                        print(f"    evolution_graph({name}): '{y_axis}' not in profiles")
                        missing = True
                    continue
                age = getattr(prof, "star_age")
                col = getattr(prof, y_axis)
                if y_axis == "mass":
                    val = col[0]
                else:
                    val = sum(col)
                raw.append((age, val))

            if not raw:
                print(f"    evolution_graph({name}): no data for '{y_axis}'")
                continue

            raw.sort(key=lambda x: x[0])
            ages, vals = zip(*raw)

            kwargs, is_winds, base_name = utils.get_line_kwargs(name, base_colors)
            (line,) = plt.plot(ages, vals, **kwargs)
            color = line.get_color()
            if not is_winds:
                base_colors[base_name] = color
            plt.plot(ages[-1], vals[-1], "o", markersize=6, color=color)
            plt.plot(ages[0], vals[0], "o", markersize=6, color=color)
            any_plotted = True

    if not any_plotted:
        plt.close()
        print(f"    evolution_graph: nothing to plot for '{y_axis}'")
        return

    plt.xlabel("Star Age [yr]")
    plt.ylabel(utils.pretty_axis_name(y_axis))
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    out_dir = f"{utils.PLOTS_ROOT}/plot{mass}m{phase_suffix}"
    os.makedirs(out_dir, exist_ok=True)
    plt.savefig(f"{out_dir}/Mass_{mass}_{y_axis}_vs_time.png", dpi=300)
    plt.close()
    print(f"    Saved: {out_dir}/Mass_{mass}_{y_axis}_vs_time.png")


def plot_percentages(mass: int, names: list[str], phase: str = ""):
    """Plot core hydrogen and stellar mass percentages vs star age.

    Only the ``ours`` run is used. Core hydrogen (center_h1) is shown as
    a percentage of the mass fraction; mass is shown as a percentage of
    its initial value. The initial mass is written on the plot.
    """
    phase_suffix = f"_{phase}" if phase else ""
    plt.figure()
    any_plotted = False
    mass_initial = None

    for name in names:
        if "ours" not in name:
            continue
        try:
            history = mr.MesaData(f"{mass}m-{name}/LOGS/history.data")
        except Exception as e:
            print(f"    percentages({name}): could not read history: {e}")
            continue
        if not (
            hasattr(history, "star_age")
            and hasattr(history, "star_mass")
            and hasattr(history, "center_h1")
        ):
            print(f"    percentages({name}): missing history columns")
            continue

        age = np.array(history.star_age, dtype=float)
        mass_star = np.array(history.star_mass, dtype=float)
        h_core = np.array(history.center_h1, dtype=float)

        order = np.argsort(age)
        age = age[order]
        mass_star = mass_star[order]
        h_core = h_core[order]

        mass_initial = mass_star[0]
        if not np.isfinite(mass_initial) or mass_initial == 0:
            print(f"    percentages({name}): invalid initial mass")
            continue

        plt.plot(
            age, 100.0 * mass_star / mass_initial, color="C0", linestyle="-",
            linewidth=1.6, label="Mass (% of $M_{\\rm init}$)",
        )
        plt.plot(
            age, 100.0 * h_core, color="C1", linestyle="--",
            linewidth=1.6, label="H in core (%)",
        )
        any_plotted = True
        break

    if not any_plotted:
        plt.close()
        print("    percentages: nothing to plot")
        return

    plt.xlabel("Star Age [yr]")
    plt.ylabel("Percentage")
    plt.ylim(0, 105)
    plt.grid(True)
    plt.legend(fontsize=9)
    plt.text(
        0.02, 0.95, f"$M_{{\\rm init}}={mass_initial:.4g}\\,M_\\odot$",
        transform=plt.gca().transAxes, va="top", fontsize=10,
        bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "none"},
    )
    plt.tight_layout()

    out_dir = f"{utils.PLOTS_ROOT}/plot{mass}m{phase_suffix}"
    os.makedirs(out_dir, exist_ok=True)
    out_path = f"{out_dir}/Mass_{mass}_H_and_mass_percent_vs_time.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"    Saved: {out_path}")


def plot_taus(mass: int, names: list[str], phase: str = ""):
    """Plot tauM, tauH, tauHt vs star age for the 'ours' model only.

    tauM = (dM/dt)/M, M = star_mass
    tauH = (dH/dt)/H, H = center_h1 (hydrogen in core)
    tauHt = (dHt/dt)/Ht, Ht = total hydrogen % in star = 100*total_mass_h1/star_mass
    """
    phase_suffix = f"_{phase}" if phase else ""
    plt.figure()
    any_plotted = False
    mass_initial = None

    for name in names:
        if "ours" not in name:
            continue
        try:
            history = mr.MesaData(f"{mass}m-{name}/LOGS/history.data")
        except Exception as e:
            print(f"    taus({name}): could not read history: {e}")
            continue
        if not (
            hasattr(history, "star_age")
            and hasattr(history, "star_mass")
            and hasattr(history, "center_h1")
            and hasattr(history, "total_mass_h1")
        ):
            print(f"    taus({name}): missing history columns")
            continue

        age = np.array(history.star_age, dtype=float)
        M = np.array(history.star_mass, dtype=float)
        Hc = np.array(history.center_h1, dtype=float)
        Ht_total = np.array(history.total_mass_h1, dtype=float)

        order = np.argsort(age)
        age = age[order]
        M = M[order]
        Hc = Hc[order]
        Ht_total = Ht_total[order]

        Ht_pct = 100.0 * Ht_total / M
        if not any_plotted:
            mass_initial = M[0]

        def _fr(x, t):
            x = np.asarray(x, dtype=float)
            t = np.asarray(t, dtype=float)
            o = np.argsort(t)
            x = x[o]
            t = t[o]
            dxdt = np.gradient(x, t)
            res = np.full_like(dxdt, np.nan, dtype=float)
            np.divide(dxdt, x, out=res, where=x != 0)
            return t[o], res

        t_m, tau_m = _fr(M, age)
        t_h, tau_h = _fr(Hc, age)
        t_ht, tau_ht = _fr(Ht_pct, age)

        t_end = np.max(age)
        x_m = np.maximum(t_end - t_m, 0.0)
        x_h = np.maximum(t_end - t_h, 0.0)
        x_ht = np.maximum(t_end - t_ht, 0.0)

        plt.plot(x_m, tau_m, color="C0", linestyle="-", linewidth=1.6, label=r"$\tau_M = \dot M/M$")
        plt.plot(x_h, tau_h, color="C1", linestyle="--", linewidth=1.6, label=r"$\tau_H = \dot H/H$ ($H=\text{center }H_1$)")
        plt.plot(x_ht, tau_ht, color="C2", linestyle="-.", linewidth=1.6, label=r"$\tau_{Ht} = \dot H_t/H_t$ ($H_t=\text{total }H/M\times100\%$)")
        any_plotted = True
        break

    if not any_plotted:
        plt.close()
        print("    taus: nothing to plot")
        return

    ax = plt.gca()
    plt.xlabel("Star Time To Live [years]")
    plt.ylabel(r"$\tau = \frac{1}{X}\frac{dX}{dt}$ [1/yr]")
    ax.invert_xaxis()
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1e-15)
    plt.grid(True, which="both")
    plt.legend(fontsize=8)
    if mass_initial is not None and np.isfinite(mass_initial):
        plt.text(
            0.02, 0.95, f"$M_{{\\rm init}}={mass_initial:.4g}\\,M_\\odot$",
            transform=ax.transAxes, va="top", fontsize=10,
            bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "none"},
        )
    plt.tight_layout()

    out_dir = f"{utils.PLOTS_ROOT}/plot{mass}m{phase_suffix}"
    os.makedirs(out_dir, exist_ok=True)
    out_path = f"{out_dir}/Mass_{mass}_tau_M_H_HT_vs_time.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"    Saved: {out_path}")
