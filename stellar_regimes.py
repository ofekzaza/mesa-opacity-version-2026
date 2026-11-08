"""Physics behind the log T - log rho regime diagram.

Everything here is plain numpy and cgs units. Conventions follow Prialnik
(ch. 3 and 7) and Kippenhahn: boundaries between EOS regimes are drawn where
the corresponding pressures are equal, so all lines meet at a common triple
point.

Contents
--------
* EOS boundaries (P_gas = P_rad, ideal/degenerate, NR/relativistic degenerate)
* Gamma_1 < 4/3 zones: pair creation (computed from an e+- EOS), H/He
  ionization (Saha), Fe photodisintegration and electron capture (schematic)
* Nuclear ignition curves  eps(rho, T) = eps_thr
"""

import os

import numpy as np

# --- constants (cgs) ---------------------------------------------------------
A_RAD = 7.5657e-15      # radiation constant
R_GAS = 8.31447e7       # gas constant
K_B = 1.380649e-16
H_PLANCK = 6.62607e-27
M_E = 9.10938e-28
M_U = 1.66054e-24
C_LIGHT = 2.99792e10
EV = 1.602177e-12
E2 = 2.30708e-19        # e^2 in erg cm
ALPHA_FS = 1 / 137.036
HBAR = H_PLANCK / (2 * np.pi)
K1_NR = 1.0036e13       # P_e = K1 (rho/mu_e)^(5/3)  (NR degenerate)
K2_ER = 1.2435e15       # P_e = K2 (rho/mu_e)^(4/3)  (ER degenerate)

# Composition used for the EOS boundary lines (ionized, H-rich envelope-like).
# Pure C/O has mu ~ 1.8 and mu_e = 2; the degenerate/relativistic lines then
# shift up by log10(mu_e / MU_E) ~ 0.23 dex.
MU = 0.6
MU_E = 1.17


# --- EOS boundaries ----------------------------------------------------------
def log_rho_prad_eq_pgas(log_t, mu=MU):
    """P_gas = P_rad:  rho R T / mu = a T^4 / 3."""
    return 3.0 * log_t + np.log10(A_RAD * mu / (3.0 * R_GAS))


def log_rho_degeneracy(log_t, mu_e=MU_E):
    """Electron ideal pressure = NR-degenerate pressure."""
    return 1.5 * log_t + 1.5 * np.log10(R_GAS / K1_NR) + np.log10(mu_e)


def log_rho_ideal_er(log_t, mu_e=MU_E):
    """Electron ideal pressure = ER-degenerate pressure."""
    return 3.0 * log_t + 3.0 * np.log10(R_GAS / K2_ER) + np.log10(mu_e)


def log_rho_nr_er(mu_e=MU_E):
    """NR-degenerate pressure = ER-degenerate pressure (horizontal line)."""
    return 3.0 * np.log10(K2_ER / K1_NR) + np.log10(mu_e)


def log_t_triple_point():
    """Where the NR-deg, ER-deg and ideal-gas regimes meet (mu_e independent)."""
    return np.log10(K2_ER ** 2 / (K1_NR * R_GAS))


LOG_T_THERMAL_RELATIVISTIC = np.log10(M_E * C_LIGHT ** 2 / K_B)   # kT = m_e c^2


# --- partial-ionization zones (Saha) ----------------------------------------
def _log_t_ionization(chi_ev, log_rho, n_per_rho, x):
    """logT at which a species with number density n_per_rho*rho has ionized
    fraction x (single-step Saha, statistical weights ~1)."""
    n = 10.0 ** np.asarray(log_rho) * n_per_rho
    target = x ** 2 / (1.0 - x ** 2)
    log_t = np.linspace(3.0, 6.5, 7001)
    t = 10.0 ** log_t
    out = []
    for ni in np.atleast_1d(n):
        s = (2 * np.pi * M_E * K_B * t / H_PLANCK ** 2) ** 1.5 \
            * np.exp(-chi_ev * EV / (K_B * t)) / ni
        out.append(log_t[np.argmax(s > target)])
    return np.array(out)


def ionization_bands(log_rho, x_lo=0.1, x_hi=0.9, h_mass_frac=0.70,
                     he_mass_frac=0.28):
    """(logT_low, logT_high) of the H, He I and He II partial-ionization zones
    as functions of log rho. Gamma_1 < 4/3 roughly where 10% < x < 90%."""
    m_h = 1.6735e-24
    species = {
        "H": (13.598, h_mass_frac / m_h),
        "He I": (24.587, he_mass_frac / (4 * m_h)),
        "He II": (54.418, he_mass_frac / (4 * m_h)),
    }
    bands = {}
    for name, (chi, n_per_rho) in species.items():
        bands[name] = (
            _log_t_ionization(chi, log_rho, n_per_rho, x_lo),
            _log_t_ionization(chi, log_rho, n_per_rho, x_hi),
        )
    return bands


# --- pair creation: Gamma_1 from a full e+- EOS -----------------------------
def _fd_moments(eta, theta, n_pts=240):
    """Dimensionless moments of the e- / e+ distribution.

    eta   : chemical potential of electrons incl. rest mass, in m_e c^2
    theta : kT / m_e c^2
    Returns net number, pressure and energy (all in units K = 1/(pi^2 lambda_C^3),
    pressure and energy additionally times m_e c^2).
    """
    x_top = np.maximum(eta, 0.0) + 45.0 * theta + 2.0
    u = np.linspace(0.0, 1.0, n_pts)
    x = x_top[:, None] * u[None, :]
    eps = np.sqrt(1.0 + x * x)
    arg_m = np.clip((eps - eta[:, None]) / theta[:, None], -200, 200)
    arg_p = np.clip((eps + eta[:, None]) / theta[:, None], -200, 200)
    f_m = 1.0 / (np.exp(arg_m) + 1.0)
    f_p = 1.0 / (np.exp(arg_p) + 1.0)

    def integ(y):
        return np.trapezoid(y, dx=1.0, axis=1) * x_top / (n_pts - 1)

    n_net = integ(x * x * (f_m - f_p))
    pres = integ(x ** 4 / eps * (f_m + f_p)) / 3.0
    ener = integ(x * x * eps * (f_m + f_p))
    return n_net, pres, ener


def _pair_eos(rho, temp, ye, a_ion):
    """Total pressure and energy density of ions + radiation + e+- gas."""
    lam = HBAR / (M_E * C_LIGHT)
    k_norm = 1.0 / (np.pi ** 2 * lam ** 3)
    theta = K_B * temp / (M_E * C_LIGHT ** 2)
    n_target = rho * ye / M_U / k_norm                      # in units of K
    x_f = (3.0 * n_target) ** (1.0 / 3.0)                   # p_F / m c
    lo = np.zeros_like(rho)
    hi = np.sqrt(1.0 + x_f ** 2) + 1e-9
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        n_net, _, _ = _fd_moments(mid, theta)
        too_big = n_net > n_target
        hi = np.where(too_big, mid, hi)
        lo = np.where(too_big, lo, mid)
    eta = 0.5 * (lo + hi)
    _, p_e, u_e = _fd_moments(eta, theta)
    mc2 = M_E * C_LIGHT ** 2
    n_ion = rho / (a_ion * M_U)
    pres = p_e * k_norm * mc2 + n_ion * K_B * temp + A_RAD * temp ** 4 / 3.0
    ener = u_e * k_norm * mc2 + 1.5 * n_ion * K_B * temp + A_RAD * temp ** 4
    return pres, ener


def gamma1_grid(log_t, log_rho, ye=0.5, a_ion=16.0):
    """Gamma_1 on a (log_t x log_rho) grid for ions + radiation + e+- gas.

    Gamma_1 = chi_rho + P chi_T^2 / (rho T c_v)  with finite-difference
    derivatives of the full EOS. Ye = 0.5, A = 16 mimics a C/O/Ne core.
    """
    tt, rr = np.meshgrid(log_t, log_rho, indexing="ij")
    shape = tt.shape
    t0 = 10.0 ** tt.ravel()
    r0 = 10.0 ** rr.ravel()
    d = 2e-3                                                 # dex
    f = 10.0 ** d

    p00, e00 = _pair_eos(r0, t0, ye, a_ion)
    p_tp, e_tp = _pair_eos(r0, t0 * f, ye, a_ion)
    p_tm, e_tm = _pair_eos(r0, t0 / f, ye, a_ion)
    p_rp, _ = _pair_eos(r0 * f, t0, ye, a_ion)
    p_rm, _ = _pair_eos(r0 / f, t0, ye, a_ion)

    chi_t = (np.log(p_tp) - np.log(p_tm)) / (2 * d * np.log(10))
    chi_r = (np.log(p_rp) - np.log(p_rm)) / (2 * d * np.log(10))
    # c_v = d(U/rho)/dT at fixed rho
    c_v = (e_tp - e_tm) / r0 / (t0 * f - t0 / f)
    g1 = chi_r + p00 * chi_t ** 2 / (r0 * t0 * c_v)
    return g1.reshape(shape)


_PAIR_CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "tmp", "pair_gamma1_v2.npz")
_PAIR_MEMO = {}


def pair_gamma1_map():
    """(log_t, log_rho, Gamma_1) for the pair-creation region, cached on disk."""
    if "map" in _PAIR_MEMO:
        return _PAIR_MEMO["map"]
    try:
        with np.load(_PAIR_CACHE) as z:
            res = (z["log_t"], z["log_rho"], z["g1"])
    except Exception:
        log_t = np.linspace(8.6, 10.3, 69)
        log_rho = np.linspace(-3.0, 7.0, 81)
        res = (log_t, log_rho, gamma1_grid(log_t, log_rho))
        try:
            os.makedirs(os.path.dirname(_PAIR_CACHE), exist_ok=True)
            np.savez(_PAIR_CACHE, log_t=res[0], log_rho=res[1], g1=res[2])
        except OSError:
            pass
    # Below logrho ~ -3 the pair gas dominates the ions completely and Gamma_1
    # no longer depends on rho (checked on the grid), so extend the first row.
    log_t, log_rho, g1 = res
    ext = np.array([-14.0, -8.0])
    log_rho = np.concatenate([ext, log_rho])
    g1 = np.concatenate([np.repeat(g1[:, :1], len(ext), axis=1), g1], axis=1)
    res = (log_t, log_rho, g1)
    _PAIR_MEMO["map"] = res
    return res


# --- Fe photodisintegration & electron capture (schematic) -------------------
def fe_photodisintegration_edge(log_rho):
    """Onset of Fe -> alpha + n (Saha: weak positive dependence on rho)."""
    return 9.0 + 0.11 * np.asarray(log_rho)


FE_PHOTO_LOG_RHO = (5.5, 10.0)
FE_PHOTO_LOG_T_MAX = 10.4        # beyond this Gamma_1 recovers (full dissociation)
ELECTRON_CAPTURE_LOG_RHO = (9.5, 11.0)   # Fermi energy above capture thresholds;
                                         # essentially independent of T


# --- nuclear ignition curves -------------------------------------------------
EPS_THRESHOLD = 10.0     # erg/g/s, ~ the rate needed to stall Kelvin-Helmholtz contraction
LOG_T_MAX_SCAN = 10.4


def _eps_pp(t, rho, x=0.70):
    t6 = t / 1e6
    return 2.4e4 * rho * x ** 2 * t6 ** (-2 / 3) * np.exp(-33.8 * t6 ** (-1 / 3))


def _eps_cno(t, rho, x=0.70, z_cno=0.01):
    t6 = t / 1e6
    return 4.4e25 * rho * x * z_cno * t6 ** (-2 / 3) * np.exp(-152.28 * t6 ** (-1 / 3))


def _eps_3alpha(t, rho, y=1.0):
    t8 = t / 1e8
    return 5.1e8 * rho ** 2 * y ** 3 * t8 ** (-3) * np.exp(-44.027 / t8)


def _gamow_tau(t, z1, z2, a_red):
    """tau = 3 (E_G / 4kT)^(1/3),  E_G = 2 mu c^2 (pi alpha Z1 Z2)^2."""
    mu_c2 = a_red * M_U * C_LIGHT ** 2
    e_g = 2.0 * mu_c2 * (np.pi * ALPHA_FS * z1 * z2) ** 2
    return 3.0 * (e_g / (4.0 * K_B * t)) ** (1.0 / 3.0)


def _screening_exponent(t, rho, z1, z2, mu_e, zeta):
    """Salpeter weak/strong electron screening, H = min(weak, strong)."""
    t6 = t / 1e6
    h_weak = 0.188 * z1 * z2 * np.sqrt(zeta * rho) * t6 ** (-1.5)
    n_e = rho / (mu_e * M_U)
    a_e = (3.0 / (4.0 * np.pi * n_e)) ** (1.0 / 3.0)
    gamma_e = E2 / (a_e * K_B * t)
    h_strong = 0.9 * gamma_e * ((z1 + z2) ** (5 / 3) - z1 ** (5 / 3) - z2 ** (5 / 3))
    return np.minimum(h_weak, h_strong)


def _make_gamow_rate(z1, z2, a_red, mu_e, zeta, anchor):
    """eps = A rho T^(-2/3) exp(-tau + H_screen), A fixed so the curve passes
    through the canonical ignition point `anchor` = (logT, logrho)."""
    def shape(t, rho):
        h = _screening_exponent(t, rho, z1, z2, mu_e, zeta)
        return rho * t ** (-2 / 3) * np.exp(-_gamow_tau(t, z1, z2, a_red) + h)

    a0 = EPS_THRESHOLD / shape(10.0 ** anchor[0], 10.0 ** anchor[1])
    return lambda t, rho: a0 * shape(t, rho)


def _make_power_rate(nu, lam, anchor):
    """eps = A rho^lam T^nu through the anchor point (photodisintegration-driven
    stages, where the effective nu ~ Q/kT)."""
    t0, r0 = 10.0 ** anchor[0], 10.0 ** anchor[1]
    return lambda t, rho: EPS_THRESHOLD * (rho / r0) ** lam * (t / t0) ** nu


# Stage definitions. `rate(T, rho)` returns eps in erg/g/s.
BURNING_STAGES = {
    "H": dict(rate=lambda t, r: _eps_pp(t, r) + _eps_cno(t, r),
              log_rho=(-3.0, 6.0), log_t_min=6.3, color="#ffb3b3",
              note="pp + CNO, X=0.7, Z_CNO=0.01"),
    "He": dict(rate=_eps_3alpha, log_rho=(0.0, 8.0), log_t_min=7.0, color="#ff6666",
               note="triple-alpha, Y=1"),
    "C": dict(rate=_make_gamow_rate(6, 6, 6.0, 2.0, 3.5, anchor=(8.8, 5.0)),
              log_rho=(2.0, 10.0), log_t_min=7.8, color="#cc0000",
              note="C+C, Gamow shape + screening, through (8.8, 5.0)"),
    "Ne": dict(rate=_make_power_rate(nu=46.0, lam=1.0, anchor=(9.15, 6.5)),
               log_rho=(1.0, 10.0), log_t_min=8.3, color="#b34700",
               note="photodisintegration, nu~Q/kT, through (9.15, 6.5)"),
    "O": dict(rate=_make_gamow_rate(8, 8, 8.0, 2.0, 4.5, anchor=(9.3, 7.0)),
              log_rho=(4.0, 10.0), log_t_min=8.5, color="#800000",
              note="O+O, Gamow shape + screening, through (9.3, 7.0)"),
    "Si": dict(rate=_make_power_rate(nu=50.0, lam=1.0, anchor=(9.55, 7.7)),
               log_rho=(5.0, 10.0), log_t_min=8.8, color="#4d0000",
               note="quasi-equilibrium Si burning, schematic, through (9.55, 7.7)"),
}


def ignition_curve(stage, n=200):
    """(logT, log rho) of the coolest temperature at which eps(rho, T) reaches
    EPS_THRESHOLD. Points where the crossing falls outside the scanned range
    [log_t_min, 10.4] (e.g. the screening/pycnonuclear regime at the dense,
    cool end) are dropped rather than extrapolated."""
    spec = BURNING_STAGES[stage]
    lo, hi = spec["log_rho"]
    log_rho = np.linspace(lo, hi, n)
    scan = np.linspace(spec["log_t_min"], LOG_T_MAX_SCAN, 2001)
    t = 10.0 ** scan
    log_t = np.full(n, np.nan)
    with np.errstate(over="ignore", invalid="ignore"):
        for i, lr in enumerate(log_rho):
            above = spec["rate"](t, 10.0 ** lr) >= EPS_THRESHOLD
            j = np.argmax(above)
            if above.any() and j > 0:
                log_t[i] = scan[j]
    ok = np.isfinite(log_t)
    return log_t[ok], log_rho[ok]
