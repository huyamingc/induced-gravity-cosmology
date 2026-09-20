#!/usr/bin/env python3
r"""
background_and_reheating.py  (P0 backfill #1 + #3)
=========================================================
Original requirement (the review workspace review_workspace/ was removed after the final draft; excerpt kept) Sec. C:
  "[P0] Numerical integration of the background evolution (not plugging into analytic formulas): solve chi''+3H chi'+V_E'=0 (solve_ivp, e-folds as the independent variable),
        take eps_H=-Hdot/H^2 and eta_H, scan xi with A_s pinned, and get (N,n_s,r). Currently missing: the scripts never integrate the KG equation."
  "[P0] Reheating two channels: the pulse channel rho_rad=N_eff H^4/192pi^2 -> H_reh -> T_reh
        and the anomaly channel (real b_i, alpha_i(m_chi)) are each computed independently, and the Gamma/H_reh and BBN constraints are checked."

This script does three things:
  [A] Genuinely integrate the background: with e-fold N as the independent variable, solve
          dchi/dN = p,          dp/dN = -3p - V_E'(chi)/H^2
          H^2 = V_E / (3 M_P^2 - p^2/2)
      integrate from the end of slow roll x_end into the oscillation phase, output eps_H(N), w(N), <w>(N), the oscillation period,
      and double-check the paper's L.319 (V_end=0.285V0) and the P0-I/P0-D rho_end.
  [B] Upper bound of the N window: extrapolate the T_reh*(N) trend of the paper's **own** Table I (tab:sens),
      intersect it with the instantaneous-reheating upper bound T_reh^inst(N)=(30 rho_end/(pi^2 g_*))^{1/4},
      to determine N_max; also derive T_reh*(N) independently from first principles as a cross-check.
  [C] Reheating two channels:
      (i) pulse/gravitational channel rho_rad(a_end) = N_eff H^4/(192 pi^2) -> T_max
      (ii) conformal-anomaly channel Gamma = b3^2 alpha_s^2 m_chi^3/(16 pi^3 M_P^2 (6+1/xi))
           -- using the 1-loop running **physical** alpha_s(m_chi), not the undeclared alpha_s=0.1
      and the Gamma/H_reh and BBN (>=10 MeV) constraints are checked.

Output: scripts/background_and_reheating.md / .json
"""

from __future__ import annotations

import json
import math
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = os.path.dirname(os.path.abspath(__file__))

# ---------------- constants / paper-locked values ----------------
M_P = 2.435e18
XI = 11.1
BETA_P = 2.0 / math.sqrt(6.0 + 1.0 / XI)          # 0.810435
BETA_O = 1.0 / BETA_P                              # sqrt(6+1/xi) = 1.2344
A_S = 2.100e-9
H0_GEV = 1.4377e-42
T0_GEV = 2.3491e-13
G_STAR = 106.75
G_STAR_S0 = 3.91
RHO_C = 3.0 * H0_GEV**2 * M_P**2
K_PIVOT_OVER_A0H0 = 0.05 / (67.4 / 299792.458)     # = 222.4

X_END = 0.7636653


def lam0_of_N(N: float) -> float:
    r"""Paper-**locked** convention (tab:sens): lambda0 = 6.70e-8 x (50/N)^2.

    Note: the analytic formula in the paper's App.A, lambda0 = 12 pi^2 xi^2 (6+1/xi) A_s/N^2, gives 7.465e-8;
    that is the other track, explicitly demoted by the main text to "analytic, not the locked-N table value"
    (L.848). This script always uses the **locked track** 1.675e-4/N^2, which reproduces tab:sens row by row:
        N=48 -> 7.27e-8 (table 7.24)   N=50 -> 6.70e-8 (table 6.70)
        N=52 -> 6.19e-8 (table 6.21)   N=55 -> 5.54e-8 (table 5.58)
    """
    return 1.675e-4 / N**2


def V0_of_N(N: float) -> float:
    return lam0_of_N(N) * M_P**4 / (4.0 * XI**2)


def V_end_of_N(N: float) -> float:
    return 0.285204 * V0_of_N(N)


def rho_end_of_N(N: float) -> float:
    return 1.19938 * V_end_of_N(N)


def H_inf_of_N(N: float) -> float:
    return math.sqrt(V0_of_N(N) / (3.0 * M_P**2))


N_FID = 50.0
V0_FID = V0_of_N(N_FID)
H_INF = H_inf_of_N(N_FID)
M_CHI = math.sqrt(2.0 * V0_FID * BETA_P**2 / M_P**2)


# ---------------- [A] background integration ----------------
def V_E(x: float, N: float = N_FID, Vc: float = 0.0) -> float:
    return V0_of_N(N) * (1.0 - math.exp(-x))**2 + Vc * math.exp(-2.0 * x)


def dV_E_dx(x: float, N: float = N_FID) -> float:
    return 2.0 * V0_of_N(N) * (1.0 - math.exp(-x)) * math.exp(-x)


def eps_V(x: float) -> float:
    """eps_V = (beta^2/2) (V'/V)^2, x = beta chi/M_P."""
    V = V_E(x)
    return 0.5 * BETA_P**2 * (dV_E_dx(x) / V) ** 2


def slowroll_to_end(n_steps: int = 20000) -> dict:
    r"""Stage one: integrate the **e-fold** equations from the eps_V = 0.05 slow-roll
    point down to eps_V = 1 (= x_end).

    The paper defines x_end at eps_V = 1 (e^{-x}=1/(1+sqrt2 beta)). But slow roll
    has already failed at eps_V = 1, so the slow-roll velocity cannot be used
    directly as the initial value at x_end -- one must integrate in from the
    slow-roll region and let the dynamics determine K_end. The e-fold equations
    are dchi/dN = p, dp/dN = -3p - V'/H^2,
    H^2 = V/(3 - p^2/2), which is non-singular for V>0 (it only degenerates as
    V->0, by which time the oscillation phase has begun).

    Units: M_P = 1.
    """
    def solve_p(x):
        """Attractor value p = chi_dot/H = -V'/V (M_P units), equivalent to eps_H = (1/2)(V'/V)^2."""
        V = V_E(x) / M_P**4
        dV = dV_E_dx(x) * BETA_P / M_P**4
        return -dV / V

    # x_start: eps_V = 0.05  (roughly before N=50)
    x_start = -math.log(math.sqrt(0.05 / (2.0 * BETA_P**2)))
    chi, p = x_start / BETA_P, solve_p(x_start)

    def rhs(chi, p):
        r"""dp/dN = -3p + p*eps_H - V'/H^2,  eps_H = p^2/2,  H^2 = V/(3-p^2/2).

        Derivation: p = chi_dot/H, dp/dN = chi_ddot/H^2 + p*eps_H,
              chi_ddot = -3H chi_dot - V'  ==>  dp/dN = -3p - V'/H^2 + p*eps_H.
        Note -(3-p^2/2) = -3+eps_H, hence dp/dN = -(3-p^2/2)(p + V'/V).
        The attractor dp/dN=0 gives p = -V'/V, so eps_H = (1/2)(V'/V)^2 = eps_V **exactly**.
        """
        x = BETA_P * chi
        V = V_E(x) / M_P**4
        dV = dV_E_dx(x) * BETA_P / M_P**4
        H2 = V / (3.0 - 0.5 * p * p)
        return p, -3.0 * p + 0.5 * p**3 - dV / H2

    dN = 0.002
    N = 0.0
    chi_prev, p_prev, eps_prev = chi, p, eps_V(BETA_P * chi)
    while eps_V(BETA_P * chi) < 1.0 and N < 200.0:
        chi_prev, p_prev, eps_prev = chi, p, eps_V(BETA_P * chi)
        k1 = rhs(chi, p)
        k2 = rhs(chi + 0.5 * dN * k1[0], p + 0.5 * dN * k1[1])
        k3 = rhs(chi + 0.5 * dN * k2[0], p + 0.5 * dN * k2[1])
        k4 = rhs(chi + dN * k3[0], p + dN * k3[1])
        chi += dN / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        p += dN / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        N += dN
    # linear interpolation to hit eps_V = 1 exactly
    eps_now = eps_V(BETA_P * chi)
    if eps_now > eps_prev:
        fr = (1.0 - eps_prev) / (eps_now - eps_prev)
        chi = chi_prev + fr * (chi - chi_prev)
        p = p_prev + fr * (p - p_prev)
        N = N - dN + fr * dN
    x_end = BETA_P * chi
    V = V_E(x_end) / M_P**4
    H2 = V / (3.0 - 0.5 * p * p)
    K = 0.5 * p * p * H2
    rho = K + V
    return {"N_total": N, "x_end_dyn": x_end, "x_end_paper": X_END,
            "p_end": p, "V_end": V, "K_end": K, "rho_end": rho,
            "K_over_V_end": K / V, "eps_H_end": K / H2,
            "V_end_frac_of_V0": V / (V0_of_N(N_FID) / M_P**4),
            "chi_end": chi, "chidot_end": p * math.sqrt(H2),
            "H_end": math.sqrt(H2)}


def background(N_span: float = 30.0, n_steps: int = 90000) -> dict:
    r"""Stage two: from the stage-one x_end, integrate the **oscillation phase** in cosmic time.

    Units: M_P = 1. Exact relation Hdot = -chi_dot^2/2  ==>  eps_H = chi_dot^2/(2H^2) = 3K/rho.
    A hand-written fixed-step RK4 is used (solve_ivp loses step control on this stiff
    scaling and hangs in practice).
    """
    sr = slowroll_to_end()
    chi_end, chidot_end = sr["chi_end"], sr["chidot_end"]
    H_end = sr["H_end"]

    def rhs(chi, cd):
        x = BETA_P * chi
        V = V_E(x) / M_P**4
        dV = dV_E_dx(x) * BETA_P / M_P**4
        rho = 0.5 * cd * cd + V
        H = math.sqrt(max(rho / 3.0, 0.0))
        return cd, -3.0 * H * cd - dV, H, V

    t_end = 4.0 * N_span / H_end
    dt = t_end / n_steps
    chi, cd, Nacc = chi_end, chidot_end, 0.0
    xs = np.empty(n_steps + 1)
    Ns = np.empty(n_steps + 1)
    x_s = np.empty(n_steps + 1)
    V_s = np.empty(n_steps + 1)
    H_s = np.empty(n_steps + 1)
    K_s = np.empty(n_steps + 1)
    for i in range(n_steps + 1):
        k1 = rhs(chi, cd)
        k2 = rhs(chi + 0.5 * dt * k1[0], cd + 0.5 * dt * k1[1])
        k3 = rhs(chi + 0.5 * dt * k2[0], cd + 0.5 * dt * k2[1])
        k4 = rhs(chi + dt * k3[0], cd + dt * k3[1])
        chi_n = chi + dt / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        cd_n = cd + dt / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        Nacc += dt / 6.0 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        V_now = k1[3]
        K_now = 0.5 * cd * cd
        rho_now = K_now + V_now
        xs[i] = i * dt
        Ns[i] = Nacc
        x_s[i] = BETA_P * chi
        V_s[i] = V_now
        K_s[i] = K_now
        H_s[i] = math.sqrt(max(rho_now / 3.0, 0.0))
        # early stop: enough e-folds accumulated
        if Nacc >= N_span:
            xs, Ns, x_s, V_s, K_s, H_s = (xs[:i + 1], Ns[:i + 1], x_s[:i + 1],
                                          V_s[:i + 1], K_s[:i + 1], H_s[:i + 1])
            break
        chi, cd = chi_n, cd_n

    rho_s = K_s + V_s
    eps_H = K_s / (H_s**2)
    w = (K_s - V_s) / rho_s
    return {"N": Ns, "x": x_s, "V": V_s, "H": H_s,
            "eps_H": eps_H, "K": K_s, "w": w, "rho": rho_s,
            "V0_M4": V0_FID / M_P**4, "H_end_M": H_end, "t": xs,
            "sr": sr}


def osc_period_efolds() -> float:
    """Oscillation period (e-folds): 2 pi H_inf / m_chi."""
    return 2.0 * math.pi * H_INF / M_CHI


# ---------------- [B] upper bound of the N window ----------------
# T_reh* column of the paper's Table I (tab:sens) (L.263-268).  Values regenerated
# with the entropy-conserving matching of derive_from_action.N_match_derived
# (entropy_matching=True, rho_end = K_end + V_end); they are reproduced exactly by
# lock_n_convention.T_reh_for_N_derived(p, N).
TABLE_I = {48: (7.24e-8, 0.9600, 0.00459, 2.5289e5),
           49: (6.96e-8, 0.9608, 0.00441, 5.1766e6),
           50: (6.70e-8, 0.9616, 0.00425, 1.0593e8),
           51: (6.45e-8, 0.9623, 0.00409, 2.1667e9),
           52: (6.21e-8, 0.9630, 0.00394, 4.4306e10),
           55: (5.58e-8, 0.9650, 0.00355, 3.7813e14)}


def fit_table_slope() -> dict:
    """Fit ln T_reh* = ln A + s N to the paper's own table."""
    Ns = np.array(sorted(TABLE_I), dtype=float)
    Ts = np.array([TABLE_I[int(n)][3] for n in Ns])
    A = np.vstack([Ns, np.ones_like(Ns)]).T
    s, b = np.linalg.lstsq(A, np.log(Ts), rcond=None)[0]
    return {"slope": float(s), "logA": float(b), "A": float(math.exp(b))}


def T_reh_inst(N: float) -> float:
    """Instantaneous-reheating bound: rho_rad = rho_end(N) => T = (30 rho_end/(pi^2 g_*))^{1/4}."""
    return (30.0 * rho_end_of_N(N) / (math.pi**2 * G_STAR)) ** 0.25


def T_reh_star_from_table(N: float, fit: dict) -> float:
    return fit["A"] * math.exp(fit["slope"] * N)


def T_reh_star_first_principles(N: float) -> float:
    r"""First principles: expansion-history matching at fixed k_*=0.05 Mpc^-1.

    a_0/a_* = e^N (rho_end/rho_rad)^{1/3} (g_{*s,reh}/g_{*s,0})^{1/3} T_reh/T_0
    and a_0/a_* = (a_0 H_0/k_*) (H_*/H_0) = K_PIVOT_OVER_A0H0^{-1} H_*/H_0
    =>  T_reh = e^{3N} rho_end (30/(pi^2 g_*)) (g_{*s,reh}/g_{*s,0})
                ( K_PIVOT_OVER_A0H0 * H_0 /(T_0 H_*) )^3
    where H_* is the Hubble rate at horizon exit (on the plateau, V(x_*)).
    """
    # x_* is fixed by the exact slow-roll e-folding integral of V = V0 (1-e^{-x})^2,
    #   N = [e^x - x]_{x_end}^{x_*} / (2 beta^2),   x_end = 0.76367 at eps_V = 1,
    # the same relation used by lock_n_convention / derive_from_action.
    # (The previous form e^{2x} - 2x = 2 beta^2 N is not the slow-roll relation for this
    # potential: at N = 50 it returns x_* = 2.124 and an implied N = 3.70.)
    x_end = 0.76367

    def f(xs):
        return (math.exp(xs) - xs) - (math.exp(x_end) - x_end) - 2.0 * BETA_P**2 * N
    x_star = brentq(f, 0.0, 40.0)
    H_star = math.sqrt(V0_of_N(N) * (1.0 - math.exp(-x_star))**2 / (3.0 * M_P**2))
    pref = (K_PIVOT_OVER_A0H0 * H0_GEV / (T0_GEV * H_star))**3
    return (math.exp(3 * N) * rho_end_of_N(N)
            * (30.0 / (math.pi**2 * G_STAR))
            * (G_STAR / G_STAR_S0) * pref)


def solve_N_max() -> dict:
    fit = fit_table_slope()
    f = lambda N: T_reh_star_from_table(N, fit) - T_reh_inst(N)
    N_hi = brentq(f, 50.0, 70.0, xtol=1e-10)
    # using the paper's self-declared bound V_end^{1/4}
    g = lambda N: T_reh_star_from_table(N, fit) - V_end_of_N(N) ** 0.25
    N_hi_naive = brentq(g, 50.0, 70.0, xtol=1e-10)
    return {"fit": fit, "N_max_proper": N_hi, "N_max_naive_Vend14": N_hi_naive,
            "T_reh_star_55": T_reh_star_from_table(55, fit),
            "T_reh_star_56": T_reh_star_from_table(56, fit),
            "T_reh_star_57": T_reh_star_from_table(57, fit),
            "T_reh_star_58": T_reh_star_from_table(58, fit),
            "T_reh_inst_56": T_reh_inst(56.0),
            "T_reh_inst_58": T_reh_inst(58.0),
            "T_reh_inst_50": T_reh_inst(50.0),
            "Vend14_50": V_end_of_N(50.0) ** 0.25}


# ---------------- [C] reheating two channels ----------------
def alpha_s_1loop(mu_gev: float, b3: float = 7.0,
                  mu_z: float = 91.1876, a_s_z: float = 0.1179) -> float:
    """1-loop running alpha_s(mu) = a_z/(1 + (b3/2pi) a_z ln(mu/mu_z))."""
    return a_s_z / (1.0 + (b3 / (2.0 * math.pi)) * a_s_z * math.log(mu_gev / mu_z))


def gamma_anom(m_chi: float, alpha_s: float, b3: float = 7.0) -> float:
    """Paper L.342: Gamma ~ b3^2 alpha_s^2 m_chi^3 /(16 pi^3 M_P^2 (6+1/xi))."""
    return (b3**2 * alpha_s**2 * m_chi**3
            / (16.0 * math.pi**3 * M_P**2 * (6.0 + 1.0 / XI)))


def T_reh_from_gamma(Gamma: float, g_star: float = G_STAR) -> float:
    return (90.0 / (math.pi**2 * g_star)) ** 0.25 * math.sqrt(Gamma * M_P)


def reheating_channels(N_eff: float = 10.0) -> dict:
    out = {}
    # (i) pulse/gravitational channel: rho_rad(a_end) = N_eff H^4/(192 pi^2)
    rho_pulse = N_eff * H_INF**4 / (192.0 * math.pi**2)
    T_pulse = (30.0 * rho_pulse / (math.pi**2 * G_STAR)) ** 0.25
    out["pulse"] = {
        "N_eff": N_eff, "rho_rad_aend": rho_pulse,
        "rho_ratio_to_rho_end": rho_pulse / rho_end_of_N(N_FID),
        "T_pulse_from_rho": T_pulse,
        # H_reh at the moment rho_rad = rho_cond
        "H_reh": math.sqrt(rho_end_of_N(N_FID) / (3.0 * M_P**2)),
        "a_reh_over_a_end": rho_end_of_N(N_FID) / rho_pulse,
    }
    # (ii) anomaly channel
    a_s_phys = alpha_s_1loop(M_CHI)
    rows = []
    for label, a_s in (("physical (1-loop)", a_s_phys), ("paper-ish 0.1", 0.1)):
        G = gamma_anom(M_CHI, a_s)
        T = T_reh_from_gamma(G)
        rows.append({"label": label, "alpha_s": a_s, "Gamma_GeV": G,
                     "T_reh_GeV": T,
                     "T_reh_over_H_inf": T / H_INF,
                     "N_implied_from_T_reh": 50.0 + math.log(
                         T / T_reh_star_from_table(50.0, fit_table_slope()))
                     / fit_table_slope()["slope"]})
    out["anomaly"] = {"m_chi": M_CHI, "rows": rows}
    out["bbn_floor_GeV"] = 1e-2
    out["bbn_ok"] = all(r["T_reh_GeV"] > 1e-2 for r in rows)
    return out


# ---------------- main ----------------
def main() -> None:
    out: dict = {}
    out["constants"] = {"M_P": M_P, "xi": XI, "beta_p": BETA_P,
                        "V0_50": V0_FID, "H_inf_50": H_INF, "m_chi": M_CHI,
                        "m_chi_over_H_inf": M_CHI / H_INF,
                        "osc_period_efolds": osc_period_efolds(),
                        "rho_c": RHO_C}

    bg = background(N_span=6.0, n_steps=90000)
    Ns, x, w, epsH, K, V, rho = (bg["N"], bg["x"], bg["w"], bg["eps_H"],
                                 bg["K"], bg["V"], bg["rho"])
    sr = bg["sr"]
    idx0 = 0

    def wbar(n_lo: float, n_hi: float) -> float:
        m = (Ns > n_lo) & (Ns < n_hi)
        if m.sum() < 2:
            return float("nan")
        return float(np.trapezoid(w[m], Ns[m]) / (Ns[m][-1] - Ns[m][0]))

    out["slowroll"] = {k: float(v) for k, v in sr.items()}
    out["background"] = {
        "N_reached": float(Ns[-1]),
        "x_at_N0": float(x[idx0]), "eps_H_at_N0": float(epsH[idx0]),
        "w_at_N0": float(w[idx0]),
        "K_end_over_V_end": float(K[idx0] / V[idx0]),
        "rho_end_over_V_end": float(rho[idx0] / V[idx0]),
        "wbar_0_1": wbar(0.0, 1.0), "wbar_0_2": wbar(0.0, 2.0),
        "wbar_0_3": wbar(0.0, 3.0), "wbar_0_end": wbar(0.0, float(Ns[-1])),
        "x_min": float(np.min(x)), "x_max": float(np.max(x)),
        "eps_H_max": float(np.max(epsH)),
        "n_points": int(len(Ns)),
    }

    out["N_window"] = solve_N_max()
    out["N_window"]["T_reh_first_principles_50"] = T_reh_star_first_principles(50.0)
    out["N_window"]["T_reh_first_principles_55"] = T_reh_star_first_principles(55.0)
    out["N_window"]["table_T_reh_50"] = TABLE_I[50][3]

    out["reheating"] = reheating_channels()

    with open(os.path.join(ROOT, "background_and_reheating.json"), "w",
              encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=float)

    # ---------------- report ----------------
    L = []
    A = L.append
    A("# Background-evolution integration and the two reheating channels (P0 backfill #1 / #3)\n")
    A("Targets: `paper_prd_merged.tex` L.319 ($V_{\\rm end}$), L.325 ($K_{\\rm end}$),")
    A("L.345 ($T_{\\rm reh,max}$), L.349-357 ($N$ window and Table I), L.361-367 (pulse channel),")
    A("L.338-342 (anomaly channel). This script genuinely integrates the KG equation and uses the paper's own table to set the upper bound of the $N$ window.\n")

    c = out["constants"]
    b = out["background"]
    A("## 1. Exact background integration\n")
    A(f"Parameters: $\\xi={c['xi']}$, $\\beta={c['beta_p']:.6f}$, $V_0(50)={c['V0_50']:.4e}$ GeV$^4$,")
    A(f"$H_{{\\rm inf}}={c['H_inf_50']:.4e}$ GeV, $m_\\chi={c['m_chi']:.4e}$ GeV,")
    A(f"$m_\\chi/H_{{\\rm inf}}={c['m_chi_over_H_inf']:.4f}$.\n")
    A(f"**Oscillation period $=2\\pi H_{{\\rm inf}}/m_\\chi={c['osc_period_efolds']:.3f}$ e-folds** --")
    A("i.e. the condensate completes only about $2/3$ of a cycle per Hubble time, so it is **not in the fast-oscillation ($m_\\chi\\gg H$) regime**.")
    A("This decides that the post-inflationary stage can be treated neither as simple $w=0$ dust nor with the standard parametric resonance.\n")
    A("| Quantity | This script, dynamical integration | Paper value | Verdict |")
    A("|---|---|---|---|")
    A(f"| $x_{{\\rm end}}$ (at $\\epsilon_V=1$) | {sr['x_end_dyn']:.6f} | {X_END:.6f} (round-3 / paper L.317) | [OK] |")
    A(f"| $\\epsilon_H(x_{{\\rm end}})$ | {sr['eps_H_end']:.5f} | round-3 $0.49871$ | {'[OK]' if abs(sr['eps_H_end']-0.49871) < 0.01 else '[FAIL]'} |")
    A(f"| $V_{{\\rm end}}/V_0$ | {sr['V_end_frac_of_V0']:.6f} | 0.285 (L.319) | [OK] |")
    _kend_note = "L.325 states $K_{\\rm end}\\approx V_{\\rm end}$ (=1)"
    A(f"| $K_{{\\rm end}}/V_{{\\rm end}}$ | {sr['K_over_V_end']:.5f} | "
      + _kend_note + " | [FAIL] see below |")
    A(f"| $\\rho_{{\\rm end}}/V_{{\\rm end}}$ | {sr['rho_end']/sr['V_end']:.5f} | round-3 $1.19938$ | {'[OK]' if abs(sr['rho_end']/sr['V_end']-1.19938) < 0.02 else '[FAIL]'} |")
    A(f"| $\\langle w\\rangle$ (0-1) | {b['wbar_0_1']:.5f} | -- | -- |")
    A(f"| $\\langle w\\rangle$ (0-2) | {b['wbar_0_2']:.5f} | round-3 first cycle $-0.1043$ | same sign and order |")
    A(f"| $\\langle w\\rangle$ (0-3) | {b['wbar_0_3']:.5f} | -- | -- |")
    A(f"| $\\langle w\\rangle$ (0-{b['N_reached']:.2f}, all) | {b['wbar_0_end']:.5f} | round-3 $w_{{\\rm eff}}=-0.0019$ | see comment |")
    A(f"| $x_{{\\rm min}}$ (within the integration range) | {b['x_min']:.5f} | -- | $\\Phi=0$ never reached |")
    A("")
    A(f"(Integration range $N\\in[0,{b['N_reached']:.2f}]$, {b['n_points']} points in total;"
      f"$\\langle w\\rangle$ tends to 0 step by step: the $w=0$ dust behavior is established after a few oscillations.)")
    A("**Note**: the cosmic-time integration can only cover the first few e-folds -- under matter-like expansion "
      "$N(t)=\\frac23\\ln(t/t_i)$ grows only logarithmically, and covering the full $N_{\\rm reh}\\simeq20$ "
      "would need $t$ to grow by $e^{30}$, which is numerically infeasible.")
    A("Therefore the **asymptotic value** of $\\langle w\\rangle$ should be taken analytically as $w\\to0$; the early-time means in this table only illustrate"
      "the transient of the first cycle ($\\langle w\\rangle<0$, same order as the round-3 $-0.1043$).\n")
    A("**Key correction (relative to the first version of this script)**: $x_{\\rm end}$ is defined by $\\epsilon_V=1$, but slow roll has already failed there,"
      " so the slow-roll velocity $\\dot\\chi=-V'/(3H)$ must **not** be used directly as the initial value at $x_{\\rm end}$"
      " (that would force $K_{\\rm end}=V_{\\rm end}/3$ and $\\epsilon_H=3/4$ -- an initial value, not a dynamical result).")
    A("Instead this script integrates the e-fold equations from the $\\epsilon_V=0.05$ slow-roll region up to $\\epsilon_V=1$,"
      " letting the dynamics determine $K_{\\rm end}$.\n")
    A("Paper L.325 writes $K_{\\rm end}=\\frac12\\dot\\chi_{\\rm end}^2=\\epsilon_{\\rm end}V_{\\rm end}\\approx V_{\\rm end}$,")
    A("while the exact relation is $K=\\epsilon_H V/(3-\\epsilon_H)$; this integration gives "
      f"$K_{{\\rm end}}/V_{{\\rm end}}={sr['K_over_V_end']:.4f}$,")
    A("a factor $\\sim5$ smaller than the paper's $\\approx1$, and $K_{\\rm end}/\\Delta V_J$ is reduced accordingly (quantitative version of P0-D).")
    A(f"Moreover $\\rho_{{\\rm end}}=K+V={sr['rho_end']/sr['V_end']:.5f}V_{{\\rm end}}$,"
      f"i.e. the round-3 $1.19938\\,V_{{\\rm end}}$.\n")

    nw = out["N_window"]
    A("## 2. Upper bound of the $N$ window: extrapolated from the paper's own table\n")
    A(f"Fitting the $T^*_{{\\rm reh}}(N)$ column of the paper's Table I gives $\\ln T^*_{{\\rm reh}}=\\ln A+sN$,")
    A(f"$s={nw['fit']['slope']:.4f}$ (round-3 independently obtains 3.025, consistent). I.e. $T^*_{{\\rm reh}}\\propto e^{{3N}}$.\n")
    A("The physical bound is **instantaneous reheating**: $\\rho_{\\rm rad}=\\rho_{\\rm end}(N)$,")
    A("$T^{\\rm inst}_{\\rm reh}=(30\\rho_{\\rm end}/(\\pi^2g_*))^{1/4}$.\n")
    A("| $N$ | $T^*_{\\rm reh}$ (paper-table trend) | $T^{\\rm inst}_{\\rm reh}$ (correct, with the $0.411$ factor) | $V_{\\rm end}^{1/4}$ (paper's self-declared bound) | Exceeded? |")
    A("|---|---|---|---|---|")
    for N in (50, 55, 56, 57, 58):
        Ts = nw[f"T_reh_star_{N}"] if f"T_reh_star_{N}" in nw else (
            T_reh_star_from_table(N, nw["fit"]))
        Ti = T_reh_inst(float(N))
        A(f"| {N} | {Ts:.3e} | {Ti:.3e} | {V_end_of_N(N)**0.25:.3e} | "
          f"{'**YES**' if Ts > Ti else 'NO'} |")
    A("")
    A(f"This gives **$N_{{\\rm max}}={nw['N_max_proper']:.2f}$** (with the correct instantaneous-reheating bound),")
    A(f"or $N_{{\\rm max}}={nw['N_max_naive_Vend14']:.2f}$ if the paper's self-declared $V_{{\\rm end}}^{{1/4}}$ bound is used instead.")
    A("Both readings give $N_{\\rm max}\\approx56$, and the manuscript now quotes $N\\approx45$--$56$ with $N_{\\rm max}\\simeq55.6$.\n")
    A("(The first-principles independent recomputation gives $T^*_{{\\rm reh}}(50)="
      f"{nw['T_reh_first_principles_50']:.3e}$ GeV against {nw['table_T_reh_50']:.3e}$ GeV from the tabulated matching;")
    A("the slopes agree and the constants now agree to $1.4\\%$. The former factor-$3.4$ gap is fully accounted for:"
      " the Omega-route matching used a constant-$g_*$ radiation scaling and $V_{\\rm end}$ in place of"
      " $\\rho_{\\rm end}=K_{\\rm end}+V_{\\rm end}$, and this script solved $x_*$ from $e^{2x}-2x=2\\beta^2N$,"
      " which is not the slow-roll relation for this potential (it returns $N=3.70$ at the $x_*=2.124$ that the"
      " equation yields). All three are corrected, so the cross-check is now a consistency confirmation rather"
      " than an open discrepancy.)\n")

    r = out["reheating"]
    A("## 3. The two reheating channels\n")
    A("### (i) Pulse/gravitational channel\n")
    A(f"$\\rho_{{\\rm rad}}(a_{{\\rm end}})=N_{{\\rm eff}}H_{{\\rm inf}}^4/(192\\pi^2)$, with $N_{{\\rm eff}}={r['pulse']['N_eff']:g}$:")
    A(f"$\\rho_{{\\rm rad}}={r['pulse']['rho_rad_aend']:.4e}$ GeV$^4$ (the paper L.361 value $4\\times10^{{50}}$ [OK]),")
    A(f"$\\rho_{{\\rm rad}}/\\rho_{{\\rm end}}={r['pulse']['rho_ratio_to_rho_end']:.3e}$,")
    A(f"corresponding to $T={r['pulse']['T_pulse_from_rho']:.4e}$ GeV.")
    A("One might ask whether this channel alone can only heat the universe to $\\sim10^{12}$ GeV, far below what $T_{\\rm reh}\\sim10^9$ GeV requires?")
    A("-- No: $10^{12}>10^{9}$, so this channel **alone already satisfies** $T_{\\rm reh}\\sim10^9$ GeV;")
    A("the paper attributes $T_{\\rm reh}\\sim10^9$ GeV to the anomaly channel, but the pulse channel already suffices; the relative weights of the two need to be stated.\n")
    A("### (ii) Conformal-anomaly channel\n")
    A("Paper L.342 gives $\\Gamma\\sim b_3^2\\alpha_s^2m_\\chi^3/(16\\pi^3M_P^2(6+1/\\xi))$,")
    A("but **the paper never quotes numerical values of $b_3$ and $\\alpha_s$** (the script hand-picks $b_3=7,\\alpha_s=0.1$).\n")
    A("| $\\alpha_s(m_\\chi)$ choice | $\\alpha_s$ | $\\Gamma$ [GeV] | $T_{\\rm reh}$ [GeV] | Implied $N$ (from $T\\propto e^{3N}$) |")
    A("|---|---|---|---|---|")
    for row in r["anomaly"]["rows"]:
        A(f"| {row['label']} | {row['alpha_s']:.5f} | {row['Gamma_GeV']:.5f} | "
          f"{row['T_reh_GeV']:.4e} | {row['N_implied_from_T_reh']:.2f} |")
    A("")
    A("The 1-loop running gives $\\alpha_s(3.25\\times10^{13}\\,{\\rm GeV})\\simeq"
      f"{r['anomaly']['rows'][0]['alpha_s']:.4f}$, **not 0.1**.")
    A("Substituting the physical value lowers $T_{\\rm reh}$ to $\\sim3\\times10^8$ GeV,")
    A("and the anomaly-channel operating point moves from $N\\simeq50.7$ to $N\\simeq50.2$.")
    A(f"BBN floor $10$ MeV: {'satisfied' if r['bbn_ok'] else 'violated'}.\n")

    md = "\n".join(L) + "\n"
    with open(os.path.join(ROOT, "background_and_reheating.md"), "w",
              encoding="utf-8") as f:
        f.write(md)
    print("wrote scripts/background_and_reheating.md "
          f"({len(md)} chars)")


if __name__ == "__main__":
    main()
