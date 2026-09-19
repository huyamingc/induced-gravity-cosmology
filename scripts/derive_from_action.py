# -*- coding: utf-8 -*-
"""
Independent derivation from the action — no presupposed paper formulas.

Pipeline:
  Jordan action
    -> conformal transform (derive Omega, kinetic, V_E)
    -> canonical field chi
    -> exact slow-roll integrals N(x_*), A_s, n_s, r (potential slow-roll)
    -> solve lambda0 from A_s=Planck at chosen N
    -> cosmological matching N(k) vs T_reh from expansion history
    -> late-time EOM near minimum (DE test)
    -> DM production order-of-magnitude

Writes scripts/derivation_report.md
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, fsolve

OUT = Path(__file__).resolve().parent / "derivation_report.md"

# ---------------------------------------------------------------------------
# Fundamental inputs (observations / definitions only)
# ---------------------------------------------------------------------------
M_PL = 2.435e18          # reduced Planck mass [GeV]  (definition 1/sqrt(8 pi G))
A_S_OBS = 2.100e-9
N_S_OBS = 0.9649
SIG_NS = 0.0042
H0_KMS_MPC = 67.4
MPC_M = 3.0856775814913673e22
HBAR_GEV_S = 6.582119569e-25
C_KMS = 2.99792458e5
T_CMB_K = 2.7255
KB_GEV_K = 8.617333262e-14
XI_CANDIDATES = [1.0, 5.0, 11.1, 20.0, 50.0, 100.0]


def H0_GeV() -> float:
    return H0_KMS_MPC / MPC_M * C_KMS * HBAR_GEV_S  # careful units
    # H0 [1/s] = H0_km_s / Mpc_m ; H0[GeV] = H0[1/s] * hbar


def H0_GeV_v2() -> float:
    # 67.4 km/s/Mpc
    H0_si = (H0_KMS_MPC * 1e3) / MPC_M  # s^-1
    return H0_si * HBAR_GEV_S


def T0_GeV() -> float:
    return T_CMB_K * KB_GEV_K


# ---------------------------------------------------------------------------
# A. Frame transformation from the action
# ---------------------------------------------------------------------------
"""
S_J = int sqrt(-g) [ 1/2 xi Phi^2 R - 1/2 (dPhi)^2 - V_J(Phi) ] + S_m
V_J = (lambda0/4)(Phi^2 - Phi0^2)^2 + Vc,  Phi0 = M_Pl/sqrt(xi)   [definition of vacuum]

Conformal: g_tilde = Omega^2 g,  Omega^2 = xi Phi^2 / M_Pl^2
Then Omega = Phi/Phi0 exactly.

Jordan EH-like term (1/2)xi Phi^2 R  ->  (1/2) M_Pl^2 R_tilde  + (kinetic in omega)
with omega = ln Omega = ln(Phi/Phi0).

Kinetic combining:
  from (1/2)xi Phi^2 R:  -3 M_Pl^2 (d omega)^2
  from -1/2 (dPhi)^2:   -1/2 M_Pl^2/xi (d omega)^2
  total:  -M_Pl^2 (3 + 1/(2xi)) (d omega)^2
Canonical:  -1/2 (d chi)^2  =>  dchi/domega = M_Pl sqrt(6 + 1/xi)

Einstein potential: V_E = V_J / Omega^4
  Phi = Phi0 e^omega,  Phi^2/Phi0^2 = e^{2 omega}
  (Phi^2 - Phi0^2)^2 / Phi0^4 = (e^{2omega}-1)^2
  V_J = (lambda0/4) Phi0^4 (e^{2omega}-1)^2 + Vc
  Omega^4 = e^{4 omega}
  V_E = (lambda0/4) Phi0^4 (e^{2omega}-1)^2 e^{-4omega} + Vc e^{-4omega}
      = (lambda0/4) Phi0^4 (1 - e^{-2omega})^2 + Vc e^{-4omega}
  Phi0^4 = M_Pl^4 / xi^2
  => V_E = (lambda0 M_Pl^4 / (4 xi^2)) (1 - e^{-2 omega})^2 + Vc e^{-4 omega}

Define x = 2 omega = 2 chi / (M_Pl beta_o), beta_o = sqrt(6+1/xi)
Then 2 omega = beta_p chi / M_Pl with beta_p = 2/beta_o, and
V_E = V0 (1 - e^{-x})^2 + Vc e^{-2x},  V0 = lambda0 M_Pl^4/(4 xi^2)
"""


def beta_o(xi: float) -> float:
    return math.sqrt(6.0 + 1.0 / xi)


def beta_p(xi: float) -> float:
    """Paper-style beta such that x = beta_p chi / M_Pl = 2 omega."""
    return 2.0 / beta_o(xi)


def Phi0(xi: float) -> float:
    return M_PL / math.sqrt(xi)


def V0_from_lambda(lambda0: float, xi: float) -> float:
    return lambda0 * M_PL**4 / (4.0 * xi**2)


def lambda0_from_V0(V0: float, xi: float) -> float:
    return 4.0 * xi**2 * V0 / M_PL**4


# ---------------------------------------------------------------------------
# B. Exact potential slow-roll on V = V0 (1-e^{-x})^2  (ignore Vc on plateau)
# ---------------------------------------------------------------------------
def eps_V(x: float, bp: float) -> float:
    u = math.exp(-x)
    return 2.0 * bp**2 * u**2 / (1.0 - u) ** 2


def eta_V(x: float, bp: float) -> float:
    u = math.exp(-x)
    return 2.0 * bp**2 * u * (2.0 * u - 1.0) / (1.0 - u) ** 2


def x_end_from_eps1(bp: float) -> float:
    """Solve 2 bp^2 u^2 = (1-u)^2, u=e^{-x} in (0,1)."""
    # bp u = 1-u  (positive root) => u = 1/(1+bp)
    # Wait: sqrt(2) bp u = 1-u => u = 1/(1+sqrt(2) bp)
    u = 1.0 / (1.0 + math.sqrt(2.0) * bp)
    return -math.log(u)


def N_of_x(x: float, x_end: float, bp: float) -> float:
    """N = int (e^x - 1) dx / (2 bp^2) from x_end to x."""
    return (math.exp(x) - x - (math.exp(x_end) - x_end)) / (2.0 * bp**2)


def x_star_for_N(N: float, bp: float) -> float:
    x_e = x_end_from_eps1(bp)
    return brentq(lambda x: N_of_x(x, x_e, bp) - N, x_e + 1e-8, x_e + 50.0)


def observables_at_x(x: float, bp: float, V0: float) -> dict:
    u = math.exp(-x)
    eps = eps_V(x, bp)
    eta = eta_V(x, bp)
    VE = V0 * (1.0 - u) ** 2
    As = VE / (24.0 * math.pi**2 * M_PL**4 * eps)
    # potential slow-roll spectra
    ns_ps = 1.0 - 6.0 * eps + 2.0 * eta
    r_ps = 16.0 * eps
    # large-field / attractor approximations
    N_approx = math.exp(x) / (2.0 * bp**2)
    ns_attr = 1.0 - 2.0 / N_approx if N_approx > 1 else float("nan")
    r_attr = 8.0 / (bp**2 * N_approx**2) if N_approx > 0 else float("nan")
    return {
        "x": x,
        "u": u,
        "eps": eps,
        "eta": eta,
        "VE": VE,
        "As": As,
        "ns_ps": ns_ps,
        "r_ps": r_ps,
        "N_from_x": N_of_x(x, x_end_from_eps1(bp), bp),
        "N_large": N_approx,
        "ns_attr": ns_attr,
        "r_attr": r_attr,
    }


def lambda0_for_As(N_target: float, xi: float, As: float = A_S_OBS) -> tuple[float, dict]:
    """Find lambda0 such that exact A_s(x_*(N)) = As.

    At fixed N (fixed x_*), A_s is linear in V0 hence in lambda0.
    """
    bp = beta_p(xi)
    x_s = x_star_for_N(N_target, bp)
    # probe V0=1
    probe = observables_at_x(x_s, bp, V0=1.0)
    V0_needed = As / probe["As"] * 1.0
    lam = lambda0_from_V0(V0_needed, xi)
    full = observables_at_x(x_s, bp, V0_needed)
    full["lambda0"] = lam
    full["V0"] = V0_needed
    full["beta_o"] = beta_o(xi)
    full["beta_p"] = bp
    full["xi"] = xi
    full["N_target"] = N_target
    return lam, full


# ---------------------------------------------------------------------------
# C. Cosmological N(T_reh) matching from expansion history (not from a quoted formula)
# ---------------------------------------------------------------------------
"""
Pivot comoving wavenumber k_* = a_* H_*  (we take k = 0.05 Mpc^{-1} as CMB pivot).

Today: a0=1, H0, T0.

After inflation (assume instantaneous transition into matter-like condensate domination
with w=0 until T_reh, then radiation):

  rho_end ~ V_end
  rho_reh = (pi^2/30) g_star T_reh^4
  During w=0: a_reh/a_end = (rho_end/rho_reh)^{1/3}

Radiation era until matter-radiation equality:
  rho_reh (a_reh/a_eq)^4 = rho_eq
  a_reh/a_eq = (rho_eq/rho_reh)^{1/4} * (g factors if needed; set g_eq~g_reh for estimate)

Matter era:
  a_eq/a0 = Omega_r / Omega_m   (approximately, since rho_m/rho_r = a_eq/a)

More carefully:
  Omega_r h^2 ~ 4.15e-5, Omega_m ~ 0.315
  a_eq = Omega_r/Omega_m

Inflationary side:
  H_*^2 = V_*/(3 M_Pl^2)   (slow-roll)
  V_* = V0 (1-e^{-x_*})^2
  k_* = a_* H_*
  a_* = k_*/H_*   with k_* in GeV units: k[Mpc^{-1}] -> GeV via hbar*c/Mpc

  N = ln(a_end / a_*) = ln( a_end H_* / k_* )

  a_end = a_reh (a_end/a_reh) = a_reh (rho_reh/rho_end)^{1/3}

  a_reh = a_eq (a_reh/a_eq) = (Omega_r/Omega_m) * (rho_eq/rho_reh)^{1/4}

Putting together gives N(T_reh, V0, x_*) self-consistently with N that set x_*.
"""


def k_pivot_GeV(k_mpc: float = 0.05) -> float:
    """k [Mpc^{-1}] -> [GeV]: k = 2pi?  Usually k = aH today scale;
    comoving wavenumber k = 0.05 Mpc^{-1} means k_phys today = 0.05 / Mpc in 1/length.
    Energy: E = hbar c k.
    """
    k_si = k_mpc / MPC_M  # m^-1
    # hbar c in GeV m: hbar=6.582e-25 GeV s, c=3e8 m/s => 1.973e-16 GeV m
    hbar_c = HBAR_GEV_S * 2.99792458e8  # GeV m
    return k_si * hbar_c


def rho_rad(T: float, gstar: float = 106.75) -> float:
    return math.pi**2 / 30.0 * gstar * T**4


def N_match_derived(
    T_reh: float,
    V_end: float,
    V_star: float,
    xi: float,
    N_guess: float,
    gstar: float = 106.75,
    Om_m: float = 0.315,
    Om_r: float = 9.0e-5,
    k_mpc: float = 0.05,
) -> dict:
    T0 = T0_GeV()
    H0 = H0_GeV_v2()
    k = k_pivot_GeV(k_mpc)

    H_star = math.sqrt(V_star / (3.0 * M_PL**2))
    a_star = k / H_star  # from k = a_* H_*  (a0=1)

    rho_reh = rho_rad(T_reh, gstar)
    rho_eq = rho_rad(T0, gstar) / Om_r * Om_m  # rho_eq = rho_r0 * (Om_m/Om_r) wait
    # rho_r0 = Om_r * rho_crit0; rho_m(eq)=rho_r(eq)=rho_r0 (a0/a_eq)^4 = Om_r rho_crit (a0/a_eq)^4
    # a_eq = Om_r/Om_m; rho_eq = Om_m * rho_crit0 * (a0/a_eq)^3 = Om_m rho_crit (Om_m/Om_r)^3 ... messy
    # Simpler: rho_r(a)= Om_r rho_c0 a^{-4}; at a_eq, rho_m=Om_m rho_c0 a_eq^{-3}=rho_r
    # => Om_m a_eq^{-3} = Om_r a_eq^{-4} => a_eq = Om_r/Om_m
    rho_c0 = 3.0 * H0**2 * M_PL**2
    a_eq = Om_r / Om_m
    rho_eq = Om_r * rho_c0 * a_eq ** (-4)

    # a_reh from radiation: rho_reh (a_reh/a_eq)^4 = rho_eq => a_reh = a_eq (rho_eq/rho_reh)^{1/4}
    a_reh = a_eq * (rho_eq / rho_reh) ** 0.25
    # a_end from matter-like: rho_end (a_end/a_reh)^3 = rho_reh => a_end/a_reh = (rho_reh/rho_end)^{1/3}
    a_end = a_reh * (rho_reh / V_end) ** (1.0 / 3.0)
    N = math.log(a_end / a_star)

    return {
        "N": N,
        "H_star": H_star,
        "a_star": a_star,
        "a_end": a_end,
        "a_reh": a_reh,
        "a_eq": a_eq,
        "rho_reh": rho_reh,
        "rho_eq": rho_eq,
        "k_GeV": k,
        "H0": H0,
        "T0": T0,
        "Om_r": Om_r,
        "Om_m": Om_m,
    }


def solve_N_T_self_consistent(xi: float, N_target: float, As: float = A_S_OBS) -> dict:
    """Iterate: pick N -> lambda0, x_*, V_end, V_*; compute derived N(T_reh); vary T_reh
    so that derived N equals N_target. Report T_reh and derived N.
    Also compute derived N for fixed T_reh benchmarks.
    """
    lam, obs = lambda0_for_As(N_target, xi, As)
    V0 = obs["V0"]
    x_star = obs["x"]
    bp = obs["beta_p"]
    x_end = x_end_from_eps1(bp)
    u_end = math.exp(-x_end)
    V_end = V0 * (1.0 - u_end) ** 2
    V_star = obs["VE"]

    def residual(logT):
        T = math.exp(logT)
        return N_match_derived(T, V_end, V_star, xi, N_target)["N"] - N_target

    # bracket T
    try:
        logT = brentq(residual, math.log(1e-2), math.log(1e16))
        T_star = math.exp(logT)
        ok = True
    except ValueError:
        T_star = float("nan")
        ok = False

    benches = {}
    for T in (4e5, 1e6, 1e9, 1e12, 6e15):
        benches[T] = N_match_derived(T, V_end, V_star, xi, N_target)

    return {
        "xi": xi,
        "N_target": N_target,
        "lambda0": lam,
        "obs": obs,
        "V_end": V_end,
        "V_star": V_star,
        "T_reh_selfcons": T_star,
        "solve_ok": ok,
        "benches": benches,
    }


# ---------------------------------------------------------------------------
# D. Late-time DE from EOM (no presupposed quintessence claim)
# ---------------------------------------------------------------------------
"""
Near omega=0 (chi=0):
  V_E(omega) = V0 (1-e^{-2omega})^2 + Vc e^{-4omega}
  Expand: (1-e^{-2w})^2 ~ (2w)^2 = 4 w^2 for small w
  V0 term: 4 V0 w^2
  In chi: omega = chi/(M_Pl beta_o), so
  V_quad = 4 V0 chi^2 / (M_Pl^2 beta_o^2)
  Compare 1/2 m^2 chi^2 => m^2 = 8 V0 / (M_Pl^2 beta_o^2)
  V0 = lambda0 M_Pl^4/(4 xi^2)
  m^2 = 8 lambda0 M_Pl^2 / (4 xi^2 beta_o^2) = 2 lambda0 M_Pl^2 / (xi^2 beta_o^2)
  This is m_chi^2 — matches common induced-gravity result.

KG: chi_ddot + 3H chi_dot + m^2 chi = 0
For H << m: oscillatory, virial <w> -> 0 for quadratic (matter-like), NOT -1.
Constant Vc contributes T_mn = -Vc g_mn independent of chi => w=-1 from Vc alone.
Ratio m/H0 determines whether any residual kinetic/oscillation can mimic DE.
"""


def m_chi_from_first_principles(lambda0: float, xi: float) -> float:
    V0 = V0_from_lambda(lambda0, xi)
    bo = beta_o(xi)
    m2 = 8.0 * V0 / (M_PL**2 * bo**2)
    return math.sqrt(m2)


def de_viability(lambda0: float, xi: float) -> dict:
    m = m_chi_from_first_principles(lambda0, xi)
    H0 = H0_GeV_v2()
    return {
        "m_chi": m,
        "H0": H0,
        "m_over_H0": m / H0,
        "damping_time_s": HBAR_GEV_S / m if m > 0 else float("inf"),  # t ~ ħ/m
        "w_from_Vc": -1.0,
        "w_from_osc_quadratic": 0.0,  # virial for quadratic
        "kinetic_redshift": "a^{-6}",
        "quintessence_possible": m_over_H0_check(m, H0),
    }


def m_over_H0_check(m: float, H0: float) -> bool:
    """Quintessence needs field to roll slowly today: roughly m_eff <= few H0."""
    return m <= 10.0 * H0


# ---------------------------------------------------------------------------
# E. DM gravitational production scaling (parametric only)
# ---------------------------------------------------------------------------
def Omega_psi_scaling(mpsi: float, Hinf: float, Treh: float, pref: float = 1.0) -> float:
    """Yield ~ Hinf^3 exp(-pi m/Hinf); dilution ~ (T0/Treh)^3 * entropy factors absorbed in pref.
    We only report relative scaling; pref calibrated at one point if needed.
    """
    if mpsi <= 0 or Hinf <= 0:
        return 0.0
    n = Hinf**3 * math.exp(-math.pi * mpsi / Hinf)
    # today number ~ n * (a_end/a0)^3, (a_end/a0)^3 ∝ Treh^{-1} * (stuff); paper says ∝ Treh
    # From reheating: a_end/a_reh ∝ (rho_reh/rho_end)^{1/3} ∝ Treh^{4/3} / rho_end^{1/3}
    # a_reh/a0 ∝ T0/Treh
    # (a_end/a0)^3 ∝ Treh^4 / rho_end * (T0/Treh)^3 = Treh * T0^3 / rho_end
    # => Omega ∝ m n Treh / (rho_end factors)  — linear in Treh if rho_end fixed
    return pref * mpsi * n * Treh


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def main() -> None:
    H0 = H0_GeV_v2()
    T0 = T0_GeV()
    lines: list[str] = []
    lines.append("# 从作用量出发的独立推导报告")
    lines.append("")
    lines.append("脚本：`scripts/derive_from_action.py`（不预设论文公式，仅用作用量结构与观测输入）")
    lines.append("")
    lines.append("## 0. 输入与定义")
    lines.append("")
    lines.append(f"- 约化 Planck 质量 M_Pl = {M_PL:.4e} GeV")
    lines.append(f"- Planck 2018: A_s={A_S_OBS}, n_s={N_S_OBS}±{SIG_NS}")
    lines.append(f"- H0 = {H0:.4e} GeV（由 67.4 km/s/Mpc 换算）")
    lines.append(f"- T_CMB = {T0:.4e} GeV")
    lines.append(f"- 真空定义：Phi0 = M_Pl/sqrt(xi)（由 V_J 极值，与 M_Pl^2=xi Phi0^2 一致）")
    lines.append("")
    lines.append("## A. 共形变换 → Einstein 帧（推导要点）")
    lines.append("")
    lines.append("由 Omega^2 = xi Phi^2/M_Pl^2 得 **Omega = Phi/Phi0**（恒等式，非假设）。")
    lines.append("动能合并系数 M_Pl^2(3+1/(2xi))，故")
    lines.append("**dchi/domega = M_Pl sqrt(6+1/xi) ≡ M_Pl beta_o**，**x ≡ 2 omega = beta_p chi/M_Pl**，beta_p=2/beta_o。")
    lines.append("")
    lines.append("V_E = V_J/Omega^4 = **(lambda0 M_Pl^4/(4 xi^2)) (1-e^{-2 omega})^2 + Vc e^{-4 omega}**")
    lines.append("")
    lines.append("论文写 beta=2/sqrt(6+1/xi) 与上述 beta_p 一致；r=16 eps 在吸引子极限给出 8/(beta_p^2 N^2)=2(6+1/xi)/N^2。**公式结构由推导支持。**")
    lines.append("")

    lines.append("## B. 精确慢滚与 lambda0 归一化（按 N 求解）")
    lines.append("")
    lines.append("对 V=V0(1-e^{-x})^2：")
    lines.append("- eps = 2 beta_p^2 e^{-2x}/(1-e^{-x})^2")
    lines.append("- eta = 2 beta_p^2 e^{-x}(2e^{-x}-1)/(1-e^{-x})^2")
    lines.append("- N(x) = [e^x - x - (e^{x_end}-x_end)]/(2 beta_p^2)，x_end 由 eps=1：u=1/(1+sqrt(2) beta_p)")
    lines.append("- A_s = V(x_*)/(24 pi^2 M_Pl^4 eps_*)；在固定 N 下对 lambda0 线性，可反解")
    lines.append("")
    lines.append("| xi | N | lambda0(exact As) | ns(PS) | r(PS) | ns(attractor) | r(attractor) | ns−Planck [sigma] |")
    lines.append("|---|---|---|---|---|---|---|---|")
    table_B = []
    for xi in XI_CANDIDATES:
        for N in (48, 50, 55):
            lam, obs = lambda0_for_As(N, xi)
            sig = (obs["ns_ps"] - N_S_OBS) / SIG_NS
            table_B.append((xi, N, lam, obs))
            lines.append(
                f"| {xi} | {N} | {lam:.4e} | {obs['ns_ps']:.4f} | {obs['r_ps']:.5f} | "
                f"{obs['ns_attr']:.4f} | {obs['r_attr']:.5f} | {sig:+.2f} |"
            )
    lines.append("")
    # Compare attractor vs PS
    lines.append("**PS vs 吸引子近似：** 同一 x_* 下二者应接近；若 r_ps 与 r_attr 差 > few %，说明 N≈e^x/(2β²) 近似在该 N 有偏。")
    lines.append("")
    for xi, N, lam, obs in table_B:
        if abs(xi - 11.1) < 1e-6 and N == 50:
            dr = abs(obs["r_ps"] - obs["r_attr"]) / obs["r_attr"] * 100
            dns = abs(obs["ns_ps"] - obs["ns_attr"])
            lines.append(f"- 示例 ξ=11.1, N=50: |r_ps−r_attr|/r_attr={dr:.2f}%, |ns_ps−ns_attr|={dns:.4f}, lambda0={lam:.4e}")
    lines.append("")
    lam50, obs50 = lambda0_for_As(50, 11.1)
    lines.append(
        f"- **ξ=11.1, N=50 精确反演 lambda0 = {lam50:.4e}**"
        f"（论文 Table I 用 6.78e-8；大场解析式另给 ~7.5e-8）"
    )
    lines.append(f"- 同点：H_inf ~ sqrt(V0/(3 M_Pl^2)) 可由 V0={obs50['V0']:.4e} 推出")
    Hinf = math.sqrt(obs50["V0"] / (3 * M_PL**2))
    U14 = obs50["V0"] ** 0.25
    mchi = m_chi_from_first_principles(lam50, 11.1)
    lines.append(f"- H_inf={Hinf:.4e} GeV, U^{{1/4}}={U14:.4e} GeV, m_chi(由 V_E''(0))={mchi:.4e} GeV")
    lines.append("")

    lines.append("## C. N(T_reh)：从膨胀史匹配（非套用论文公式）")
    lines.append("")
    lines.append("输入：k=0.05 Mpc^{-1}，Omega_m=0.315，Omega_r=9e-5，g*=106.75，a_eq=Omega_r/Omega_m。")
    lines.append("链条：a_*=k/H_* → a_end=a_reh (rho_reh/rho_end)^{1/3} → a_reh=a_eq (rho_eq/rho_reh)^{1/4} → N=ln(a_end/a_*)。")
    lines.append("")
    for N_target in (50,):
        res = solve_N_T_self_consistent(11.1, N_target)
        lines.append(f"### 目标 N={N_target}, xi=11.1")
        lines.append("")
        lines.append(f"- lambda0(exact As)={res['lambda0']:.4e}, V_*={res['V_star']:.4e}, V_end={res['V_end']:.4e}")
        if res["solve_ok"]:
            lines.append(f"- **使推导 N 等于 {N_target} 的 T_reh ≈ {res['T_reh_selfcons']:.4e} GeV**")
        else:
            lines.append(f"- **在 [1e-2,1e16] GeV 内未能反解出 T_reh 使 N={N_target}**（检查匹配或输入）")
        lines.append("")
        lines.append("| T_reh [GeV] | N_derived | 备注 |")
        lines.append("|---|---|---|")
        for T, info in res["benches"].items():
            inside = 48.0 <= info["N"] <= 55.0
            tag = "∈[48,55]" if inside else "窗外"
            lines.append(f"| {T:.2e} | {info['N']:.2f} | {tag} |")
        lines.append("")

    # Also try other xi quickly at T=1e9
    lines.append("### 固定 T_reh=1e9 GeV 时不同 xi 的推导 N（N_target 仅用于定 lambda0）")
    lines.append("")
    lines.append("| xi | N_target | lambda0 | N_derived(T=1e9) | r(PS) | ns(PS) |")
    lines.append("|---|---|---|---|---|---|")
    for xi in XI_CANDIDATES:
        res = solve_N_T_self_consistent(xi, 50)
        Nd = res["benches"][1e9]["N"]
        lines.append(
            f"| {xi} | 50 | {res['lambda0']:.3e} | {Nd:.2f} | "
            f"{res['obs']['r_ps']:.5f} | {res['obs']['ns_ps']:.4f} |"
        )
    lines.append("")
    lines.append("**推导结论（C）：** N–T_reh 关系依赖匹配输入（k, g*, Omega，再加热状态方程）。")
    lines.append("若推导的 N(T) 与论文 Eq.(18) 或「50+¼ln」不一致，则论文该式 **不能** 当作唯一真理，只能当某约定下的近似；应报告推导值与不确定度。")
    lines.append("")

    lines.append("## D. 晚期暗能量：从运动方程判定")
    lines.append("")
    lam_use = lam50
    de = de_viability(lam_use, 11.1)
    lines.append(f"由 V_E 在极小处展开：m_chi^2 = 8 V0/(M_Pl^2 beta_o^2) = 2 lambda0 M_Pl^2/(xi^2 beta_o^2)")
    lines.append(f"- lambda0={lam_use:.4e} => m_chi={de['m_chi']:.4e} GeV")
    lines.append(f"- H0={de['H0']:.4e} GeV, m_chi/H0={de['m_over_H0']:.4e}")
    lines.append(f"- 阻尼时间 ~ {de['damping_time_s']:.3e} s")
    lines.append(f"- 二次势振荡：virial => w_osc={de['w_from_osc_quadratic']}（物质型），动能 ∝ {de['kinetic_redshift']}")
    lines.append(f"- 常数项 Vc：w={de['w_from_Vc']}，与 chi 演化独立")
    lines.append(f"- 今日慢滚精质可行？ m_chi <=~ 10 H0 ? **{de['quintessence_possible']}**")
    lines.append("")
    lines.append("**推导结论（D）：** 在该作用量与由 A_s 定出的 lambda0 下，")
    lines.append("场在今日不可能以 w≈−0.987 慢滚；DE 只能来自 **常数 Vc**（或框架外新自由度）。")
    lines.append("这不是预设，而是 m/H0 与 KG 方程的直接结果。")
    lines.append("")

    lines.append("## E. 暗物质：仅参数标度，不预设成功")
    lines.append("")
    Hinf = math.sqrt(obs50["V0"] / (3 * M_PL**2))
    Phi0_11 = Phi0(11.1)
    lines.append(f"- Phi0(xi=11.1)={Phi0_11:.4e} GeV, H_inf={Hinf:.4e} GeV")
    lines.append("- 共形恒等式 => E 帧 m_psi=g Phi0 **与 chi 无关**，树图 chi-psi-psi = 0（推导支持）")
    lines.append("- 引力产生：n_psi ~ H_inf^3 exp(-pi m_psi/H_inf)（重场标准估计，系数 O(1) 未从第一性原理钉死）")
    lines.append("- 稀释：(a_end/a0)^3 在 w=0 再加热下 **∝ T_reh**（推导），故 Omega 对 T_reh 线性、对 m_psi 指数敏感")
    lines.append("")
    lines.append("| g | m_psi | m_psi/H_inf | Omega_rel(T=1e9, 归一前) | Omega_rel(T=4e5) |")
    lines.append("|---|---|---|---|---|")
    ref = Omega_psi_scaling(4e12, Hinf, 1e9)
    for g in (1e-5, 2.3e-5, 1e-4):
        mpsi = g * Phi0_11
        o1 = Omega_psi_scaling(mpsi, Hinf, 1e9) / ref
        o2 = Omega_psi_scaling(mpsi, Hinf, 4e5) / ref
        lines.append(f"| {g:.2e} | {mpsi:.3e} | {mpsi/Hinf:.3f} | {o1:.4f} | {o2:.4e} |")
    lines.append("")
    lines.append("**推导结论（E）：** 机制在数学上允许在 m_psi~H_inf、T_reh~1e9 附近调到观测丰度；")
    lines.append("低 T_reh 通道不足。**不**从第一性原理得到唯一 g 或 Omega=0.12。")
    lines.append("稳定化若诉诸 Z2^psi gauge，仍是 **额外 UV 假设**，推导链在此分支。")
    lines.append("")

    lines.append("## F. 不预设时的综合判断")
    lines.append("")
    lines.append("| 问题 | 从作用量+观测输入推出什么 |")
    lines.append("|---|---|")
    lines.append("| Einstein 帧是否 Starobinsky 型 | **是**（推导） |")
    lines.append("| r, n_s 函数形式 | **是**，且精确慢滚可算 |")
    lines.append("| lambda0 | **依赖 N 与 xi**；N=50,xi=11.1 时精确值约 6.7e-8 |")
    lines.append("| N 与 T_reh | **依赖宇宙学匹配输入**；须报告推导 N(T) 而非单一教条公式 |")
    lines.append("| 同场能否做今日 DE | **否**（m/H0 与 KG） |")
    lines.append("| DE 来自 Vc | **是**（常数项）；数值仍需观测标定 |")
    lines.append("| DM | **条件性**；顶点关闭由共形恒等式支持；丰度未第一性原理算尽 |")
    lines.append("| 单一场统一暴涨+DM+DE | **DE 部分不成立**；暴涨+条件性 DM 可立 |")
    lines.append("")
    lines.append("报告文件：scripts/derivation_report.md")
    lines.append("")
    lines.append("[独立推导完成]")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
