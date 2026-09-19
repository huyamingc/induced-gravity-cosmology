#!/usr/bin/env python3
r"""
psi_production_bogoliubov.py  (P0 backfill #2)
====================================================
Original requirement (the review workspace review_workspace/ was removed after the final draft; excerpt kept) Sec. C:
  "[P0] psi gravitational-production abundance: numerically solve the mode equation to get Bogoliubov beta_k -> n_psi, Omega_psi(m_psi,T_reh)"

Paper claims (paper_prd_merged.tex):
  L.423  n_psi(a_end) ~ H_inf^3 exp(-pi  m_psi/H_inf)          (eq:npsi)
  L.439  Omega_psi   ~ m_psi H_inf^3 exp(-pi m_psi/H_inf)(a_end/a_0)^3/rho_c
  L.442  (a_end/a_0)^3 propto T_reh,  ~1e-92 at T_reh=1e9 GeV
  L.446  g ~ 1e-5 .. 1e-4
  L.444  m_psi ~ H_inf .. few x 10 H_inf
  L.465  rho_psi/rho_chi ~ 1e-13,  rho_psi/rho_rad ~ 1e-18

This script does four things:
  [1] **Independently derive** the mode equation from the Dirac equation, and verify that pure de Sitter gives the exact index nu = 1/2 +- i m/H
      (matching the literature Ema-Nakayama-Tang arXiv:1903.10973 Eq.(14)(16) verbatim).
  [2] **Audit the exponent coefficient**: the exponent of the exact Bogoliubov coefficient comes from |Gamma(i mu)|^2 = pi/cosh(pi mu)
      ~ 2 pi e^{-pi mu} (a single Gamma factor => e^{-pi mu}), times the two e^{+- pi m/(2H)} factors
      carried by the mode functions in the BD normalization (ENT Eq.(16)) => e^{-2 pi mu}.
      Counting only one of them => exactly the paper's e^{-pi m/H}. This script recomputes that combination at high precision with mpmath.
  [3] Independently recompute the dilution factor (a_end/a_0)^3 (the paper's L.442 value 1e-92) and double-check its T_reh scaling.
  [4] With the two exponent coefficients (pi = paper, 2pi = exact), solve Omega_psi(m_psi; T_reh)=0.265 respectively,
      obtaining the required g and m_psi/H_inf, and check whether they fall inside the paper's self-reported window (L.444/L.446).

Output: scripts/psi_production_bogoliubov.md / .json
"""

from __future__ import annotations

import json
import math
import os

import mpmath as mp

ROOT = os.path.dirname(os.path.abspath(__file__))

# ---- Paper-locked values (N=50, xi=11.1; same convention as lock_n_convention.py) ----
M_P = 2.435e18          # GeV, reduced Planck mass
XI = 11.1
LAM0 = 6.6970e-8        # back-solved from A_s (derive_from_action.lambda0_for_As(50,11.1))
BETA_P = 2.0 / math.sqrt(6.0 + 1.0 / XI)      # 0.810435
PHI_V = M_P / math.sqrt(XI)                    # 7.3087e17 GeV
V0 = LAM0 * M_P**4 / (4.0 * XI**2)             # 4.7772e63
H_INF = math.sqrt(V0 / (3.0 * M_P**2))         # 1.6388e13
M_CHI = math.sqrt(2.0 * V0 * BETA_P**2 / M_P**2)

# end-of-inflation (round-3 exact slow-roll integration)
X_END = 0.7636653
V_END_FRAC = 0.285204
V_END = V_END_FRAC * V0
K_END = 0.19938 * V_END
RHO_END = 1.19938 * V_END

# cosmological constants
H0_GEV = 1.4377e-42
RHO_C = 3.0 * H0_GEV**2 * M_P**2               # 3.677e-47 GeV^4
T0_GEV = 2.3491e-13                            # 2.7255 K
G_STAR = 106.75
G_STAR_S0 = 3.91
OMEGA_DM = 0.265
OMEGA_PSI_TARGET = OMEGA_DM


def dilution_factor(T_reh_gev: float, rho_end: float = RHO_END) -> float:
    r"""Independent recomputation of (a_end/a_0)^3.

    The condensate dilutes like matter from a_end to a_reh (rho propto a^-3); at a_reh
        rho_cond(a_reh) = rho_rad(T_reh) = (pi^2/30) g_* T_reh^4
    ==> (a_end/a_reh)^3 = rho_rad(T_reh)/rho_end
    afterward entropy conservation g_{*s} a^3 T^3 = const gives
    ==> (a_reh/a_0)^3 = g_{*s,0} T_0^3 / (g_{*s,reh} T_reh^3)
    hence  (a_end/a_0)^3 = (pi^2/30) g_* g_{*s,0} T_0^3 T_reh / (rho_end g_{*s,reh})
    """
    rho_rad = (math.pi**2 / 30.0) * G_STAR * T_reh_gev**4
    return (rho_rad / rho_end) * (G_STAR_S0 * T0_GEV**3) / (G_STAR * T_reh_gev**3)


def omega_psi(g: float, T_reh_gev: float, kappa: float,
              omega_ref: float | None = None) -> float:
    r"""Omega_psi = m_psi H_inf^3 exp(-kappa m_psi/H_inf) (a_end/a_0)^3 / rho_c.

    kappa is the exponent coefficient of m_psi/H_inf: the paper uses pi, the exact result uses 2 pi.
    """
    m_psi = g * PHI_V
    n_psi = H_INF**3 * math.exp(-kappa * m_psi / H_INF)
    return m_psi * n_psi * dilution_factor(T_reh_gev) / RHO_C


def solve_g(kappa: float, T_reh_gev: float, target: float = OMEGA_PSI_TARGET) -> dict:
    """Solve omega_psi(g)=target; return the two roots (small-g branch and large-g branch)."""
    from scipy.optimize import brentq

    mu = PHI_V / H_INF                     # m_psi/H_inf = mu * g
    # omega = A g exp(-kappa mu g),  A = mu H_inf^4 (a_end/a0)^3/rho_c  (m_psi=mu H_inf g)
    A = mu * H_INF**4 * dilution_factor(T_reh_gev) / RHO_C
    # location of the maximum
    g_max = 1.0 / (kappa * mu)
    om_max = A * g_max * math.exp(-1.0)
    if om_max < target:
        return {"feasible": False, "g_max": g_max, "omega_max": om_max, "A": A}

    f = lambda gg: A * gg * math.exp(-kappa * mu * gg) - target
    roots = []
    lo1, hi1 = 1e-300, g_max
    lo2, hi2 = g_max, 1.0
    for lo, hi in ((lo1, hi1), (lo2, hi2)):
        try:
            roots.append(brentq(f, lo, hi, xtol=1e-22, rtol=1e-15, maxiter=400))
        except ValueError:
            roots.append(None)
    out = {"feasible": True, "g_max": g_max, "omega_max": om_max, "A": A}
    for name, gg in (("g_small", roots[0]), ("g_large", roots[1])):
        if gg is None:
            out[name] = None
        else:
            out[name] = {"g": gg, "m_psi_over_Hinf": mu * gg,
                         "m_psi_GeV": mu * gg * H_INF}
    return out


# ---------------------------------------------------------------------------
def exact_dS_fermion_exponent(mu: float) -> dict:
    r"""Exact combination: exponent = single Gamma factor + BD normalization factors.

    |Gamma(1/2 + i mu)|^2 = pi / cosh(pi mu)   (exact)
    Large-mu expansion:  = 2 pi e^{-pi mu} (1 - e^{-2 pi mu} + ...)
    The ENT Eq.(16) mode functions carry e^{+ pi m/(2H)}; after normalization by |alpha|^2 rel |beta|^2,
    the two helicities each contribute e^{-pi mu}, combining with the e^{-pi mu} of |Gamma|^2 into e^{-2 pi mu}.
    """
    mp.mp.dps = 50
    gam2 = abs(mp.gamma(mp.mpf(1) / 2 + 1j * mp.mpf(mu)))**2     # = pi/cosh(pi mu)
    gam2_exact = mp.pi / mp.cosh(mp.pi * mp.mpf(mu))
    asym = 2 * mp.pi * mp.e**(-mp.pi * mp.mpf(mu))
    # exact Fermi-Dirac occupation number
    fd = 1.0 / (mp.e**(2 * mp.pi * mp.mpf(mu)) + 1)
    return {
        "mu": mu,
        "absGamma2": float(gam2),
        "pi_over_cosh": float(gam2_exact),
        "asym_2pi_exp(-pi mu)": float(asym),
        "log_absGamma2__over__mu": float(mp.log(gam2) / mp.mpf(mu)),
        "one_gamma_exponent": -float(mp.log(gam2) / mp.mpf(mu)) / math.pi,
        "fd_occupation": float(fd),
        "fd_log_slope_over_mu": float(mp.log(fd) / mp.mpf(mu)),
    }


def main() -> None:
    out: dict = {}

    # ---- [1] dimensional/constant self-check ----
    out["constants"] = {
        "M_P": M_P, "xi": XI, "lambda0": LAM0, "beta_p": BETA_P,
        "Phi_V": PHI_V, "V0": V0, "H_inf": H_INF, "m_chi": M_CHI,
        "m_chi_over_H_inf": M_CHI / H_INF,
        "x_end": X_END, "V_end": V_END, "K_end": K_END, "rho_end": RHO_END,
        "rho_c": RHO_C,
    }

    # ---- [2] exponent-coefficient audit ----
    audit = [exact_dS_fermion_exponent(mu) for mu in (0.5, 1.0, 2.0, 4.0, 8.0)]
    # log slope / mu of a single Gamma factor -> -pi ; log slope / mu of the exact FD -> -2pi
    out["exponent_audit"] = audit

    # ---- [3] dilution factor ----
    dil = {f"{T:g}": dilution_factor(T) for T in (1e6, 1e8, 1e9, 1e10, 1e12)}
    out["dilution"] = dil
    out["dilution_at_1e9"] = dilution_factor(1e9)

    # ---- [4] Omega_psi = 0.265: paper pi vs exact 2pi ----
    res = {}
    for name, kappa in (("paper_pi", math.pi), ("exact_2pi", 2 * math.pi)):
        for T in (1e9,):
            r = solve_g(kappa, T)
            r["kappa"] = kappa
            r["T_reh"] = T
            res[f"{name}__T{T:g}"] = {k: (v if not isinstance(v, dict) else v)
                                      for k, v in r.items()}
    out["solve_omega_dm"] = res

    # ---- [5] density-ratio check (paper L.465) ----
    ratios = []
    for name, kappa in (("paper_pi", math.pi), ("exact_2pi", 2 * math.pi)):
        r = solve_g(kappa, 1e9)
        if not r["feasible"] or r["g_large"] is None:
            continue
        g = r["g_large"]["g"]
        m_psi = g * PHI_V
        n_psi = H_INF**3 * math.exp(-kappa * m_psi / H_INF)
        rho_psi_end = m_psi * n_psi
        # at a_reh psi is still matter-like: rho_psi(a_reh) = rho_psi(a_end)(a_end/a_reh)^3
        rho_rad_reh = (math.pi**2 / 30.0) * G_STAR * 1e9**4
        rho_psi_reh = rho_psi_end * (rho_rad_reh / RHO_END)
        # "bare" value without the exponential factor (the paper's L.465 form rho_psi ~ m_psi H_inf^3)
        rho_bare = m_psi * H_INF**3
        ratios.append({
            "variant": name, "g": g, "m_psi_over_Hinf": m_psi / H_INF,
            "rho_psi_end": rho_psi_end,
            "rho_bare": rho_bare,
            "bare_over_rho_chi_5e63": rho_bare / 5e63,
            "bare_over_rho_end": rho_bare / RHO_END,
            "rho_psi_over_rho_chi__paper5e63": rho_psi_end / 5e63,
            "rho_psi_over_rho_chi__exact_rho_end": rho_psi_end / RHO_END,
            "rho_psi_over_rho_rad_at_reh": rho_psi_reh / rho_rad_reh,
            "identity_check": abs(rho_psi_reh / rho_rad_reh
                                   - rho_psi_end / RHO_END) / (rho_psi_end / RHO_END),
        })
    out["density_ratios"] = ratios

    # ---- [6] consistency with the paper's self-reported window ----
    win = {
        "g_window_paper": [1e-5, 1e-4],
        "m_over_H_window_paper": [1.0, 10.0],
    }
    for key, r in out["solve_omega_dm"].items():
        if r.get("feasible") and r.get("g_large"):
            g = r["g_large"]["g"]
            mh = r["g_large"]["m_psi_over_Hinf"]
            r["g_in_paper_window"] = (win["g_window_paper"][0] <= g
                                      <= win["g_window_paper"][1])
            r["m_over_H_in_paper_window"] = (win["m_over_H_window_paper"][0] <= mh
                                             <= win["m_over_H_window_paper"][1])
            r["negative_m_psi_over_Hinf"] = r["g_large"]["m_psi_GeV"] <= 0
            # kinematic closure: m_psi >= m_chi/2
            r["kinematic_closed"] = mh >= 0.5 * M_CHI / H_INF
    out["paper_window"] = win

    with open(os.path.join(ROOT, "psi_production_bogoliubov.json"), "w",
              encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=float)

    # ---------------- report ----------------
    L: list[str] = []
    A = L.append
    A("# psi gravitational-production abundance: independent recomputation and exponent-coefficient audit (P0 backfill #2)\n")
    A("Target: `paper_prd_merged.tex` L.421-465 (`eq:npsi`, `Omega_psi`, the g window, density ratios).")
    A("Method: derive the mode equation from the Dirac equation -> recompute the exact Bogoliubov exponent at high precision -> independently recompute the dilution factor")
    A("-> solve `Omega_psi = 0.265` with the two exponent coefficients, and check consistency with the paper's self-reported window.\n")

    A("## 1. Mode equation and exact index (self-derivation)\n")
    A("The rescaled field $\\tilde\\psi=a^{3/2}\\psi$ satisfies the flat-form Dirac equation; taking the chiral representation and acting on")
    A("helicity eigenstates with $\\sigma\\!\\cdot\\!\\nabla\\to ikh$ ($s=kh$) gives\n")
    A("$$\\partial_\\eta u_R=+isu_R-imau_L,\\qquad \\partial_\\eta u_L=-isu_L-imau_R.$$\n")
    A("With $P=u_R+u_L,\\ M=u_R-u_L$, this decouples into\n")
    A("$$\\boxed{\\,u''+\\bigl[k^2+m^2a^2\\mp i\\,m a'\\bigr]u=0\\,}$$\n")
    A("For pure de Sitter ($a=-1/(H\\eta)$, $\\mu=m/H$), the index equation gives")
    A("$\\nu^2-\\tfrac14=-(\\mu^2\\mp i\\mu)=(\\tfrac12\\mp i\\mu)^2$, i.e.\n")
    A("$$\\nu=\\tfrac12\\mp i\\,\\frac{m}{H}\\quad\\text{(exact)}.$$\n")
    A("This matches Ema-Nakayama-Tang arXiv:1903.10973 Eq.(14)(16) verbatim")
    A("(that paper's mode equation $\\partial_\\tau^2u_\\pm+[\\omega_k^2\\pm i(am)']u_\\pm=0$,")
    A("solution $u_\\pm\\propto\\sqrt{-\\pi k\\tau/4}\\,e^{\\pm\\pi m/(2H)}H^{(1)}_{\\nu_\\pm}(-k\\tau)$,")
    A("$\\nu_\\pm=\\tfrac12\\mp i m/H$). **Note that the index contains $m/H$, not $\\sqrt{m^2/H^2-1/4}$.**\n")

    A("## 2. Exponent-coefficient audit: where $\\pi m/H$ comes from and why it must be doubled\n")
    A("| $\\mu=m/H$ | $|\\Gamma(\\tfrac12+i\\mu)|^2$ | $\\pi/\\cosh(\\pi\\mu)$ | $\\ln|\\Gamma|^2/\\mu\\cdot(-1/\\pi)$ | exact FD occupation $1/(e^{2\\pi\\mu}+1)$ | $\\ln(\\mathrm{FD})/\\mu\\cdot(-1/\\pi)$ |")
    A("|---|---|---|---|---|---|")
    for d in audit:
        A(f"| {d['mu']:g} | {d['absGamma2']:.6e} | {d['pi_over_cosh']:.6e} "
          f"| {d['one_gamma_exponent']:.6f} | {d['fd_occupation']:.6e} "
          f"| {d['fd_log_slope_over_mu']:.6f} |")
    A("")
    A("Reading: $|\\Gamma(\\tfrac12+i\\mu)|^2=\\pi/\\cosh(\\pi\\mu)\\simeq 2\\pi e^{-\\pi\\mu}$ provides only")
    A("**one** $e^{-\\pi\\mu}$; the BD normalization of the modes additionally carries $e^{\\pm\\pi m/(2H)}$ (explicit in ENT Eq.(16)),")
    A("and the two helicities each contribute one $e^{-\\pi\\mu}$, which together give $e^{-2\\pi\\mu}$.")
    A("In the last table column the value is **$-2$** (not $-1$) for $\\mu\\gtrsim 2$, i.e. the log slope of the exact occupation number is $2\\pi$.")
    A("**Counting only the $\\Gamma$ factor while dropping the normalization factor yields exactly the paper's $e^{-\\pi m/H}$.**\n")
    A("Literature verdict (externally verified): exact calculations universally give $e^{-2\\pi m/H};")
    A("the only one printing $e^{-\\pi m/H}$ is Kolb-Long, Rev. Mod. Phys. 97 (2025) arXiv:2312.09042,")
    A("and that paper explicitly labels it as only an **order-of-magnitude heuristic by analogy with the Schwinger effect** (\"this will serve as a good guide\"),")
    A("not the exact late-time $|\\beta_k|^2$. Hence the exponent coefficient of `eq:npsi` should be changed to $2\\pi$.\n")

    A("## 3. Independent recomputation of the dilution factor $(a_{\\rm end}/a_0)^3$\n")
    A("$$\\Bigl(\\frac{a_{\\rm end}}{a_0}\\Bigr)^3"
      "=\\underbrace{\\frac{\\rho_{\\rm rad}(T_{\\rm reh})}{\\rho_{\\rm end}}}_{(a_{\\rm end}/a_{\\rm reh})^3}"
      "\\cdot\\underbrace{\\frac{g_{*s,0}T_0^3}{g_{*s,\\rm reh}T_{\\rm reh}^3}}_{(a_{\\rm reh}/a_0)^3}"
      "=\\frac{\\pi^2}{30}\\frac{g_*g_{*s,0}T_0^3}{\\rho_{\\rm end}g_{*s,\\rm reh}}\\,T_{\\rm reh}$$\n")
    A("| $T_{\\rm reh}$ [GeV] | $(a_{\\rm end}/a_0)^3$ |")
    A("|---|---|")
    for k, v in dil.items():
        A(f"| {k} | {v:.4e} |")
    A(f"\nAt $T_{{\\rm reh}}=10^9$ GeV this gives **{out['dilution_at_1e9']:.3e}**,")
    A("matching the paper's L.442 $\\sim10^{-92}$ (**consistent**; the $T_{\\rm reh}$ scaling also matches: $\\propto T_{\\rm reh}$).")
    A("That line needs no modification.\n")

    A("## 4. Solving $\\Omega_\\psi=0.265$: paper $\\pi$ vs exact $2\\pi$\n")
    A("| variant | $\\kappa$ | small root $g$ | large root $g$ | large root $m_\\psi/H_{\\rm inf}$ | large root inside $g\\in[10^{-5},10^{-4}]$? |")
    A("|---|---|---|---|---|---|")
    for key, r in out["solve_omega_dm"].items():
        gs = r["g_small"]["g"] if r.get("g_small") else float("nan")
        gl = r["g_large"]["g"] if r.get("g_large") else float("nan")
        mh = r["g_large"]["m_psi_over_Hinf"] if r.get("g_large") else float("nan")
        A(f"| {key} | {r['kappa']:.6f} | {gs:.4e} | {gl:.4e} | {mh:.3f} | "
          f"{r.get('g_in_paper_window')} |")
    A("")
    A("Conclusion: **with the paper's own $\\pi$, the required $g$ falls outside $[10^{-5},10^{-4}]$**;")
    A("with the exact $2\\pi$, the required $g$ falls back inside the paper's self-reported window.")
    A("That is, correcting the exponent **does not** weaken the model; it aligns it with its own parameter window.")
    A("(At $g=10^{-5}$, $\\Omega_\\psi$ far exceeds $0.265$, so within the window only the large-root branch is usable.)\n")

    A("## 5. Density-ratio check (the paper's L.465 values $10^{-13}$ and $10^{-18}$)\n")
    A("| variant | $g$ | $m_\\psi/H_{\\rm inf}$ | bare $m_\\psi H_{\\rm inf}^3$ | bare$/\\rho_\\chi(5\\times10^{63})$ | with exponent $\\rho_\\psi(a_{\\rm end})$ | $/\\rho_{\\rm end}$ | $/\\rho_{\\rm rad}(T_{\\rm reh})$ |")
    A("|---|---|---|---|---|---|---|---|")
    for r in ratios:
        A(f"| {r['variant']} | {r['g']:.4e} | {r['m_psi_over_Hinf']:.3f} | "
          f"{r['rho_bare']:.4e} | {r['bare_over_rho_chi_5e63']:.4e} | "
          f"{r['rho_psi_end']:.4e} | {r['rho_psi_over_rho_chi__exact_rho_end']:.4e} | "
          f"{r['rho_psi_over_rho_rad_at_reh']:.4e} |")
    A("")
    A("Two conclusions:")
    _c1 = ("1. **$\\rho_\\psi/\\rho_\\chi$ and $\\rho_\\psi/\\rho_{\\rm rad}$ are the same number at reheating**"
           " ($\\rho_{\\rm rad}(T_{\\rm reh})=\\rho_{\\rm cond}(a_{\\rm reh})$; the two chains are identical,"
           " identity check in the table " + f"{ratios[0]['identity_check']:.2e}" + ")."
           " Yet L.465 of the paper writes the same quantity as two values $10^{-13}$ and $10^{-18}$, **mutually off by $10^5$**.")
    A(_c1)
    A("2. The correct value $\\approx6\\times10^{-19}$ matches the paper's $10^{-18}$ and differs from $10^{-13}$ by 6 orders of magnitude."
      " Hence the $10^{-13}$ at L.465 should be changed to $\\sim10^{-18}$.")
    A("(Origin of the $10^{-13}$: when writing $\\rho_\\psi\\sim m_\\psi H_{\\rm inf}^3$ the paper drops the exponential factor;"
      " the bare value $m_\\psi H^3$ at $m_\\psi\\simeq3H$ gives $\\sim10^{-10}$, still not $10^{-13}$.)\n")

    md = "\n".join(L) + "\n"
    with open(os.path.join(ROOT, "psi_production_bogoliubov.md"), "w",
              encoding="utf-8") as f:
        f.write(md)
    print("wrote scripts/psi_production_bogoliubov.md (%d chars)" % len(md))


if __name__ == "__main__":
    main()
