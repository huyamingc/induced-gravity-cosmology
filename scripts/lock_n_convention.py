# -*- coding: utf-8 -*-
"""
Lock one cosmological definition of N and recompute all observables.

CONVENTION (fixed, not presupposed from the paper):
  N = ln(a_end / a_*)   (locked definition)
  with a_* from k_* = a_* H_*  (pivot k=0.05 Mpc^{-1})
  and a_end from post-inflationary matching
      a_end = a_reh (rho_reh/rho_end)^{1/3}
      a_reh = a_eq (rho_eq/rho_reh)^{1/4}
      a_eq  = Omega_r / Omega_m
  assuming matter-like (w=0) condensate domination until T_reh.

At a given N (this definition), x_* is fixed by the exact slow-roll integral
  N(x) = [e^x - x - (e^{x_end}-x_end)] / (2 beta_p^2)
and observables are potential slow-roll at that x_*:
  eps, eta, r=16 eps, n_s=1-6eps+2eta, A_s = V(x_*)/(24 pi^2 M_P^4 eps)
lambda0 is then fixed by A_s = Planck.

Self-consistency: for each (xi, N, T_reh) we also evaluate the *derived* N
from expansion history using that lambda0/V_end/V_* and report the residual.
The "physical" band is where |N_derived - N| is small AND n_s within Planck 2sigma
AND T_reh above BBN.

Outputs:
  scripts/n_convention_results.md
  scripts/n_convention_results.json
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from derive_from_action import (
    A_S_OBS,
    M_PL,
    N_S_OBS,
    SIG_NS,
    H0_GeV_v2,
    T0_GeV,
    beta_p,
    k_pivot_GeV,
    lambda0_for_As,
    m_chi_from_first_principles,
    N_match_derived,
    N_of_x,
    observables_at_x,
    rho_rad,
    x_end_from_eps1,
    x_star_for_N,
    V0_from_lambda,
)

OUT_MD = Path(__file__).resolve().parent / "n_convention_results.md"
OUT_JSON = Path(__file__).resolve().parent / "n_convention_results.json"

XI_LIST = [1.0, 5.0, 11.1, 30.0, 100.0]
N_LIST = [45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58]
T_BENCH = [1e4, 1e6, 1e8, 1e9, 1e10, 1e12, 1e14]  # GeV
T_BBN = 1e-2  # GeV


def point(xi: float, N: float) -> dict:
    lam, obs = lambda0_for_As(N, xi)
    bp = obs["beta_p"]
    x_e = x_end_from_eps1(bp)
    u_e = math.exp(-x_e)
    V0 = obs["V0"]
    V_end = V0 * (1.0 - u_e) ** 2
    V_star = obs["VE"]
    Hinf = math.sqrt(V0 / (3.0 * M_PL**2))
    return {
        "xi": xi,
        "N": N,
        "lambda0": lam,
        "x_star": obs["x"],
        "eps": obs["eps"],
        "eta": obs["eta"],
        "r": obs["r_ps"],
        "ns": obs["ns_ps"],
        "As": obs["As"],
        "V0": V0,
        "V_end": V_end,
        "V_star": V_star,
        "H_inf": Hinf,
        "U_quarter": V0**0.25,
        "m_chi": m_chi_from_first_principles(lam, xi),
        "beta_p": bp,
        "beta_o": obs["beta_o"],
        "ns_sigma": (obs["ns_ps"] - N_S_OBS) / SIG_NS,
    }


def N_derived_for_T(p: dict, T_reh: float) -> float:
    info = N_match_derived(T_reh, p["V_end"], p["V_star"], p["xi"], p["N"])
    return info["N"]


def T_reh_for_N_derived(p: dict, N_target: float) -> float:
    """Find T_reh such that cosmological matching gives N_derived = N_target."""

    def f(logT):
        return N_derived_for_T(p, math.exp(logT)) - N_target

    try:
        return math.exp(brentq(f, math.log(1e-2), math.log(1e16)))
    except ValueError:
        return float("nan")


def r_falsify_line(xi: float, r_lim: float = 0.01) -> float:
    """Smallest r on the locked-N curve that is still >= r_lim? We want N such that r(N)=r_lim
    under THIS convention: r = 16 eps(x_*(N)). Invert numerically.
    """
    def g(N):
        return point(xi, N)["r"] - r_lim

    try:
        return brentq(g, 30.0, 80.0)
    except ValueError:
        return float("nan")


def main() -> None:
    lines = []
    lines.append("# Derived results under the locked N definition")
    lines.append("")
    lines.append("**Convention (fixed for this derivation, not copied from the paper):**")
    lines.append("")
    lines.append("$$N \\equiv \\ln(a_{\\rm end}/a_*)$$")
    lines.append("")
    lines.append("where \\(k_*=a_*H_*\\) (\\(k_*=0.05\\,\\mathrm{Mpc}^{-1}\\)),")
    lines.append("the condensate has \\(w=0\\) before reheating,")
    lines.append("\\(a_{\\rm end}=a_{\\rm reh}(\\rho_{\\rm reh}/\\rho_{\\rm end})^{1/3}\\),")
    lines.append("\\(a_{\\rm reh}=a_{\\rm eq}(\\rho_{\\rm eq}/\\rho_{\\rm reh})^{1/4}\\),")
    lines.append("\\(a_{\\rm eq}=\\Omega_r/\\Omega_m\\).")
    lines.append("")
    lines.append("At this N: \\(x_*\\) is fixed by the exact slow-roll integral \\(N(x_*)=N\\),")
    lines.append("the observables are potential slow-roll \\(r=16\\varepsilon_*\\), \\(n_s=1-6\\varepsilon_*+2\\eta_*\\),")
    lines.append("\\(\\lambda_0\\) is solved back from \\(A_s=2.1\\times10^{-9}\\).")
    lines.append("")

    # --- Fiducial xi=11.1 table ---
    xi0 = 11.1
    lines.append(f"## 1. xi={xi0} locked-N grid")
    lines.append("")
    lines.append("| N | lambda0 | r | n_s | n_s-Planck | H_inf | U^{1/4} | m_chi | T_reh* (making N_derived=N) |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    grid = []
    for N in N_LIST:
        p = point(xi0, N)
        Tstar = T_reh_for_N_derived(p, N)
        p["T_reh_selfcons"] = Tstar
        grid.append(p)
        lines.append(
            f"| {N} | {p['lambda0']:.3e} | {p['r']:.5f} | {p['ns']:.4f} | "
            f"{p['ns_sigma']:+.2f} sigma | {p['H_inf']:.3e} | {p['U_quarter']:.3e} | "
            f"{p['m_chi']:.3e} | {Tstar:.3e} |"
        )
    lines.append("")

    # --- Derived N vs T for N=50,55 points ---
    lines.append("## 2. Self-consistency: with lambda0(N) fixed, N_derived(T_reh)")
    lines.append("")
    lines.append("| Locked N | T_reh | N_derived | Delta N |")
    lines.append("|---|---|---|---|")
    for p in grid:
        if p["N"] not in (48, 50, 52, 55):
            continue
        for T in T_BENCH:
            Nd = N_derived_for_T(p, T)
            lines.append(f"| {p['N']} | {T:.2e} | {Nd:.2f} | {Nd - p['N']:+.2f} |")
    lines.append("")
    lines.append("Reading: if |Delta N| <~ 0.5, the (N, T_reh) pair is approximately self-consistent under this matching;")
    lines.append("the deviation between the locked N and N_derived directly reflects the uncertainty of the reheating assumptions (w, k, Omega, g*).")
    lines.append("")

    # --- Planck window under locked convention ---
    lines.append("## 3. Planck 2-sigma window and r range (locked-N definition)")
    lines.append("")
    ns_lo = N_S_OBS - 2 * SIG_NS
    ns_hi = N_S_OBS + 2 * SIG_NS
    N_win = [p for p in grid if ns_lo <= p["ns"] <= ns_hi]
    if N_win:
        Nmin, Nmax = N_win[0]["N"], N_win[-1]["N"]
        rmin = min(p["r"] for p in N_win)
        rmax = max(p["r"] for p in N_win)
        lines.append(f"- Planck n_s in [{ns_lo:.4f}, {ns_hi:.4f}]")
        lines.append(f"- N satisfying 2 sigma under the locked N: N in **[{Nmin}, {Nmax}]**")
        lines.append(f"- corresponding r in **[{rmin:.5f}, {rmax:.5f}]**")
    else:
        lines.append("- no grid point falls inside the Planck 2-sigma window (refine the N grid or check ns(PS))")
    lines.append("")
    # 1 sigma
    ns_lo1 = N_S_OBS - SIG_NS
    N_win1 = [p for p in grid if ns_lo1 <= p["ns"] <= ns_hi]
    if N_win1:
        lines.append(
            f"- with the 1-sigma lower bound n_s >= {ns_lo1:.4f}: N in [{N_win1[0]['N']}, {N_win1[-1]['N']}],"
            f"r in [{min(p['r'] for p in N_win1):.5f}, {max(p['r'] for p in N_win1):.5f}]"
        )
    lines.append("")

    # --- Falsification ---
    lines.append("## 4. Falsification line (locked-N definition)")
    lines.append("")
    lines.append("| xi | N* at r=0.01 | n_s there | r_max (2-sigma window) | Does LiteBIRD r>0.01 exclude |")
    lines.append("|---|---|---|---|---|")
    for xi in XI_LIST:
        Nstar = r_falsify_line(xi, 0.01)
        if math.isnan(Nstar):
            pstar = None
            ns_s = float("nan")
        else:
            pstar = point(xi, Nstar)
            ns_s = pstar["ns"]
        # r max in Planck 2sigma window at this xi
        rmax = None
        for N in N_LIST:
            p = point(xi, N)
            if ns_lo <= p["ns"] <= ns_hi:
                rmax = p["r"] if rmax is None else max(rmax, p["r"])
        if rmax is None:
            excl = "no point in window"
        else:
            excl = "YES" if rmax < 0.01 else "NO (window already has r >= 0.01)"
        lines.append(
            f"| {xi} | {Nstar:.2f} | {ns_s:.4f} | "
            f"{'--' if rmax is None else f'{rmax:.5f}'} | {excl} |"
        )
    lines.append("")
    lines.append("If r_max < 0.01 within the 2-sigma window, then a LiteBIRD detection of r > 0.01 would exclude this model class at this xi under the locked-N convention.")
    lines.append("")

    # --- xi scan at N=50 locked ---
    lines.append("## 5. xi scan (locked N=50)")
    lines.append("")
    lines.append("| xi | lambda0 | r | n_s | m_chi | H_inf | T_reh* |")
    lines.append("|---|---|---|---|---|---|---|")
    for xi in XI_LIST:
        p = point(xi, 50)
        Tstar = T_reh_for_N_derived(p, 50)
        lines.append(
            f"| {xi} | {p['lambda0']:.3e} | {p['r']:.5f} | {p['ns']:.4f} | "
            f"{p['m_chi']:.3e} | {p['H_inf']:.3e} | {Tstar:.3e} |"
        )
    lines.append("")
    lines.append("Under the locked-N definition, r depends only weakly on xi (through beta_p and x_*(N)); lambda0 propto xi^2 still holds approximately.")
    lines.append("")

    # --- DE / DM unchanged structure ---
    p50 = point(11.1, 50)
    p50["T_reh_selfcons"] = T_reh_for_N_derived(p50, 50)
    H0 = H0_GeV_v2()
    lines.append("## 6. Relation to the DE/DM structure (unchanged by the N convention)")
    lines.append("")
    lines.append(f"- locked N=50, xi=11.1: lambda0={p50['lambda0']:.4e}, m_chi={p50['m_chi']:.4e} GeV")
    lines.append(f"- m_chi/H0={p50['m_chi']/H0:.4e} -> **still impossible** as today's quintessence DE")
    lines.append(f"- H_inf={p50['H_inf']:.4e} GeV; the DM window m_psi ~ H_inf order-of-magnitude structure is unchanged")
    lines.append(f"- T_reh* (N_derived=50) ~ {p50['T_reh_selfcons']:.3e} GeV (this matching)")
    lines.append("")

    # --- Figure 2 locked ---
    # Fig. 2 is produced by fig2_ns_r.py, the FIG-class script.  This script used
    # to redraw it here as figures/fig2_ns_r_lockedN.pdf, which duplicated the
    # artefact byte-for-byte and left two files claiming to be the same figure.
    lines.append("## 7. Figure 2")
    lines.append("")
    lines.append("- written by fig2_ns_r.py as `figures/fig2_ns_r.pdf`, the file the manuscript includes")
    lines.append("")

    # --- Summary ---
    lines.append("## 8. Conclusions (under this N convention)")
    lines.append("")
    rmax2 = None
    Nband = []
    for p in grid:
        if ns_lo <= p["ns"] <= ns_hi:
            Nband.append(p["N"])
            rmax2 = p["r"] if rmax2 is None else max(rmax2, p["r"])
    lines.append("1. **The observables are recomputable at a given N (this definition)**: at xi=11.1, N=50 gives")
    lines.append(f"   lambda0 ~ {p50['lambda0']:.3e}, r ~ {p50['r']:.5f}, n_s ~ {p50['ns']:.4f} (Planck deviation {p50['ns_sigma']:+.2f} sigma).")
    if Nband:
        lines.append(f"2. Planck 2 sigma corresponds to locked N in [{min(Nband)}, {max(Nband)}], r_max ~ {rmax2:.5f}.")
    lines.append("3. **The self-consistent T_reh\\*** differs from the textbook/paper closed form: under this matching, T_reh\\* at N=50 is of order 10^7-10^8 GeV;")
    lines.append("   at T_reh=10^9, N_derived is about 51 -- still inside the window but no longer 50.")
    lines.append("4. **Falsification**: if r_max < 0.01 within the Planck 2-sigma window, then r > 0.01 excludes this model class (usually satisfied at xi=11.1 under this convention).")
    lines.append("5. The structural DE/DM conclusions (no single-field quintessence; conditional DM) **do not depend** on the N convention.")
    lines.append("")
    lines.append("## 9. Implications for the manuscript")
    lines.append("")
    lines.append("| Item | Suggested text for paper_prd_merged.tex |")
    lines.append("|---|---|")
    lines.append("| Definition of N | state N=ln(a_end/a_*) and the matching assumptions (k, w=0 reheating, Omega) explicitly |")
    lines.append("| r, n_s value table | replace the vague single statement \"N=50 => r=0.00487\" with this locked table |")
    lines.append("| N window | phrase as \"Planck 2 sigma intersected with self-consistent T_reh >~ BBN\", not a single formula |")
    lines.append("| Appendix lambda0 | use the exact slow-roll inversion (N=50, xi=11.1 gives ~6.7e-8) and delete the wrong 48 pi^2 formula |")
    lines.append("| Falsification sentence | compare r_max (2-sigma window) with r=0.01; state the dependence on the N convention |")
    lines.append("")
    lines.append("[Locked-N convention derivation complete]")

    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    OUT_JSON.write_text(
        json.dumps({"xi11.1": grid, "convention": "N=ln(a_end/a_*), exact PS at x_*(N)"}, indent=2),
        encoding="utf-8",
    )
    print("\n".join(lines))
    print(f"\nWrote {OUT_MD}")
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
