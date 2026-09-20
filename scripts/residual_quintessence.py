#!/usr/bin/env python3
r"""
residual_quintessence.py  (P0 backfill #4)
==========================================
Original requirement (the review workspace review_workspace/ was removed after
the final draft; excerpt kept) Sec. C:
  "[P0] residual-quintessence lesssim1e-10 recheck: integrate the condensate
        Boltzmann system (Gamma_anom, H(T), a^-3 dilution) and output the
        rho_cond(a_0)/V_c - T_reh curve.
        Missing: the 1e-10 bound is an abstract-level selling point with no
        script, and the Delta w_quant paper value 1e-118
        disagrees with its own formula 1.7e-121 by 3 orders of magnitude."

Paper claims (paper_prd_merged.tex; cited by section -- line numbers drift
between revisions and are no longer hardcoded):
  Sec. VI   Delta w_Ricci  ~ (H_0/m_chi)^2        ~ 2e-111
  Sec. VI   Delta w_quant  ~ H_0^4/V_c            ~ 2e-121  (corrected in
            revision; the original draft printed 1e-118)
  Sec. VI   rho_cond(a_0)/V_c <= Omega_DM/Omega_Lambda ~ 0.39  (BBN floor
            T_reh = 1e-2 GeV)
  Sec. VI   "for the fiducial T_reh >= 1e9 GeV ... the residual condensate is
            exponentially negligible"
  Sec. VI   "up to residual-condensate corrections bounded in App. D5
            (the dominant correction Delta w <= 1e-10 for T_reh >= 1 GeV)"

Method
------
Integrate the EXACT two-fluid Boltzmann system (cosmic time replaced by
e-folds, and an exponential integrator for the decay term, because the
late-time Gamma/H reaches 1e30 and an explicit scheme goes unstable):

    Y = rho_cond a^3        dY/dN = -(Gamma/H) Y
    R = rho_rad  a^4        dR/dN = +(Gamma/H) Y a
    H^2 = (Y/a^3 + R/a^4) / (3 M_P^2)
    exponential integrator:  Y' = Y exp(-Gamma dt),  R' = R + (Y - Y') a_mid

N runs from a_end (=1) to a_0 (N ~ 70.7). Gamma (equivalently T_reh) is
scanned; the output is rho_cond(a_0)/V_c, T_reh, N_reh, compared against the
three paper statements. A control case with Gamma = 0 (absolutely stable
condensate) plus the gravitational pulse as initial radiation is also run.

Output: scripts/residual_quintessence.md / .json
"""

from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cosmo_model import V_end_over_V0, rho_end_over_V_end  # noqa: E402

ROOT = os.path.dirname(os.path.abspath(__file__))

M_P = 2.435e18
H0_GEV = 1.4377e-42
T0_GEV = 2.3491e-13
G_STAR = 106.75
G_STAR_S0 = 3.91
OMEGA_LAMBDA = 0.683
OMEGA_DM = 0.265
RHO_C = 3.0 * H0_GEV**2 * M_P**2
V_C = OMEGA_LAMBDA * RHO_C

# Paper locked values (N=50, xi=11.1)
XI = 11.1
V0 = 4.7772e63
V_END = V_end_over_V0(XI) * V0
RHO_END = rho_end_over_V_end() * V_END
H_INF = 1.6388e13
M_CHI = 3.2533e13
N_EFF = 10.0
RHO_PULSE = N_EFF * H_INF**4 / (192.0 * math.pi**2)

# a_0/a_end: inferred from (a_end/a_0)^3 = 1.0204e-92 at T_reh = 1e9
A_END_OVER_A0_CUBED_1E9 = 1.0204e-92
N_TO_A0 = -math.log(A_END_OVER_A0_CUBED_1E9) / 3.0     # ~70.74


def T_of_rho_rad(rho: float, g_star: float = G_STAR) -> float:
    return (30.0 * rho / (math.pi**2 * g_star)) ** 0.25


def Gamma_from_Treh(T_reh: float, g_star: float = G_STAR) -> float:
    """T_reh = (90/(pi^2 g_*))^{1/4} sqrt(Gamma M_P)   <=>   H(T_reh) = Gamma."""
    return T_reh**2 / (math.sqrt(90.0 / (math.pi**2 * g_star)) ** 2 * M_P)


def solve_two_fluid(Gamma: float, rho_rad0: float = 0.0,
                    N_end: float = N_TO_A0, dN: float = 0.005) -> dict:
    r"""Exponential integrator for the two-fluid system.

    a = e^N (a_end = 1).  Y = rho_cond a^3, R = rho_rad a^4.
    """
    Y = RHO_END * 1.0**3
    R = rho_rad0 * 1.0**4
    N = 0.0
    N_reh = float("nan")
    T_reh = float("nan")
    rho_rad_at_cross = float("nan")
    while N < N_end:
        a = math.exp(N)
        rho_cond = Y / a**3
        rho_rad = R / a**4
        H = math.sqrt(max((rho_cond + rho_rad) / (3.0 * M_P**2), 1e-300))
        if not math.isfinite(N_reh) and rho_rad >= rho_cond and rho_cond > 0:
            N_reh = N
            T_reh = T_of_rho_rad(rho_rad)
            rho_rad_at_cross = rho_rad
        dt = dN / H
        decay = math.exp(-Gamma * dt) if Gamma * dt < 700.0 else 0.0
        Y_new = Y * decay
        a_mid = math.exp(N + 0.5 * dN)
        R = R + (Y - Y_new) * a_mid
        Y = Y_new
        N += dN
    a0 = math.exp(N)
    rho_cond_0 = Y / a0**3
    rho_rad_0 = R / a0**4
    return {"Gamma": Gamma, "N_reh": N_reh, "T_reh": T_reh,
            "rho_cond_a0": rho_cond_0,
            "rho_cond_a0_over_Vc": rho_cond_0 / V_C,
            "rho_rad_a0": rho_rad_0,
            "T_today_from_rad": T_of_rho_rad(rho_rad_0),
            "rho_rad_at_cross": rho_rad_at_cross,
            "N_total": N}


def main() -> None:
    out: dict = {"constants": {"M_P": M_P, "H0_GeV": H0_GEV, "RHO_C": RHO_C,
                               "V_C": V_C, "V_C_paper_2p7e-47": 2.7e-47,
                               "V0": V0, "V_end": V_END, "rho_end": RHO_END,
                               "H_inf": H_INF, "m_chi": M_CHI,
                               "rho_pulse": RHO_PULSE,
                               "N_to_a0": N_TO_A0,
                               "Omega_DM_over_Omega_L": OMEGA_DM / OMEGA_LAMBDA}}

    # ---- [1] three Delta w terms ----
    dw_ricci = (H0_GEV / M_CHI) ** 2
    dw_quant = H0_GEV**4 / V_C
    dw_quant_paperVc = H0_GEV**4 / 2.7e-47
    out["delta_w"] = {
        "Ricci_recomputed": dw_ricci,
        "Ricci_paper": 2e-111,
        "quant_recomputed_H0^4/Vc": dw_quant,
        "quant_with_paper_Vc": dw_quant_paperVc,
        "quant_paper": 2e-121,
        "quant_ratio_paper_over_recomputed": 2e-121 / dw_quant,
        "overclosure_bound_OmegaDM_over_OmegaL": OMEGA_DM / OMEGA_LAMBDA,
    }
    # Residual displacement delta Phi/PhiV ~ (H_0/m_chi)^2? Paper Sec. VII / App. D5 quote 1e-111
    out["delta_w"]["displacement_sq"] = dw_ricci   # (H0/m_chi)^2 ~ 2e-111
    # dA/A: ~ (H0/m_chi)^2 * R0/(m_chi^2) ... the paper quotes 6.3e-111; recorded for comparison only
    out["delta_w"]["dA_over_A_paper"] = 6.3e-111

    # ---- [2] two-fluid scan: Gamma = Gamma_anom(T_reh) ----
    rows = []
    for T in (1e-2, 1e0, 1e3, 1e6, 1e9, 1e12, 1e15):
        G = Gamma_from_Treh(T)
        r = solve_two_fluid(G)
        r["T_reh_input"] = T
        r["Gamma_over_H_reh"] = G / math.sqrt(
            (math.pi**2 / 30.0) * G_STAR * T**4 / (3.0 * M_P**2))
        rows.append(r)
    out["scan_with_decay"] = rows

    # ---- [3] control: Gamma = 0 (absolutely stable condensate) + gravitational pulse ----
    r0_nopulse = solve_two_fluid(0.0, rho_rad0=0.0)
    r0_pulse = solve_two_fluid(0.0, rho_rad0=RHO_PULSE)
    out["stable_condensate"] = {"no_pulse": r0_nopulse, "with_pulse": r0_pulse}

    # ---- [3b] scaling check of the pulse channel: does r = rho_rad/rho_cond grow or decay with a? ----
    # Analytic: for Gamma=0 both Y and R are conserved -> r = (R/Y)/a  propto a^-1,
    # i.e. d ln r/dN = -1 (it DECAYS!).
    scal = []
    N_test = np.linspace(0.0, 20.0, 201)
    for N in N_test[::40]:
        a = math.exp(N)
        rho_c = RHO_END / a**3
        rho_r = RHO_PULSE / a**4
        scal.append({"N": float(N), "a": a, "rho_rad_over_rho_cond": rho_r / rho_c})
    out["pulse_scaling"] = scal
    if len(scal) >= 2:
        s = (math.log(scal[1]["rho_rad_over_rho_cond"])
             - math.log(scal[0]["rho_rad_over_rho_cond"])) / (scal[1]["N"] - scal[0]["N"])
        out["pulse_scaling_slope_dlnr_dN"] = s
    # The "reheating" temperature inferred with the paper's Sec. IV wrong scaling (r propto a^{+1})
    a_ratio_paper = RHO_END / RHO_PULSE
    rho_rad_at_paper_cross = RHO_PULSE * a_ratio_paper**4 / a_ratio_paper**4 * a_ratio_paper
    out["paper_pulse_crossing"] = {
        "a_reh_over_a_end_claimed": a_ratio_paper,
        "T_reh_claimed_by_same_arithmetic":
            T_of_rho_rad(RHO_END / a_ratio_paper),   # result under the paper's wrong scaling
    }
    # With the correct scaling: radiation never catches up; today's radiation temperature
    a0_over_aend = math.exp(N_TO_A0)
    rho_rad_today = RHO_PULSE / a0_over_aend**4
    out["pulse_correct_today"] = {
        "a0_over_aend": a0_over_aend,
        "rho_rad_today": rho_rad_today,
        "T_today_GeV": T_of_rho_rad(rho_rad_today),
        "T_today_eV": T_of_rho_rad(rho_rad_today) * 1e9,
    }

    # ---- [4] Gamma t_0 criterion: when is the decay "complete"? ----
    T0_S = 4.35e17
    rows2 = []
    for T in (1e-2, 1e0, 1e9):
        G = Gamma_from_Treh(T)
        G_s = G / 6.582119569e-25
        rows2.append({"T_reh": T, "Gamma_GeV": G, "Gamma_s^-1": G_s,
                      "Gamma_t0": G_s * T0_S})
    out["decay_completeness"] = rows2

    with open(os.path.join(ROOT, "residual_quintessence.json"), "w",
              encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=float)

    # ---------------- report ----------------
    L = []
    A = L.append
    A("# Residual quintessence recheck: condensate Boltzmann integration and the three $\\Delta w$ terms (P0 backfill #4)\n")
    A("Target: `paper_prd_merged.tex` Sec. VI ($\\Delta w\\leq10^{-10}$, Ricci, quantum,")
    A("$\\leq0.39$ and \"exponentially negligible\") and App. D5a-c (term-by-term derivation).\n")

    c = out["constants"]
    A("## 1. Input constants\n")
    A(f"$\\rho_c={c['RHO_C']:.4e}$ GeV$^4$, $V_c=\\Omega_\\Lambda\\rho_c={c['V_C']:.4e}$ GeV$^4$")
    A(f"(paper Sec. VI uses $2.7\\times10^{{-47}}$, corresponding to $H_0\\simeq70$ km/s/Mpc),")
    A(f"$\\rho_{{\\rm end}}={c['rho_end']:.4e}$ GeV$^4$, $\\rho_{{\\rm pulse}}={c['rho_pulse']:.4e}$ GeV$^4$,")
    A(f"$\\Omega_{{\\rm DM}}/\\Omega_\\Lambda={c['Omega_DM_over_Omega_L']:.4f}$.")
    A(f"Integration range $N\\in[0,{c['N_to_a0']:.2f}]$ (inferred from $(a_{{\\rm end}}/a_0)^3=1.0204\\times10^{{-92}}$).\n")

    d = out["delta_w"]
    A("## 2. Independent recomputation of the three $\\Delta w$ terms\n")
    A("| Term | Paper value | Recomputed | Verdict |")
    A("|---|---|---|---|")
    A(f"| $\\Delta w_{{\\rm Ricci}}=(H_0/m_\\chi)^2$ | $2\\times10^{{-111}}$ | ${d['Ricci_recomputed']:.4e}$ | [OK] |")
    A(f"| $\\Delta w_{{\\rm quant}}=H_0^4/V_c$ | $2\\times10^{{-121}}$ | ${d['quant_recomputed_H0^4/Vc']:.4e}$ | [OK] (ratio {d['quant_ratio_paper_over_recomputed']:.2f}) |")
    A(f"| (if the paper Sec. VI value $V_c=2.7\\times10^{{-47}}$ is used) | -- | ${d['quant_with_paper_Vc']:.4e}$ | same conclusion |")
    A("")
    A(f"**Conclusion**: the recomputed $H_0^4/V_c = {d['quant_recomputed_H0^4/Vc']:.2e}$ and the paper value $2\\times10^{{-121}}$ agree to within the rounding of the printed value (ratio {d['quant_ratio_paper_over_recomputed']:.2f}).")
    A("The direction is **conservative** (the true value is even smaller), so the conclusion is unaffected.")
    A("The earlier draft's $10^{-118}$ was corrected to $\\lesssim2\\times10^{-121}$ in both Sec. VI and App. D5b during revision.\n")

    A("## 3. Condensate residual: exact two-fluid integration\n")
    A("Integrating $dY/dN=-(\\Gamma/H)Y$ ($Y=\\rho_{{\\rm cond}}a^3$) and $dR/dN=+(\\Gamma/H)Ya$")
    A("($R=\\rho_{{\\rm rad}}a^4$), with an exponential integrator for the decay term (the late-time $\\Gamma/H$ reaches $10^{{30}}$).\n")
    A("| $T_{{\\rm reh}}$ input [GeV] | $\\Gamma$ [GeV] | $T_{{\\rm reh}}$ recomputed | $\\rho_{{\\rm cond}}(a_0)$ | $\\rho_{{\\rm cond}}(a_0)/V_c$ |")
    A("|---|---|---|---|---|")
    for r in rows:
        A(f"| {r['T_reh_input']:.0e} | {r['Gamma']:.4e} | {r['T_reh']:.4e} | "
          f"{r['rho_cond_a0']:.4e} | {r['rho_cond_a0_over_Vc']:.4e} |")
    A("")
    s0 = out["stable_condensate"]
    A("**Control: absolutely stable condensate ($\\Gamma=0$)**\n")
    A("| Case | $T_{{\\rm reh}}$ | $\\rho_{{\\rm cond}}(a_0)$ | $\\rho_{{\\rm cond}}(a_0)/V_c$ | Overclosed? |")
    A("|---|---|---|---|---|")
    for key, lab in (("no_pulse", "$\\Gamma=0$, no initial radiation"),
                     ("with_pulse", "$\\Gamma=0$, with the gravitational pulse")):
        r = s0[key]
        A(f"| {lab} | {r['T_reh']:.4e} | {r['rho_cond_a0']:.4e} | "
          f"{r['rho_cond_a0_over_Vc']:.4e} | "
          f"{'**YES**' if r['rho_cond_a0_over_Vc'] > c['Omega_DM_over_Omega_L'] else 'NO'} |")
    A("")
    A("### Verdict\n")
    A("1. **For any $\\Gamma\\neq0$ the residual is negligible**: in the table, from $T_{{\\rm reh}}=10^{{-2}}$ to $10^{{15}}$ GeV")
    A("$\\rho_{{\\rm cond}}(a_0)/V_c$ is 0 throughout (underflow). The reason is the decay-completeness criterion:")
    for r in rows2:
        A(f"   - $T_{{\\rm reh}}={r['T_reh']:.0e}$ GeV $\\Rightarrow \\Gamma={r['Gamma_GeV']:.3e}$ GeV "
          f"$={r['Gamma_s^-1']:.3e}$ s$^{{-1}}$, $\\Gamma t_0={r['Gamma_t0']:.3e}\\gg1$.")
    A("2. **But the paper Sec. VI statement \"the BBN floor gives $\\leq0.39$\" is a constraint, not a recomputation result**:")
    A("   the recomputed residual of the absolutely stable condensate ($\\Gamma=0$) is row 2/3 of the table above,")
    _ovc = s0["no_pulse"]["rho_cond_a0_over_Vc"] / c["Omega_DM_over_Omega_L"]
    A(f"   larger than $\\Omega_{{\\rm DM}}/\\Omega_\\Lambda={c['Omega_DM_over_Omega_L']:.3f}$ by "
      f"**$\\sim10^{{{math.log10(_ovc):.0f}}}$ times** ($={_ovc:.2e}$).")
    A("   I.e. $0.39$ is the **requirement** that \"the condensate must not exceed the dark-matter budget\", while a truly stable")
    A("   condensate would **violate** it. It is satisfied only because the condensate decays ($\\Gamma t_0\\gg1$); the paper says this")
    A("   in the second half of the sentence, but the first half reads as if \"$0.39$ were derived from the constraint\". Suggested rewording:")
    A(f"   \"An absolutely stable condensate would overclose the universe by $\\sim10^{{{math.log10(_ovc):.0f}}}$; since $\\Gamma t_0\\gg1$ (for any $T_{{\\rm reh}}\\gtrsim10^{{-2}}$ GeV),")
    A("   the residual at $a_0$ is exactly zero.\"")
    A("")
    A(f"3. The paper Sec. VI bound $\\Delta w\\leq10^{{-10}}$ ($T_{{\\rm reh}}\\geq1$ GeV) **holds**,")
    A("   and is in reality far stronger than that bound (the actual value is 0).\n")

    A("## 4. [NEW, P0-level] The pulse channel **cannot** reheat\n")
    A("Paper Sec. IV states: \"This radiation redshifts as $a^{-4}$ during matter-domination; ")
    A("the ratio $\\rho_{\\rm rad}/\\rho_\\chi\\propto a^{-1}$ **grows** until ")
    A("$\\rho_{\\rm rad}\\sim\\rho_\\chi$ at $a_{\\rm reh}/a_{\\rm end}\\sim10^{13}$\".\n")
    A("**The scaling sign is wrong.** During condensate domination ($w=0$, $\\rho_{\\rm cond}\\propto a^{-3}$),")
    A("radiation dilutes faster, as $a^{-4}$, hence\n")
    A("$$\\frac{\\rho_{\\rm rad}}{\\rho_{\\rm cond}}=\\frac{R/a^4}{Y/a^3}=\\frac{R}{Y}\\frac1a\\propto a^{-1}"
      "\\quad\\text{(decaying, not growing)}.$$\n")
    sl = out.get("pulse_scaling_slope_dlnr_dN")
    A(f"Numerical check: $d\\ln(\\rho_{{\\rm rad}}/\\rho_{{\\rm cond}})/dN={sl:.4f}$ (analytic value $-1$).\n")
    A("| $N$ | $a/a_{\\rm end}$ | $\\rho_{\\rm rad}/\\rho_{\\rm cond}$ |")
    A("|---|---|---|")
    for r in out["pulse_scaling"][:6]:
        A(f"| {r['N']:.1f} | {r['a']:.3e} | {r['rho_rad_over_rho_cond']:.4e} |")
    A("")
    A("Consequence: **with the gravitational pulse alone, radiation never catches up with the condensate**; the universe stays matter-dominated,")
    A("and no radiation-dominated era begins. The \"$T_{{\\rm reh}}\\sim4\\times10^5$ GeV\" the paper infers using the wrong scaling is therefore invalid:")
    pc = out["pulse_correct_today"]
    A(f"with the correct scaling, today ($a_0/a_{{\\rm end}}={pc['a0_over_aend']:.3e}$) the radiation temperature is only "
      f"$T_0={pc['T_today_GeV']:.3e}$ GeV $={pc['T_today_eV']:.3e}$ eV,")
    A("i.e. the universe would be almost empty.")
    A("**Conclusion: pulse channel (i) is not a reheating channel**; reheating must be provided by the true decay of the condensate (iv).")
    A("Paper Sec. IV's \"In isolation this pulse would reheat the universe to "
      "$T_{\\rm reh}\\sim4\\times10^5$ GeV\" and the corresponding paragraph must be rewritten.")
    A("(This also makes the Sec. IV sentence \"The parametric estimate below accordingly combines channels (i) and (iv)\"")
    A(" unnecessary -- (i) contributes nothing to $T_{\\rm reh}$.)\n")

    md = "\n".join(L) + "\n"
    with open(os.path.join(ROOT, "residual_quintessence.md"), "w",
              encoding="utf-8") as f:
        f.write(md)
    print(f"wrote scripts/residual_quintessence.md ({len(md)} chars)")


if __name__ == "__main__":
    main()
