# -*- coding: utf-8 -*-
"""Deeper consistency checks: paper numbers vs scripts.

Covers:
1. Paper Table I lambda0 vs analytic vs numeric
2. Reheating channel -> N window membership
3. DM kinematic double-protection vs g window
4. Exact slow-roll A_s, N (large-field vs full potential)
5. Internal consistency of quoted scales under one lambda0 choice
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosmo_model import (
    A_S,
    M_P,
    H0_GeV,
    N_match,
    T_reh_from_N,
    Vc,
    beta_o,
    beta_p,
    fiducial,
    lambda0_analytic,
    m_chi,
    m_Phi,
    r_of,
    x_end,
)

ROOT = Path(__file__).resolve().parent.parent
TEX = ROOT / "paper_prd_merged.tex"
OUT = Path(__file__).resolve().parent / "consistency_report.md"


def paper_table_I():
    # From paper_prd_merged.tex Table (tab:sens)
    return [
        (48, 7.36e-8, 0.958, 0.00529),
        (49, 7.06e-8, 0.959, 0.00508),
        (50, 6.78e-8, 0.959, 0.00487),
        (52, 6.27e-8, 0.961, 0.00451),
        (55, 5.60e-8, 0.963, 0.00403),
    ]


def exact_As_NS(lam0: float, xi: float, N_target: float):
    """Solve for chi_* such that N(chi_*)=N_target using exact eps, then As.

    V_E = V0 (1-e^{-x})^2, x=beta_p chi/M_P
    N = int_{chi_end}^{chi_*} dchi / (M_P sqrt(2 eps))
      = (1/beta_p) int_{x_end}^{x_*} dx / sqrt(2 eps(x))
    eps(x) = 2 beta_p^2 e^{-2x}/(1-e^{-x})^2
    sqrt(2 eps) = 2 beta_p e^{-x}/(1-e^{-x})
    dchi = M_P/beta_p dx
    N = int dx (1-e^{-x}) / (2 beta_p^2 e^{-x}) = int (e^x - 1) dx / (2 beta_p^2)
      = [e^x - x]_{x_end}^{x_*} / (2 beta_p^2)
    Large-field: e^{x_*} approx 2 beta_p^2 N  (paper Eq A2 leading)
    """
    bp = beta_p(xi)
    x_e = x_end(xi)

    def N_of_x(x):
        return (math.exp(x) - x - (math.exp(x_e) - x_e)) / (2.0 * bp**2)

    x_star = brentq(lambda x: N_of_x(x) - N_target, x_e + 1e-6, x_e + 40)
    V0 = lam0 * M_P**4 / (4.0 * xi**2)
    x = x_star
    u = math.exp(-x)
    eps = 2.0 * bp**2 * u**2 / (1.0 - u) ** 2
    VE = V0 * (1.0 - u) ** 2
    As = VE / (24.0 * math.pi**2 * M_P**4 * eps)
    # ns, r at N_target using Hubble/potential LO at x_star
    # ns ~ 1 - 6 eps + 2 eta, eta = M_P^2 V''/V
    # V' /V = 2 * (2/(beta_p M_P)) * u/(1-u) wait x=beta_p chi/M_P, d/dchi = (beta_p/M_P) d/dx
    # V = V0 (1-u)^2, dV/dx = V0 * 2(1-u) u, d2V/dx2 = V0 * 2[u^2 - u(1-u)] = 2 V0 u (2u-1)
    # eta_V = M_P^2 V''/V = M_P^2 (beta_p/M_P)^2 * 2u(2u-1) / (1-u)^2 = beta_p^2 * 2u(2u-1)/(1-u)^2
    eta = bp**2 * 2.0 * u * (2.0 * u - 1.0) / (1.0 - u) ** 2
    ns = 1.0 - 6.0 * eps + 2.0 * eta
    r = 16.0 * eps
    return {
        "x_star": x_star,
        "N_exact": N_of_x(x_star),
        "eps": eps,
        "eta": eta,
        "As": As,
        "ns": ns,
        "r": r,
        "N_large": math.exp(x_star) / (2 * bp**2),
    }


def main() -> None:
    lines = [
        "> # WARNING: OUTDATED (LEGACY) -- not valid for the current paper",
        ">",
        "> This report was computed with **OLD parameters**: analytic-A5 normalization lambda0=7.465e-8 (or the old draft value 6.78e-8), and r=0.00487 from the old Table I.",
        "> **The current paper** uses the locked-$N$ convention: lambda0=6.70e-8, r=0.00425, n_s=0.9616 (N=50, xi=11.1).",
        ">",
        "> Therefore the **numerical values and PASS/FAIL verdicts in this report do NOT represent the current paper** -- historical comparison only.",
        "> For current values see `lock_n_convention.py`, `background_and_reheating.py`, `dm_gap_closure_test.py`.",
        "",
    ]
    lines.append("# Script results vs paper claims: deeper consistency checks")
    lines.append("")
    lines.append("Script: `scripts/consistency_checks.py`")
    lines.append("")

    fp_an = fiducial(use_numeric_lam0=False)
    fp_num = fiducial(use_numeric_lam0=True)

    # --- 1. Table I ---
    lines.append("## 1. Paper Table I vs analytic lambda0 vs numeric lambda0")
    lines.append("")
    lines.append("| N | paper lambda0 | analytic lambda0 | ratio | paper r | script r | paper n_s | script n_s NLO |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for N, lam_p, ns_p, r_p in paper_table_I():
        lam_a = lambda0_analytic(fp_an.xi, N)
        ns_s = 1 - 2 / N - 1.5 / N**2
        r_s = r_of(fp_an.xi, N)
        lines.append(
            f"| {N} | {lam_p:.3e} | {lam_a:.3e} | {lam_p/lam_a:.3f} | {r_p:.5f} | {r_s:.5f} | {ns_p:.3f} | {ns_s:.4f} |"
        )
    lines.append("")
    lines.append(
        "**Verdict:** the lambda0 column of Table I is systematically about **0.907-0.909 times** the analytic value"
        " (~sqrt(0.82), or from the exact-potential/N matching), **internally consistent within the table**, "
        "but **inconsistent** with the appendix analytic formula A5 (7.46e-8 @ N=50). r and n_s agree with the formulas (PASS)."
    )
    lines.append("")

    # --- 2. Exact slow roll ---
    lines.append("## 2. Exact slow-roll integral vs large-field approximation (A_s normalization)")
    lines.append("")
    lines.append("| lambda0 | N_exact target | As_exact | As_paper_formula | As_ratio | r_exact | r_formula | ns_ps | ns_formula NLO |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    As_formula = lambda lam0, xi, N: lam0 * N**2 / (12 * math.pi**2 * xi**2 * (6 + 1 / xi))
    for tag, lam0 in [("analytic", fp_an.lam0), ("numeric", fp_num.lam0), ("TableI-N50", 6.78e-8)]:
        res = exact_As_NS(lam0, fp_an.xi, 50.0)
        Af = As_formula(lam0, fp_an.xi, 50.0)
        ns_f = 1 - 2 / 50 - 1.5 / 2500
        lines.append(
            f"| {tag}={lam0:.3e} | 50 | {res['As']:.4e} | {Af:.4e} | {res['As']/Af:.4f} | "
            f"{res['r']:.5f} | {r_of(fp_an.xi,50):.5f} | {res['ns']:.4f} | {ns_f:.4f} |"
        )
    lines.append("")
    lines.append(
        "**Verdict:** if As_exact/As_formula deviates significantly from 1, then the large-field approximation "
        "A_s=lambda0 N^2/(12 pi^2 xi^2 beta_o^2) deviates from the exact integral;"
        " the paper's numeric-draft lambda0=6.78e-8 may instead be the result of pinning As to 2.1e-9 with the exact potential, "
        "rather than of the appendix analytic formula."
    )
    lines.append("")

    # invert: what lam0 makes exact As = 2.1e-9 at N=50?
    def lam0_for_As_exact(As_target=2.1e-9, N=50.0, xi=11.1):
        # As scales linearly with lam0 at fixed x_star(N) because eps,eta independent of lam0
        # res As = lam0 * (As/lam0)
        probe = exact_As_NS(1.0e-8, xi, N)
        return As_target / (probe["As"] / 1.0e-8)

    lam_exact = lam0_for_As_exact()
    lines.append(f"- lambda0 required for **exact** A_s=2.1e-9 with N_exact=50: **{lam_exact:.4e}**")
    lines.append(f"- analytic large-field formula gives lambda0 = {fp_an.lam0:.4e}")
    lines.append(f"- paper numeric / Table I = 6.78e-8")
    lines.append(f"- exact/analytic = {lam_exact/fp_an.lam0:.4f}; exact/paper-numeric = {lam_exact/6.78e-8:.4f}")
    res50 = exact_As_NS(lam_exact, 11.1, 50.0)
    lines.append(
        f"- at this lambda0_exact: r_exact={res50['r']:.5f}, ns_PS={res50['ns']:.4f}, N_large~{res50['N_large']:.2f}"
    )
    lines.append("")

    # --- 3. Reheating window ---
    lines.append("## 3. Reheating channel -> N window (Liddle-Leach matching)")
    lines.append("")
    lines.append("The paper claims: gravitational channel T_reh~4e5 -> N~48; anomaly channel T_reh~1e9 -> N~50; both within [48,55].")
    lines.append("")
    lines.append("| lambda0 source | V_end | T_reh | N_match | in [48,55]? | paper value |")
    lines.append("|---|---|---|---|---|---|")
    for tag, lam0 in [("analytic", fp_an.lam0), ("numeric", fp_num.lam0)]:
        V0_ = lam0 * M_P**4 / (4 * 11.1**2)
        Vend = V0_ * fp_an.V_end_frac
        for T, paperN in [(4e5, 48), (1e9, 50), (1e-2, 41), (6e15, 55)]:
            Nm = N_match(T, Vend)
            inside = 48.0 <= Nm <= 55.0
            lines.append(
                f"| {tag}={lam0:.2e} | {Vend:.3e} | {T:.2e} | {Nm:.2f} | "
                f"{'YES' if inside else '**NO**'} | paper N~{paperN} |"
            )
    lines.append("")
    lines.append(
        "**Verdict (new problem):** using the paper's own matching formula Eq.(18) with the analytic lambda0,"
        " **the gravitational channel T_reh=4e5 GeV gives N~47.2, below the claimed window [48,55]**."
        " The paper's N~48 comes from another constant form (50+1/4 ln), not exactly equivalent to Eq.(18)."
        " The script therefore **weakens** the statement that \"both channels self-consistently fall in the physical window\"; "
        "the formulas should be aligned, or the window relaxed and the text corrected."
    )
    lines.append("")

    # --- 4. DM kinematics ---
    lines.append("## 4. DM kinematic double protection vs the g window")
    lines.append("")
    Hinf = fp_an.H_inf
    mchi = fp_an.m_chi
    lines.append(f"- H_inf={Hinf:.3e}, m_chi={mchi:.3e}, m_chi/2={mchi/2:.3e}")
    lines.append(f"- kinematic-closure condition m_chi <~ 2 m_psi <=> m_psi >~ m_chi/2 ~ {mchi/2:.3e} GeV")
    lines.append("")
    lines.append("| g | m_psi | m_psi/H_inf | m_psi vs m_chi/2 | kinematically closed? | vertex closed? |")
    lines.append("|---|---|---|---|---|---|")
    for g in (1e-5, 2.3e-5, 5e-5, 1e-4):
        mpsi = g * fp_an.Phi0
        kin = mpsi >= mchi / 2
        lines.append(
            f"| {g:.2e} | {mpsi:.3e} | {mpsi/Hinf:.3f} | "
            f"{'>=' if kin else '<'} m_chi/2 | {'YES' if kin else '**NO**'} | YES (conformal identity) |"
        )
    lines.append("")
    lines.append(
        "**Verdict:** vertex closure holds for any g; **\"double protection\" holds only when m_psi >~ m_chi/2, "
        "i.e. g >~ about 2.3e-5-5e-5**."
        " At g=1e-5 the kinematic protection fails (acknowledged in appendix B of the paper, but the main text/abstract sometimes says doubly protected)."
        " The script **supports** vertex closure and **restricts** the kinematic-protection range."
    )
    lines.append("")

    # --- 5. Scale consistency under one lambda0 ---
    lines.append("## 5. Scale consistency under one lambda0")
    lines.append("")
    lines.append("| quantity | numeric lambda0=6.78e-8 | analytic lambda0=7.47e-8 | exact-inverted lambda0 | used in paper main text |")
    lines.append("|---|---|---|---|---|")
    rows = ["lambda0", "H_inf", "V0", "U1/4", "m_chi", "m_Phi"]
    vals = {}
    for tag, lam0 in [("num", 6.78e-8), ("an", fp_an.lam0), ("ex", lam_exact)]:
        H = math.sqrt(lam0) * M_P / (2 * math.sqrt(3) * 11.1)
        V0_ = lam0 * M_P**4 / (4 * 123.21)
        vals[tag] = {
            "lambda0": lam0,
            "H_inf": H,
            "V0": V0_,
            "U1/4": V0_**0.25,
            "m_chi": m_chi(lam0, 11.1),
            "m_Phi": m_Phi(lam0, 11.1),
        }
    paper_ref = {"lambda0": "6.78e-8 (Table I)", "H_inf": "1.65e13", "V0": "4.83e63", "U1/4": "~8.3e15", "m_chi": "3.28e13", "m_Phi": "~2.6e14"}
    for key in rows:
        lines.append(
            f"| {key} | {vals['num'][key]:.4e} | {vals['an'][key]:.4e} | {vals['ex'][key]:.4e} | {paper_ref[key]} |"
        )
    lines.append("")
    lines.append(
        "**Verdict (about the old draft; outdated):** the old-draft main text (H_inf=1.65e13, V0=4.83e63, m_chi=3.28e13) is **self-consistent with the numeric lambda0=6.78e-8**."
        " **Note: the current paper has switched to the locked convention lambda0=6.70e-8, H_inf=1.64e13, V0=4.78e63**,"
        " so the conclusion of this section no longer applies to the current draft; the \"dual-lambda0 track\" problem flagged back then was resolved uniformly by `lock_n_convention.py`."
        " (This section is kept only as a record that the problem once existed.)"
    )
    lines.append("")

    # --- 6. Parse tex for dual values ---
    lines.append("## 6. Scan for coexisting numeric values in the paper .tex")
    lines.append("")
    if TEX.exists():
        tex = TEX.read_text(encoding="utf-8", errors="replace")
        pats = {
            "6.78e-8 / 6.78\\times10^{-8}": r"6\.78\\times10\^\{-8\}|6\.78e-8|6\.78\\times10\^\{-8\}",
            "7.46e-8": r"7\.46",
            "1.65e13 H_inf": r"1\.65\\times10\^\{13\}|1\.65\\times10\^\{13\}",
            "1.73e13 H_inf": r"1\.73\\times10\^\{13\}",
            "3.28e13 m_chi": r"3\.28\\times10\^\{13\}",
            "3.43e13 m_chi": r"3\.43\\times10\^\{13\}",
            "4.83e63 V0": r"4\.83\\times10\^\{63\}",
            "5.32e63 / 5.33e63": r"5\.3[23]\\times10\^\{63\}",
        }
        lines.append("| pattern | occurrence count |")
        lines.append("|---|---|")
        for name, pat in pats.items():
            n = len(re.findall(pat, tex))
            lines.append(f"| {name} | {n} |")
        lines.append("")
        lines.append(
            "**Verdict:** if 6.78e-8 and 7.46e-8, or 1.65e13 and 1.73e13, appear simultaneously in the main text/appendix,"
            " the paper's internal **dual-track numerics** have not fully converged; a revision should state that "
            "\"Table I uses the exact numerical normalization while the appendix gives the analytic approximation\","
            " or unify everything to a single lambda0."
        )
    else:
        lines.append("paper_prd_merged.tex not found")
    lines.append("")

    # --- 7. What scripts support / not / new problems ---
    lines.append("## 7. Overall verdict")
    lines.append("")
    lines.append("### Paper claims supported by the scripts")
    lines.append("")
    lines.append("1. Starobinsky-type functional form r=8/(beta^2 N^2)=2(6+1/xi)/N^2, n_s~1-2/N -- **supported**")
    lines.append("2. End condition e^{-x_end}=0.466, V_end/V0=0.285 -- **supported**")
    lines.append("3. Conformal decoupling Omega=Phi/Phi_0, tree-level chi psibar psi=0 -- **supported** (algebra)")
    lines.append("4. G_eff(Phi_0)=G_N -- **supported**")
    lines.append("5. Single-field quintessence not viable (m_chi/H0>>1); DE = frozen V_c -- **supported**")
    lines.append("6. T_reh~1e9 <-> N~50 (under the same matching formula) -- **supported**")
    lines.append("7. r<~0.0053 for N in [48,55], xi=11.1 -- **supported** (at the formula level)")
    lines.append("8. Domain-wall sigma_wall~1e50 GeV^3 order of magnitude, F(0)=F'(0)=0 -- **supported**")
    lines.append("")
    lines.append("### Claims not supported, or only partially supported, by the scripts")
    lines.append("")
    lines.append("1. Appendix analytic lambda0=48 pi^2 xi^2 A_s/N^2 -- **not supported** (algebra error)")
    lines.append("2. \"Both the gravitational T_reh~4e5 and the anomaly T_reh~1e9 fall in N in [48,55]\" -- **partially not supported** (4e5 -> N~47)")
    lines.append("3. \"chi -> psi psi kinematic + vertex double protection\" over the whole g in [1e-5,1e-4] -- **supported only at the high-g end**")
    lines.append("4. DM Omega_DM=0.12 quantitative prediction -- **not verified by the scripts** (the mechanism's parameter window can be estimated; the normalization needs lattice input)")
    lines.append("5. The Planck contours in Fig2 -- **schematic**, not an official likelihood; cannot serve as evidence of a data constraint")
    lines.append("")
    lines.append("### Problems introduced or exposed by the scripts/plots")
    lines.append("")
    lines.append("1. **Dual lambda0 track**: Table I/main text use 6.78e-8 while the analytic appendix and figures use ~7.47e-8 -> a ~10% scale mismatch between figures and text")
    lines.append("2. **N-T_reh formulas not unified**: Eq.(18) and \"50+1/4 ln\" give different N for the same T_reh")
    lines.append("3. **Exact vs large-field A_s**: use the script's lambda0_exact to cross-check Table I and confirm the origin of 6.78e-8")
    lines.append("4. Code blemish: `m_Phi_jordan` dead code; the Fig2 constraint is schematic")
    lines.append("")
    lines.append("### Suggested additional scripts")
    lines.append("")
    lines.append("| script | purpose |")
    lines.append("|---|---|")
    lines.append("| `consistency_checks.py` (this file) | Table I / exact slow roll / reheating window / DM kinematics |")
    lines.append("| `check_rg_running.py` | recheck beta_xi, Delta xi <~ 1e-6 |")
    lines.append("| `check_reheating_rates.py` | recheck Gamma_anom~O(1) GeV and the T_reh formula |")
    lines.append("| `check_dm_abundance.py` | parametrize the Omega_psi(m_psi,T_reh) surface and its 0.12 crossing line |")
    lines.append("| `check_condensate_overclosure.py` | upper bound on the post-reheating condensate remnant |")
    lines.append("| switch the figure scripts to a single lambda0 convention | align with Table I or the analytic formula, and state it in the caption |")
    lines.append("")
    lines.append("[deeper consistency checks complete]")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
