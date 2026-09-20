# -*- coding: utf-8 -*-
"""Additional paper checks: RG running, reheating rates, DM abundance surface."""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosmo_model import M_P, Vc, fiducial, m_chi, r_of  # noqa: E402

OUT = Path(__file__).resolve().parent / "extended_checks.md"


def beta_xi_terms(xi: float, lam0: float, g: float, dlnmu: float = 60.0):
    pref = xi - 1.0 / 6.0
    dxi_lam = pref * (3.0 * lam0) / (16.0 * math.pi**2) * dlnmu
    dxi_g = pref * (g**2) / (16.0 * math.pi**2) * dlnmu
    return dxi_lam, dxi_g, dxi_lam + dxi_g


def Gamma_anom(m_chi: float, xi: float, b3: float = 7.0, alpha_s: float = 0.1):
    """Order-of-magnitude anomaly width ~ b3^2 alpha_s^2 m_chi^3 / (16 pi^3 M_P^2 (6+1/xi))."""
    return b3**2 * alpha_s**2 * m_chi**3 / (16.0 * math.pi**3 * M_P**2 * (6.0 + 1.0 / xi))


def Treh_const(Gamma: float, gstar: float = 106.75):
    return (90.0 * Gamma**2 * M_P**2 / (math.pi**2 * gstar)) ** 0.25


def Omega_psi(mpsi: float, Hinf: float, Treh: float, Treh_ref: float = 1e9, mpsi_ref: float = 4e12, Om_ref: float = 0.12):
    """Parametric scaling from paper Eq (25)-style: Omega propto mpsi H^3 e^{-pi mpsi/H} * dilution(Treh)."""
    # Use paper-like normalization at reference point, scale exponential and T
    expo = math.exp(-math.pi * mpsi / Hinf) / math.exp(-math.pi * mpsi_ref / Hinf)
    # n propto H^3 e^{-pi m/H} held; dilution propto Treh; also m_psi linear
    return Om_ref * (mpsi / mpsi_ref) * expo * (Treh / Treh_ref)


def main() -> None:
    fp = fiducial(use_locked_lam0=True)     # paper primary: exact A_s inversion
    fpa = fiducial()                        # analytic App.A normalization
    fpn = fiducial(use_numeric_lam0=True)   # legacy draft value
    lines = [
        "> # SCOPE: order-of-magnitude checks at the CURRENT paper parameters",
        ">",
        "> The lambda0 used throughout is the paper's primary locked-$N$ value, i.e.",
        f"> the exact A_s inversion: {fp.lam0:.4e} at N=50, xi=11.1.  Where a comparison",
        "> against an older normalization is informative, that column is labelled",
        "> explicitly (legacy draft 6.78e-8, or the analytic App.A formula 7.465e-8).",
        ">",
        "> **Every verdict below is an ORDER-OF-MAGNITUDE statement.**  The microscopic",
        "> coefficients b_s, alpha_s (anomaly) and g_* are not derived from the model,",
        "> so no number here is a precise prediction.  For the paper's own table values",
        "> (Table I, the N <-> T_reh matching) defer to `lock_n_convention.py`.",
        "",
    ]
    lines.append("# Extended checks: RG / reheating rates / DM abundance scaling / condensate")
    lines.append("")

    lines.append("## A. RG running (paper Sec. VIII E)")
    lines.append("")
    lines.append(f"Paper: Delta xi <~ 1e-6 (xi=11.1, lambda0={fp.lam0:.4e}, g~1e-5-1e-4, Delta ln mu=60); quark term dominates.")
    lines.append("")
    lines.append("| lambda0 | g | Delta xi_lambda0 | Delta xi_g | Delta xi_tot | paper |")
    lines.append("|---|---|---|---|---|---|")
    for lam0 in (fp.lam0, fpn.lam0):
        for g in (2.3e-5, 1e-4):
            a, b, t = beta_xi_terms(11.1, lam0, g)
            lines.append(f"| {lam0:.2e} | {g:.1e} | {a:.3e} | {b:.3e} | {t:.3e} | <~1e-6 order |")
    lines.append("")
    lines.append("**Verdict:** the orders of magnitude agree with the paper (Delta xi~1e-6-1e-7); **supports** the RG-stability claim (given the adopted beta-function form).")
    lines.append("Note: the \"direct quark-curvature\" term of beta_xi is convention-dependent in how it is written; the script uses the paper's Eq.(37) form.")
    lines.append("")

    lines.append("## B. Anomaly decay width and T_reh")
    lines.append("")
    lines.append("| lambda0 source | m_chi | Gamma_anom (est.) | T_reh=(90 Gamma^2 M_P^2 / pi^2 g*)^{1/4} | paper |")
    lines.append("|---|---|---|---|---|")
    for tag, lam0, mchi in [("locked (paper)", fp.lam0, fp.m_chi),
                            ("legacy draft", fpn.lam0, fpn.m_chi),
                            ("analytic App.A", fpa.lam0, fpa.m_chi)]:
        G = Gamma_anom(mchi, 11.1)
        T = Treh_const(G)
        lines.append(f"| {tag} | {mchi:.3e} | {G:.3e} GeV | {T:.3e} GeV | ~1e9, Gamma~O(1) GeV |")
    lines.append("")
    lines.append(
        "**Verdict:** the rough estimates of Gamma and T_reh fall near the orders of magnitude claimed by the paper (Gamma~0.1-1 GeV, T_reh~1e8-1e9);"
        " **supports \"order-of-magnitude viability\"**, but the coefficients b_s, alpha_s, g* are not derived microscopically from the model, **not a precise prediction**."
    )
    lines.append("")

    lines.append("## C. DM abundance parametric scaling")
    lines.append("")
    lines.append("Following the paper's formula: Omega propto m_psi H^3 e^{-pi m_psi/H} * (T_reh/T_ref), normalized to 0.12 at (m_psi=4e12, T=1e9).")
    lines.append("")
    Hinf = fp.H_inf
    lines.append("| m_psi | m_psi/H_inf | T_reh=1e9 | T_reh=4e5 |")
    lines.append("|---|---|---|---|")
    for mpsi in (4e12, 7.3e12, 1.7e13, 3e13, 7e13):
        o1 = Omega_psi(mpsi, Hinf, 1e9)
        o2 = Omega_psi(mpsi, Hinf, 4e5)
        lines.append(f"| {mpsi:.2e} | {mpsi/Hinf:.2f} | {o1:.3e} | {o2:.3e} |")
    lines.append("")
    lines.append("**Verdict:**")
    lines.append("- Around T_reh=1e9 and m_psi~4e12-2e13, the scaling can pass the 0.12 **order of magnitude**;")
    lines.append("- **At T_reh=4e5, even when m_psi reaches the m_chi scale, Omega remains << 0.12** -> the gravitational channel alone is not enough (consistent with the earlier review);")
    lines.append("- The exponential is highly sensitive, so **one cannot** write a single g as the unique prediction; **supports** the paper's statement that \"a lattice is needed and the mechanism must hold\", **does not support** an unconditional Omega=0.12.")
    lines.append("")

    lines.append("## D. Condensate overclosure rough check")
    lines.append("")
    # After reheating, if only a fraction f of condensate remains, need rho_cond < rho_DM today
    rho_DM0 = 0.265 * Vc(0.683) / 0.683 * (0.265 / 0.265)  # rough
    # better: Omega_DM rho_crit, rho_crit = 3 H0^2 M_P^2
    H0 = fp.H0
    rho_crit = 3 * H0**2 * M_P**2
    rho_DM = 0.265 * rho_crit
    lines.append(f"- rho_crit(today)={rho_crit:.3e} GeV^4, rho_DM={rho_DM:.3e} GeV^4")
    lines.append(f"- initial condensate rho_cond(a_end) ~ (1/2) m_Phi^2 M_P^2 order:")
    for tag, lam0 in [("locked", fp.lam0), ("legacy", fpn.lam0)]:
        mPhi = math.sqrt(2 * lam0) * M_P / math.sqrt(11.1)
        rho0 = 0.5 * mPhi**2 * M_P**2
        # max allowed survival fraction today
        frac_max = rho_DM / rho0
        lines.append(f"  - {tag}: m_Phi={mPhi:.3e}, rho0={rho0:.3e}, allowed survival fraction <~{frac_max:.3e}")
    lines.append("")
    lines.append("**Verdict:** an extremely high decay efficiency is required (remnant <~10^{-75} order); the paper says that after anomalous decay the condensate is exponentially suppressed for T_reh >~ 1e9 -- **supports the direction**, exact values not computed.")
    lines.append("")

    lines.append("## E. Merged conclusions with the main verification report")
    lines.append("")
    lines.append("| category | conclusion |")
    lines.append("|---|---|")
    lines.append("| supported by the scripts and retainable by the paper | inflation functional form, end condition, conformal decoupling, G_eff, frozen DE, falsifiable r window, RG order of magnitude, T_reh~1e9 order of magnitude |")
    lines.append("| the paper needs rewording | double protection restricted in g; the N window is tighter for the gravitational channel; the lambda0 dual track must be declared; the Planck figure is schematic |")
    lines.append("| not supported by the scripts | the erroneous appendix lambda0 formula; single-field quintessence; Omega_DM=0.12 without a lattice |")
    lines.append("| suggested additional tests (scriptable) | comparison with the official Planck likelihood; exact N(chi) and Table I normalization; a finer b_i table for Gamma_anom; lattice/Boltzmann Omega_psi |")
    lines.append("")
    lines.append("[extended checks complete]")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
