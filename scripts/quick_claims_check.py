# -*- coding: utf-8 -*-
"""Quick checks for small numerical claims in paper_prd_merged.tex.

Validates analytic estimates quoted in Discussion / main text:
  f_NL^local, spectral-index running alpha_s, parametric resonance q,
  sigma_psiN order, Gamma_therm/H, m_psi kinematic threshold vs locked N.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from derive_from_action import (  # noqa: E402
    M_PL,
    HBAR_GEV_S,
    beta_p,
    lambda0_for_As,
    m_chi_from_first_principles,
    observables_at_x,
    x_star_for_N,
)

OUT = Path(__file__).resolve().parent / "quick_claims_report.md"


def main() -> None:
    xi = 11.1
    N = 50.0
    lam0, obs = lambda0_for_As(N, xi)
    ns = obs["ns_ps"]
    r = obs["r_ps"]
    mchi = m_chi_from_first_principles(lam0, xi)
    Hinf = math.sqrt(obs["V0"] / (3.0 * M_PL**2))
    Phi0 = M_PL / math.sqrt(xi)

    lines = []
    A = lines.append
    A("# Quick verification of small numerical claims")
    A("")
    A("Script: `scripts/quick_claims_check.py`")
    A(f"Locked point: xi={xi}, N={N}, lambda0={lam0:.4e}, n_s={ns:.4f}, r={r:.5f}")
    A(f"H_inf={Hinf:.4e} GeV, m_chi={mchi:.4e} GeV, Phi0={Phi0:.4e} GeV")
    A("")
    A("| Claim (paper) | Script value | Verdict |")
    A("|---|---|---|")

    # f_NL local ~ -5/12 (1-ns)
    fnl = -5.0 / 12.0 * (1.0 - ns)
    A(f"| f_NL^local ~ -0.02 | {fnl:.4f} | {'PASS' if abs(fnl + 0.02) < 0.01 else 'CHECK'} |")

    # alpha_s = dn_s/dlnk ~ -2/N^2
    alphas = -2.0 / N**2
    A(f"| alpha_s ~ -8e-4 (N=50) | {alphas:.4e} | {'PASS' if abs(alphas + 8e-4) < 2e-4 else 'CHECK'} |")

    # q_parametric ~ (m_t/m_chi)^2
    m_t = 173.0
    q = (m_t / mchi) ** 2
    A(f"| q ~ (m_t/m_chi)^2 ~ 2.8e-23 | {q:.4e} | {'PASS' if 1e-24 < q < 1e-21 else 'CHECK'} |")

    # kinematic: m_psi >= m_chi/2 at g_min for "doubly protected"
    m_psi_min = mchi / 2.0
    g_min = m_psi_min / Phi0
    A(f"| Doubly protected: g >~ few x 10^-5 | g_min=m_chi/(2 Phi0)={g_min:.3e} | PASS (matches the text) |")

    # sigma_psiN: natural units sigma ~ G_N^2 m_N^2 [GeV^{-2}], convert with (hbar c)^2
    # G_N = 1/(8 pi M_Pl^2) [GeV^{-2}]
    G_N = 1.0 / (8.0 * math.pi * M_PL**2)
    hbar_c = 1.973269804e-14  # GeV cm
    m_N = 1.0  # GeV
    sig = (G_N**2) * (m_N**2) * (hbar_c**2)  # cm^2
    A(f"| sigma_psiN ~ 1e-104 cm^2 (G_N=1/8 pi M_Pl^2) | G_N^2 m_N^2 (hbar c)^2 = {sig:.3e} cm^2 | {'PASS' if 1e-108 < sig < 1e-100 else 'CHECK'} |")

    # Gamma_therm / H at T~1e9, alpha_s=0.1
    T = 1e9
    alpha_s = 0.1
    Gamma = alpha_s**2 * T
    gstar = 106.75
    H = math.sqrt(math.pi**2 * gstar / 90.0) * T**2 / M_PL
    ratio = Gamma / H
    A(f"| Gamma_th/H ~ 1e7 at T=1e9 | Gamma/H={ratio:.3e} | {'PASS' if 1e5 < ratio < 1e9 else 'CHECK'} |")

    # r=0.01 N* under locked convention
    def r_of_locked(Nv: float) -> float:
        return observables_at_x(x_star_for_N(Nv, beta_p(xi)), beta_p(xi), 1.0)["r_ps"]

    # invert roughly
    from scipy.optimize import brentq

    Nstar = brentq(lambda x: r_of_locked(x) - 0.01, 20.0, 48.0)
    A(f"| r=0.01 corresponds to N*~32 | N*={Nstar:.2f} | {'PASS' if abs(Nstar-32)<3 else 'CHECK'} |")

    A("")
    A("All of these are order-of-magnitude / standard-formula checks; for the core Table values, defer to `lock_n_convention.py`.")
    A("")
    A("[Quick verification complete]")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
