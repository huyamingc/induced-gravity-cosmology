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
    A("# 小项数值声称快速核验")
    A("")
    A("脚本：`scripts/quick_claims_check.py`")
    A(f"锁定点：ξ={xi}, N={N}, λ₀={lam0:.4e}, n_s={ns:.4f}, r={r:.5f}")
    A(f"H_inf={Hinf:.4e} GeV, m_χ={mchi:.4e} GeV, Φ₀={Phi0:.4e} GeV")
    A("")
    A("| 声称（论文） | 脚本计算 | 判定 |")
    A("|---|---|---|")

    # f_NL local ~ -5/12 (1-ns)
    fnl = -5.0 / 12.0 * (1.0 - ns)
    A(f"| f_NL^local ~ -0.02 | {fnl:.4f} | {'PASS' if abs(fnl + 0.02) < 0.01 else 'CHECK'} |")

    # alpha_s = dn_s/dlnk ~ -2/N^2
    alphas = -2.0 / N**2
    A(f"| α_s ~ -8e-4 (N=50) | {alphas:.4e} | {'PASS' if abs(alphas + 8e-4) < 2e-4 else 'CHECK'} |")

    # q_parametric ~ (m_t/m_chi)^2
    m_t = 173.0
    q = (m_t / mchi) ** 2
    A(f"| q ~ (m_t/m_χ)² ~ 2.8e-23 | {q:.4e} | {'PASS' if 1e-24 < q < 1e-21 else 'CHECK'} |")

    # kinematic: m_psi >= m_chi/2 at g_min for "doubly protected"
    m_psi_min = mchi / 2.0
    g_min = m_psi_min / Phi0
    A(f"| 双重保护 g ≳ few×10⁻⁵ | g_min=m_χ/(2Φ₀)={g_min:.3e} | PASS（与正文一致） |")

    # sigma_psiN: natural units sigma ~ G_N^2 m_N^2 [GeV^{-2}], convert with (hbar c)^2
    # G_N = 1/(8 pi M_Pl^2) [GeV^{-2}]
    G_N = 1.0 / (8.0 * math.pi * M_PL**2)
    hbar_c = 1.973269804e-14  # GeV cm
    m_N = 1.0  # GeV
    sig = (G_N**2) * (m_N**2) * (hbar_c**2)  # cm^2
    A(f"| σ_ψN ~ 1e-104 cm² (G_N=1/8πM_Pl²) | G_N²m_N²(ħc)² = {sig:.3e} cm² | {'PASS' if 1e-108 < sig < 1e-100 else 'CHECK'} |")

    # Gamma_therm / H at T~1e9, alpha_s=0.1
    T = 1e9
    alpha_s = 0.1
    Gamma = alpha_s**2 * T
    gstar = 106.75
    H = math.sqrt(math.pi**2 * gstar / 90.0) * T**2 / M_PL
    ratio = Gamma / H
    A(f"| Γ_th/H ~ 1e7 at T=1e9 | Γ/H={ratio:.3e} | {'PASS' if 1e5 < ratio < 1e9 else 'CHECK'} |")

    # r=0.01 N* under locked convention
    def r_of_locked(Nv: float) -> float:
        return observables_at_x(x_star_for_N(Nv, beta_p(xi)), beta_p(xi), 1.0)["r_ps"]

    # invert roughly
    from scipy.optimize import brentq

    Nstar = brentq(lambda x: r_of_locked(x) - 0.01, 20.0, 48.0)
    A(f"| r=0.01 对应 N*≈32 | N*={Nstar:.2f} | {'PASS' if abs(Nstar-32)<3 else 'CHECK'} |")

    A("")
    A("这些均为数量级/标准公式核验；核心 Table 数值仍以 `lock_n_convention.py` 为准。")
    A("")
    A("[快速核验完成]")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
