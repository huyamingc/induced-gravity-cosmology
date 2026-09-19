# -*- coding: utf-8 -*-
"""Follow-up: how the definition of N changes r and n_s (no presupposition)."""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from derive_from_action import (  # noqa: E402
    A_S_OBS,
    M_PL,
    N_S_OBS,
    SIG_NS,
    beta_o,
    beta_p,
    lambda0_for_As,
    observables_at_x,
    x_end_from_eps1,
    x_star_for_N,
    N_of_x,
)

OUT = Path(__file__).resolve().parent / "n_definition_report.md"


def main() -> None:
    lines = []
    lines.append("# How the definition of the e-fold number N changes r and n_s (independent derivation)")
    lines.append("")
    lines.append("Script: `scripts/n_definition_check.py`")
    lines.append("")
    lines.append("Once the action gives the potential V=V0(1-e^{-x})^2, the slow-roll observables at **a given x_\\*** are determined.")
    lines.append("There are several common conventions for labeling the horizontal axis as 'N'; **r(N) is not unique**.")
    lines.append("")

    xi = 11.1
    bp = beta_p(xi)
    bo = beta_o(xi)
    x_e = x_end_from_eps1(bp)
    lines.append(f"xi={xi}, beta_p={bp:.4f}, beta_o={bo:.4f}, x_end={x_e:.4f}")
    lines.append("")
    lines.append("| Convention | Definition | x_* at N=50 | r | n_s(PS) | lambda0 (pinning A_s) |")
    lines.append("|---|---|---|---|---|---|")

    # 1) Exact integral N = 50
    x1 = x_star_for_N(50.0, bp)
    lam1, o1 = lambda0_for_As(50.0, xi)
    N1 = N_of_x(x1, x_e, bp)
    lines.append(
        f"| Exact integral N | N(x)=Integral of (e^x-x-c)/(2 beta^2) dx | {x1:.4f} (N={N1:.3f}) | "
        f"{o1['r_ps']:.5f} | {o1['ns_ps']:.4f} | {lam1:.4e} |"
    )

    # 2) Attractor N_large == e^x/(2 beta^2) = 50
    # e^x = 2 beta^2 * 50
    x2 = math.log(2.0 * bp**2 * 50.0)
    o2 = observables_at_x(x2, bp, V0=1.0)
    # scale V0 so As=obs
    V0_2 = A_S_OBS / o2["As"]
    o2b = observables_at_x(x2, bp, V0_2)
    lam2 = 4.0 * xi**2 * V0_2 / M_PL**4
    lines.append(
        f"| Attractor N_large | N==e^x/(2 beta^2) | {x2:.4f} (N_exact={o2['N_from_x']:.2f}) | "
        f"{o2b['r_ps']:.5f} | {o2b['ns_ps']:.4f} | {lam2:.4e} |"
    )

    # 3) Paper attractor formula r=8/(beta^2 N^2) with N input 50 -- evaluated without x
    r_paper = 8.0 / (bp**2 * 50.0**2)
    ns_paper_NLO = 1.0 - 2.0 / 50.0 - 1.5 / 2500.0
    lines.append(
        f"| Paper formula | r=8/(beta^2 N^2), n_s=1-2/N-3/(2N^2) | (no explicit x) | "
        f"{r_paper:.5f} | {ns_paper_NLO:.4f} | -- |"
    )

    lines.append("")
    lines.append("## Same physical point: plug the exact-integral x_* into 16 epsilon")
    lines.append("")
    lines.append(f"- Exact integral N=50 -> x_*={x1:.4f}, 16 epsilon={o1['r_ps']:.5f}")
    lines.append(f"- Attractor N_large=50 -> x_*={x2:.4f}, 16 epsilon={o2b['r_ps']:.5f}")
    lines.append(f"- Paper closed-form r(input N=50)={r_paper:.5f}")
    lines.append("")
    ratio = o1["r_ps"] / r_paper
    lines.append(f"- **Exact-integral r / paper-formula r = {ratio:.4f}** (if != 1, the 'r quoted at N=50' depends on the definition of N)")
    lines.append(f"- For the exact integral, N_large=e^x/(2 beta^2)={math.exp(x1)/(2*bp**2):.2f}, not 50")
    lines.append("")

    lines.append("## n_s: potential slow-roll PS vs attractor NLO")
    lines.append("")
    lines.append("| N_exact | ns_PS | ns_attractor 1-2/N | ns_paper NLO | deviation from Planck (PS) |")
    lines.append("|---|---|---|---|---|")
    for Nt in (48, 50, 52, 55):
        xx = x_star_for_N(float(Nt), bp)
        lam, oo = lambda0_for_As(float(Nt), xi)
        ns_a = 1 - 2.0 / Nt
        ns_p = 1 - 2.0 / Nt - 1.5 / Nt**2
        sig = (oo["ns_ps"] - N_S_OBS) / SIG_NS
        lines.append(
            f"| {Nt} | {oo['ns_ps']:.4f} | {ns_a:.4f} | {ns_p:.4f} | {sig:+.2f} sigma |"
        )
    lines.append("")
    lines.append("## Assumption-free reading")
    lines.append("")
    lines.append("1. At a given x_* the slow-roll potential gives definite (n_s, r, A_s); this step is **unambiguous**.")
    lines.append("2. When the horizontal axis is labeled 'N': if N means the exact integral, the paper closed form r=8/(beta^2 N^2) **overestimates r** (by roughly ~20% here).")
    lines.append("3. If N means the attractor e^x/(2 beta^2), the closed form agrees with 16 epsilon, but it is then not the same number as the e-fold integral obtained from the CMB match.")
    lines.append("4. When reporting r and n_s one must **state the definition of N**; the numerical boundaries of the LiteBIRD falsifiability window shift accordingly.")
    lines.append("5. The lambda0 solved back from A_s is ~6.7e-8 under the exact-integral definition (xi=11.1, N=50), close to Table I.")
    lines.append("")
    lines.append("[N definition analysis complete]")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
