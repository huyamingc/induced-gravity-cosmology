# -*- coding: utf-8 -*-
"""Fig.3 -- Domain-wall two-frame comparison: F(varphi) and V_E / sigma_E artifact.

Type: FIG
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosmo_model import M_P, F_of_varphi, VE_of_varphi, fiducial  # noqa: E402

FIGDIR = Path(__file__).resolve().parent.parent / "figures"
FIGDIR.mkdir(exist_ok=True)


def main() -> None:
    # Locked-N convention (paper Sec. III / Table I): exact slow-roll lambda0 at
    # N=50, xi=11.1 -> 6.70e-8.  Previously used use_numeric_lam0=True (legacy
    # draft value 6.78e-8), which was inconsistent with the paper and with Fig. 1.
    fp = fiducial(use_locked_lam0=True)
    varphi = np.linspace(-3.0, 2.0, 600)
    F = F_of_varphi(varphi, fp.xi)
    VE = VE_of_varphi(varphi, fp.lam0, fp.xi, fp.Vc)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8), dpi=200)

    ax = axes[0]
    ax.plot(varphi, F, color="#1f4e79", lw=2)
    ax.axhline(0.0, color="k", lw=0.8)
    ax.axhline(M_P**2, color="#548235", ls="--", lw=1, label=r"$F(\Phi_0)=M_{\rm Pl}^2$")
    ax.set_xlabel(r"$\varphi=\ln(\Phi/\Phi_0)$")
    ax.set_ylabel(r"$F=\xi\Phi^2=M_{\rm Pl}^2 e^{2\varphi}$")
    ax.set_title(r"Jordan: $F(\varphi)$, $F(0)=F'(0)=0$")
    ax.set_ylim(-0.05 * M_P**2, 1.15 * M_P**2)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    # mark Phi=0 as varphi -> -inf
    ax.annotate(r"$\Phi\to0\Rightarrow\varphi\to-\infty$", xy=(-2.8, 0.02 * M_P**2), fontsize=8)

    ax = axes[1]
    ax.semilogy(varphi, np.maximum(VE, 1e20), color="#833c0c", lw=2)
    ax.axhline(fp.V0, color="#548235", ls=":", lw=1.2, label=r"$V_0$ plateau")
    ax.axhline(fp.Vc, color="#c45911", ls="--", lw=1.2, label=r"$V_c$ DE")
    ax.set_xlabel(r"$\varphi$")
    ax.set_ylabel(r"$V_E$ [GeV$^4$]")
    ax.set_title(r"Einstein: $V_E\to\infty$ as $\varphi\to-\infty$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    ax.annotate("conformal artifact", xy=(-2.5, 1e65), fontsize=8, color="#833c0c")

    fig.suptitle("Fig.3  Domain-wall / EFT boundary at $\\Phi=0$ (two frames)", fontsize=11)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"fig3_domain_wall.{ext}")
    print(f"Wrote {FIGDIR/'fig3_domain_wall.pdf'}")
    print(f"F at varphi=-2: {F_of_varphi(-2.0, fp.xi):.4e} (vs M_P^2={M_P**2:.4e})")


if __name__ == "__main__":
    main()
