# -*- coding: utf-8 -*-
"""Fig.1 -- Einstein-frame potential V_E(varphi): plateau + frozen V_c.

Type: FIG
"""
from __future__ import annotations

import sys
from pathlib import Path
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosmo_model import M_P, Vc, fiducial  # noqa: E402

FIGDIR = Path(__file__).resolve().parent.parent / "figures"
FIGDIR.mkdir(exist_ok=True)


def main() -> None:
    # Locked-N convention (Sec. III of paper): exact slow-roll A_s inversion at
    # N=50, xi=11.1 -> lambda0 = 6.70e-8.  Use the shared helper instead of
    # hand-patching the three derived fields (lam0 / V0 / H_inf).
    fp = fiducial(use_locked_lam0=True)
    Vc_val = fp.Vc
    varphi = np.linspace(-2.0, 6.0, 800)
    V = (fp.lam0 * M_P**4 / (4 * fp.xi**2)) * (1 - np.exp(-2 * varphi)) ** 2 + Vc_val * np.exp(-4 * varphi)
    V_plateau = np.full_like(varphi, fp.V0)

    fig, ax = plt.subplots(figsize=(6.2, 4.2), dpi=200)
    ax.semilogy(varphi, np.maximum(V, 1e-60), color="#1f4e79", lw=2, label=r"$V_E(\varphi)$ exact")
    ax.axhline(Vc_val, color="#c45911", ls="--", lw=1.5, label=rf"$V_c\approx{Vc_val:.1e}\,\mathrm{{GeV^4}}$ (DE)")
    ax.axhline(fp.V0, color="#548235", ls=":", lw=1.5, label=rf"$V_0\approx{fp.V0:.2e}\,\mathrm{{GeV^4}}$ (plateau)")
    ax.axvline(1.0, color="gray", ls="-.", lw=0.8)
    # Keep the inflation marker away from the upper-left legend.
    ax.text(1.05, 1e40, r"$\varphi=1$: inflation", fontsize=9, color="gray")
    ax.text(-1.7, Vc_val * 30, "0: DE", fontsize=9, color="#c45911")
    ax.set_xlabel(r"$\varphi=\ln(\Phi/\Phi_0)$")
    ax.set_ylabel(r"$V_E$  [GeV$^4$]")
    ax.set_title(r"Fig.1  Einstein-frame potential $V_E(\varphi)$")
    ax.set_ylim(1e-50, 1e67)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"fig1_VE_potential.{ext}")
    print(f"Wrote {FIGDIR/'fig1_VE_potential.pdf'}")
    print(f"V0={fp.V0:.4e}, Vc={Vc_val:.4e}, hierarchy={fp.V0/Vc_val:.4e}")


if __name__ == "__main__":
    main()
