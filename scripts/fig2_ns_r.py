# -*- coding: utf-8 -*-
"""Fig.2 -- n_s-r under LOCKED cosmological N (exact potential slow-roll).

MUST match paper Table tab:sens / scripts/lock_n_convention.py.
Writes figures/fig2_ns_r.pdf, the file the manuscript includes.  It used to
also write a second, byte-identical copy under an alias name "so tex
includegraphics stays valid regardless of filename choice"; the manuscript never
used that name and lock_n_convention.py wrote the same bytes again, so the alias
was removed.

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
from derive_from_action import (  # noqa: E402
    N_S_OBS,
    SIG_NS,
    lambda0_for_As,
    observables_at_x,
    x_star_for_N,
)

FIGDIR = Path(__file__).resolve().parent.parent / "figures"
FIGDIR.mkdir(exist_ok=True)


def locked_point(N: float, xi: float = 11.1) -> dict:
    _, obs = lambda0_for_As(N, xi)
    return obs


def main() -> None:
    xi = 11.1
    Ns = np.array([48, 49, 50, 51, 52, 55], dtype=float)
    pts = [locked_point(float(N), xi) for N in Ns]
    ns = np.array([p["ns_ps"] for p in pts])
    r = np.array([p["r_ps"] for p in pts])

    Nfine = np.linspace(44, 60, 80)
    ns_f, r_f = [], []
    for N in Nfine:
        p = locked_point(float(N), xi)
        ns_f.append(p["ns_ps"])
        r_f.append(p["r_ps"])
    ns_f = np.array(ns_f)
    r_f = np.array(r_f)

    fig, ax = plt.subplots(figsize=(6.2, 4.6), dpi=200)
    # schematic Planck 1sigma (not official likelihood)
    th = np.linspace(0, 2 * np.pi, 300)
    ax.fill(
        N_S_OBS + SIG_NS * np.cos(th),
        0.004 * np.sin(th),
        color="#d6dce4",
        alpha=0.7,
        label="Planck 1 sigma (schematic)",
    )
    ax.plot(
        ns_f,
        r_f,
        color="#1f4e79",
        lw=2,
        label=r"locked $N=\ln(a_{\rm end}/a_*)$, exact PS",
    )
    for N, x, y in zip(Ns, ns, r):
        ax.scatter([x], [y], c="#c00000", s=40, zorder=5)
        # Stagger labels so the crowded N=48..52 cluster stays readable.
        label = f"N={int(N)}"
        if N <= 52:
            idx = int(N) - 48
            # Spread the crowded N=48..52 cluster in a diagonal fan.
            offsets = [(-18, 12), (-6, 18), (4, 8), (12, -4), (2, -14)]
            xytext = offsets[idx]
        else:
            xytext = (10, 2)
        ax.annotate(
            label,
            (x, y),
            textcoords="offset points",
            xytext=xytext,
            fontsize=8,
            ha="left" if N <= 52 else "left",
        )
    ax.axhline(0.036, color="k", ls="--", lw=1, label=r"$r<0.036$ BICEP/Keck")
    ax.axhline(0.01, color="#7030a0", ls=":", lw=1.2, label=r"$r=0.01$ LiteBIRD")
    ax.set_xlabel(r"$n_s$")
    ax.set_ylabel(r"$r$")
    ax.set_title(r"Fig.2 (locked $N$)  $\xi=11.1$, exact potential slow-roll")
    ax.set_xlim(0.954, 0.972)
    ax.set_ylim(-0.002, 0.040)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, loc="upper right")
    ax.text(
        0.955,
        0.032,
        "Planck ellipse schematic; N = cosmological e-folds",
        fontsize=7,
        color="#555555",
    )
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"fig2_ns_r.{ext}")
    print("Wrote locked-N Fig.2 as fig2_ns_r.pdf")
    for N, p in zip(Ns, pts):
        print(f"N={int(N)}: ns={p['ns_ps']:.4f}, r={p['r_ps']:.5f}")


if __name__ == "__main__":
    main()
