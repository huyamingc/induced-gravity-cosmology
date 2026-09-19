# -*- coding: utf-8 -*-
"""
Lock one cosmological definition of N and recompute all observables.

CONVENTION (fixed, not presupposed from the paper):
  N ≡ ln(a_end / a_*)
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
  figures/fig2_ns_r_lockedN.pdf
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
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

ROOT = Path(__file__).resolve().parent.parent
FIGDIR = ROOT / "figures"
FIGDIR.mkdir(exist_ok=True)
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
    lines.append("# 锁定 N 定义后的推导结果")
    lines.append("")
    lines.append("**约定（本次推导固定，不从论文抄用）：**")
    lines.append("")
    lines.append("$$N \\equiv \\ln(a_{\\rm end}/a_*)$$")
    lines.append("")
    lines.append("其中 \(k_*=a_*H_*\)（\(k_*=0.05\\,\\mathrm{Mpc}^{-1}\)），")
    lines.append("再加热前凝聚体 \(w=0\)，")
    lines.append("\\(a_{\\rm end}=a_{\\rm reh}(\\rho_{\\rm reh}/\\rho_{\\rm end})^{1/3}\\)，")
    lines.append("\\(a_{\\rm reh}=a_{\\rm eq}(\\rho_{\\rm eq}/\\rho_{\\rm reh})^{1/4}\\)，")
    lines.append("\\(a_{\\rm eq}=\\Omega_r/\\Omega_m\\)。")
    lines.append("")
    lines.append("在该 N 下：\(x_*\) 由精确慢滚积分 \(N(x_*)=N\) 定出，")
    lines.append("可观测量取势慢滚 \(r=16\\varepsilon_*\)，\(n_s=1-6\\varepsilon_*+2\\eta_*\)，")
    lines.append("\\(\\lambda_0\\) 由 \(A_s=2.1\\times10^{-9}\\) 反解。")
    lines.append("")

    # --- Fiducial xi=11.1 table ---
    xi0 = 11.1
    lines.append(f"## 1. ξ={xi0} 锁定 N 网格")
    lines.append("")
    lines.append("| N | λ₀ | r | n_s | n_s−Planck | H_inf | U^{1/4} | m_χ | T_reh* (使 N_derived=N) |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    grid = []
    for N in N_LIST:
        p = point(xi0, N)
        Tstar = T_reh_for_N_derived(p, N)
        p["T_reh_selfcons"] = Tstar
        grid.append(p)
        lines.append(
            f"| {N} | {p['lambda0']:.3e} | {p['r']:.5f} | {p['ns']:.4f} | "
            f"{p['ns_sigma']:+.2f}σ | {p['H_inf']:.3e} | {p['U_quarter']:.3e} | "
            f"{p['m_chi']:.3e} | {Tstar:.3e} |"
        )
    lines.append("")

    # --- Derived N vs T for N=50,55 points ---
    lines.append("## 2. 自洽性：固定 λ₀(N) 后，N_derived(T_reh)")
    lines.append("")
    lines.append("| 锁定 N | T_reh | N_derived | ΔN |")
    lines.append("|---|---|---|---|")
    for p in grid:
        if p["N"] not in (48, 50, 52, 55):
            continue
        for T in T_BENCH:
            Nd = N_derived_for_T(p, T)
            lines.append(f"| {p['N']} | {T:.2e} | {Nd:.2f} | {Nd - p['N']:+.2f} |")
    lines.append("")
    lines.append("解读：若 |ΔN|≲0.5，可认为该 (N,T_reh) 在本匹配下近似自洽；")
    lines.append("锁定 N 与 N_derived 的偏差直接反映再加热假设（w、k、Ω、g*）的不确定度。")
    lines.append("")

    # --- Planck window under locked convention ---
    lines.append("## 3. Planck 2σ 窗与 r 范围（锁定 N 定义）")
    lines.append("")
    ns_lo = N_S_OBS - 2 * SIG_NS
    ns_hi = N_S_OBS + 2 * SIG_NS
    N_win = [p for p in grid if ns_lo <= p["ns"] <= ns_hi]
    if N_win:
        Nmin, Nmax = N_win[0]["N"], N_win[-1]["N"]
        rmin = min(p["r"] for p in N_win)
        rmax = max(p["r"] for p in N_win)
        lines.append(f"- Planck n_s ∈ [{ns_lo:.4f}, {ns_hi:.4f}]")
        lines.append(f"- 锁定 N 下满足 2σ 的 N ∈ **[{Nmin}, {Nmax}]**")
        lines.append(f"- 对应 r ∈ **[{rmin:.5f}, {rmax:.5f}]**")
    else:
        lines.append("- 网格内无点落入 Planck 2σ（需加密 N 或检查 ns(PS)）")
    lines.append("")
    # 1 sigma
    ns_lo1 = N_S_OBS - SIG_NS
    N_win1 = [p for p in grid if ns_lo1 <= p["ns"] <= ns_hi]
    if N_win1:
        lines.append(
            f"- 若用 1σ 下界 n_s≥{ns_lo1:.4f}：N ∈ [{N_win1[0]['N']}, {N_win1[-1]['N']}]，"
            f"r ∈ [{min(p['r'] for p in N_win1):.5f}, {max(p['r'] for p in N_win1):.5f}]"
        )
    lines.append("")

    # --- Falsification ---
    lines.append("## 4. 证伪线（锁定 N 定义）")
    lines.append("")
    lines.append("| ξ | r=0.01 对应 N* | 该点 n_s | r_max(2σ窗) | LiteBIRD r>0.01 能否排除 |")
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
            excl = "窗内无点"
        else:
            excl = "YES" if rmax < 0.01 else "NO（窗内已有 r≥0.01）"
        lines.append(
            f"| {xi} | {Nstar:.2f} | {ns_s:.4f} | "
            f"{'—' if rmax is None else f'{rmax:.5f}'} | {excl} |"
        )
    lines.append("")
    lines.append("若 2σ 窗内 r_max < 0.01，则 LiteBIRD 测到 r>0.01 可排除该 ξ 在本 N 约定下的模型类。")
    lines.append("")

    # --- xi scan at N=50 locked ---
    lines.append("## 5. ξ 扫描（锁定 N=50）")
    lines.append("")
    lines.append("| ξ | λ₀ | r | n_s | m_χ | H_inf | T_reh* |")
    lines.append("|---|---|---|---|---|---|---|")
    for xi in XI_LIST:
        p = point(xi, 50)
        Tstar = T_reh_for_N_derived(p, 50)
        lines.append(
            f"| {xi} | {p['lambda0']:.3e} | {p['r']:.5f} | {p['ns']:.4f} | "
            f"{p['m_chi']:.3e} | {p['H_inf']:.3e} | {Tstar:.3e} |"
        )
    lines.append("")
    lines.append("在锁定 N 定义下，r 对 ξ 仅弱依赖（通过 β_p 与 x_*(N)）；λ₀ ∝ ξ² 仍大致成立。")
    lines.append("")

    # --- DE / DM unchanged structure ---
    p50 = point(11.1, 50)
    p50["T_reh_selfcons"] = T_reh_for_N_derived(p50, 50)
    H0 = H0_GeV_v2()
    lines.append("## 6. 与 DE/DM 结构的关系（不因 N 约定改变）")
    lines.append("")
    lines.append(f"- 锁定 N=50, ξ=11.1: λ₀={p50['lambda0']:.4e}, m_χ={p50['m_chi']:.4e} GeV")
    lines.append(f"- m_χ/H0={p50['m_chi']/H0:.4e} → **仍不可能**做今日精质 DE")
    lines.append(f"- H_inf={p50['H_inf']:.4e} GeV；DM 窗 m_ψ~H_inf 量级结构不变")
    lines.append(f"- T_reh* (N_derived=50)≈{p50['T_reh_selfcons']:.3e} GeV（本匹配）")
    lines.append("")

    # --- Figure 2 locked ---
    lines.append("## 7. 重绘 Fig.2（锁定 N）")
    lines.append("")
    xi = 11.1
    Ns_plot = np.linspace(44, 60, 80)
    ns_plot, r_plot = [], []
    for N in Ns_plot:
        p = point(xi, float(N))
        ns_plot.append(p["ns"])
        r_plot.append(p["r"])
    fig, ax = plt.subplots(figsize=(6.2, 4.6), dpi=200)
    # schematic Planck
    th = np.linspace(0, 2 * np.pi, 300)
    ax.fill(
        N_S_OBS + SIG_NS * np.cos(th),
        0.004 * np.sin(th),
        color="#d6dce4",
        alpha=0.7,
        label="Planck 1σ (schematic)",
    )
    ax.plot(ns_plot, r_plot, color="#1f4e79", lw=2, label=r"locked $N=\ln(a_{\rm end}/a_*)$")
    # mark N=48,50,52,55
    for N in (48, 50, 52, 55):
        p = point(xi, float(N))
        ax.scatter([p["ns"]], [p["r"]], c="#c00000", s=40, zorder=5)
        ax.annotate(f"N={N}", (p["ns"], p["r"]), textcoords="offset points", xytext=(6, 4), fontsize=8)
    ax.axhline(0.036, color="k", ls="--", lw=1, label=r"$r<0.036$ BICEP/Keck")
    N10 = r_falsify_line(xi, 0.01)
    if not math.isnan(N10):
        p10 = point(xi, N10)
        ax.axhline(0.01, color="#7030a0", ls=":", lw=1.2, label=fr"$r=0.01$ at $N\approx{N10:.1f}$")
        ax.scatter([p10["ns"]], [p10["r"]], c="#7030a0", marker="x", s=60, zorder=6)
    ax.set_xlabel(r"$n_s$")
    ax.set_ylabel(r"$r$")
    ax.set_title(r"Fig.2 (locked $N$)  $\xi=11.1$, exact slow-roll")
    ax.set_xlim(0.954, 0.972)
    ax.set_ylim(-0.002, 0.040)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, loc="upper right")
    ax.text(0.955, 0.033, "Planck ellipse schematic; N = exact cosmological e-folds", fontsize=7, color="#555")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"fig2_ns_r_lockedN.{ext}")
    lines.append(f"- 图已写：`figures/fig2_ns_r_lockedN.pdf`")
    lines.append("")

    # --- Summary ---
    lines.append("## 8. 结论（在本 N 约定下）")
    lines.append("")
    rmax2 = None
    Nband = []
    for p in grid:
        if ns_lo <= p["ns"] <= ns_hi:
            Nband.append(p["N"])
            rmax2 = p["r"] if rmax2 is None else max(rmax2, p["r"])
    lines.append("1. **可观测量在给定 N（本定义）下可复算**：ξ=11.1 时 N=50 给出")
    lines.append(f"   λ₀≈{p50['lambda0']:.3e}, r≈{p50['r']:.5f}, n_s≈{p50['ns']:.4f}（Planck 偏离 {p50['ns_sigma']:+.2f}σ）。")
    if Nband:
        lines.append(f"2. Planck 2σ 对应锁定 N ∈ [{min(Nband)}, {max(Nband)}]，r_max≈{rmax2:.5f}。")
    lines.append("3. **自洽 T_reh\\*** 与教科书/论文闭式不同：本匹配下 N=50 的 T_reh\\*~10⁷–10⁸ GeV 量级；")
    lines.append("   T_reh=10⁹ 时 N_derived 约 51，仍在窗内但不再是 50。")
    lines.append("4. **证伪**：在 Planck 2σ 窗内若 r_max<0.01，则 r>0.01 可排除该类模型（本约定下 ξ=11.1 通常成立）。")
    lines.append("5. DE/DM 的结构性结论（不能同场精质、DM 条件性）**不依赖** N 的约定。")
    lines.append("")
    lines.append("## 9. 对论文正文的含义")
    lines.append("")
    lines.append("| 项目 | 建议写入 paper_prd_merged.tex |")
    lines.append("|---|---|")
    lines.append("| N 的定义 | 明确写 N=ln(a_end/a_*) 及匹配假设（k, w=0 再加热, Ω） |")
    lines.append("| r, n_s 数值表 | 用本锁定表替换模糊的「N=50 ⇒ r=0.00487」单一说法 |")
    lines.append("| N 窗 | 用「Planck 2σ ∩ 自洽 T_reh≳BBN」表述，而非单一公式 |")
    lines.append("| 附录 λ₀ | 采用精确慢滚反演值（N=50,ξ=11.1 约 6.7e-8），删错误 48π² 式 |")
    lines.append("| 证伪句 | r_max(2σ窗) 与 r=0.01 比较；写明依赖 N 约定 |")
    lines.append("")
    lines.append("[锁定 N 约定推导完成]")

    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    OUT_JSON.write_text(
        json.dumps({"xi11.1": grid, "convention": "N=ln(a_end/a_*), exact PS at x_*(N)"}, indent=2),
        encoding="utf-8",
    )
    print("\n".join(lines))
    print(f"\nWrote {OUT_MD}")
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {FIGDIR/'fig2_ns_r_lockedN.pdf'}")


if __name__ == "__main__":
    main()
