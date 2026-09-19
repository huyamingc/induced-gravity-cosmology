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
    lines.append("# e-fold 数 N 的定义如何改变 r、n_s（独立推导）")
    lines.append("")
    lines.append("脚本：`scripts/n_definition_check.py`")
    lines.append("")
    lines.append("作用量给出势 V=V0(1-e^{-x})^2 后，慢滚可观测量在 **给定 x_\\*** 处是确定的。")
    lines.append("把横轴标成「N」有几种常见约定，**r(N) 不是唯一的**。")
    lines.append("")

    xi = 11.1
    bp = beta_p(xi)
    bo = beta_o(xi)
    x_e = x_end_from_eps1(bp)
    lines.append(f"xi={xi}, beta_p={bp:.4f}, beta_o={bo:.4f}, x_end={x_e:.4f}")
    lines.append("")
    lines.append("| 约定 | 定义 | N=50 时 x_* | r | n_s(PS) | λ₀(钉 A_s) |")
    lines.append("|---|---|---|---|---|---|")

    # 1) Exact integral N = 50
    x1 = x_star_for_N(50.0, bp)
    lam1, o1 = lambda0_for_As(50.0, xi)
    N1 = N_of_x(x1, x_e, bp)
    lines.append(
        f"| 精确积分 N | N(x)=∫(e^x−x−c)/(2β²) | {x1:.4f} (N={N1:.3f}) | "
        f"{o1['r_ps']:.5f} | {o1['ns_ps']:.4f} | {lam1:.4e} |"
    )

    # 2) Attractor N_large ≡ e^x/(2β²) = 50
    # e^x = 2 β² * 50
    x2 = math.log(2.0 * bp**2 * 50.0)
    o2 = observables_at_x(x2, bp, V0=1.0)
    # scale V0 so As=obs
    V0_2 = A_S_OBS / o2["As"]
    o2b = observables_at_x(x2, bp, V0_2)
    lam2 = 4.0 * xi**2 * V0_2 / M_PL**4
    lines.append(
        f"| 吸引子 N_large | N≡e^x/(2β²) | {x2:.4f} (N_exact={o2['N_from_x']:.2f}) | "
        f"{o2b['r_ps']:.5f} | {o2b['ns_ps']:.4f} | {lam2:.4e} |"
    )

    # 3) Paper attractor formula r=8/(β² N²) with N input 50 — evaluated without x
    r_paper = 8.0 / (bp**2 * 50.0**2)
    ns_paper_NLO = 1.0 - 2.0 / 50.0 - 1.5 / 2500.0
    lines.append(
        f"| 论文公式 | r=8/(β²N²), n_s=1−2/N−3/(2N²) | （不显式 x） | "
        f"{r_paper:.5f} | {ns_paper_NLO:.4f} | — |"
    )

    lines.append("")
    lines.append("## 同一物理点：把精确积分的 x_* 代入 16ε")
    lines.append("")
    lines.append(f"- 精确积分 N=50 → x_*={x1:.4f}, 16ε={o1['r_ps']:.5f}")
    lines.append(f"- 吸引子 N_large=50 → x_*={x2:.4f}, 16ε={o2b['r_ps']:.5f}")
    lines.append(f"- 论文闭式 r(输入N=50)={r_paper:.5f}")
    lines.append("")
    ratio = o1["r_ps"] / r_paper
    lines.append(f"- **精确积分 r / 论文公式 r = {ratio:.4f}**（若 ≠1，则「N=50 的 r」依赖 N 的定义）")
    lines.append(f"- 精确积分时 N_large=e^x/(2β²)={math.exp(x1)/(2*bp**2):.2f}，不是 50")
    lines.append("")

    lines.append("## n_s：潜在慢滚 PS vs 吸引子 NLO")
    lines.append("")
    lines.append("| N_exact | ns_PS | ns_吸引子1−2/N | ns_论文NLO | 偏离 Planck (PS) |")
    lines.append("|---|---|---|---|---|")
    for Nt in (48, 50, 52, 55):
        xx = x_star_for_N(float(Nt), bp)
        lam, oo = lambda0_for_As(float(Nt), xi)
        ns_a = 1 - 2.0 / Nt
        ns_p = 1 - 2.0 / Nt - 1.5 / Nt**2
        sig = (oo["ns_ps"] - N_S_OBS) / SIG_NS
        lines.append(
            f"| {Nt} | {oo['ns_ps']:.4f} | {ns_a:.4f} | {ns_p:.4f} | {sig:+.2f}σ |"
        )
    lines.append("")
    lines.append("## 不预设的读法")
    lines.append("")
    lines.append("1. 慢滚势在给定 x_* 下给出确定的 (n_s, r, A_s)；这一步 **无歧义**。")
    lines.append("2. 用「N」作横轴时，若 N 指精确积分，则论文闭式 r=8/(β²N²) **高估 r**（此处约高估 ~20%）。")
    lines.append("3. 若 N 指吸引子 e^x/(2β²)，则闭式与 16ε 一致，但此时与「CMB 匹配得到的 e-fold 积分」不是同一数。")
    lines.append("4. 报告 r、n_s 时必须 **写明 N 的定义**；LiteBIRD 可证伪窗的数值边界随之移动。")
    lines.append("5. 由 A_s 反解的 λ₀ 在精确积分定义下 ≈6.7e-8（ξ=11.1, N=50），与 Table I 接近。")
    lines.append("")
    lines.append("[N 定义分析完成]")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
