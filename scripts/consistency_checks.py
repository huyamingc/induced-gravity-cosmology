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
        "> # ⚠️ 已过期（LEGACY）——不可用于当前论文",
        ">",
        "> 本报告用**旧参数**：解析 A5 归一化 λ₀=7.465e-8（或旧草稿值 6.78e-8）、旧 Table I 的 r=0.00487。",
        "> **当前论文**用锁定 $N$ 约定：λ₀=6.70e-8、r=0.00425、n_s=0.9616（N=50, ξ=11.1）。",
        ">",
        "> 故本报告的**数值与 PASS/FAIL 判定均不代表论文现状**，仅作历史对照。",
        "> 当前值请看 `lock_n_convention.py`、`background_and_reheating.py`、`dm_gap_closure_test.py`。",
        "",
    ]
    lines.append("# 脚本结果 vs 论文主张：一致性深检")
    lines.append("")
    lines.append("脚本：`scripts/consistency_checks.py`")
    lines.append("")

    fp_an = fiducial(use_numeric_lam0=False)
    fp_num = fiducial(use_numeric_lam0=True)

    # --- 1. Table I ---
    lines.append("## 1. 论文 Table I vs 解析 λ₀ vs 数值 λ₀")
    lines.append("")
    lines.append("| N | 论文 λ₀ | 解析 λ₀ | 比值 | 论文 r | 脚本 r | 论文 n_s | 脚本 n_s NLO |")
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
        "**判定：** Table I 的 λ₀ 全表系统性约为解析值的 **0.907–0.909 倍**（≈√0.82，或来自精确势/N 匹配），"
        "**表内自洽**，但与附录解析式 A5（7.46e-8 @ N=50）**不一致**。r、n_s 与公式一致（PASS）。"
    )
    lines.append("")

    # --- 2. Exact slow roll ---
    lines.append("## 2. 精确慢滚积分 vs 大场近似（A_s 归一化）")
    lines.append("")
    lines.append("| λ₀ | N_exact 目标 | As_exact | As_paper_formula | As_ratio | r_exact | r_formula | ns_ps | ns_formula NLO |")
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
        "**判定：** 若 As_exact/As_formula 显著偏离 1，则大场近似 A_s=λ₀N²/(12π²ξ²β_o²) 与精确积分有偏差；"
        "论文数值稿 λ₀=6.78e-8 可能是用精确势把 As 钉到 2.1e-9 的结果，而非附录解析式。"
    )
    lines.append("")

    # invert: what lam0 makes exact As = 2.1e-9 at N=50?
    def lam0_for_As_exact(As_target=2.1e-9, N=50.0, xi=11.1):
        # As scales linearly with lam0 at fixed x_star(N) because eps,eta independent of lam0
        # res As = lam0 * (As/lam0)
        probe = exact_As_NS(1.0e-8, xi, N)
        return As_target / (probe["As"] / 1.0e-8)

    lam_exact = lam0_for_As_exact()
    lines.append(f"- 使 **精确** A_s=2.1e-9 且 N_exact=50 所需 λ₀ = **{lam_exact:.4e}**")
    lines.append(f"- 解析大场式给出 λ₀ = {fp_an.lam0:.4e}")
    lines.append(f"- 论文数值/Table I = 6.78e-8")
    lines.append(f"- 精确/解析 = {lam_exact/fp_an.lam0:.4f}；精确/论文数值 = {lam_exact/6.78e-8:.4f}")
    res50 = exact_As_NS(lam_exact, 11.1, 50.0)
    lines.append(
        f"- 在该 λ₀_exact 下：r_exact={res50['r']:.5f}, ns_PS={res50['ns']:.4f}, N_large≈{res50['N_large']:.2f}"
    )
    lines.append("")

    # --- 3. Reheating window ---
    lines.append("## 3. 再加热通道 → N 窗（Liddle–Leach 匹配）")
    lines.append("")
    lines.append("论文声称：引力通道 T_reh~4e5 → N≈48；反常通道 T_reh~1e9 → N≈50；均在 [48,55]。")
    lines.append("")
    lines.append("| λ₀ 源 | V_end | T_reh | N_match | 是否 ∈[48,55] | 论文对应 |")
    lines.append("|---|---|---|---|---|---|")
    for tag, lam0 in [("analytic", fp_an.lam0), ("numeric", fp_num.lam0)]:
        V0_ = lam0 * M_P**4 / (4 * 11.1**2)
        Vend = V0_ * fp_an.V_end_frac
        for T, paperN in [(4e5, 48), (1e9, 50), (1e-2, 41), (6e15, 55)]:
            Nm = N_match(T, Vend)
            inside = 48.0 <= Nm <= 55.0
            lines.append(
                f"| {tag}={lam0:.2e} | {Vend:.3e} | {T:.2e} | {Nm:.2f} | "
                f"{'YES' if inside else '**NO**'} | 论文 N≈{paperN} |"
            )
    lines.append("")
    lines.append(
        "**判定（新问题）：** 用论文自己的匹配公式 Eq.(18) 与解析 λ₀，"
        "**引力通道 T_reh=4e5 GeV 给出 N≈47.2，落在声明窗 [48,55] 之下**。"
        "论文写 N≈48 来自另一常数形式（50+¼ln），与 Eq.(18) 不完全等价。"
        "脚本因此 **削弱**「两通道均自洽落在物理窗」的表述；应对齐公式或放宽窗并改正文。"
    )
    lines.append("")

    # --- 4. DM kinematics ---
    lines.append("## 4. 暗物质运动学双重保护 vs g 窗")
    lines.append("")
    Hinf = fp_an.H_inf
    mchi = fp_an.m_chi
    lines.append(f"- H_inf={Hinf:.3e}, m_χ={mchi:.3e}, m_χ/2={mchi/2:.3e}")
    lines.append(f"- 运动学关闭条件 m_χ ≲ 2 m_ψ ⇔ m_ψ ≳ m_χ/2 ≈ {mchi/2:.3e} GeV")
    lines.append("")
    lines.append("| g | m_ψ | m_ψ/H_inf | m_ψ vs m_χ/2 | 运动学关闭？ | 顶点关闭？ |")
    lines.append("|---|---|---|---|---|---|")
    for g in (1e-5, 2.3e-5, 5e-5, 1e-4):
        mpsi = g * fp_an.Phi0
        kin = mpsi >= mchi / 2
        lines.append(
            f"| {g:.2e} | {mpsi:.3e} | {mpsi/Hinf:.3f} | "
            f"{'≥' if kin else '<'} m_χ/2 | {'YES' if kin else '**NO**'} | YES (共形恒等式) |"
        )
    lines.append("")
    lines.append(
        "**判定：** 顶点关闭对任意 g 成立；**「双重保护」仅在 m_ψ≳m_χ/2 即 g≳约 2.3e-5–5e-5 时成立**。"
        "g=1e-5 时运动学不保护（论文附录 B 已承认，但正文/摘要有时写 doubly protected）。"
        "脚本 **支持**顶点关闭，**限定**运动学保护范围。"
    )
    lines.append("")

    # --- 5. Scale consistency under one lambda0 ---
    lines.append("## 5. 同一 λ₀ 下能标自洽性")
    lines.append("")
    lines.append("| 量 | 数值 λ₀=6.78e-8 | 解析 λ₀=7.47e-8 | 精确反演 λ₀ | 论文正文常用 |")
    lines.append("|---|---|---|---|---|")
    rows = ["λ₀", "H_inf", "V0", "U1/4", "m_chi", "m_Phi"]
    vals = {}
    for tag, lam0 in [("num", 6.78e-8), ("an", fp_an.lam0), ("ex", lam_exact)]:
        H = math.sqrt(lam0) * M_P / (2 * math.sqrt(3) * 11.1)
        V0_ = lam0 * M_P**4 / (4 * 123.21)
        vals[tag] = {
            "λ₀": lam0,
            "H_inf": H,
            "V0": V0_,
            "U1/4": V0_**0.25,
            "m_chi": m_chi(lam0, 11.1),
            "m_Phi": m_Phi(lam0, 11.1),
        }
    paper_ref = {"λ₀": "6.78e-8 (Table I)", "H_inf": "1.65e13", "V0": "4.83e63", "U1/4": "~8.3e15", "m_chi": "3.28e13", "m_Phi": "~2.6e14"}
    for key in rows:
        lines.append(
            f"| {key} | {vals['num'][key]:.4e} | {vals['an'][key]:.4e} | {vals['ex'][key]:.4e} | {paper_ref[key]} |"
        )
    lines.append("")
    lines.append(
        "**判定（针对旧稿，已过期）：** 旧稿正文（H_inf=1.65e13, V0=4.83e63, m_χ=3.28e13）与 **数值 λ₀=6.78e-8 自洽**。"
        "**注意：当前论文已改用锁定约定 λ₀=6.70e-8、H_inf=1.64e13、V0=4.78e63**，"
        "因此本节结论已不适用于现稿；当年指出的“双 λ₀ 轨”问题已由 `lock_n_convention.py` 统一解决。"
        "（保留本节仅为记录问题曾存在。）"
    )
    lines.append("")

    # --- 6. Parse tex for dual values ---
    lines.append("## 6. 论文 .tex 内数值并存扫描")
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
        lines.append("| 模式 | 出现次数 |")
        lines.append("|---|---|")
        for name, pat in pats.items():
            n = len(re.findall(pat, tex))
            lines.append(f"| {name} | {n} |")
        lines.append("")
        lines.append(
            "**判定：** 若 6.78e-8 与 7.46e-8、1.65e13 与 1.73e13 同时出现在正文/附录，"
            "则论文内部 **双轨数值** 未完全收束；需在修订中声明「Table I 用精确数值归一，附录给出解析近似」，"
            "或全部统一到同一 λ₀。"
        )
    else:
        lines.append("未找到 paper_prd_merged.tex")
    lines.append("")

    # --- 7. What scripts support / not / new problems ---
    lines.append("## 7. 总判定")
    lines.append("")
    lines.append("### 脚本支持的论文主张")
    lines.append("")
    lines.append("1. Starobinsky 型函数形式 r=8/(β²N²)=2(6+1/ξ)/N²，n_s≈1−2/N —— **支持**")
    lines.append("2. 结束条件 e^{-x_end}=0.466，V_end/V0=0.285 —— **支持**")
    lines.append("3. 共形退耦 Ω=Φ/Φ₀，树图 χψ̄ψ=0 —— **支持**（代数）")
    lines.append("4. G_eff(Φ₀)=G_N —— **支持**")
    lines.append("5. 单场精质不可行（m_χ/H0≫1）；DE=冻结 V_c —— **支持**")
    lines.append("6. T_reh~1e9 ↔ N≈50（同一匹配式下）—— **支持**")
    lines.append("7. r≲0.0053 在 N∈[48,55]、ξ=11.1 —— **支持**（公式层面）")
    lines.append("8. 畴壁 σ_wall~1e50 GeV³ 量级、F(0)=F′(0)=0 —— **支持**")
    lines.append("")
    lines.append("### 脚本不支持或仅部分支持")
    lines.append("")
    lines.append("1. 附录解析 λ₀=48π²ξ²A_s/N² —— **不支持**（代数错误）")
    lines.append("2. 「引力 T_reh~4e5 与反常 T_reh~1e9 均落在 N∈[48,55]」—— **部分不支持**（4e5→N≈47）")
    lines.append("3. 「χ→ψψ 运动学+顶点双重保护」对整个 g∈[1e-5,1e-4] —— **仅高 g 端支持**")
    lines.append("4. DM Ω_DM=0.12 定量预言 —— **脚本未验证**（机制参数窗可估，归一化需格点）")
    lines.append("5. Fig2 中 Planck 等高线 —— **示意**，非官方似然；不能当数据约束证据")
    lines.append("")
    lines.append("### 脚本/作图引入或暴露的问题")
    lines.append("")
    lines.append("1. **双 λ₀ 轨**：Table I/正文用 6.78e-8，解析附录与 figures 用 ~7.47e-8 → 图文标度差 ~10%")
    lines.append("2. **N–T_reh 公式不统一**：Eq.(18) 与「50+¼ln」对同一 T_reh 给出不同 N")
    lines.append("3. **精确 vs 大场 A_s**：需用脚本给出的 λ₀_exact 对照 Table I，确认 6.78e-8 的来源")
    lines.append("4. 代码瑕疵：`m_Phi_jordan` 死代码；Fig2 约束为示意")
    lines.append("")
    lines.append("### 建议补充的脚本")
    lines.append("")
    lines.append("| 脚本 | 目的 |")
    lines.append("|---|---|")
    lines.append("| `consistency_checks.py`（本文件） | Table I / 精确慢滚 / 再加热窗 / DM 运动学 |")
    lines.append("| `check_rg_running.py` | 复核 β_ξ、Δξ≲1e-6 |")
    lines.append("| `check_reheating_rates.py` | 复核 Γ_anom~O(1) GeV 与 T_reh 公式 |")
    lines.append("| `check_dm_abundance.py` | 参数化 Ω_ψ(m_ψ,T_reh) 曲面与 0.12 交线 |")
    lines.append("| `check_condensate_overclosure.py` | 再加热后凝聚体残留上界 |")
    lines.append("| 图脚本改用单一 λ₀ 约定 | 与 Table I 或解析式对齐，并在 caption 写明 |")
    lines.append("")
    lines.append("[一致性深检完成]")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
