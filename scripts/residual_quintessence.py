#!/usr/bin/env python3
r"""
residual_quintessence.py  (缺失脚本补全 #4  [P0])
================================================
需求原文（评审工作区 review_workspace/ 已在定稿后清理；此处保留摘录）§C:
  "[P0] 残余精质 lesssim1e-10 复算：积分凝聚体 Boltzmann 方程
        （Gamma_anom、H(T)、a^-3 稀释），输出 rho_cond(a_0)/V_c – T_reh 曲线。
        现缺：1e-10 是摘要级卖点却无脚本，且 Delta w_quant 论文值 1e-118
        与自身公式 1.7e-121 差 3 个量级。"

论文声称 (paper_prd_merged.tex；按节引用——行号随改稿漂移，不再硬编码):
  §VI    Delta w_Ricci  ~ (H_0/m_chi)^2        ~ 2e-111
  §VI    Delta w_quant  ~ H_0^4/V_c            ~ 1e-118
  §VI    rho_cond(a_0)/V_c <= Omega_DM/Omega_Lambda ~ 0.39  (BBN 地板 T_reh=1e-2 GeV)
  §VI    "for the fiducial T_reh >= 1e9 GeV ... the residual condensate is
          exponentially negligible"
  §VI    "up to residual-condensate corrections bounded in App. D5
          (the dominant correction Delta w <= 1e-10 for T_reh >= 1 GeV)"

方法
----
积分**精确**的两流体 Boltzmann 系统（宇宙时 -> 改用 e-fold，且对衰变项用指数积分器，
因晚期 Gamma/H 可达 1e30，显式格式会失稳）:

    Y = rho_cond a^3        dY/dN = -(Gamma/H) Y
    R = rho_rad  a^4        dR/dN = +(Gamma/H) Y a
    H^2 = (Y/a^3 + R/a^4) / (3 M_P^2)
    指数积分器:  Y' = Y exp(-Gamma dt),  R' = R + (Y - Y') a_mid

N 从 a_end (=1) 起算到 a_0 (N ~ 70.7)。扫描 Gamma（等价于 T_reh），输出
rho_cond(a_0)/V_c、T_reh、N_reh，并对照论文的三个说法。
另算 Gamma=0（绝对稳定凝聚体）+ 引力脉冲作为初始辐射的对照情形。

输出: scripts/residual_quintessence.md / .json
"""

from __future__ import annotations

import json
import math
import os

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))

M_P = 2.435e18
H0_GEV = 1.4377e-42
T0_GEV = 2.3491e-13
G_STAR = 106.75
G_STAR_S0 = 3.91
OMEGA_LAMBDA = 0.683
OMEGA_DM = 0.265
RHO_C = 3.0 * H0_GEV**2 * M_P**2
V_C = OMEGA_LAMBDA * RHO_C

# 论文锁定值 (N=50, xi=11.1)
V0 = 4.7772e63
V_END = 0.285204 * V0
RHO_END = 1.19938 * V_END
H_INF = 1.6388e13
M_CHI = 3.2533e13
N_EFF = 10.0
RHO_PULSE = N_EFF * H_INF**4 / (192.0 * math.pi**2)

# a_0/a_end : 由 (a_end/a_0)^3 = 1.0204e-92 at T_reh=1e9 反推
A_END_OVER_A0_CUBED_1E9 = 1.0204e-92
N_TO_A0 = -math.log(A_END_OVER_A0_CUBED_1E9) / 3.0     # ~70.74


def T_of_rho_rad(rho: float, g_star: float = G_STAR) -> float:
    return (30.0 * rho / (math.pi**2 * g_star)) ** 0.25


def Gamma_from_Treh(T_reh: float, g_star: float = G_STAR) -> float:
    """T_reh = (90/(pi^2 g_*))^{1/4} sqrt(Gamma M_P)   <=>   H(T_reh) = Gamma."""
    return T_reh**2 / (math.sqrt(90.0 / (math.pi**2 * g_star)) ** 2 * M_P)


def solve_two_fluid(Gamma: float, rho_rad0: float = 0.0,
                    N_end: float = N_TO_A0, dN: float = 0.005) -> dict:
    r"""指数积分器解两流体系统.

    a = e^N (a_end = 1).  Y = rho_cond a^3, R = rho_rad a^4.
    """
    Y = RHO_END * 1.0**3
    R = rho_rad0 * 1.0**4
    N = 0.0
    N_reh = float("nan")
    T_reh = float("nan")
    rho_rad_at_cross = float("nan")
    while N < N_end:
        a = math.exp(N)
        rho_cond = Y / a**3
        rho_rad = R / a**4
        H = math.sqrt(max((rho_cond + rho_rad) / (3.0 * M_P**2), 1e-300))
        if not math.isfinite(N_reh) and rho_rad >= rho_cond and rho_cond > 0:
            N_reh = N
            T_reh = T_of_rho_rad(rho_rad)
            rho_rad_at_cross = rho_rad
        dt = dN / H
        decay = math.exp(-Gamma * dt) if Gamma * dt < 700.0 else 0.0
        Y_new = Y * decay
        a_mid = math.exp(N + 0.5 * dN)
        R = R + (Y - Y_new) * a_mid
        Y = Y_new
        N += dN
    a0 = math.exp(N)
    rho_cond_0 = Y / a0**3
    rho_rad_0 = R / a0**4
    return {"Gamma": Gamma, "N_reh": N_reh, "T_reh": T_reh,
            "rho_cond_a0": rho_cond_0,
            "rho_cond_a0_over_Vc": rho_cond_0 / V_C,
            "rho_rad_a0": rho_rad_0,
            "T_today_from_rad": T_of_rho_rad(rho_rad_0),
            "rho_rad_at_cross": rho_rad_at_cross,
            "N_total": N}


def main() -> None:
    out: dict = {"constants": {"M_P": M_P, "H0_GeV": H0_GEV, "RHO_C": RHO_C,
                               "V_C": V_C, "V_C_paper_2p7e-47": 2.7e-47,
                               "V0": V0, "V_end": V_END, "rho_end": RHO_END,
                               "H_inf": H_INF, "m_chi": M_CHI,
                               "rho_pulse": RHO_PULSE,
                               "N_to_a0": N_TO_A0,
                               "Omega_DM_over_Omega_L": OMEGA_DM / OMEGA_LAMBDA}}

    # ---- [1] Delta w 三项 ----
    dw_ricci = (H0_GEV / M_CHI) ** 2
    dw_quant = H0_GEV**4 / V_C
    dw_quant_paperVc = H0_GEV**4 / 2.7e-47
    out["delta_w"] = {
        "Ricci_recomputed": dw_ricci,
        "Ricci_paper": 2e-111,
        "quant_recomputed_H0^4/Vc": dw_quant,
        "quant_with_paper_Vc": dw_quant_paperVc,
        "quant_paper": 1e-118,
        "quant_ratio_paper_over_recomputed": 1e-118 / dw_quant,
        "overclosure_bound_OmegaDM_over_OmegaL": OMEGA_DM / OMEGA_LAMBDA,
    }
    # 残余位移 delta Phi/PhiV ~ (H_0/m_chi)^2 ? 论文 §VII / App. D5 说 1e-111
    out["delta_w"]["displacement_sq"] = dw_ricci   # (H0/m_chi)^2 ~ 2e-111
    # A/A: ~ (H0/m_chi)^2 * R0/(m_chi^2) ... 论文给 6.3e-111；仅记录论文值供对照
    out["delta_w"]["dA_over_A_paper"] = 6.3e-111

    # ---- [2] 两流体扫描: Gamma = Gamma_anom(T_reh) ----
    rows = []
    for T in (1e-2, 1e0, 1e3, 1e6, 1e9, 1e12, 1e15):
        G = Gamma_from_Treh(T)
        r = solve_two_fluid(G)
        r["T_reh_input"] = T
        r["Gamma_over_H_reh"] = G / math.sqrt(
            (math.pi**2 / 30.0) * G_STAR * T**4 / (3.0 * M_P**2))
        rows.append(r)
    out["scan_with_decay"] = rows

    # ---- [3] 对照: Gamma = 0 (绝对稳定凝聚体) + 引力脉冲 ----
    r0_nopulse = solve_two_fluid(0.0, rho_rad0=0.0)
    r0_pulse = solve_two_fluid(0.0, rho_rad0=RHO_PULSE)
    out["stable_condensate"] = {"no_pulse": r0_nopulse, "with_pulse": r0_pulse}

    # ---- [3b] 脉冲通道的标度检验: r = rho_rad/rho_cond 随 a 是增长还是衰减 ----
    # 解析: Gamma=0 时 Y,R 均守恒 -> r = (R/Y)/a  ∝ a^{-1}, 即 d ln r/dN = -1 (衰减!)
    scal = []
    N_test = np.linspace(0.0, 20.0, 201)
    for N in N_test[::40]:
        a = math.exp(N)
        rho_c = RHO_END / a**3
        rho_r = RHO_PULSE / a**4
        scal.append({"N": float(N), "a": a, "rho_rad_over_rho_cond": rho_r / rho_c})
    out["pulse_scaling"] = scal
    if len(scal) >= 2:
        s = (math.log(scal[1]["rho_rad_over_rho_cond"])
             - math.log(scal[0]["rho_rad_over_rho_cond"])) / (scal[1]["N"] - scal[0]["N"])
        out["pulse_scaling_slope_dlnr_dN"] = s
    # 若按论文 §IV 的错误标度 (r ∝ a^{+1}) 反推的"再加热"温度
    a_ratio_paper = RHO_END / RHO_PULSE
    rho_rad_at_paper_cross = RHO_PULSE * a_ratio_paper**4 / a_ratio_paper**4 * a_ratio_paper
    out["paper_pulse_crossing"] = {
        "a_reh_over_a_end_claimed": a_ratio_paper,
        "T_reh_claimed_by_same_arithmetic":
            T_of_rho_rad(RHO_END / a_ratio_paper),   # 论文用的算错标度下的结果
    }
    # 正确标度下: 辐射永远追不上; 今天的辐射温度
    a0_over_aend = math.exp(N_TO_A0)
    rho_rad_today = RHO_PULSE / a0_over_aend**4
    out["pulse_correct_today"] = {
        "a0_over_aend": a0_over_aend,
        "rho_rad_today": rho_rad_today,
        "T_today_GeV": T_of_rho_rad(rho_rad_today),
        "T_today_eV": T_of_rho_rad(rho_rad_today) * 1e9,
    }

    # ---- [4] Gamma t_0 判据: 何时衰变才算"完全" ----
    T0_S = 4.35e17
    rows2 = []
    for T in (1e-2, 1e0, 1e9):
        G = Gamma_from_Treh(T)
        G_s = G / 6.582119569e-25
        rows2.append({"T_reh": T, "Gamma_GeV": G, "Gamma_s^-1": G_s,
                      "Gamma_t0": G_s * T0_S})
    out["decay_completeness"] = rows2

    with open(os.path.join(ROOT, "residual_quintessence.json"), "w",
              encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=float)

    # ---------------- 报告 ----------------
    L = []
    A = L.append
    A("# 残余精质复算：凝聚体 Boltzmann 积分与 $\\Delta w$ 三项（缺失脚本补全 #4）\n")
    A("对象：`paper_prd_merged.tex` §VI（$\\Delta w\\leq10^{-10}$、Ricci、quantum、")
    A("$\\leq0.39$ 与\"指数可忽略\"）与 App. D5a–c（逐项推导）。\n")

    c = out["constants"]
    A("## 1. 输入常量\n")
    A(f"$\\rho_c={c['RHO_C']:.4e}$ GeV$^4$，$V_c=\\Omega_\\Lambda\\rho_c={c['V_C']:.4e}$ GeV$^4$")
    A(f"（论文 §VI 用 $2.7\\times10^{{-47}}$，对应 $H_0\\simeq70$ km/s/Mpc），")
    A(f"$\\rho_{{\\rm end}}={c['rho_end']:.4e}$ GeV$^4$，$\\rho_{{\\rm pulse}}={c['rho_pulse']:.4e}$ GeV$^4$，")
    A(f"$\\Omega_{{\\rm DM}}/\\Omega_\\Lambda={c['Omega_DM_over_Omega_L']:.4f}$。")
    A(f"积分区间 $N\\in[0,{c['N_to_a0']:.2f}]$（由 $(a_{{\\rm end}}/a_0)^3=1.0204\\times10^{{-92}}$ 反推）。\n")

    d = out["delta_w"]
    A("## 2. $\\Delta w$ 三项独立复算\n")
    A("| 项 | 论文值 | 复算值 | 判定 |")
    A("|---|---|---|---|")
    A(f"| $\\Delta w_{{\\rm Ricci}}=(H_0/m_\\chi)^2$ | $2\\times10^{{-111}}$ | ${d['Ricci_recomputed']:.4e}$ | ✓ |")
    A(f"| $\\Delta w_{{\\rm quant}}=H_0^4/V_c$ | $10^{{-118}}$ | ${d['quant_recomputed_H0^4/Vc']:.4e}$ | ✗ 差 {d['quant_ratio_paper_over_recomputed']:.0f} 倍 |")
    A(f"| （若用论文 §VI 的 $V_c=2.7\\times10^{{-47}}$） | — | ${d['quant_with_paper_Vc']:.4e}$ | 仍非 $10^{{-118}}$ |")
    A("")
    A("**结论**：$\\Delta w_{{\\rm quant}}$ 的论文值 $10^{-118}$ 比自身公式 $H_0^4/V_c$ 的结果")
    A(f"大 $\\sim{d['quant_ratio_paper_over_recomputed']:.0f}$ 倍。论文方向是**保守**的")
    A("（真值更小），故结论不受影响，但数字须改。")
    A(f"论文 §VI 与 App. D5b 的 $10^{{-118}}$ 应改为 $\\lesssim2\\times10^{{-121}}$。\n")

    A("## 3. 凝聚体残余：精确两流体积分\n")
    A("积分 $dY/dN=-(\\Gamma/H)Y$（$Y=\\rho_{{\\rm cond}}a^3$）与 $dR/dN=+(\\Gamma/H)Ya$")
    A("（$R=\\rho_{{\\rm rad}}a^4$），衰变项用指数积分器（晚期 $\\Gamma/H$ 可达 $10^{{30}}$）。\n")
    A("| $T_{{\\rm reh}}$ 输入 [GeV] | $\\Gamma$ [GeV] | $T_{{\\rm reh}}$ 实算 | $\\rho_{{\\rm cond}}(a_0)$ | $\\rho_{{\\rm cond}}(a_0)/V_c$ |")
    A("|---|---|---|---|---|")
    for r in rows:
        A(f"| {r['T_reh_input']:.0e} | {r['Gamma']:.4e} | {r['T_reh']:.4e} | "
          f"{r['rho_cond_a0']:.4e} | {r['rho_cond_a0_over_Vc']:.4e} |")
    A("")
    s0 = out["stable_condensate"]
    A("**对照：绝对稳定凝聚体（$\\Gamma=0$）**\n")
    A("| 情形 | $T_{{\\rm reh}}$ | $\\rho_{{\\rm cond}}(a_0)$ | $\\rho_{{\\rm cond}}(a_0)/V_c$ | 超闭合？ |")
    A("|---|---|---|---|---|")
    for key, lab in (("no_pulse", "$\\Gamma=0$，无初始辐射"),
                     ("with_pulse", "$\\Gamma=0$，含引力脉冲")):
        r = s0[key]
        A(f"| {lab} | {r['T_reh']:.4e} | {r['rho_cond_a0']:.4e} | "
          f"{r['rho_cond_a0_over_Vc']:.4e} | "
          f"{'**是**' if r['rho_cond_a0_over_Vc'] > c['Omega_DM_over_Omega_L'] else '否'} |")
    A("")
    A("### 判定\n")
    A("1. **只要 $\\Gamma\\neq0$，残余就是可忽略的**：表中 $T_{{\\rm reh}}$ 从 $10^{-2}$ 到 $10^{15}$ GeV")
    A("一路 $\\rho_{{\\rm cond}}(a_0)/V_c$ 都是 0（下溢）。原因是衰变完整性判据：")
    for r in rows2:
        A(f"   - $T_{{\\rm reh}}={r['T_reh']:.0e}$ GeV $\\Rightarrow \\Gamma={r['Gamma_GeV']:.3e}$ GeV "
          f"$={r['Gamma_s^-1']:.3e}$ s$^{{-1}}$，$\\Gamma t_0={r['Gamma_t0']:.3e}\\gg1$。")
    A("2. **但论文 §VI 的\"BBN 地板给 $\\leq0.39$\"是约束、不是复算结果**：")
    A("   绝对稳定凝聚体（$\\Gamma=0$）的实算残余是上表第二/三行，")
    _ovc = s0["no_pulse"]["rho_cond_a0_over_Vc"] / c["Omega_DM_over_Omega_L"]
    A(f"   比 $\\Omega_{{\\rm DM}}/\\Omega_\\Lambda={c['Omega_DM_over_Omega_L']:.3f}$ 大 "
      f"**$\\sim10^{{{math.log10(_ovc):.0f}}}$ 倍**（$={_ovc:.2e}$）。")
    A("   即：$0.39$ 是\"凝聚体不得超过暗物质预算\"这一**要求**，而凝聚体若真稳定则**违反**该要求。")
    A("   满足它的唯一原因是凝聚体衰变了（$\\Gamma t_0\\gg1$），这一点论文在 §VI 后半句说了，")
    A("   但前半句的写法读起来像\"由约束导出了 $0.39$\"。建议改写为：")
    A(f"   \"若凝聚体绝对稳定则会超闭合 $\\sim10^{{{math.log10(_ovc):.0f}}}$ 倍；由于 $\\Gamma t_0\\gg1$（对任何 $T_{{\\rm reh}}\\gtrsim10^{{-2}}$ GeV），")
    A("   残余在 $a_0$ 时严格为零。\"")
    A("")
    A(f"3. 论文 §VI 的 $\\Delta w\\leq10^{{-10}}$（$T_{{\\rm reh}}\\geq1$ GeV）**成立**，")
    A("   且实际远强于该界（实际为 0）。\n")

    A("## 4. 【新发现，P0 级】脉冲通道**不能**再加热\n")
    A("论文 §IV 写：\"This radiation redshifts as $a^{-4}$ during matter-domination; ")
    A("the ratio $\\rho_{\\rm rad}/\\rho_\\chi\\propto a^{-1}$ **grows** until ")
    A("$\\rho_{\\rm rad}\\sim\\rho_\\chi$ at $a_{\\rm reh}/a_{\\rm end}\\sim10^{13}$\"。\n")
    A("**标度符号错了。** 在凝聚体主导（$w=0$，$\\rho_{\\rm cond}\\propto a^{-3}$）期间，")
    A("辐射以 $a^{-4}$ 稀释得更快，故\n")
    A("$$\\frac{\\rho_{\\rm rad}}{\\rho_{\\rm cond}}=\\frac{R/a^4}{Y/a^3}=\\frac{R}{Y}\\frac1a\\propto a^{-1}"
      "\\quad\\text{（衰减，不是增长）}.$$\n")
    sl = out.get("pulse_scaling_slope_dlnr_dN")
    A(f"数值检验：$d\\ln(\\rho_{{\\rm rad}}/\\rho_{{\\rm cond}})/dN={sl:.4f}$（解析值 $-1$）。\n")
    A("| $N$ | $a/a_{\\rm end}$ | $\\rho_{\\rm rad}/\\rho_{\\rm cond}$ |")
    A("|---|---|---|")
    for r in out["pulse_scaling"][:6]:
        A(f"| {r['N']:.1f} | {r['a']:.3e} | {r['rho_rad_over_rho_cond']:.4e} |")
    A("")
    A("后果：**引力脉冲单独存在时，辐射永远追不上凝聚体**，宇宙保持物质主导，")
    A("不存在辐射主导期。论文按错误标度反推的 \"$T_{{\\rm reh}}\\sim4\\times10^5$ GeV\" 因此不成立：")
    pc = out["pulse_correct_today"]
    A(f"按正确标度，今天（$a_0/a_{{\\rm end}}={pc['a0_over_aend']:.3e}$）的辐射温度只有 "
      f"$T_0={pc['T_today_GeV']:.3e}$ GeV $={pc['T_today_eV']:.3e}$ eV，")
    A("即宇宙几乎是空的。")
    A("**结论：脉冲通道（i）不是再加热通道**；再加热必须由凝聚体的真实衰变（iv）提供。")
    A("论文 §IV 的 \"In isolation this pulse would reheat the universe to "
      "$T_{\\rm reh}\\sim4\\times10^5$ GeV\" 与相应整段须重写。")
    A("（这同时使 §IV \"The parametric estimate below accordingly combines channels (i) and (iv)\"")
    A(" 的措辞变得不必要——(i) 对 $T_{\\rm reh}$ 没有贡献。）\n")

    md = "\n".join(L) + "\n"
    with open(os.path.join(ROOT, "residual_quintessence.md"), "w",
              encoding="utf-8") as f:
        f.write(md)
    print(f"wrote scripts/residual_quintessence.md ({len(md)} chars)")


if __name__ == "__main__":
    main()
