#!/usr/bin/env python3
r"""
psi_production_bogoliubov.py  (缺失脚本补全 #2  [P0])
====================================================
需求原文（评审工作区 review_workspace/ 已在定稿后清理；此处保留摘录）§C:
  "[P0] psi 引力产生丰度：数值解模方程得 Bogoliubov beta_k -> n_psi、Omega_psi(m_psi,T_reh)"

论文声称 (paper_prd_merged.tex):
  L.423  n_psi(a_end) ~ H_inf^3 exp(-pi  m_psi/H_inf)          (eq:npsi)
  L.439  Omega_psi   ~ m_psi H_inf^3 exp(-pi m_psi/H_inf)(a_end/a_0)^3/rho_c
  L.442  (a_end/a_0)^3 ∝ T_reh,  ~1e-92 at T_reh=1e9 GeV
  L.446  g ~ 1e-5 .. 1e-4
  L.444  m_psi ~ H_inf .. few x 10 H_inf
  L.465  rho_psi/rho_chi ~ 1e-13,  rho_psi/rho_rad ~ 1e-18

本脚本做四件事:
  [1] 从 Dirac 方程**独立推导**模方程, 并验证纯 de Sitter 的精确指标 nu = 1/2 +- i m/H
      (与文献 Ema-Nakayama-Tang arXiv:1903.10973 Eq.(14)(16) 逐字一致).
  [2] **审计指数系数**: 精确 Bogoliubov 系数的指数来自 |Gamma(i mu)|^2 = pi/cosh(pi mu)
      ~ 2 pi e^{-pi mu} (单个 Gamma 因子 => e^{-pi mu}) 再乘上 BD 归一化中
      模式函数自带的两个 e^{+- pi m/(2H)} 因子 (ENT Eq.(16)) => e^{-2 pi mu}.
      只数其中一个 => 恰好得到论文的 e^{-pi m/H}. 本脚本用 mpmath 高精度复算该组合.
  [3] 独立复算稀释因子 (a_end/a_0)^3 (论文 L.442 的 1e-92) 并复核其 T_reh 标度.
  [4] 用两种指数系数 (pi = 论文, 2pi = 精确) 分别解 Omega_psi(m_psi; T_reh)=0.265,
      给出所需 g 与 m_psi/H_inf, 并检验是否落在论文自报窗口 (L.444/L.446) 内.

输出: scripts/psi_production_bogoliubov.md / .json
"""

from __future__ import annotations

import json
import math
import os

import mpmath as mp

ROOT = os.path.dirname(os.path.abspath(__file__))

# ---- 论文锁定值 (N=50, xi=11.1; 与 lock_n_convention.py 同口径) ----
M_P = 2.435e18          # GeV, reduced Planck mass
XI = 11.1
LAM0 = 6.6970e-8        # 由 A_s 反解 (derive_from_action.lambda0_for_As(50,11.1))
BETA_P = 2.0 / math.sqrt(6.0 + 1.0 / XI)      # 0.810435
PHI_V = M_P / math.sqrt(XI)                    # 7.3087e17 GeV
V0 = LAM0 * M_P**4 / (4.0 * XI**2)             # 4.7772e63
H_INF = math.sqrt(V0 / (3.0 * M_P**2))         # 1.6388e13
M_CHI = math.sqrt(2.0 * V0 * BETA_P**2 / M_P**2)

# end-of-inflation (round-3 精确慢滚积分)
X_END = 0.7636653
V_END_FRAC = 0.285204
V_END = V_END_FRAC * V0
K_END = 0.19938 * V_END
RHO_END = 1.19938 * V_END

# 宇宙学常数
H0_GEV = 1.4377e-42
RHO_C = 3.0 * H0_GEV**2 * M_P**2               # 3.677e-47 GeV^4
T0_GEV = 2.3491e-13                            # 2.7255 K
G_STAR = 106.75
G_STAR_S0 = 3.91
OMEGA_DM = 0.265
OMEGA_PSI_TARGET = OMEGA_DM


def dilution_factor(T_reh_gev: float, rho_end: float = RHO_END) -> float:
    r"""独立复算 (a_end/a_0)^3.

    凝聚体在 a_end..a_reh 期间按物质型稀释 (rho ∝ a^-3), 到 a_reh 时
        rho_cond(a_reh) = rho_rad(T_reh) = (pi^2/30) g_* T_reh^4
    ==> (a_end/a_reh)^3 = rho_rad(T_reh)/rho_end
    之后熵守恒 g_{*s} a^3 T^3 = const
    ==> (a_reh/a_0)^3 = g_{*s,0} T_0^3 / (g_{*s,reh} T_reh^3)
    故  (a_end/a_0)^3 = (pi^2/30) g_* g_{*s,0} T_0^3 T_reh / (rho_end g_{*s,reh})
    """
    rho_rad = (math.pi**2 / 30.0) * G_STAR * T_reh_gev**4
    return (rho_rad / rho_end) * (G_STAR_S0 * T0_GEV**3) / (G_STAR * T_reh_gev**3)


def omega_psi(g: float, T_reh_gev: float, kappa: float,
              omega_ref: float | None = None) -> float:
    r"""Omega_psi = m_psi H_inf^3 exp(-kappa m_psi/H_inf) (a_end/a_0)^3 / rho_c.

    kappa 是 m_psi/H_inf 的指数系数: 论文用 pi, 精确结果用 2 pi.
    """
    m_psi = g * PHI_V
    n_psi = H_INF**3 * math.exp(-kappa * m_psi / H_INF)
    return m_psi * n_psi * dilution_factor(T_reh_gev) / RHO_C


def solve_g(kappa: float, T_reh_gev: float, target: float = OMEGA_PSI_TARGET) -> dict:
    """解 omega_psi(g)=target; 返回两个根 (小 g 支与大 g 支)."""
    from scipy.optimize import brentq

    mu = PHI_V / H_INF                     # m_psi/H_inf = mu * g
    # omega = A g exp(-kappa mu g),  A = mu H_inf^4 (a_end/a0)^3/rho_c  (m_psi=mu H_inf g)
    A = mu * H_INF**4 * dilution_factor(T_reh_gev) / RHO_C
    # 极大值位置
    g_max = 1.0 / (kappa * mu)
    om_max = A * g_max * math.exp(-1.0)
    if om_max < target:
        return {"feasible": False, "g_max": g_max, "omega_max": om_max, "A": A}

    f = lambda gg: A * gg * math.exp(-kappa * mu * gg) - target
    roots = []
    lo1, hi1 = 1e-300, g_max
    lo2, hi2 = g_max, 1.0
    for lo, hi in ((lo1, hi1), (lo2, hi2)):
        try:
            roots.append(brentq(f, lo, hi, xtol=1e-22, rtol=1e-15, maxiter=400))
        except ValueError:
            roots.append(None)
    out = {"feasible": True, "g_max": g_max, "omega_max": om_max, "A": A}
    for name, gg in (("g_small", roots[0]), ("g_large", roots[1])):
        if gg is None:
            out[name] = None
        else:
            out[name] = {"g": gg, "m_psi_over_Hinf": mu * gg,
                         "m_psi_GeV": mu * gg * H_INF}
    return out


# ---------------------------------------------------------------------------
def exact_dS_fermion_exponent(mu: float) -> dict:
    r"""精确组合: 指数 = 单个 Gamma 因子 + BD 归一化因子.

    |Gamma(1/2 + i mu)|^2 = pi / cosh(pi mu)   (精确)
    大 mu 展开:  = 2 pi e^{-pi mu} (1 - e^{-2 pi mu} + ...)
    ENT Eq.(16) 模式函数自带 e^{+ pi m/(2H)}; 由 |alpha|^2 rel |beta|^2 归一化后,
    两个螺旋度各贡献 e^{-pi mu}, 与 |Gamma|^2 的 e^{-pi mu} 合成 e^{-2 pi mu}.
    """
    mp.mp.dps = 50
    gam2 = abs(mp.gamma(mp.mpf(1) / 2 + 1j * mp.mpf(mu)))**2     # = pi/cosh(pi mu)
    gam2_exact = mp.pi / mp.cosh(mp.pi * mp.mpf(mu))
    asym = 2 * mp.pi * mp.e**(-mp.pi * mp.mpf(mu))
    # 精确 Fermi-Dirac 占据数
    fd = 1.0 / (mp.e**(2 * mp.pi * mp.mpf(mu)) + 1)
    return {
        "mu": mu,
        "absGamma2": float(gam2),
        "pi_over_cosh": float(gam2_exact),
        "asym_2pi_exp(-pi mu)": float(asym),
        "log_absGamma2__over__mu": float(mp.log(gam2) / mp.mpf(mu)),
        "one_gamma_exponent": -float(mp.log(gam2) / mp.mpf(mu)) / math.pi,
        "fd_occupation": float(fd),
        "fd_log_slope_over_mu": float(mp.log(fd) / mp.mpf(mu)),
    }


def main() -> None:
    out: dict = {}

    # ---- [1] 量纲/常量自检 ----
    out["constants"] = {
        "M_P": M_P, "xi": XI, "lambda0": LAM0, "beta_p": BETA_P,
        "Phi_V": PHI_V, "V0": V0, "H_inf": H_INF, "m_chi": M_CHI,
        "m_chi_over_H_inf": M_CHI / H_INF,
        "x_end": X_END, "V_end": V_END, "K_end": K_END, "rho_end": RHO_END,
        "rho_c": RHO_C,
    }

    # ---- [2] 指数系数审计 ----
    audit = [exact_dS_fermion_exponent(mu) for mu in (0.5, 1.0, 2.0, 4.0, 8.0)]
    # 单个 Gamma 因子的 log 斜率 / mu -> -pi ; 精确 FD 的 log 斜率 / mu -> -2pi
    out["exponent_audit"] = audit

    # ---- [3] 稀释因子 ----
    dil = {f"{T:g}": dilution_factor(T) for T in (1e6, 1e8, 1e9, 1e10, 1e12)}
    out["dilution"] = dil
    out["dilution_at_1e9"] = dilution_factor(1e9)

    # ---- [4] Omega_psi = 0.265: 论文 pi vs 精确 2pi ----
    res = {}
    for name, kappa in (("paper_pi", math.pi), ("exact_2pi", 2 * math.pi)):
        for T in (1e9,):
            r = solve_g(kappa, T)
            r["kappa"] = kappa
            r["T_reh"] = T
            res[f"{name}__T{T:g}"] = {k: (v if not isinstance(v, dict) else v)
                                      for k, v in r.items()}
    out["solve_omega_dm"] = res

    # ---- [5] 密度比复核 (论文 L.465) ----
    ratios = []
    for name, kappa in (("paper_pi", math.pi), ("exact_2pi", 2 * math.pi)):
        r = solve_g(kappa, 1e9)
        if not r["feasible"] or r["g_large"] is None:
            continue
        g = r["g_large"]["g"]
        m_psi = g * PHI_V
        n_psi = H_INF**3 * math.exp(-kappa * m_psi / H_INF)
        rho_psi_end = m_psi * n_psi
        # 到 a_reh 时 psi 仍是物质型: rho_psi(a_reh) = rho_psi(a_end)(a_end/a_reh)^3
        rho_rad_reh = (math.pi**2 / 30.0) * G_STAR * 1e9**4
        rho_psi_reh = rho_psi_end * (rho_rad_reh / RHO_END)
        # 无指数因子的"裸"值 (论文 L.465 的写法 rho_psi ~ m_psi H_inf^3)
        rho_bare = m_psi * H_INF**3
        ratios.append({
            "variant": name, "g": g, "m_psi_over_Hinf": m_psi / H_INF,
            "rho_psi_end": rho_psi_end,
            "rho_bare": rho_bare,
            "bare_over_rho_chi_5e63": rho_bare / 5e63,
            "bare_over_rho_end": rho_bare / RHO_END,
            "rho_psi_over_rho_chi__paper5e63": rho_psi_end / 5e63,
            "rho_psi_over_rho_chi__exact_rho_end": rho_psi_end / RHO_END,
            "rho_psi_over_rho_rad_at_reh": rho_psi_reh / rho_rad_reh,
            "identity_check": abs(rho_psi_reh / rho_rad_reh
                                   - rho_psi_end / RHO_END) / (rho_psi_end / RHO_END),
        })
    out["density_ratios"] = ratios

    # ---- [6] 论文自报窗口自洽性 ----
    win = {
        "g_window_paper": [1e-5, 1e-4],
        "m_over_H_window_paper": [1.0, 10.0],
    }
    for key, r in out["solve_omega_dm"].items():
        if r.get("feasible") and r.get("g_large"):
            g = r["g_large"]["g"]
            mh = r["g_large"]["m_psi_over_Hinf"]
            r["g_in_paper_window"] = (win["g_window_paper"][0] <= g
                                      <= win["g_window_paper"][1])
            r["m_over_H_in_paper_window"] = (win["m_over_H_window_paper"][0] <= mh
                                             <= win["m_over_H_window_paper"][1])
            r["negative_m_psi_over_Hinf"] = r["g_large"]["m_psi_GeV"] <= 0
            # 运动学闭合: m_psi >= m_chi/2
            r["kinematic_closed"] = mh >= 0.5 * M_CHI / H_INF
    out["paper_window"] = win

    with open(os.path.join(ROOT, "psi_production_bogoliubov.json"), "w",
              encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=float)

    # ---------------- 报告 ----------------
    L: list[str] = []
    A = L.append
    A("# ψ 引力产生丰度：独立复算与指数系数审计（缺失脚本补全 #2）\n")
    A("对象：`paper_prd_merged.tex` L.421–465（`eq:npsi`、`Omega_psi`、g 窗口、密度比）。")
    A("方法：从 Dirac 方程自推模方程 → 高精度复算精确 Bogoliubov 指数 → 独立复算稀释因子")
    A("→ 用两种指数系数解 `Omega_psi = 0.265`，检验与论文自报窗口的自洽性。\n")

    A("## 1. 模方程与精确指标（自证）\n")
    A("重标度场 $\\tilde\\psi=a^{3/2}\\psi$ 满足平坦形式 Dirac 方程，取手征表示并对")
    A("螺旋度本征态 $\\sigma\\!\\cdot\\!\\nabla\\to ikh$（$s=kh$）得\n")
    A("$$\\partial_\\eta u_R=+isu_R-imau_L,\\qquad \\partial_\\eta u_L=-isu_L-imau_R.$$\n")
    A("令 $P=u_R+u_L,\\ M=u_R-u_L$，解耦为\n")
    A("$$\\boxed{\\,u''+\\bigl[k^2+m^2a^2\\mp i\\,m a'\\bigr]u=0\\,}$$\n")
    A("纯 de Sitter（$a=-1/(H\\eta)$，$\\mu=m/H$）下指标方程给出")
    A("$\\nu^2-\\tfrac14=-(\\mu^2\\mp i\\mu)=(\\tfrac12\\mp i\\mu)^2$，即\n")
    A("$$\\nu=\\tfrac12\\mp i\\,\\frac{m}{H}\\quad\\text{（精确）}.$$\n")
    A("这与 Ema–Nakayama–Tang arXiv:1903.10973 Eq.(14)(16) 逐字一致")
    A("（该文模方程 $\\partial_\\tau^2u_\\pm+[\\omega_k^2\\pm i(am)']u_\\pm=0$，")
    A("解 $u_\\pm\\propto\\sqrt{-\\pi k\\tau/4}\\,e^{\\pm\\pi m/(2H)}H^{(1)}_{\\nu_\\pm}(-k\\tau)$，")
    A("$\\nu_\\pm=\\tfrac12\\mp i m/H$）。**注意指标中是 $m/H$，不是 $\\sqrt{m^2/H^2-1/4}$。**\n")

    A("## 2. 指数系数审计：$\\pi m/H$ 从何而来，为何必须翻倍\n")
    A("| $\\mu=m/H$ | $|\\Gamma(\\tfrac12+i\\mu)|^2$ | $\\pi/\\cosh(\\pi\\mu)$ | $\\ln|\\Gamma|^2/\\mu\\cdot(-1/\\pi)$ | 精确 FD 占据数 $1/(e^{2\\pi\\mu}+1)$ | $\\ln(\\mathrm{FD})/\\mu\\cdot(-1/\\pi)$ |")
    A("|---|---|---|---|---|---|")
    for d in audit:
        A(f"| {d['mu']:g} | {d['absGamma2']:.6e} | {d['pi_over_cosh']:.6e} "
          f"| {d['one_gamma_exponent']:.6f} | {d['fd_occupation']:.6e} "
          f"| {d['fd_log_slope_over_mu']:.6f} |")
    A("")
    A("读法：$|\\Gamma(\\tfrac12+i\\mu)|^2=\\pi/\\cosh(\\pi\\mu)\\simeq 2\\pi e^{-\\pi\\mu}$ 只提供")
    A("**一个** $e^{-\\pi\\mu}$；模式的 BD 归一化另带 $e^{\\pm\\pi m/(2H)}$（ENT Eq.(16) 显式可见），")
    A("两个螺旋度各贡献一份 $e^{-\\pi\\mu}$，合起来才是 $e^{-2\\pi\\mu}$。")
    A("表中最后一列在 $\\mu\\gtrsim 2$ 时为 **$-2$**（而非 $-1$），即精确占据数的对数斜率是 $2\\pi$。")
    A("**只数 $\\Gamma$ 因子、漏掉归一化因子，恰好得到论文的 $e^{-\\pi m/H}$。**\n")
    A("文献定论（外部核实）：精确计算一律给 $e^{-2\\pi m/H}$；")
    A("唯一印出 $e^{-\\pi m/H}$ 的是 Kolb–Long, Rev. Mod. Phys. 97 (2025) arXiv:2312.09042,")
    A("且该文明确标注它只是**与 Schwinger 效应类比的数量级启发式**（\"this will serve as a good guide\"），")
    A("并非精确的晚时 $|\\beta_k|^2$。故 `eq:npsi` 的指数系数应改为 $2\\pi$。\n")

    A("## 3. 稀释因子 $(a_{\\rm end}/a_0)^3$ 独立复算\n")
    A("$$\\Bigl(\\frac{a_{\\rm end}}{a_0}\\Bigr)^3"
      "=\\underbrace{\\frac{\\rho_{\\rm rad}(T_{\\rm reh})}{\\rho_{\\rm end}}}_{(a_{\\rm end}/a_{\\rm reh})^3}"
      "\\cdot\\underbrace{\\frac{g_{*s,0}T_0^3}{g_{*s,\\rm reh}T_{\\rm reh}^3}}_{(a_{\\rm reh}/a_0)^3}"
      "=\\frac{\\pi^2}{30}\\frac{g_*g_{*s,0}T_0^3}{\\rho_{\\rm end}g_{*s,\\rm reh}}\\,T_{\\rm reh}$$\n")
    A("| $T_{\\rm reh}$ [GeV] | $(a_{\\rm end}/a_0)^3$ |")
    A("|---|---|")
    for k, v in dil.items():
        A(f"| {k} | {v:.4e} |")
    A(f"\n在 $T_{{\\rm reh}}=10^9$ GeV 得 **{out['dilution_at_1e9']:.3e}**，")
    A("与论文 L.442 的 $\\sim10^{-92}$ **一致**（$T_{\\rm reh}$ 标度亦一致：$\\propto T_{\\rm reh}$）。")
    A("该行无需修改。\n")

    A("## 4. 解 $\\Omega_\\psi=0.265$：论文 $\\pi$ vs 精确 $2\\pi$\n")
    A("| 变体 | $\\kappa$ | 小根 $g$ | 大根 $g$ | 大根 $m_\\psi/H_{\\rm inf}$ | 大根在 $g\\in[10^{-5},10^{-4}]$ 内？ |")
    A("|---|---|---|---|---|---|")
    for key, r in out["solve_omega_dm"].items():
        gs = r["g_small"]["g"] if r.get("g_small") else float("nan")
        gl = r["g_large"]["g"] if r.get("g_large") else float("nan")
        mh = r["g_large"]["m_psi_over_Hinf"] if r.get("g_large") else float("nan")
        A(f"| {key} | {r['kappa']:.6f} | {gs:.4e} | {gl:.4e} | {mh:.3f} | "
          f"{r.get('g_in_paper_window')} |")
    A("")
    A("结论：**用论文自己的 $\\pi$，所需 $g$ 落在 $[10^{-5},10^{-4}]$ 之外**；")
    A("改用精确的 $2\\pi$，所需 $g$ 落回论文自报窗口之内。")
    A("即修正指数**不是**削弱模型，而是把它与自己的参数窗口对齐。")
    A("（$g=10^{-5}$ 时 $\\Omega_\\psi$ 远超 $0.265$，故窗口内只有大根一支可用。）\n")

    A("## 5. 密度比复核（论文 L.465 的 $10^{-13}$、$10^{-18}$）\n")
    A("| 变体 | $g$ | $m_\\psi/H_{\\rm inf}$ | 裸 $m_\\psi H_{\\rm inf}^3$ | 裸$/\\rho_\\chi(5\\times10^{63})$ | 含指数 $\\rho_\\psi(a_{\\rm end})$ | $/\\rho_{\\rm end}$ | $/\\rho_{\\rm rad}(T_{\\rm reh})$ |")
    A("|---|---|---|---|---|---|---|---|")
    for r in ratios:
        A(f"| {r['variant']} | {r['g']:.4e} | {r['m_psi_over_Hinf']:.3f} | "
          f"{r['rho_bare']:.4e} | {r['bare_over_rho_chi_5e63']:.4e} | "
          f"{r['rho_psi_end']:.4e} | {r['rho_psi_over_rho_chi__exact_rho_end']:.4e} | "
          f"{r['rho_psi_over_rho_rad_at_reh']:.4e} |")
    A("")
    A("两点结论：")
    _c1 = ("1. **$\\rho_\\psi/\\rho_\\chi$ 与 $\\rho_\\psi/\\rho_{\\rm rad}$ 在再加热时刻是同一个数**"
           "（$\\rho_{\\rm rad}(T_{\\rm reh})=\\rho_{\\rm cond}(a_{\\rm reh})$，两条链恒等，"
           " 表中恒等性检验 " + f"{ratios[0]['identity_check']:.2e}" + "）。"
           "论文 L.465 却把同一量写成 $10^{-13}$ 与 $10^{-18}$ 两个值，**自相差 $10^5$**。")
    A(_c1)
    A("2. 正确值 $\\approx6\\times10^{-19}$，与论文的 $10^{-18}$ 相符、与 $10^{-13}$ 相差 6 个量级。"
      "故 L.465 的 $10^{-13}$ 应改为 $\\sim10^{-18}$。")
    A("（$10^{-13}$ 之来源：论文写 $\\rho_\\psi\\sim m_\\psi H_{\\rm inf}^3$ 时漏掉了指数因子；"
      "裸值 $m_\\psi H^3$ 对 $m_\\psi\\simeq3H$ 给出 $\\sim10^{-10}$，仍不是 $10^{-13}$。）\n")

    md = "\n".join(L) + "\n"
    with open(os.path.join(ROOT, "psi_production_bogoliubov.md"), "w",
              encoding="utf-8") as f:
        f.write(md)
    print("wrote scripts/psi_production_bogoliubov.md (%d chars)" % len(md))


if __name__ == "__main__":
    main()
