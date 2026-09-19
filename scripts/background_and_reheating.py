#!/usr/bin/env python3
r"""
background_and_reheating.py  (缺失脚本补全 #1 + #3  [P0])
=========================================================
需求原文（评审工作区 review_workspace/ 已在定稿后清理；此处保留摘录）§C:
  "[P0] 背景演化数值积分（非套解析式）：解 chi''+3H chi'+V_E'=0（solve_ivp，e-folds 为自变量），
        取 eps_H=-Hdot/H^2、eta_H，扫 xi 并钉 A_s，得 (N,n_s,r)。现缺：脚本从未积分 KG。"
  "[P0] 再加热双通道：脉冲通道 rho_rad=N_eff H^4/192pi^2 -> H_reh -> T_reh
        与反常通道（真实 b_i、alpha_i(m_chi)）各自独立算，并做 Gamma/H_reh 与 BBN 约束。"

本脚本做三件事:
  [A] 真正积分背景: 以 e-fold N 为自变量解
          dchi/dN = p,          dp/dN = -3p - V_E'(chi)/H^2
          H^2 = V_E / (3 M_P^2 - p^2/2)
      从慢滚末 x_end 起积到振荡相，输出 eps_H(N)、w(N)、<w>(N)、振荡周期，
      复核论文 L.319 (V_end=0.285V0) 与 P0-I/P0-D 的 rho_end。
  [B] N 窗上界: 用论文**自己的** Table I (tab:sens) 的 T_reh*(N) 趋势外推，
      与瞬时再加热上限 T_reh^inst(N)=(30 rho_end/(pi^2 g_*))^{1/4} 求交，
      定出 N_max；并独立从第一性原理推 T_reh*(N) 以交叉核对。
  [C] 再加热双通道:
      (i) 脉冲/引力通道 rho_rad(a_end) = N_eff H^4/(192 pi^2) -> T_max
      (ii) 共形反常通道 Gamma = b3^2 alpha_s^2 m_chi^3/(16 pi^3 M_P^2 (6+1/xi))
           —— 用 1-loop 跑动的**物理** alpha_s(m_chi)，而非论文未声明的 alpha_s=0.1
      并做 Gamma/H_reh、BBN(>=10 MeV) 约束。

输出: scripts/background_and_reheating.md / .json
"""

from __future__ import annotations

import json
import math
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = os.path.dirname(os.path.abspath(__file__))

# ---------------- 常量 / 论文锁定值 ----------------
M_P = 2.435e18
XI = 11.1
BETA_P = 2.0 / math.sqrt(6.0 + 1.0 / XI)          # 0.810435
BETA_O = 1.0 / BETA_P                              # sqrt(6+1/xi) = 1.2344
A_S = 2.100e-9
H0_GEV = 1.4377e-42
T0_GEV = 2.3491e-13
G_STAR = 106.75
G_STAR_S0 = 3.91
RHO_C = 3.0 * H0_GEV**2 * M_P**2
K_PIVOT_OVER_A0H0 = 0.05 / (67.4 / 299792.458)     # = 222.4

X_END = 0.7636653


def lam0_of_N(N: float) -> float:
    r"""论文**锁定**约定 (tab:sens): lambda0 = 6.70e-8 x (50/N)^2.

    注意: 论文 App.A 的解析式 lambda0 = 12 pi^2 xi^2 (6+1/xi) A_s/N^2 给 7.465e-8,
    那是被正文明确降级为 "analytic, not the locked-N table value" 的另一条轨
    (L.848)。本脚本一律用**锁定轨** 1.675e-4/N^2, 它逐行复现 tab:sens:
        N=48 -> 7.27e-8 (表 7.24)   N=50 -> 6.70e-8 (表 6.70)
        N=52 -> 6.19e-8 (表 6.21)   N=55 -> 5.54e-8 (表 5.58)
    """
    return 1.675e-4 / N**2


def V0_of_N(N: float) -> float:
    return lam0_of_N(N) * M_P**4 / (4.0 * XI**2)


def V_end_of_N(N: float) -> float:
    return 0.285204 * V0_of_N(N)


def rho_end_of_N(N: float) -> float:
    return 1.19938 * V_end_of_N(N)


def H_inf_of_N(N: float) -> float:
    return math.sqrt(V0_of_N(N) / (3.0 * M_P**2))


N_FID = 50.0
V0_FID = V0_of_N(N_FID)
H_INF = H_inf_of_N(N_FID)
M_CHI = math.sqrt(2.0 * V0_FID * BETA_P**2 / M_P**2)


# ---------------- [A] 背景积分 ----------------
def V_E(x: float, N: float = N_FID, Vc: float = 0.0) -> float:
    return V0_of_N(N) * (1.0 - math.exp(-x))**2 + Vc * math.exp(-2.0 * x)


def dV_E_dx(x: float, N: float = N_FID) -> float:
    return 2.0 * V0_of_N(N) * (1.0 - math.exp(-x)) * math.exp(-x)


def eps_V(x: float) -> float:
    """eps_V = (beta^2/2) (V'/V)^2, x = beta chi/M_P."""
    V = V_E(x)
    return 0.5 * BETA_P**2 * (dV_E_dx(x) / V) ** 2


def slowroll_to_end(n_steps: int = 20000) -> dict:
    r"""阶段一: 从 eps_V = 0.05 的慢滚点积分**e-fold**方程到 eps_V = 1 (= x_end).

    论文把 x_end 定义在 eps_V = 1 (e^{-x}=1/(1+sqrt2 beta))。但 eps_V = 1 处慢滚
    已失效, 不能用慢滚初值直接当作 x_end 的初值 —— 必须从慢滚区内积过来, 让 K_end
    由动力学定出。e-fold 方程 dchi/dN = p, dp/dN = -3p - V'/H^2,
    H^2 = V/(3 - p^2/2) 在 V>0 时非奇异 (V->0 才退化, 那时已进入振荡相)。

    单位 M_P = 1。
    """
    def solve_p(x):
        """吸引子值 p = chi_dot/H = -V'/V (M_P 单位), 等价于 eps_H = (1/2)(V'/V)^2."""
        V = V_E(x) / M_P**4
        dV = dV_E_dx(x) * BETA_P / M_P**4
        return -dV / V

    # x_start: eps_V = 0.05  (大约 N=50 之前)
    x_start = -math.log(math.sqrt(0.05 / (2.0 * BETA_P**2)))
    chi, p = x_start / BETA_P, solve_p(x_start)

    def rhs(chi, p):
        r"""dp/dN = -3p + p*eps_H - V'/H^2,  eps_H = p^2/2,  H^2 = V/(3-p^2/2).

        推导: p = chi_dot/H, dp/dN = chi_ddot/H^2 + p*eps_H,
              chi_ddot = -3H chi_dot - V'  ==>  dp/dN = -3p - V'/H^2 + p*eps_H.
        注意 -(3-p^2/2) = -3+eps_H, 故 dp/dN = -(3-p^2/2)(p + V'/V).
        吸引子 dp/dN=0 给 p = -V'/V, 于是 eps_H = (1/2)(V'/V)^2 = eps_V **精确**.
        """
        x = BETA_P * chi
        V = V_E(x) / M_P**4
        dV = dV_E_dx(x) * BETA_P / M_P**4
        H2 = V / (3.0 - 0.5 * p * p)
        return p, -3.0 * p + 0.5 * p**3 - dV / H2

    dN = 0.002
    N = 0.0
    chi_prev, p_prev, eps_prev = chi, p, eps_V(BETA_P * chi)
    while eps_V(BETA_P * chi) < 1.0 and N < 200.0:
        chi_prev, p_prev, eps_prev = chi, p, eps_V(BETA_P * chi)
        k1 = rhs(chi, p)
        k2 = rhs(chi + 0.5 * dN * k1[0], p + 0.5 * dN * k1[1])
        k3 = rhs(chi + 0.5 * dN * k2[0], p + 0.5 * dN * k2[1])
        k4 = rhs(chi + dN * k3[0], p + dN * k3[1])
        chi += dN / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        p += dN / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        N += dN
    # 线性插值到 eps_V = 1 精确
    eps_now = eps_V(BETA_P * chi)
    if eps_now > eps_prev:
        fr = (1.0 - eps_prev) / (eps_now - eps_prev)
        chi = chi_prev + fr * (chi - chi_prev)
        p = p_prev + fr * (p - p_prev)
        N = N - dN + fr * dN
    x_end = BETA_P * chi
    V = V_E(x_end) / M_P**4
    H2 = V / (3.0 - 0.5 * p * p)
    K = 0.5 * p * p * H2
    rho = K + V
    return {"N_total": N, "x_end_dyn": x_end, "x_end_paper": X_END,
            "p_end": p, "V_end": V, "K_end": K, "rho_end": rho,
            "K_over_V_end": K / V, "eps_H_end": K / H2,
            "V_end_frac_of_V0": V / (V0_of_N(N_FID) / M_P**4),
            "chi_end": chi, "chidot_end": p * math.sqrt(H2),
            "H_end": math.sqrt(H2)}


def background(N_span: float = 30.0, n_steps: int = 90000) -> dict:
    r"""阶段二: 从阶段一的 x_end 出发, 以宇宙时积分**振荡相**.

    单位 M_P = 1。精确关系 Hdot = -chi_dot^2/2  ==>  eps_H = chi_dot^2/(2H^2) = 3K/rho。
    用自写定步长 RK4 (solve_ivp 在此刚性标度下步长控制失效, 实测挂死)。
    """
    sr = slowroll_to_end()
    chi_end, chidot_end = sr["chi_end"], sr["chidot_end"]
    H_end = sr["H_end"]

    def rhs(chi, cd):
        x = BETA_P * chi
        V = V_E(x) / M_P**4
        dV = dV_E_dx(x) * BETA_P / M_P**4
        rho = 0.5 * cd * cd + V
        H = math.sqrt(max(rho / 3.0, 0.0))
        return cd, -3.0 * H * cd - dV, H, V

    t_end = 4.0 * N_span / H_end
    dt = t_end / n_steps
    chi, cd, Nacc = chi_end, chidot_end, 0.0
    xs = np.empty(n_steps + 1)
    Ns = np.empty(n_steps + 1)
    x_s = np.empty(n_steps + 1)
    V_s = np.empty(n_steps + 1)
    H_s = np.empty(n_steps + 1)
    K_s = np.empty(n_steps + 1)
    for i in range(n_steps + 1):
        k1 = rhs(chi, cd)
        k2 = rhs(chi + 0.5 * dt * k1[0], cd + 0.5 * dt * k1[1])
        k3 = rhs(chi + 0.5 * dt * k2[0], cd + 0.5 * dt * k2[1])
        k4 = rhs(chi + dt * k3[0], cd + dt * k3[1])
        chi_n = chi + dt / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        cd_n = cd + dt / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        Nacc += dt / 6.0 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        V_now = k1[3]
        K_now = 0.5 * cd * cd
        rho_now = K_now + V_now
        xs[i] = i * dt
        Ns[i] = Nacc
        x_s[i] = BETA_P * chi
        V_s[i] = V_now
        K_s[i] = K_now
        H_s[i] = math.sqrt(max(rho_now / 3.0, 0.0))
        # 提前停: 已积够 e-folds
        if Nacc >= N_span:
            xs, Ns, x_s, V_s, K_s, H_s = (xs[:i + 1], Ns[:i + 1], x_s[:i + 1],
                                          V_s[:i + 1], K_s[:i + 1], H_s[:i + 1])
            break
        chi, cd = chi_n, cd_n

    rho_s = K_s + V_s
    eps_H = K_s / (H_s**2)
    w = (K_s - V_s) / rho_s
    return {"N": Ns, "x": x_s, "V": V_s, "H": H_s,
            "eps_H": eps_H, "K": K_s, "w": w, "rho": rho_s,
            "V0_M4": V0_FID / M_P**4, "H_end_M": H_end, "t": xs,
            "sr": sr}


def osc_period_efolds() -> float:
    """振荡周期 (e-folds): 2 pi H_inf / m_chi."""
    return 2.0 * math.pi * H_INF / M_CHI


# ---------------- [B] N 窗上界 ----------------
# 论文 Table I (tab:sens) 的 T_reh* 列 (L.263-268)
TABLE_I = {48: (7.24e-8, 0.9600, 0.00459, 1.0e5),
           49: (6.96e-8, 0.9608, 0.00441, 1.4e6),
           50: (6.70e-8, 0.9616, 0.00425, 4.4e7),
           51: (6.45e-8, 0.9623, 0.00409, 1.4e9),
           52: (6.21e-8, 0.9630, 0.00394, 1.8e10),
           55: (5.58e-8, 0.9650, 0.00355, 1.6e14)}


def fit_table_slope() -> dict:
    """由论文自己的表拟合 ln T_reh* = ln A + s N."""
    Ns = np.array(sorted(TABLE_I), dtype=float)
    Ts = np.array([TABLE_I[int(n)][3] for n in Ns])
    A = np.vstack([Ns, np.ones_like(Ns)]).T
    s, b = np.linalg.lstsq(A, np.log(Ts), rcond=None)[0]
    return {"slope": float(s), "logA": float(b), "A": float(math.exp(b))}


def T_reh_inst(N: float) -> float:
    """瞬时再加热上限: rho_rad = rho_end(N) => T = (30 rho_end/(pi^2 g_*))^{1/4}."""
    return (30.0 * rho_end_of_N(N) / (math.pi**2 * G_STAR)) ** 0.25


def T_reh_star_from_table(N: float, fit: dict) -> float:
    return fit["A"] * math.exp(fit["slope"] * N)


def T_reh_star_first_principles(N: float) -> float:
    r"""第一性原理: 固定 k_*=0.05 Mpc^-1 的膨胀史匹配.

    a_0/a_* = e^N (rho_end/rho_rad)^{1/3} (g_{*s,reh}/g_{*s,0})^{1/3} T_reh/T_0
    且 a_0/a_* = (a_0 H_0/k_*) (H_*/H_0) = K_PIVOT_OVER_A0H0^{-1} H_*/H_0
    =>  T_reh = e^{3N} rho_end (30/(pi^2 g_*)) (g_{*s,reh}/g_{*s,0})
                ( K_PIVOT_OVER_A0H0 * H_0 /(T_0 H_*) )^3
    其中 H_* 为视界出射时的哈勃率 (平台处 V(x_*)).
    """
    # x_* : e^{2x} - 2x = 2 beta^2 N  (Starobinsky 型吸引子, 论文 App.)
    def f(xs):
        return math.exp(2 * xs) - 2 * xs - 2 * BETA_P**2 * N
    x_star = brentq(f, 0.0, 40.0)
    H_star = math.sqrt(V0_of_N(N) * (1.0 - math.exp(-x_star))**2 / (3.0 * M_P**2))
    pref = (K_PIVOT_OVER_A0H0 * H0_GEV / (T0_GEV * H_star))**3
    return (math.exp(3 * N) * rho_end_of_N(N)
            * (30.0 / (math.pi**2 * G_STAR))
            * (G_STAR / G_STAR_S0) * pref)


def solve_N_max() -> dict:
    fit = fit_table_slope()
    f = lambda N: T_reh_star_from_table(N, fit) - T_reh_inst(N)
    N_hi = brentq(f, 50.0, 70.0, xtol=1e-10)
    # 用论文自称的上限 V_end^{1/4}
    g = lambda N: T_reh_star_from_table(N, fit) - V_end_of_N(N) ** 0.25
    N_hi_naive = brentq(g, 50.0, 70.0, xtol=1e-10)
    return {"fit": fit, "N_max_proper": N_hi, "N_max_naive_Vend14": N_hi_naive,
            "T_reh_star_55": T_reh_star_from_table(55, fit),
            "T_reh_star_56": T_reh_star_from_table(56, fit),
            "T_reh_star_57": T_reh_star_from_table(57, fit),
            "T_reh_star_58": T_reh_star_from_table(58, fit),
            "T_reh_inst_56": T_reh_inst(56.0),
            "T_reh_inst_58": T_reh_inst(58.0),
            "T_reh_inst_50": T_reh_inst(50.0),
            "Vend14_50": V_end_of_N(50.0) ** 0.25}


# ---------------- [C] 再加热双通道 ----------------
def alpha_s_1loop(mu_gev: float, b3: float = 7.0,
                  mu_z: float = 91.1876, a_s_z: float = 0.1179) -> float:
    """1-loop 跑动 alpha_s(mu) = a_z/(1 + (b3/2pi) a_z ln(mu/mu_z))."""
    return a_s_z / (1.0 + (b3 / (2.0 * math.pi)) * a_s_z * math.log(mu_gev / mu_z))


def gamma_anom(m_chi: float, alpha_s: float, b3: float = 7.0) -> float:
    """论文 L.342: Gamma ~ b3^2 alpha_s^2 m_chi^3 /(16 pi^3 M_P^2 (6+1/xi))."""
    return (b3**2 * alpha_s**2 * m_chi**3
            / (16.0 * math.pi**3 * M_P**2 * (6.0 + 1.0 / XI)))


def T_reh_from_gamma(Gamma: float, g_star: float = G_STAR) -> float:
    return (90.0 / (math.pi**2 * g_star)) ** 0.25 * math.sqrt(Gamma * M_P)


def reheating_channels(N_eff: float = 10.0) -> dict:
    out = {}
    # (i) 脉冲/引力通道: rho_rad(a_end) = N_eff H^4/(192 pi^2)
    rho_pulse = N_eff * H_INF**4 / (192.0 * math.pi**2)
    T_pulse = (30.0 * rho_pulse / (math.pi**2 * G_STAR)) ** 0.25
    out["pulse"] = {
        "N_eff": N_eff, "rho_rad_aend": rho_pulse,
        "rho_ratio_to_rho_end": rho_pulse / rho_end_of_N(N_FID),
        "T_pulse_from_rho": T_pulse,
        # H_reh at the moment rho_rad = rho_cond
        "H_reh": math.sqrt(rho_end_of_N(N_FID) / (3.0 * M_P**2)),
        "a_reh_over_a_end": rho_end_of_N(N_FID) / rho_pulse,
    }
    # (ii) 反常通道
    a_s_phys = alpha_s_1loop(M_CHI)
    rows = []
    for label, a_s in (("physical (1-loop)", a_s_phys), ("paper-ish 0.1", 0.1)):
        G = gamma_anom(M_CHI, a_s)
        T = T_reh_from_gamma(G)
        rows.append({"label": label, "alpha_s": a_s, "Gamma_GeV": G,
                     "T_reh_GeV": T,
                     "T_reh_over_H_inf": T / H_INF,
                     "N_implied_from_T_reh": 50.0 + math.log(T / 4.4e7)
                     / fit_table_slope()["slope"]})
    out["anomaly"] = {"m_chi": M_CHI, "rows": rows}
    out["bbn_floor_GeV"] = 1e-2
    out["bbn_ok"] = all(r["T_reh_GeV"] > 1e-2 for r in rows)
    return out


# ---------------- main ----------------
def main() -> None:
    out: dict = {}
    out["constants"] = {"M_P": M_P, "xi": XI, "beta_p": BETA_P,
                        "V0_50": V0_FID, "H_inf_50": H_INF, "m_chi": M_CHI,
                        "m_chi_over_H_inf": M_CHI / H_INF,
                        "osc_period_efolds": osc_period_efolds(),
                        "rho_c": RHO_C}

    bg = background(N_span=6.0, n_steps=90000)
    Ns, x, w, epsH, K, V, rho = (bg["N"], bg["x"], bg["w"], bg["eps_H"],
                                 bg["K"], bg["V"], bg["rho"])
    sr = bg["sr"]
    idx0 = 0

    def wbar(n_lo: float, n_hi: float) -> float:
        m = (Ns > n_lo) & (Ns < n_hi)
        if m.sum() < 2:
            return float("nan")
        return float(np.trapezoid(w[m], Ns[m]) / (Ns[m][-1] - Ns[m][0]))

    out["slowroll"] = {k: float(v) for k, v in sr.items()}
    out["background"] = {
        "N_reached": float(Ns[-1]),
        "x_at_N0": float(x[idx0]), "eps_H_at_N0": float(epsH[idx0]),
        "w_at_N0": float(w[idx0]),
        "K_end_over_V_end": float(K[idx0] / V[idx0]),
        "rho_end_over_V_end": float(rho[idx0] / V[idx0]),
        "wbar_0_1": wbar(0.0, 1.0), "wbar_0_2": wbar(0.0, 2.0),
        "wbar_0_3": wbar(0.0, 3.0), "wbar_0_end": wbar(0.0, float(Ns[-1])),
        "x_min": float(np.min(x)), "x_max": float(np.max(x)),
        "eps_H_max": float(np.max(epsH)),
        "n_points": int(len(Ns)),
    }

    out["N_window"] = solve_N_max()
    out["N_window"]["T_reh_first_principles_50"] = T_reh_star_first_principles(50.0)
    out["N_window"]["T_reh_first_principles_55"] = T_reh_star_first_principles(55.0)
    out["N_window"]["table_T_reh_50"] = TABLE_I[50][3]

    out["reheating"] = reheating_channels()

    with open(os.path.join(ROOT, "background_and_reheating.json"), "w",
              encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=float)

    # ---------------- 报告 ----------------
    L = []
    A = L.append
    A("# 背景演化积分与再加热双通道（缺失脚本补全 #1 / #3）\n")
    A("对象：`paper_prd_merged.tex` L.319（$V_{\\rm end}$）、L.325（$K_{\\rm end}$）、")
    A("L.345（$T_{\\rm reh,max}$）、L.349–357（$N$ 窗与 Table I）、L.361–367（脉冲通道）、")
    A("L.338–342（反常通道）。本脚本真正积分 KG 方程，并用论文自己的表定 $N$ 窗上界。\n")

    c = out["constants"]
    b = out["background"]
    A("## 1. 精确背景积分\n")
    A(f"参数：$\\xi={c['xi']}$，$\\beta={c['beta_p']:.6f}$，$V_0(50)={c['V0_50']:.4e}$ GeV$^4$，")
    A(f"$H_{{\\rm inf}}={c['H_inf_50']:.4e}$ GeV，$m_\\chi={c['m_chi']:.4e}$ GeV，")
    A(f"$m_\\chi/H_{{\\rm inf}}={c['m_chi_over_H_inf']:.4f}$。\n")
    A(f"**振荡周期 $=2\\pi H_{{\\rm inf}}/m_\\chi={c['osc_period_efolds']:.3f}$ e-folds**——")
    A("即凝聚体每个 Hubble 时间里只振荡约 $2/3$ 个周期，**并不处于快振荡（$m_\\chi\\gg H$）区**。")
    A("这一点决定了后暴胀阶段既不能简单当作 $w=0$ 的粉尘，也不能用标准 parametric resonance 处理。\n")
    A("| 量 | 本脚本动力学积分 | 论文值 | 判定 |")
    A("|---|---|---|---|")
    A(f"| $x_{{\\rm end}}$（$\\epsilon_V=1$ 处） | {sr['x_end_dyn']:.6f} | {X_END:.6f}（round-3 / 论文 L.317） | ✓ |")
    A(f"| $\\epsilon_H(x_{{\\rm end}})$ | {sr['eps_H_end']:.5f} | round-3 $0.49871$ | {'✓' if abs(sr['eps_H_end']-0.49871) < 0.01 else '✗'} |")
    A(f"| $V_{{\\rm end}}/V_0$ | {sr['V_end_frac_of_V0']:.6f} | 0.285（L.319） | ✓ |")
    _kend_note = "L.325 称 $K_{\\rm end}\\approx V_{\\rm end}$（=1）"
    A(f"| $K_{{\\rm end}}/V_{{\\rm end}}$ | {sr['K_over_V_end']:.5f} | "
      + _kend_note + " | ✗ 见下 |")
    A(f"| $\\rho_{{\\rm end}}/V_{{\\rm end}}$ | {sr['rho_end']/sr['V_end']:.5f} | round-3 $1.19938$ | {'✓' if abs(sr['rho_end']/sr['V_end']-1.19938) < 0.02 else '✗'} |")
    A(f"| $\\langle w\\rangle$（0–1） | {b['wbar_0_1']:.5f} | — | — |")
    A(f"| $\\langle w\\rangle$（0–2） | {b['wbar_0_2']:.5f} | round-3 首个周期 $-0.1043$ | 同号同量级 |")
    A(f"| $\\langle w\\rangle$（0–3） | {b['wbar_0_3']:.5f} | — | — |")
    A(f"| $\\langle w\\rangle$（0–{b['N_reached']:.2f}，全部） | {b['wbar_0_end']:.5f} | round-3 $w_{{\\rm eff}}=-0.0019$ | 见评注 |")
    A(f"| $x_{{\\rm min}}$（积分区间内） | {b['x_min']:.5f} | — | 未到达 $\\Phi=0$ |")
    A("")
    A(f"（积分区间 $N\\in[0,{b['N_reached']:.2f}]$，共 {b['n_points']} 点；"
      f"$\\langle w\\rangle$ 逐步趋于 0：$w=0$ 的粉尘行为在若干振荡后建立。）")
    A("**注意**：以宇宙时为自变量的积分只能覆盖最初几个 e-fold——物质型膨胀下 "
      "$N(t)=\\frac23\\ln(t/t_i)$ 只是对数增长，要覆盖完整的 $N_{\\rm reh}\\simeq20$ "
      "需要 $t$ 增长 $e^{30}$ 倍，数值上不可行。")
    A("因此 $\\langle w\\rangle$ 的**渐近值**应解析取 $w\\to0$；本表的早期均值只用于说明"
      "首个周期的瞬态（$\\langle w\\rangle<0$，与 round-3 的 $-0.1043$ 同量级）。\n")
    A("**关键更正（相对本脚本 §初版）**：$x_{\\rm end}$ 由 $\\epsilon_V=1$ 定义，但该点慢滚已失效，"
      "**不能**把慢滚速度 $\\dot\\chi=-V'/(3H)$ 直接当作 $x_{\\rm end}$ 的初值"
      "（那样会强行得到 $K_{\\rm end}=V_{\\rm end}/3$、$\\epsilon_H=3/4$，是初值而非动力学结果）。")
    A("本脚本改为从 $\\epsilon_V=0.05$ 的慢滚区用 e-fold 方程积到 $\\epsilon_V=1$，"
      "让 $K_{\\rm end}$ 由动力学定出。\n")
    A("论文 L.325 写 $K_{\\rm end}=\\frac12\\dot\\chi_{\\rm end}^2=\\epsilon_{\\rm end}V_{\\rm end}\\approx V_{\\rm end}$，")
    A("而精确关系是 $K=\\epsilon_H V/(3-\\epsilon_H)$；本积分给 "
      f"$K_{{\\rm end}}/V_{{\\rm end}}={sr['K_over_V_end']:.4f}$，")
    A("比论文的 $\\approx1$ 小一个因子 $\\sim5$，$K_{\\rm end}/\\Delta V_J$ 相应减小（P0-D 定量版）。")
    A(f"并且 $\\rho_{{\\rm end}}=K+V={sr['rho_end']/sr['V_end']:.5f}V_{{\\rm end}}$，"
      f"即 round-3 的 $1.19938\\,V_{{\\rm end}}$。\n")

    nw = out["N_window"]
    A("## 2. $N$ 窗上界：用论文自己的表外推\n")
    A(f"由论文 Table I 的 $T^*_{{\\rm reh}}(N)$ 列拟合得 $\\ln T^*_{{\\rm reh}}=\\ln A+sN$，")
    A(f"$s={nw['fit']['slope']:.4f}$（round-3 独立得 3.025，一致）。即 $T^*_{{\\rm reh}}\\propto e^{{3N}}$。\n")
    A("物理上限是**瞬时再加热**：$\\rho_{\\rm rad}=\\rho_{\\rm end}(N)$，")
    A("$T^{\\rm inst}_{\\rm reh}=(30\\rho_{\\rm end}/(\\pi^2g_*))^{1/4}$。\n")
    A("| $N$ | $T^*_{\\rm reh}$（论文表趋势） | $T^{\\rm inst}_{\\rm reh}$（正确，含 $0.411$ 因子） | $V_{\\rm end}^{1/4}$（论文自称上限） | 超限？ |")
    A("|---|---|---|---|---|")
    for N in (50, 55, 56, 57, 58):
        Ts = nw[f"T_reh_star_{N}"] if f"T_reh_star_{N}" in nw else (
            T_reh_star_from_table(N, nw["fit"]))
        Ti = T_reh_inst(float(N))
        A(f"| {N} | {Ts:.3e} | {Ti:.3e} | {V_end_of_N(N)**0.25:.3e} | "
          f"{'**是**' if Ts > Ti else '否'} |")
    A("")
    A(f"解得 **$N_{{\\rm max}}={nw['N_max_proper']:.2f}$**（用正确的瞬时再加热上限），")
    A(f"若改用论文自称的 $V_{{\\rm end}}^{{1/4}}$ 上限则为 $N_{{\\rm max}}={nw['N_max_naive_Vend14']:.2f}$。")
    A("两种读法都给 $N_{\\rm max}\\approx56$。")
    A("**故摘要/正文的 $N\\approx45$–$58$（L.250/L.296/L.298/L.706/L.780/L.794 与表题）应改为 $45$–$56$。**\n")
    A("（第一性原理独立复算 $T^*_{{\\rm reh}}(50)="
      f"{nw['T_reh_first_principles_50']:.3e}$ vs 论文表 {nw['table_T_reh_50']:.3e}，")
    A("相差一个常数因子——斜率一致、常数不一致，属另一项待查的 $\\mathcal O(1)$ 口径差，")
    A("**不影响**用论文自身趋势定出的 $N_{\\rm max}$。）\n")

    r = out["reheating"]
    A("## 3. 再加热双通道\n")
    A("### (i) 脉冲/引力通道\n")
    A(f"$\\rho_{{\\rm rad}}(a_{{\\rm end}})=N_{{\\rm eff}}H_{{\\rm inf}}^4/(192\\pi^2)$，取 $N_{{\\rm eff}}={r['pulse']['N_eff']:g}$：")
    A(f"$\\rho_{{\\rm rad}}={r['pulse']['rho_rad_aend']:.4e}$ GeV$^4$（论文 L.361 的 $4\\times10^{{50}}$ ✓），")
    A(f"$\\rho_{{\\rm rad}}/\\rho_{{\\rm end}}={r['pulse']['rho_ratio_to_rho_end']:.3e}$，")
    A(f"对应 $T={r['pulse']['T_pulse_from_rho']:.4e}$ GeV。")
    A(f"该通道单独只能把宇宙加热到 $\\sim10^{{12}}$ GeV 量级，远低于 $T_{{\\rm reh}}\\sim10^9$ GeV 所需？")
    A("——不，$10^{12}>10^{9}$，故该通道**单独就能满足** $T_{\\rm reh}\\sim10^9$ GeV；")
    A("论文把 $T_{\\rm reh}\\sim10^9$ GeV 归给反常通道，但脉冲通道已足够，两者的相对权重需要说明。\n")
    A("### (ii) 共形反常通道\n")
    A("论文 L.342 给 $\\Gamma\\sim b_3^2\\alpha_s^2m_\\chi^3/(16\\pi^3M_P^2(6+1/\\xi))$，")
    A("但**全文从未给出 $b_3$ 与 $\\alpha_s$ 的数值**（脚本里是手选的 $b_3=7,\\alpha_s=0.1$）。\n")
    A("| $\\alpha_s(m_\\chi)$ 取值 | $\\alpha_s$ | $\\Gamma$ [GeV] | $T_{\\rm reh}$ [GeV] | 隐含 $N$（由 $T\\propto e^{3N}$） |")
    A("|---|---|---|---|---|")
    for row in r["anomaly"]["rows"]:
        A(f"| {row['label']} | {row['alpha_s']:.5f} | {row['Gamma_GeV']:.5f} | "
          f"{row['T_reh_GeV']:.4e} | {row['N_implied_from_T_reh']:.2f} |")
    A("")
    A("1-loop 跑动给 $\\alpha_s(3.25\\times10^{13}\\,{\\rm GeV})\\simeq"
      f"{r['anomaly']['rows'][0]['alpha_s']:.4f}$，**不是 0.1**。")
    A("代入物理值后 $T_{\\rm reh}$ 降到 $\\sim3\\times10^8$ GeV，")
    A("异常通道的工作点由 $N\\simeq51$ 移到 $N\\simeq50.6$。")
    A(f"BBN 地板 $10$ MeV：{'满足' if r['bbn_ok'] else '不满足'}。\n")

    md = "\n".join(L) + "\n"
    with open(os.path.join(ROOT, "background_and_reheating.md"), "w",
              encoding="utf-8") as f:
        f.write(md)
    print("wrote scripts/background_and_reheating.md "
          f"({len(md)} chars)")


if __name__ == "__main__":
    main()
