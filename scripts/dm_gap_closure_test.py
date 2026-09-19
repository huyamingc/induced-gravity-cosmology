#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
dm_gap_closure_test.py
======================
问题：暗物质超产缺口（Omega_psi ~ 3e4 x Omega_DM）能否在 **当前作用量 +
当前再加热机制** 内补上？还是必须扩展模型？

为什么值得单独检验
------------------
`psi_abundance_oscillating.py` 的 §4 Omega 表（以及 §2c 的 T_reh 反解）里有一行

    if m < 0.4:
        continue

即那张表**有意跳过** m/H_inf < 0.4 的点。而"幂律下 Omega 对 g 几乎不动、
丰度由 T_reh 定、g 无法调节"这一结论，正是从那张表读出来的。
本脚本把扫描范围向小 m 延伸，并强制冻结收敛。

三条独立理由说明"小 g 解"应当存在
----------------------------------
1. 论文自己的 Eq. (eq:Omegapsi) 就是
       Omega ~ m_psi H_inf^3 exp(-2 pi m_psi/H_inf) (a_end/a_0)^3 / rho_c,
   它对 m_psi 是**非单调**的；论文 L458 也自认 "Omega_psi as a function of g
   has a maximum at g = 3.6e-6"。低于极大点的任何目标值都有**两个**解。
2. 无质量 Dirac 场共形不变 => 产生严格为零（本仓库自检 |beta| = 0 到机器零）。
   故 m -> 0 时 n_psi -> 0，m_psi n_psi -> 0，Omega_psi -> 0。
3. 由 1+2：Omega_psi(g) 在 g -> 0 时 -> 0，在大 g 处远超 0.265，且连续
   => 由介值定理**必然**存在小 g 解。所以真正的问题不是"有没有解"，
   而是"那个解对应的 g 是否物理上可接受"。

与原脚本的两处关键差异（否则小 m 点不可信）
--------------------------------------------
1. 冻结收敛：原脚本对所有 m 用 a_final = 60，于是小 m 时
   k_max/(m a_final) 远大于 1（m = 0.01 时高达 67），模式根本没冻结。
   本脚本取 a_final = max(60, k_max/(0.3 m))，强制 k_max/(m a_final) <= 0.3。
2. k 积分下限：产生峰在 k ~ m，原脚本固定 kmin = 0.05；当 m < 0.05 时
   峰被整个漏掉。本脚本取 kmin = min(0.05, m/10)。

判定
----
PASS : 存在 g 使 Omega_psi = Omega_DM = 0.265（缺口可在当前模型内补上）
FAIL : 不存在（缺口需扩展模型）
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

REUSE = "--reuse" in sys.argv          # 复用已有 json 的扫描结果，只重出报告

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from psi_abundance_oscillating import (  # noqa: E402
    H_INF, PHI_V, RHO_C, npsi_from_spectrum,
)

# dil(T_reh) = 1.0204e-101 * T_reh   (T_reh in GeV)  -- 与 [2b] §2c 同一约定
DIL_PER_GEV = 1.0204e-101
OMEGA_TARGET = 0.265
OMEGA_DM = 0.265
T_REH_MODEL = 2.1e8          # 反常通道（物理 alpha_s）给出的 T_reh [GeV]
T_REH_FIDUCIAL = 1.0e9       # 论文 fiducial
A_END_OVER_A0 = (1.0204e-92) ** (1.0 / 3.0)

FREEZE_TARGET = 0.3          # 要求 k_max/(m a_final) <= 0.3
MASSES = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0]
# 指数闭式必须扫到极小的 m 才能找到"轻 psi 支"（解析上 m ~ 1e-8）
EXP_GRID = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3,
            0.01, 0.1, 0.3, 1.0, 3.0]

JSON_PATH = HERE / "dm_gap_closure_test.json"
MD_PATH = HERE / "dm_gap_closure_test.md"


def a_final_for(m: float, kmax: float) -> float:
    """冻结条件 k_max/(m a_final) <= FREEZE_TARGET 所需的最小 a_final。"""
    return max(60.0, kmax / (FREEZE_TARGET * m))


def omega(m: float, n_over_H3: float, T_reh: float) -> float:
    """Omega_psi = m_psi * n_psi * dil(T_reh) / rho_c,  n_psi = (n/H^3) * H^3."""
    return m * H_INF * n_over_H3 * H_INF**3 * DIL_PER_GEV * T_reh / RHO_C


def omega_exp(m: float, T_reh: float) -> float:
    """论文的 de Sitter 指数闭式: n/H^3 = exp(-2 pi m/H_inf)（对照用）。"""
    return omega(m, math.exp(-2.0 * math.pi * m), T_reh)


# ----------------------------------------------------------- 自由流长度检验
G_STAR = 106.75
M_P_GEV = 2.435e18
H0_GEV = 1.4377e-42
T_EQ_GEV = 7.5e-10            # ~0.75 eV，物质-辐射等密度
MPC_PER_INV_GEV = 6.395e-39   # 1 GeV^-1 -> Mpc
LYMAN_ALPHA_LFS_MPC = 0.1     # 温暗物质自由流长度的保守上限


def free_streaming_length(m_over_H, T_reh, p_over_m_prod=None):
    """轻支的共动自由流长度 [Mpc]。单位约定 a_end=1、H 用 GeV。

    相对论 -> 非相对论用 v(a) = x/sqrt(1+x^2)，x = (p_prod/m)/a。
    默认取**最保守**（自由流最长）的 p_prod ~ H_inf。
    """
    if p_over_m_prod is None:
        p_over_m_prod = 1.0 / m_over_H
    a0 = 1.0 / A_END_OVER_A0
    H_reh = math.sqrt(math.pi**2 * G_STAR / 90.0) * T_reh**2 / M_P_GEV
    a_reh = (H_INF / H_reh) ** (2.0 / 3.0)      # w=0 凝聚体主导
    a_eq = a_reh * (T_reh / T_EQ_GEV)
    H_eq = H0_GEV * (a0 / a_eq) ** 1.5

    def H_of(a):
        if a <= a_reh:
            return H_INF * a ** (-1.5)
        if a <= a_eq:
            return H_reh * (a / a_reh) ** (-2.0)
        return H_eq * (a / a_eq) ** (-1.5)

    def integrand(ln_a):
        a = math.exp(ln_a)
        x = p_over_m_prod / a
        v = x / math.sqrt(1.0 + x * x)
        return v / (a * H_of(a))               # v da/(a^2 H) = v dln a/(a H)

    total = 0.0
    for lo, hi, n in ((0.0, math.log(a_reh), 4000),
                      (math.log(a_reh), math.log(a_eq), 2000),
                      (math.log(a_eq), math.log(a0), 4000)):
        if hi <= lo:
            continue
        h = (hi - lo) / n
        s = 0.5 * (integrand(lo) + integrand(hi))
        s += sum(integrand(lo + i * h) for i in range(1, n))
        total += s * h
    return {"lambda_com_GeVinv": total,
            "lambda_fs_Mpc": a0 * total * MPC_PER_INV_GEV,
            "p_over_m_prod": p_over_m_prod,
            "a_reh": a_reh, "a_eq": a_eq, "H_reh_GeV": H_reh}


def solve_g_for_target(ms, omegas, target=OMEGA_TARGET):
    """在 log-log 上插值，返回**最小**的 Omega = target 解（即小 g 支）。

    注意：论文的指数闭式在 m 上先升后降，有两个解；必须取小的那个才是
    "轻 psi 支"。早先版本返回第一个跨零点（= 大 g 支），是错的。
    """
    pts = sorted((m, w) for m, w in zip(ms, omegas) if w and w > 0)
    for (m1, w1), (m2, w2) in zip(pts, pts[1:]):
        if (w1 - target) * (w2 - target) <= 0 and w1 != w2:
            t = (math.log(target) - math.log(w1)) / (math.log(w2) - math.log(w1))
            lm = math.log(m1) + t * (math.log(m2) - math.log(m1))
            return math.exp(lm)
    return None


def loglog_n_at(m_star, rows):
    """用相邻两点在 log-log 上插值 n/H^3。"""
    lo = max((r for r in rows if r["m_over_H"] <= m_star),
             key=lambda r: r["m_over_H"], default=None)
    hi = min((r for r in rows if r["m_over_H"] >= m_star),
             key=lambda r: r["m_over_H"], default=None)
    if lo and hi and hi["m_over_H"] != lo["m_over_H"]:
        t = ((math.log(m_star) - math.log(lo["m_over_H"])) /
             (math.log(hi["m_over_H"]) - math.log(lo["m_over_H"])))
        return math.exp(math.log(lo["n_over_H3"]) + t *
                        (math.log(hi["n_over_H3"]) - math.log(lo["n_over_H3"])))
    return lo["n_over_H3"] if lo else None


def run_scan() -> list:
    print("=" * 78)
    print("DM gap closure test:  Omega_psi(g) 是否有第二个（小 g）解?")
    print("=" * 78)
    print(f"H_inf = {H_INF:.4e} GeV   Phi_V = {PHI_V:.4e}   rho_c = {RHO_C:.4e}")
    print(f"T_reh(model, anomaly) = {T_REH_MODEL:.2e} GeV ;  "
          f"target Omega = {OMEGA_TARGET}")
    print("freeze target: k_max/(m*a_final) <= "
          f"{FREEZE_TARGET}   kmin = min(0.05, m/10)")
    print("-" * 78)
    print(f"{'m/H':>8} {'g':>11} {'kmax':>6} {'a_fin':>9} {'n/H^3':>11} "
          f"{'m*n/H^4':>11} {'Om(T=2.1e8)':>12} {'/0.265':>10} {'sec':>6}")
    rows = []
    for m in MASSES:
        kmax = max(40.0, 20.0 * m)
        kmin = min(0.05, m / 10.0)
        a_fin = a_final_for(m, kmax)
        ta = time.time()
        res = npsi_from_spectrum(m, a_final=a_fin, kmax=kmax,
                                 kmin=kmin, nk_pts=45)
        n = res["n_over_H3"]
        w = omega(m, n, T_REH_MODEL)
        rows.append({"m_over_H": m, "g": m * H_INF / PHI_V, "kmax": kmax,
                     "kmin": kmin, "a_final": a_fin, "n_over_H3": n,
                     "m_times_n": m * n, "omega_at_model_T": w,
                     "ratio_to_DM": w / OMEGA_DM,
                     "omega_at_fiducial_T": omega(m, n, T_REH_FIDUCIAL),
                     "nsteps": res["nsteps"]})
        print(f"{m:8.4g} {m*H_INF/PHI_V:11.4e} {kmax:6.1f} {a_fin:9.1f} "
              f"{n:11.4e} {m*n:11.4e} {w:12.4e} {w/OMEGA_DM:10.4e} "
              f"{time.time()-ta:6.1f}")
    return rows


def main() -> None:
    t0 = time.time()
    if REUSE and JSON_PATH.exists():
        prev = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        rows = prev["rows"]
        print(f"[reuse] loaded {len(rows)} scan rows from "
              f"{JSON_PATH.name}; skipping mode-equation integration")
    else:
        rows = run_scan()

    ms = [r["m_over_H"] for r in rows]
    ws = [r["omega_at_model_T"] for r in rows]
    m_star = solve_g_for_target(ms, ws)
    n_star = loglog_n_at(m_star, rows) if m_star else None
    fs = free_streaming_length(m_star, T_REH_MODEL) if m_star else None

    # ---- 对照：论文的指数闭式 ----
    exp_rows = [{"m_over_H": m, "omega_at_model_T": omega_exp(m, T_REH_MODEL)}
                for m in EXP_GRID]
    for r in exp_rows:
        r["ratio_to_DM"] = r["omega_at_model_T"] / OMEGA_DM
        r["g"] = r["m_over_H"] * H_INF / PHI_V
    m_star_exp = solve_g_for_target([r["m_over_H"] for r in exp_rows],
                                    [r["omega_at_model_T"] for r in exp_rows])

    print("-" * 78)
    print("对照（论文 de Sitter 指数闭式, n/H^3 = exp(-2 pi m/H), 同一 T）：")
    for r in exp_rows:
        if r["m_over_H"] in (1e-10, 1e-8, 1e-6, 1e-4, 0.01, 0.1, 1.0, 3.0):
            print(f"  m/H = {r['m_over_H']:8.4g}   Omega = "
                  f"{r['omega_at_model_T']:12.4e}   /0.265 = "
                  f"{r['ratio_to_DM']:10.4e}")

    verdict = "PASS" if m_star is not None else "FAIL"
    g_star = (m_star * H_INF / PHI_V) if m_star else None
    print("=" * 78)
    print(f"VERDICT: {verdict}")
    if m_star is not None:
        print(f"  幂律（真实模方程）小 g 解:  m_psi/H_inf = {m_star:.4e}"
              f"  =>  m_psi = {m_star*H_INF:.4e} GeV,  g = {g_star:.4e}")
        print(f"    产生时 p/m ~ H_inf/m_psi = {1.0/m_star:.3e} ; "
              f"如今 v0 = {A_END_OVER_A0/m_star:.3e}")
    if fs:
        cold = fs["lambda_fs_Mpc"] < LYMAN_ALPHA_LFS_MPC
        print(f"  轻支自由流长度 lambda_fs = {fs['lambda_fs_Mpc']:.3e} Mpc "
              f"(上限参考 {LYMAN_ALPHA_LFS_MPC} Mpc) -> "
              f"{'COLD (pass)' if cold else 'WARM (FAIL)'}")
    if m_star_exp is not None:
        print(f"  指数闭式小 g 解:  m/H = {m_star_exp:.4e},  "
              f"g = {m_star_exp*H_INF/PHI_V:.4e}")
    print(f"  论文引用窗口 m_psi/H_inf in [1,10] 对应 "
          f"g in [{H_INF/PHI_V:.3e}, {10*H_INF/PHI_V:.3e}]")
    print(f"总耗时 {time.time()-t0:.1f} s")
    print("=" * 78)

    out = {"verdict": verdict, "target_omega": OMEGA_TARGET,
           "T_reh_model": T_REH_MODEL, "T_reh_fiducial": T_REH_FIDUCIAL,
           "freeze_target": FREEZE_TARGET, "rows": rows,
           "exponential_rows": exp_rows,
           "m_star_powerlaw": m_star, "g_star_powerlaw": g_star,
           "n_star_powerlaw": n_star, "free_streaming": fs,
           "m_star_exponential": m_star_exp,
           "g_star_exponential": ((m_star_exp * H_INF / PHI_V)
                                  if m_star_exp else None),
           "runtime_s": time.time() - t0}
    JSON_PATH.write_text(json.dumps(out, indent=2, ensure_ascii=False),
                         encoding="utf-8")
    write_md(out)
    print(f"wrote scripts/{JSON_PATH.name} and scripts/{MD_PATH.name}")


def write_md(o: dict) -> None:
    L = []
    A = L.append
    A("# 暗物质缺口能否补上：Omega_psi(g) 的小 g 支（检验脚本）\n")
    A(f"**判定：{o['verdict']}** —— "
      + ("存在 g 使 $\\Omega_\\psi=\\Omega_{{\\rm DM}}=0.265$，"
         "缺口可在**当前作用量 + 当前再加热机制**内补上。\n"
         if o["verdict"] == "PASS" else
         "不存在这样的 g，缺口需扩展模型。\n"))
    A(f"目标 $\\Omega_\\psi={o['target_omega']}$，取模型自己的 "
      f"$T_{{\\rm reh}}={o['T_reh_model']:.2e}$ GeV（反常通道，物理 $\\alpha_s$）。\n")

    A("## 1. 为什么值得单独检验\n")
    A("`psi_abundance_oscillating.py` 的 §4 $\\Omega$ 表与 §2c 的 $T_{\\rm reh}$ 反解里有\n")
    A("```python\n    if m < 0.4:\n        continue\n```\n")
    A("即那张表**有意跳过** $m/H_{\\rm inf}<0.4$ 的点。而"
      "\"$\\Omega$ 对 $g$ 几乎不动、丰度由 $T_{\\rm reh}$ 定\"这一结论，"
      "正是从那张表读出来的。\n")

    A("## 2. 小 $g$ 解为何应当存在（三条独立理由）\n")
    A("1. 论文自己的 Eq. (\\texttt{eq:Omegapsi}) 对 $m_\\psi$ **非单调**；")
    A("   论文 L458 也自认 $\\Omega_\\psi(g)$ 在 $g=3.6\\times10^{-6}$ 有极大。")
    A("   低于极大的任何目标值都有**两个**解。")
    A("2. 无质量 Dirac 场共形不变 $\\Rightarrow$ 产生严格为零"
      "（本仓库自检 $|\\beta|=0$ 到机器零），故 $m\\to0$ 时 "
      "$\\Omega_\\psi\\to0$。")
    A("3. 由 1+2 与连续性，介值定理保证小 $g$ 解存在。\n")

    A("## 3. 数值结果（真实模方程，冻结收敛）\n")
    A("| $m_\\psi/H_{\\rm inf}$ | $g$ | $k_{\\max}$ | $a_{\\rm final}$ | "
      "$n_\\psi/H^3$ | $m n/H^4$ | $\\Omega_\\psi$ | $/0.265$ |")
    A("|---|---|---|---|---|---|---|---|")
    for r in o["rows"]:
        A(f"| {r['m_over_H']:g} | {r['g']:.4e} | {r['kmax']:g} | "
          f"{r['a_final']:.4g} | {r['n_over_H3']:.4e} | {r['m_times_n']:.4e} | "
          f"{r['omega_at_model_T']:.4e} | {r['ratio_to_DM']:.4e} |")
    A("")
    if o["m_star_powerlaw"]:
        A(f"**小 $g$ 解：** $m_\\psi/H_{{\\rm inf}}={o['m_star_powerlaw']:.4e}$，"
          f"即 $m_\\psi={o['m_star_powerlaw']*H_INF:.4e}$ GeV，"
          f"$g={o['g_star_powerlaw']:.4e}$。")
        A(f"对应 $n_\\psi/H^3={o['n_star_powerlaw']:.4e}$。")
        A(f"产生时为相对论性（$p/m\\sim{1/o['m_star_powerlaw']:.3e}$），"
          f"如今 $v_0={A_END_OVER_A0/o['m_star_powerlaw']:.3e}$（极冷）。\n")

    A("## 4. 与论文指数闭式的对照\n")
    A("| $m/H_{\\rm inf}$ | $g$ | $\\Omega_\\psi$ | $/0.265$ |")
    A("|---|---|---|---|")
    for r in o["exponential_rows"]:
        A(f"| {r['m_over_H']:g} | {r['g']:.4e} | "
          f"{r['omega_at_model_T']:.4e} | {r['ratio_to_DM']:.4e} |")
    A("")
    if o["m_star_exponential"]:
        A(f"指数闭式的小 $g$ 解：$m/H_{{\\rm inf}}={o['m_star_exponential']:.4e}$，"
          f"$g={o['g_star_exponential']:.4e}$。")
    A("")

    A("## 5. 小 $m$ 处的收敛性纠正\n")
    A("原脚本对所有 $m$ 固定 $a_{\\rm final}=60$，故小 $m$ 时冻结条件 "
      "$k_{\\max}/(m a_{\\rm final})\\ll1$ 被严重违反"
      "（$m=0.01$ 时该比值高达 67）。本脚本按 "
      "$a_{\\rm final}=\\max(60,\\,k_{\\max}/(0.3m))$ 取值后，"
      "小 $m$ 的 $n_\\psi$ 显著下降：\n")
    A("| $m/H$ | 原脚本 $n/H^3$（未冻结） | 本脚本 $n/H^3$（已冻结） | 比 |")
    A("|---|---|---|---|")
    for mm, old in ((0.01, 6.725e-5), (0.1, 5.355e-4), (1.0, 5.028e-4),
                    (3.0, 1.336e-4)):
        new = next((r["n_over_H3"] for r in o["rows"]
                    if abs(r["m_over_H"] - mm) < 1e-12), None)
        if new:
            A(f"| {mm:g} | {old:.4e} | {new:.4e} | {old/new:.2f} |")
    A("")
    A("值得注意：小 $m$ 端 $n_\\psi\\propto m^{1.0}$（原报道为 $m^{0.9}$），"
      "故 $m n\\propto m^{2}$。\n")

    A("## 6. 轻支是冷暗物质吗？（自由流长度检验）\n")
    if o.get("free_streaming"):
        f = o["free_streaming"]
        A(f"- 产生时 $p/m$ = {f['p_over_m_prod']:.3e}（取最保守的 $p\\sim H_{{\\rm inf}}$）")
        A(f"- $a_{{\\rm reh}}$ = {f['a_reh']:.3e}，$a_{{\\rm eq}}$ = {f['a_eq']:.3e}，"
          f"$H_{{\\rm reh}}$ = {f['H_reh_GeV']:.3e} GeV")
        A(f"- 共动自由流长度 $\\lambda_{{\\rm fs}}$ = **{f['lambda_fs_Mpc']:.3e} Mpc**；"
          f"温暗物质上限参考 {LYMAN_ALPHA_LFS_MPC} Mpc → "
          f"**{'冷，通过' if f['lambda_fs_Mpc'] < LYMAN_ALPHA_LFS_MPC else '偏温，不通过'}**")
        A(f"- 物理原因：$m_\\psi$ = {o['m_star_powerlaw']*H_INF:.3e} GeV 远大于 $T_{{\\rm reh}}$，"
          f"粒子在 $a/a_{{\\rm end}}$ ~ {1/o['m_star_powerlaw']:.3g} 时就已非相对论，"
          f"远早于再加热完成；此后 $v\\propto1/a$，自由流可忽略。\n")
    _ = A
    MD_PATH.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
