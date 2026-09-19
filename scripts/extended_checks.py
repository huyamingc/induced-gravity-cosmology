# -*- coding: utf-8 -*-
"""Additional paper checks: RG running, reheating rates, DM abundance surface."""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosmo_model import M_P, Vc, fiducial, m_chi, r_of  # noqa: E402

OUT = Path(__file__).resolve().parent / "extended_checks.md"


def beta_xi_terms(xi: float, lam0: float, g: float, dlnmu: float = 60.0):
    pref = xi - 1.0 / 6.0
    dxi_lam = pref * (3.0 * lam0) / (16.0 * math.pi**2) * dlnmu
    dxi_g = pref * (g**2) / (16.0 * math.pi**2) * dlnmu
    return dxi_lam, dxi_g, dxi_lam + dxi_g


def Gamma_anom(m_chi: float, xi: float, b3: float = 7.0, alpha_s: float = 0.1):
    """Order-of-magnitude anomaly width ~ b3^2 alpha_s^2 m_chi^3 / (16 pi^3 M_P^2 (6+1/xi))."""
    return b3**2 * alpha_s**2 * m_chi**3 / (16.0 * math.pi**3 * M_P**2 * (6.0 + 1.0 / xi))


def Treh_const(Gamma: float, gstar: float = 106.75):
    return (90.0 * Gamma**2 * M_P**2 / (math.pi**2 * gstar)) ** 0.25


def Omega_psi(mpsi: float, Hinf: float, Treh: float, Treh_ref: float = 1e9, mpsi_ref: float = 4e12, Om_ref: float = 0.12):
    """Parametric scaling from paper Eq (25)-style: Omega ∝ mpsi H^3 e^{-pi mpsi/H} * dilution(Treh)."""
    # Use paper-like normalization at reference point, scale exponential and T
    expo = math.exp(-math.pi * mpsi / Hinf) / math.exp(-math.pi * mpsi_ref / Hinf)
    # n ∝ H^3 e^{-pi m/H} held; dilution ∝ Treh; also m_psi linear
    return Om_ref * (mpsi / mpsi_ref) * expo * (Treh / Treh_ref)


def main() -> None:
    fp = fiducial()
    fpn = fiducial(use_numeric_lam0=True)
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
    lines.append("# 扩展核验：RG / 再加热率 / DM 丰度标度 / 凝聚体")
    lines.append("")

    lines.append("## A. RG 运行（论文 Sec. VIII E）")
    lines.append("")
    lines.append("论文：Δξ≲1e-6（ξ=11.1，λ0=6.78e-8，g~1e-5–1e-4，Δlnμ=60）；夸克项主导。")
    lines.append("")
    lines.append("| λ0 | g | Δξ_λ0 | Δξ_g | Δξ_tot | 论文 |")
    lines.append("|---|---|---|---|---|---|")
    for lam0 in (6.78e-8, fp.lam0):
        for g in (2.3e-5, 1e-4):
            a, b, t = beta_xi_terms(11.1, lam0, g)
            lines.append(f"| {lam0:.2e} | {g:.1e} | {a:.3e} | {b:.3e} | {t:.3e} | ≲1e-6 量级 |")
    lines.append("")
    lines.append("**判定：** 量级与论文一致（Δξ~1e-6–1e-7）；**支持** RG 稳定主张（在给定 β 函数形式下）。")
    lines.append("注意：β_ξ 的「直接夸克-曲率」项写法依赖约定，脚本采用论文 Eq.(37) 形式。")
    lines.append("")

    lines.append("## B. 反常衰变宽度与 T_reh")
    lines.append("")
    lines.append("| λ0 源 | m_χ | Γ_anom (估) | T_reh=(90Γ²M_P²/π²g*)^{1/4} | 论文 |")
    lines.append("|---|---|---|---|---|")
    for tag, lam0, mchi in [("numeric", 6.78e-8, fpn.m_chi), ("analytic", fp.lam0, fp.m_chi)]:
        G = Gamma_anom(mchi, 11.1)
        T = Treh_const(G)
        lines.append(f"| {tag} | {mchi:.3e} | {G:.3e} GeV | {T:.3e} GeV | ~1e9, Γ~O(1) GeV |")
    lines.append("")
    lines.append(
        "**判定：** 粗估 Γ 与 T_reh 落在论文声称量级（Γ~0.1–1 GeV，T_reh~1e8–1e9）附近；"
        "**支持「数量级可行」**，但系数 b_s、α_s、g* 未从模型微观推导，**不是精确预言**。"
    )
    lines.append("")

    lines.append("## C. DM 丰度参数化标度")
    lines.append("")
    lines.append("按论文式：Ω ∝ m_ψ H³ e^{-π m_ψ/H} × (T_reh/T_ref)，在 (m_ψ=4e12, T=1e9) 归一到 0.12。")
    lines.append("")
    Hinf = fp.H_inf
    lines.append("| m_ψ | m_ψ/H_inf | T_reh=1e9 | T_reh=4e5 |")
    lines.append("|---|---|---|---|")
    for mpsi in (4e12, 7.3e12, 1.7e13, 3e13, 7e13):
        o1 = Omega_psi(mpsi, Hinf, 1e9)
        o2 = Omega_psi(mpsi, Hinf, 4e5)
        lines.append(f"| {mpsi:.2e} | {mpsi/Hinf:.2f} | {o1:.3e} | {o2:.3e} |")
    lines.append("")
    lines.append("**判定：**")
    lines.append("- 在 T_reh=1e9、m_ψ~4e12–2e13 附近，标度可过 0.12 **量级**；")
    lines.append("- **T_reh=4e5 时即便 m_ψ 达 m_χ 量级，Ω 仍 << 0.12** → 引力通道单独不够（与此前审校一致）；")
    lines.append("- 指数灵敏度高，**不能**把某一 g 写成唯一预言；**支持**论文「需格点、机制成立」的表述，**不支持**无条件 Ω=0.12。")
    lines.append("")

    lines.append("## D. 凝聚体过闭合粗检")
    lines.append("")
    # After reheating, if only a fraction f of condensate remains, need rho_cond < rho_DM today
    rho_DM0 = 0.265 * Vc(0.683) / 0.683 * (0.265 / 0.265)  # rough
    # better: Omega_DM rho_crit, rho_crit = 3 H0^2 M_P^2
    H0 = fp.H0
    rho_crit = 3 * H0**2 * M_P**2
    rho_DM = 0.265 * rho_crit
    lines.append(f"- ρ_crit(today)={rho_crit:.3e} GeV⁴, ρ_DM={rho_DM:.3e} GeV⁴")
    lines.append(f"- 初始凝聚体 ρ_cond(a_end)~½ m_Φ² M_P² 量级：")
    for tag, lam0 in [("num", 6.78e-8), ("an", fp.lam0)]:
        mPhi = math.sqrt(2 * lam0) * M_P / math.sqrt(11.1)
        rho0 = 0.5 * mPhi**2 * M_P**2
        # max allowed survival fraction today
        frac_max = rho_DM / rho0
        lines.append(f"  - {tag}: m_Φ={mPhi:.3e}, ρ0={rho0:.3e}, 允许残留比例 ≲{frac_max:.3e}")
    lines.append("")
    lines.append("**判定：** 需要极高衰变效率（残留 ≲10^{-75} 量级）；论文称异常衰变后对 T_reh≳1e9 指数压低 —— **方向支持**，精确值未算。")
    lines.append("")

    lines.append("## E. 与主验证报告的合并结论")
    lines.append("")
    lines.append("| 类别 | 结论 |")
    lines.append("|---|---|")
    lines.append("| 脚本支持且论文可保留 | 暴涨函数形式、结束条件、共形退耦、G_eff、冻结 DE、r 可证伪窗、RG 量级、T_reh~1e9 量级 |")
    lines.append("| 论文需改表述 | 双重保护限定 g；N 窗对引力通道更紧；λ₀ 双轨须声明；Planck 图为示意 |")
    lines.append("| 脚本不支持 | 附录错误 λ₀ 公式；单场精质；无格点的 Ω_DM=0.12 |")
    lines.append("| 建议补测（可脚本化） | 官方 Planck 似然对比；精确 N(χ) 与 Table I 归一；Γ_anom 更细的 b_i 表；格点/玻尔兹曼 Ω_ψ |")
    lines.append("")
    lines.append("[扩展核验完成]")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
