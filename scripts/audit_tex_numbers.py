# -*- coding: utf-8 -*-
"""Audit paper_prd_merged.tex scientific numbers against locked-N convention.

Primary source: scripts/n_convention_results.json (if present) + locked table values.
Flags:
  OK_EXACT   - locked-N primary values
  OK_NOTE    - allowed only as attractor/historical comparison (context window)
  FLAG       - old-convention value used as if current prediction
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEX = ROOT / "paper_prd_merged.tex"
JSON = Path(__file__).resolve().parent / "n_convention_results.json"
OUT = Path(__file__).resolve().parent / "tex_number_audit.md"

# Locked-N primary at xi=11.1 (from lock_n_convention / paper table)
PRIMARY = {
    "lam0_N50": 6.70e-8,
    "r_N50": 0.00425,
    "ns_N50": 0.9616,
    "r_max_2s": 0.0052,
    "m_chi_N50": 3.25e13,
    "H_N50": 1.64e13,
    "N_anomaly": 51,
    "N_band": (45, 58),
}

# Values that are ONLY OK in comparison/historical notes
OLD_NOTE = [
    (r"0\.00487", "attractor r at N=50"),
    (r"6\.78", "old numeric lambda0 / Table I draft"),
    (r"7\.46", "large-field analytic lambda0"),
    (r"4\.90", "wrong A5 formula value"),
    (r"3\.28", "old m_chi"),
    (r"1\.65", "old H_inf"),
    (r"4\.83", "old V0"),
    (r"5\.33", "analytic V0"),
    (r"3\.43", "analytic m_chi"),
    (r"1\.73\\times10", "analytic H_inf"),
]

ALLOW_CONTEXT = [
    "attractor",
    "Attractor",
    "closed form",
    "earlier draft",
    "Earlier draft",
    "analytic",
    "Analytic",
    "must not",
    "not mixed",
    "large-field",
    "Large-field",
    "cross-check",
    "Note.",
    "printed",
    "formula",
    "guide",
]

FLAG_CONTEXT_OK = False


def context_ok(line: str) -> bool:
    return any(k in line for k in ALLOW_CONTEXT)


def main() -> None:
    t = TEX.read_text(encoding="utf-8")
    lines = t.splitlines()
    out = []
    A = out.append
    A("# tex 数值审计（锁定 N 约定）")
    A("")
    A(f"目标：`{TEX.name}`")
    A("")

    A("## 主数值出现情况")
    A("")
    A("| 锁定值 | 模式 | 次数 |")
    A("|---|---|---|")
    checks = [
        ("λ0 N=50", r"6\.70"),
        ("r N=50", r"0\.00425"),
        ("r_max 2σ", r"0\.0052"),
        ("N anomaly", r"N\\simeq51|N\simeq51"),
        ("m_chi", r"3\.25"),
        ("H_inf", r"1\.64"),
        ("N band low", r"45"),
        ("N band high", r"58"),
    ]
    for name, pat in checks:
        A(f"| {name} | `{pat}` | {len(re.findall(pat, t))} |")
    A("")

    A("## 旧值扫描")
    A("")
    A("| 模式 | 行 | 上下文判断 | 判定 |")
    A("|---|---|---|---|")
    n_flag = n_ok = 0
    flags = []
    for pat, meaning in OLD_NOTE:
        for i, ln in enumerate(lines, 1):
            if not re.search(pat, ln):
                continue
            ok = context_ok(ln)
            if ok:
                n_ok += 1
                verdict = "OK_NOTE"
            else:
                n_flag += 1
                verdict = "**FLAG**"
                flags.append((i, meaning, ln[:160]))
            ctx = "allow" if ok else "no-allow-keyword"
            A(f"| {meaning} | L{i} | {ctx} | {verdict} |")
    A("")
    A(f"统计：OK_NOTE={n_ok}, FLAG={n_flag}")
    A("")
    if flags:
        A("### 需人工确认/修改的 FLAG")
        A("")
        for i, meaning, s in flags:
            A(f"- L{i} ({meaning}): {s}")
        A("")
    else:
        A("无 FLAG（旧值均出现在对比/历史语境，或已清除）。")
        A("")

    A("## 主张级核对")
    A("")
    A("| 检查项 | 结果 |")
    A("|---|---|")
    abs_txt = t.split("\\end{abstract}")[0] if "\\end{abstract}" in t else t[:3000]
    abs_bad = ("[48,55]" in abs_txt) or ("0.0053" in abs_txt) or ("self-consistently selecting the fiducial" in abs_txt)
    A(f"| 摘要无旧窗/旧 r/旧 fiducial | {'PASS' if not abs_bad else 'FAIL'} |")
    A(f"| 定义 eq:Ndef 存在 | {'PASS' if 'eq:Ndef' in t else 'FAIL'} |")
    A(f"| Table tab:sens 存在 | {'PASS' if 'tab:sens' in t else 'FAIL'} |")
    A(f"| DE calibration 声明 | {'PASS' if 'calibrated' in t else 'FAIL'} |")
    A(f"| 拒绝残余精质 | {'PASS' if 'quintessence' in t.lower() and ('excluded' in t.lower() or 'excluded' in t) else 'FAIL'} |")
    A("")

    A("## 结论")
    A("")
    if n_flag == 0:
        A("- 数值审计：**通过**（旧值仅出现在约定对比/附录历史说明中）。")
    else:
        A(f"- 数值审计：**发现 {n_flag} 处需处理**（见上 FLAG 列表）。")
    A("- 建议：将本脚本并入 `run_all.py`，每次改稿后重跑。")
    A("")
    A("[审计完成]")

    OUT.write_text("\n".join(out), encoding="utf-8")
    print("\n".join(out))
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
