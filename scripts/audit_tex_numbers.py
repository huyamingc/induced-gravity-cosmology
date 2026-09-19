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
    # Appendix A splits the disclaimer onto the line after the aligned
    # equation (7.46e-8), and App. D5c introduces the hypothetical stable-
    # condensate value (1.65e-29) in the sentence before the equation, so a
    # small neighbor window plus these keywords is required.
    "hypothetical",
    "stable condensate",
    "stable-condensate",
    "not the locked",
]

# The disclaimer often sits on the line above/below the number (aligned
# equation environments split lines), so the context check also looks at
# neighboring lines.
CONTEXT_WINDOW = 2


def context_ok(lines: list, idx: int) -> bool:
    """idx is the 0-based index of the matched line; also scan neighbors."""
    lo = max(0, idx - CONTEXT_WINDOW)
    hi = min(len(lines), idx + CONTEXT_WINDOW + 1)
    return any(k in lines[j] for j in range(lo, hi) for k in ALLOW_CONTEXT)


def main() -> None:
    t = TEX.read_text(encoding="utf-8")
    lines = t.splitlines()
    out = []
    A = out.append
    A("# tex number audit (locked-N convention)")
    A("")
    A(f"Target: `{TEX.name}`")
    A("")

    A("## Locked primary occurrences")
    A("")
    A("| Locked value | Pattern | Count |")
    A("|---|---|---|")
    checks = [
        ("lambda0 N=50", r"6\.70"),
        ("r N=50", r"0\.00425"),
        ("r_max 2sigma", r"0\.0052"),
        ("N anomaly", r"N\\simeq51|N\simeq51"),
        ("m_chi", r"3\.25"),
        ("H_inf", r"1\.64"),
        ("N band low", r"45"),
        ("N band high", r"58"),
    ]
    for name, pat in checks:
        A(f"| {name} | `{pat}` | {len(re.findall(pat, t))} |")
    A("")

    A("## Old-value scan")
    A("")
    A("| Pattern | Line | Context check | Verdict |")
    A("|---|---|---|---|")
    n_flag = n_ok = 0
    flags = []
    for pat, meaning in OLD_NOTE:
        for i, ln in enumerate(lines, 1):
            if not re.search(pat, ln):
                continue
            ok = context_ok(lines, i - 1)
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
    A(f"Statistics: OK_NOTE={n_ok}, FLAG={n_flag}")
    A("")
    if flags:
        A("### FLAGs needing manual review/fix")
        A("")
        for i, meaning, s in flags:
            A(f"- L{i} ({meaning}): {s}")
        A("")
    else:
        A("No FLAG (all old values appear in comparison/historical context, or have been removed).")
        A("")

    A("## Claim-level checks")
    A("")
    A("| Check | Result |")
    A("|---|---|")
    abs_txt = t.split("\\end{abstract}")[0] if "\\end{abstract}" in t else t[:3000]
    abs_bad = ("[48,55]" in abs_txt) or ("0.0053" in abs_txt) or ("self-consistently selecting the fiducial" in abs_txt)
    A(f"| Abstract has no old window / old r / old fiducial | {'PASS' if not abs_bad else 'FAIL'} |")
    A(f"| Definition eq:Ndef exists | {'PASS' if 'eq:Ndef' in t else 'FAIL'} |")
    A(f"| Table tab:sens exists | {'PASS' if 'tab:sens' in t else 'FAIL'} |")
    A(f"| DE calibration statement | {'PASS' if 'calibrated' in t else 'FAIL'} |")
    A(f"| Residual quintessence rejected | {'PASS' if 'quintessence' in t.lower() and ('excluded' in t.lower() or 'excluded' in t) else 'FAIL'} |")
    A("")

    A("## Conclusion")
    A("")
    if n_flag == 0:
        A("- Number audit: **PASS** (old values appear only in convention-comparison / appendix historical notes).")
    else:
        A(f"- Number audit: **{n_flag} issue(s) found** (see the FLAG list above).")
    A("- Recommendation: merge this script into `run_all.py` and re-run it after every manuscript edit.")
    A("")
    A("[Audit complete]")

    OUT.write_text("\n".join(out), encoding="utf-8")
    print("\n".join(out))
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
