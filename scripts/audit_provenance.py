#!/usr/bin/env python3
"""
Reverse provenance audit -- every scientific magnitude in the manuscript must
have a source or an explicit exemption
=============================================================================
Type:           AUDIT
Role:           STALE (reverse coverage); complements verify_numerics (CLAIM),
                consistency_checks (TABLE) and audit_tex_numbers (STALE)
Paper Sec.:     n/a (repository traceability)
Manuscript:     paper_prd_merged.tex
Produces:       scripts/provenance_coverage.json, scripts/provenance_coverage.md
Reads:          paper_prd_merged.tex, scripts/verify_numerics.py

What it does:   verify_numerics.py runs the FORWARD direction -- each
                registered claim recomputes a printed number -- but nothing
                forced every printed number to BE registered.  This script
                runs the REVERSE direction: it extracts every scientific
                magnitude from the manuscript (a x 10^b in any notation, plus
                bare 10^b with |b| >= 3, after masking \ref/\cite/\label) and
                demands that each one falls into one of

                    CLAIM      matched (by modulus, 2%) to a registered
                               claim(...)/claim_order(...) paper mirror
                    TABLE      sits inside a tabular checked cell-by-cell by
                               consistency_checks.py
                    EXTERNAL   declared observational / definitional /
                               parameter-choice input
                    HISTORICAL comparison or superseded-value wording (the
                               same +-2-line window audit_tex_numbers uses)
                    EXEMPT     explicitly listed below WITH a reason

                Anything else is an UNTRACED magnitude and fails the audit.
                The exemption table is the ONLY door for sourceless numbers.

Maintenance:    when you add a magnitude to the manuscript, either register a
                claim for it (preferred) or add an EXEMPT entry with a reason.
                When you correct a value, update the matching claim mirror.

History:        the first sweep (review_workspace/provenance_audit_report.md)
                found 177 untraced occurrences / 77 distinct values, including
                two genuine numerical errors; the registrations and exemptions
                below close that gap.
=============================================================================
"""
from __future__ import annotations

import ast
import io
import json
import os
import re
import sys
import time
from pathlib import Path

os.environ.setdefault("PYTHONUNBUFFERED", "1")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TEX = ROOT / "paper_prd_merged.tex"
VERIFY = HERE / "verify_numerics.py"
OUT_JSON = HERE / "provenance_coverage.json"
OUT_MD = HERE / "provenance_coverage.md"

MASK = re.compile(
    r"\\(?:ref|eqref|cite|citep|citet|label|input|include|usepackage|documentclass|"
    r"includegraphics|bibliography|bibliographystyle)\s*(?:\[[^\]]*\])?\{[^{}]*\}",
    re.S,
)
SCI = re.compile(r"(\d+(?:\.\d+)?)\s*\\times\s*10\s*\^\{?([+-]?\d+)\}?")
POW10 = re.compile(r"10\s*\^\{?([+-]?\d+)\}?")

EXTERNAL = {
    2.1e-9: "A_s (Planck 2018 input)",
    2.3e-5: "Cassini |gamma-1| bound (Bertotti 2003)",
    0.1179: "alpha_s(m_Z), PDG input",
    1e-4: "CODATA Newton-constant precision",
    2.36e-13: "T_0 (definition)",
    1e9: "fiducial T_reh (stated input)",
    1e8: "T_reh window edge (stated input)",
    1e-22: "g_chi_ff lower end (electron mass input)",
    1e-17: "g_chi_ff upper end (top mass input)",
    9.15e-5: "Omega_r (Planck input)",
    4.15e-5: "Omega_r h^2 (Planck input)",
    1e15: "T_reh scan upper edge (stated input)",
    2.435e18: "M_Pl reduced (CODATA definition)",
}
HISTORY = ("earlier", "Earlier", "attractor", "Analytic", "analytic",
           "Note.", "hypothetical", "stable-condensate", "stable condensate",
           "printed", "quoted in earlier", "large-field")

# The only door for sourceless magnitudes.  Each entry needs a reason; where a
# value is guarded by another auditor, the guard is named.
EXEMPT = {
    1e-3: "bound target 'Delta w << 1e-3' (App. D5) -- a threshold, not a result",
    1e-5: "bound targets in Discussion items (PPN gamma_J, Delta N_eff, beta_PPN)",
    1e-6: "bound target Delta xi <= 1e-6, enforced by the registered boolean claim",
    1e-8: "order of the graviton-exchange width Gamma^(ii) (Sec. IV, parametric)",
    1e-10: "bound target 'residual-condensate corrections <= 1e-10' (abstract/VI/D5)",
    1e-11: "loop-induced couplings g_eff ~ 1e-11 g and PPN beta (Sec. XII, parametric)",
    1e-12: "bound target of the unitarity claim (measured 7.7e-13) + preheating partition",
    1e-13: "parametric smallness (rho_part, curvature displacement, thermalization)",
    1e-14: "bound target, Higgs-Phi portal coupling (Sec. VIII)",
    1e-15: "bound target, gravitational-wave speed c_T (Sec. XII)",
    1e-20: "bound target, thermal mass correction delta m_Phi^2 (Sec. VI)",
    1e-25: "primordial magnetic seed B_0 (Sec. XII vi, parametric chain)",
    1e-27: "chi-nu coupling m_nu/M_Pl and Yukawa range in cm (Sec. II/IV/XII)",
    1e-29: "psi decay-width order (Sec. II G, parametric)",
    1e-33: "condensate inhomogeneity parameter (Sec. IV, parametric)",
    1e-38: "residual-oscillation amplitude chain (App. D5, parametric)",
    2e-14: "residual-amplitude precision (App. D5, parametric)",
    1e-42: "H_0 order of magnitude (definition context)",
    1e-46: "indirect-detection cross-section order (Sec. II G, parametric)",
    1e-55: "scalar-photon mixing angle (Sec. XII x, parametric)",
    1e-60: "CMB distortion bound (Sec. XII x, parametric)",
    1e-80: "CMB distortion / thermalization smallness (Sec. XII, parametric)",
    4e-11: "loop-induced Yukawa-reheating coupling (Sec. II E/IV, parametric)",
    3e-11: "HSR-PSR conversion difference (Sec. III, parametric)",
    3e-4: "NLO coefficient beta-dependence (Sec. III, parametric)",
    2.3e-13: "seed-bath parameter of the Ford-1982 estimate (Sec. IV a)",
    7.46e-8: "analytic large-field lambda0 comparison -- guarded by audit_tex_numbers",
    1.02e-92: "dilution factor -- guarded by consistency_checks shared constants",
    1e55: "decade rounding of the registered m_chi/H_0 = 2.3e55",
    1e10: "narrative T ~ 1e10 GeV in the combined reheating picture (Sec. IV)",
    1e4: "narrative scale (Sec. I related work; Sec. II E distances)",
    3e-42: "preheating energy-partition fraction (Sec. XII, parametric)",
    2e-13: "radiation-pulse temperature parameter (Sec. IV a, parametric)",
    1e6: "decade rounding of the registered Gamma_therm/H = 1.4e6",
    1e63: "decade rounding of the registered V_0 = 4.78e63",
    1e17: "decade rounding of the registered Phi_V = 7.31e17",
    7e19: "suppression ratio Gamma_anom/Gamma_tt, quoted from two registered values",
    3e-21: "psi decay-lifetime estimate tau_psi (Sec. II G, parametric)",
}


def claims():
    """Printed values registered in verify_numerics.py: args[3] of every
    claim()/claim_order() call, as literals or sub-constants."""
    tree = ast.parse(io.open(VERIFY, encoding="utf-8").read())
    out = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        f = node.func
        name = f.id if isinstance(f, ast.Name) else (
            f.attr if isinstance(f, ast.Attribute) else None)
        if name not in ("claim", "claim_order") or len(node.args) < 5:
            continue
        q = node.args[1].value if isinstance(node.args[1], ast.Constant) else "?"
        vals = []
        try:
            v = ast.literal_eval(node.args[3])
            if isinstance(v, (int, float)):
                vals.append(float(v))
        except Exception:
            for sub in ast.walk(node.args[3]):
                if isinstance(sub, ast.Constant) and isinstance(sub.value, (int, float)):
                    vals.append(float(sub.value))
        out.append((node.lineno, q, vals))
    return out


def tabular_lines(lines):
    inside, on = set(), False
    for i, ln in enumerate(lines, 1):
        if "\\begin{tabular" in ln:
            on = True
        if on:
            inside.add(i)
        if "\\end{tabular}" in ln:
            on = False
    return inside


def main() -> int:
    t0 = time.time()
    print("[%s] START: reverse provenance audit" % time.strftime("%H:%M:%S"))
    raw = io.open(TEX, encoding="utf-8").read()
    mlines = MASK.sub(lambda m: " " * len(m.group(0)), raw).splitlines()
    tlines = raw.splitlines()
    intab = tabular_lines(tlines)
    n = len(mlines)
    cl = claims()

    tokens = []
    for i, mline in enumerate(mlines, 1):
        taken = [False] * len(mline)
        for rx in (SCI, POW10):
            for m in rx.finditer(mline):
                if any(taken[m.start():m.end()]):
                    continue
                if rx is SCI:
                    v = float(m.group(1)) * 10.0 ** int(m.group(2))
                else:
                    e = int(m.group(1))
                    if abs(e) < 3:
                        continue
                    v = 10.0 ** e
                taken[m.start():m.end()] = [True] * (m.end() - m.start())
                ctx = " ".join(tlines[max(0, i - 3):min(n, i + 2)])
                tokens.append((i, m.group(0).strip(), v, ctx))

    stats = {"CLAIM": 0, "TABLE": 0, "EXTERNAL": 0, "HISTORICAL": 0,
             "EXEMPT": 0, "UNTRACED": 0}
    untraced = []
    for i, tok, v, ctx in tokens:
        hits = []
        for cl_no, q, vals in cl:
            for cv in vals:
                if cv and abs(abs(v) - abs(cv)) / abs(cv) < 0.02:
                    hits.append("verify L%d '%s' (%g)" % (cl_no, q[:40], cv))
        if i in intab:
            cat, ev = "TABLE", "consistency_checks.py tab-cell recompute"
        elif hits:
            cat, ev = "CLAIM", "; ".join(hits[:2])
        else:
            ex = next((lab for e, lab in EXTERNAL.items()
                       if e and abs(abs(v) - abs(e)) / abs(e) < 0.02), None)
            xp = next((r for e, r in EXEMPT.items()
                       if e and abs(abs(v) - abs(e)) / abs(e) < 0.02), None)
            if ex:
                cat, ev = "EXTERNAL", ex
            elif xp:
                cat, ev = "EXEMPT", xp
            elif any(k in ctx for k in HISTORY):
                cat, ev = "HISTORICAL", "comparison/superseded wording"
            else:
                cat, ev = "UNTRACED", "-- no mechanism --"
                untraced.append({"line": i, "token": tok, "value": v})
        stats[cat] += 1

    distinct = len({round(v, 300) for _, _, v, _ in tokens})
    total = len(tokens)
    covered = total - stats["UNTRACED"]

    out = []
    A = out.append
    A("# Reverse provenance coverage")
    A("")
    A("Every scientific magnitude of `paper_prd_merged.tex` (a x 10^b plus bare")
    A("10^b with |b| >= 3, cross-references masked) classified by source:")
    A("")
    A("| class | occurrences |")
    A("|---|---|")
    for k in ("CLAIM", "TABLE", "EXTERNAL", "HISTORICAL", "EXEMPT", "UNTRACED"):
        A("| %s | %d |" % (k, stats[k]))
    A("| **total** | **%d** |" % total)
    A("")
    A("- distinct values: **%d**" % distinct)
    A("- coverage: **%.1f%%** of occurrences outside the exemption class"
      % (100.0 * (total - stats["UNTRACED"] - stats["EXEMPT"]) / total))
    A("- untraced (failing): **%d**" % stats["UNTRACED"])
    A("")
    if untraced:
        A("## Untraced magnitudes")
        A("")
        for u in untraced:
            A("- L%d `%s` (%.3g)" % (u["line"], u["token"], u["value"]))
        A("")
    A("Exemption table lives in `audit_provenance.py` (EXEMPT), one reason each;")
    A("external inputs, tabular cells and historical wording are classified in")
    A("code.  Registered claims: %d call sites parsed from verify_numerics.py."
      % len(cl))
    A("")
    A("Runtime: %.2f s" % (time.time() - t0))

    OUT_MD.write_text("\n".join(out), encoding="utf-8")
    with io.open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump({"totals": stats, "total": total, "distinct": distinct,
                   "untraced": untraced, "claims_parsed": len(cl)}, fh, indent=1)
    print("\n".join(out))
    print("[%s] DONE: %d magnitudes, %d untraced"
          % (time.strftime("%H:%M:%S"), total, stats["UNTRACED"]))
    return 1 if untraced else 0


if __name__ == "__main__":
    raise SystemExit(main())
