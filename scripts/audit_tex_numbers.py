# -*- coding: utf-8 -*-
"""Audit paper_prd_merged.tex scientific numbers against locked-N convention.

Type:           PAPER
Primary source: scripts/n_convention_results.json (if present) + locked table values.
Flags:
  OK_EXACT   - locked-N primary values
  OK_NOTE    - allowed only as attractor/historical comparison (context window)
  FLAG       - old-convention value used as if current prediction
Exit code:      0 when no FLAG, 1 when any FLAG -- so run_all.py's failure
                detection covers this audit the same way as the other auditors.
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
    "N_anomaly": 50.7,
    "N_band": (45, 55.6),
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
    # Stale hand-written value in the App. A verbatim listing: Table I and
    # n_convention_results.json (T_reh_selfcons at N=51) give 2.2e9 GeV, and
    # no script in scripts/ has ever produced 1.4e9.  A hand-copied literal
    # inside a code listing is exactly how the 0.285204 incident happened.
    (r"T_reh\*\(N=51\)\s*=\s*1\.4",
     "stale verbatim-listing T_reh*(N=51); Table I gives 2.2e9"),
    # Superseded route-change factor: the tabulated T_reh* ratio across the
    # 0fe113c matching-route change is 4.4e7 -> 1.1e8 GeV (about 2.5), not
    # 2.04; the e-fold equivalent follows from the e^{3.02 N} grid scaling.
    (r"factor \$2\.04\$",
     "superseded route-change factor; tabulated ratio is about 2.5"),
    (r"0\.23\$ \$e\$-folds",
     "stale e-fold equivalent; the 2.5 ratio gives about 0.30"),
    # Section pointers: after the reorganization the SM-embedding chapter is
    # Sec. VII and the domain-wall chapter is Sec. VIII.
    (r"Secs\.~IV and VIII",
     "stale section pointer; the chi-gauge discussion is in Secs. IV and VII"),
    # --- Anchoring round: hand-written magnitudes now have executing sources
    # in order_estimates.py; these patterns guard against regression. ---
    # Z2^psi is an ASSUMED anomaly-free UV remnant; only the vanishing vertex
    # follows from F = xi Phi^2.  The abstract once attributed both to F.
    (r"both consequences of \$F=\\xi\\Phi\^2\$",
     "Z2^psi is an assumed UV remnant, not a consequence of F=xi Phi^2"),
    # g-window lower edge must equal the tab:trehband T=1e9 row (4.7e-8).
    (r"0\.5\$--\$1\.5\\times10\^\{-7\}",
     "stale g-window rounding; the T=1e9 table row gives 4.7e-8"),
    # lambda0 in the lattice-parameters paragraph drifted to 6.6e-8.
    (r"lambda_0\\simeq6\.6\\times10\^\{-8\}",
     "stale lambda0 rounding; the locked N=50 value is 6.70e-8"),
    # Stale Higgs-portal magnitude (one-loop formula gives ~6e-17).
    (r"\\lambda_\{\\Phi H\}\\sim10\^\{-15\}",
     "stale portal magnitude; the one-loop formula gives ~6e-17"),
    (r"\\MP\^2\)\\sim10\^\{-15\}",
     "stale portal magnitude; the one-loop formula gives ~6e-17"),
    # Stale heavy-field threshold bound: (m_chi/Lambda_J)^2 = 2.2e-8.
    (r"\\Lambda_J\^2\)\\leq10\^\{-10\}",
     "stale threshold bound; (m_chi/Lambda_J)^2 = 2.2e-8"),
    # Discussion #14 stale quartic chain.
    (r"10\^\{-6\}\\varphi\^2",
     "stale quartic coefficient; lambda0/xi^2 = 5.4e-10"),
    (r"H_0/m_\\chi\\sim10\^\{-46\}",
     "stale frozen-field value; Ricci driving gives (H0/m_chi)^2 ~ 2e-111"),
    (r"correction is \$\\sim10\^\{-98\}",
     "stale quartic correction; the anchored chain gives ~2e-231"),
    # Thermalization: alpha_s ~ 0.1 was a hand estimate; the one-loop running
    # value is 0.0377 and Gamma_therm/H ~ 1e6.
    (r"Gamma_\{\\rm therm\}/H\\sim10\^7",
     "stale thermalization ratio; one-loop running alpha_s gives ~1e6"),
    (r"alpha_s\^2 T\\sim10\^7",
     "stale Gamma_therm magnitude; alpha_s(1e9)=0.038 gives ~1.4e6 GeV"),
    # lock_n_convention.py no longer plots; it provides the locked-N grid.
    (r"same as \\texttt\{lock\\_n\\_convention\.py\}",
     "stale caption wording; the script provides the locked-N grid, not plots"),
    # Stale single-cause attribution of the 2.5 route factor.
    (r"because \$g_\{\*s\}\$ falls from \$106\.75\$",
     "stale attribution; the 2.5 factor is the combined entropy-matching/"
     "rho_end/Omega_r effect"),
    # 'selects' next to the raw band (the L835 conclusion occurrence keeps
    # 'selects' legitimately, with its statistical-band parenthesis).
    (r"selects \$N\\approx45\$--\$56\$ with \$r\\in\[0\.0034",
     "stale 'selects' wording; Planck allows the band on n_s alone"),
    (r"this selects \$N\\simeq45\$--\$56\$ by Planck",
     "stale 'selects' wording; use 'allows' with the N_max cap"),
    # Over-strong observational-handle wording (lattice still pending).
    (r"only non-degenerate observational handle",
     "over-strong wording; observational connection, normalization pending"),
    (r"a measured \$g\$ would translate into a prediction",
     "over-strong wording; g is not directly measurable (gravitational only)"),
    (r"one non-degenerate observable handle of the framework",
     "over-strong wording; relation pending the lattice normalization"),
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
        ("N anomaly", r"N\\simeq50\.7"),
        ("m_chi", r"3\.25"),
        ("H_inf", r"1\.64"),
        ("N band low", r"45"),
        ("N band high", r"56"),
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

    A("## Window co-occurrence (N band vs N_max cap)")
    A("")
    A("Every quote of the $n_s$ band `45--56` must sit next to the physical cap")
    A("`N_max ~= 55.6` (Appendix D), otherwise a reader of the abstract or a")
    A("table caption alone concludes that N = 56 is viable, contradicting")
    A("`N_max = 55.6` and the exclusion of N = 56, 57, 58.  For the abstract and")
    A("the tab:sens caption the cap is REQUIRED on the same block; elsewhere it")
    A("is reported as info (those lines are covered by the Appendix D argument).")
    A("")
    A("| Line | Quote | Cap nearby | Requirement | Verdict |")
    A("|---|---|---|---|---|")
    band_re = re.compile(r"45\$?--\$?56|45\\text\{--\}56")
    cap_re = re.compile(r"55\.6|55\.9|N\\in\[45")
    window_flags = []
    n_window = 0
    for i, ln in enumerate(lines, 1):
        if not band_re.search(ln):
            continue
        n_window += 1
        lo = max(0, i - 1 - CONTEXT_WINDOW)
        hi = min(len(lines), i - 1 + CONTEXT_WINDOW + 1)
        near_sens_label = any("tab:sens" in lines[j] for j in range(lo, hi))
        is_mandatory = ("\\abstract{" in ln) or ("\\caption{" in ln and near_sens_label)
        capped = any(cap_re.search(lines[j]) for j in range(lo, hi))
        verdict = "OK" if capped else ("**FLAG**" if is_mandatory else "info")
        if is_mandatory and not capped:
            window_flags.append(i)
        A(f"| L{i} | `45--56` | {'yes' if capped else 'no'} "
          f"| {'required' if is_mandatory else 'optional'} | {verdict} |")
    A("")
    A(f"Statistics: band quotes={n_window}, mandatory-missing-cap={len(window_flags)}")
    A("")
    if window_flags:
        A("### Window FLAGs needing manual review/fix")
        A("")
        for i in window_flags:
            A(f"- L{i}: band quote `45--56` without an `N_max ~= 55.6` cap in the "
              "abstract/caption block")
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
    total_issues = n_flag + len(window_flags)
    if total_issues == 0:
        A("- Number audit: **PASS** (old values appear only in convention-comparison / appendix historical notes; every abstract/caption band quote carries the N_max cap).")
    else:
        A(f"- Number audit: **{total_issues} issue(s) found** ({n_flag} old-value FLAG(s), {len(window_flags)} window FLAG(s); see the FLAG lists above).")
    A("- Recommendation: merge this script into `run_all.py` and re-run it after every manuscript edit.")
    A("")
    A("[Audit complete]")

    OUT.write_text("\n".join(out), encoding="utf-8")
    print("\n".join(out))
    print("Wrote", OUT)
    return 0 if total_issues == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
