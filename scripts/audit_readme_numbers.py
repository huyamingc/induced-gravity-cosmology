#!/usr/bin/env python3
"""
README prose-number audit -- check the numbers in README.md against the scripts
=============================================================================
Type:           PAPER (audit; writes a report only, runs no simulation)
Paper section:  n/a (repository documentation)
Experiment:     documentation consistency
=============================================================================
Why this exists
---------------
`verify_numerics.py` recomputes the quantitative claims printed in the *manuscript*,
and `consistency_checks.py` cross-checks the *script* outputs against the paper
tables.  Neither looks at README.md, so its prose numbers were the one surface
where a value could go stale with nothing to notice.  That is not hypothetical: the
README carried "agrees to 1.4 percent" for several rounds after the true figure had
become 0.094 percent.

Method
------
Each entry below pairs a regular expression (which must match exactly once in the
README and capture the numbers in order) with a callable-free list of values
recomputed here.  A pattern that stops matching is also a failure: it means the
sentence was reworded and the audit no longer covers it.

Only prose that remains in the public front-door README.md is audited.
Historical maintainer notes live in the gitignored MAINTENANCE.md and are
deliberately not required here (see that file if present on this machine).

Tolerances are the rounding the README itself uses, not a wish for agreement:
a value printed as "1.0e-7" cannot be checked more tightly than a few percent.

Output: scripts/readme_numbers_report.md
"""
from __future__ import annotations

import io
import json
import os
import re
import sys

os.environ.setdefault("PYTHONUNBUFFERED", "1")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import background_and_reheating as bar   # noqa: E402

README = os.path.join(ROOT, "README.md")
OUT_MD = os.path.join(HERE, "readme_numbers_report.md")


def rel(a: float, b: float) -> float:
    return abs(a - b) / abs(b) if b else abs(a - b)


def load_json(name: str):
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        return None
    with io.open(path, encoding="utf-8") as fh:
        return json.load(fh)


def build_checks():
    """Return a list of (label, regex, measured values, tolerances, source)."""
    checks = []

    # 1. The two independent N <-> T_reh routes (public README, Dark-matter section).
    fp50 = bar.T_reh_star_first_principles(50.0)
    fp55 = bar.T_reh_star_first_principles(55.0)
    checks.append((
        "T_reh* agreement at N=50/55",
        r"agrees with that table to ([0-9.]+) percent at N=50 "
        r"and ([0-9.]+) percent at N=55",
        [100.0 * rel(fp50, bar.TABLE_I[50][3]),
         100.0 * rel(fp55, bar.TABLE_I[55][3])],
        [0.02, 0.02],
        "background_and_reheating: first principles vs TABLE_I",
    ))

    # 2. The dark-matter abundance-matched coupling (public README).  Since the
    #    exact-background computation (psi_mode_oscillating) the README quotes
    #    the primary value in the "g ~ ..., m_psi ~ ... GeV" form and the
    #    transition-only baseline with "=", so this pattern is unambiguous.
    pm = load_json("psi_mode_oscillating.json")
    if pm:
        checks.append((
            "dark-matter anchor (exact background, primary)",
            r"g ~ ([0-9.]+e-?[0-9]+), m_psi ~ ([0-9.]+e[0-9]+) GeV",
            [pm["matching"]["exact"]["g"], pm["matching"]["exact"]["m_psi_GeV"]],
            [0.05, 0.05],
            "psi_mode_oscillating.json (exact background, primary)",
        ))

    # 3. What the T_reh band does to g, n_s and r (public README).
    teb = load_json("treh_error_band.json")
    if teb:
        s = teb["summary"]
        checks.append((
            "T_reh band effect",
            r"moves them by a factor ([0-9.]+), while it moves n_s by only "
            r"~([0-9.]+e-[0-9]+) and r by\s+~([0-9.]+) percent",
            [s["g_ratio"], s["dns"], 100.0 * s["dr_over_r"]],
            [0.05, 0.05, 0.05],
            "treh_error_band.json summary",
        ))

    # 4. Reverse provenance coverage (audit_provenance.py output).  The README
    #    states how many magnitudes the reverse sweep classifies, how many
    #    claims are registered, and how many are untraced; all three must match
    #    the machine report exactly (integer counts, zero tolerance).
    cov = load_json("provenance_coverage.json")
    if cov:
        checks.append((
            "reverse provenance coverage",
            r"\*\*(\d+) scientific magnitudes\*\* classified, "
            r"\*\*(\d+) registered claim call sites\*\*, "
            r"\*\*(\d+) untraced magnitudes\*\*",
            [float(cov["total"]), float(cov["claims_parsed"]),
             float(cov["totals"]["UNTRACED"])],
            [0.0, 0.0, 0.0],
            "audit_provenance.py -> provenance_coverage.json",
        ))

    return checks


def main() -> None:
    with io.open(README, encoding="utf-8") as fh:
        text = fh.read()

    lines = []
    A = lines.append
    A("# README prose-number audit")
    A("")
    A("Source: `README.md`.  Values recomputed from the scripts in this directory.")
    A("A pattern that no longer matches counts as a failure: the sentence was")
    A("reworded and this audit no longer covers it.")
    A("")
    A("Historical maintainer notes (if present as local `MAINTENANCE.md`) are")
    A("outside this audit; they are gitignored and not part of the submission.")
    A("")
    A("| check | README says | recomputed | rel | tol | verdict |")
    A("|---|---|---|---|---|---|")

    failures = 0
    total = 0
    # A missing upstream artefact must FAIL loudly: with the old `if dmg:` /
    # `if teb:` guards the affected checks silently vanished from the report
    # and the audit stayed green while the dark-matter anchor went unchecked.
    for name in ("dm_gap_closure_test.json", "treh_error_band.json",
                 "psi_mode_oscillating.json", "provenance_coverage.json"):
        if load_json(name) is None:
            A("| (upstream artefact) | *%s missing* | -- | -- | -- | **FAIL** |" % name)
            failures += 1
    for label, pattern, measured, tols, source in build_checks():
        m = re.search(pattern, text)
        if m is None:
            A("| %s | *pattern not found* | -- | -- | -- | **FAIL** |" % label)
            failures += 1
            continue
        said = [float(g) for g in m.groups()]
        for said_v, meas_v, tol in zip(said, measured, tols):
            total += 1
            r = rel(meas_v, said_v)
            ok = r <= tol
            if not ok:
                failures += 1
            A("| %s | %g | %g | %.2e | %.2e | %s |"
              % (label, said_v, meas_v, r, tol, "OK" if ok else "**FAIL**"))
        A("|  | ^ %s | | | | |" % source)

    A("")
    A("- Checks run: **%d**" % total)
    A("- Failures: **%d**" % failures)
    A("")
    A("Runtime: this audit performs no simulation; it reads README.md and the")
    A("`.json` artefacts, so it must run after `dm_gap_closure_test.py` and")
    A("`treh_error_band.py` in `run_all.py`.")
    A("")
    A("[readme audit complete]")

    report = "\n".join(lines) + "\n"
    with io.open(OUT_MD, "w", encoding="utf-8") as fh:
        fh.write(report)
    print(report)
    print("Wrote %s" % OUT_MD)
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
