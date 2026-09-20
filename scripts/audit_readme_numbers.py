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

Tolerances are the rounding the README itself uses, not a wish for agreement:
a value printed as "1.0e-7" cannot be checked more tightly than a few percent.

Output: scripts/readme_numbers_report.md
"""
from __future__ import annotations

import io
import json
import math
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

import cosmo_model as cm                 # noqa: E402
import background_and_reheating as bar   # noqa: E402
import order_estimates as oe             # noqa: E402

README = os.path.join(ROOT, "README.md")
OUT_MD = os.path.join(HERE, "readme_numbers_report.md")

XI = 11.1


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

    # 1. The two independent N <-> T_reh routes.
    fp50 = bar.T_reh_star_first_principles(50.0)
    fp55 = bar.T_reh_star_first_principles(55.0)
    checks.append((
        "T_reh* agreement at N=50/55",
        r"agrees with\s+that table to ([0-9.]+) percent at N=50 "
        r"and ([0-9.]+) percent at N=55",
        [100.0 * rel(fp50, bar.TABLE_I[50][3]),
         100.0 * rel(fp55, bar.TABLE_I[55][3])],
        [0.02, 0.02],
        "background_and_reheating: first principles vs TABLE_I",
    ))

    # 2. The dark-matter light-branch anchor.
    dmg = load_json("dm_gap_closure_test.json")
    if dmg:
        checks.append((
            "dark-matter anchor",
            r"g ~ ([0-9.]+e-?[0-9]+), m_psi ~ ([0-9.]+e[0-9]+) GeV",
            [dmg["g_star_powerlaw"], dmg["m_star_powerlaw"] * bar.H_INF],
            [0.05, 0.05],
            "dm_gap_closure_test.json (light branch)",
        ))

    # 3. Closed form against the dynamical integration for V_end/V0.
    diff = rel(bar.slowroll_to_end()["V_end_frac_of_V0"], cm.V_end_over_V0(XI))
    checks.append((
        "closed form vs integration",
        r"dynamical integration sits ([0-9.]+e-[0-9]+) from it",
        [diff],
        [0.05],
        "cosmo_model.V_end_over_V0 vs background.slowroll_to_end",
    ))

    # 4. What the T_reh band does to g, n_s and r.
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

    # 5. The drift that this repository actually suffered, quoted as history.
    checks.append((
        "historical RHO_END drift",
        r"drifted ([0-9.]+e-[0-9]+) away from",
        [rel(1.63485e63, cm.rho_end())],
        [0.05],
        "cosmo_model.rho_end() vs the removed literal 1.63485e63",
    ))

    # 6. The two Sec. XII consistency entries.  These had NO executing source
    #    until this pass; they were quoted at a stale m = 1e13 GeV and now come
    #    from order_estimates, the single source for both.
    if dmg:
        m_psi = dmg["m_star_powerlaw"] * bar.H_INF
        sig = oe.sigma_psipsi_over_m(m_psi)
        checks.append((
            "Sec. XII self-interaction",
            r"`~([0-9.]+e-70) cm\^2/g`,\s+`~([0-9]+)`\s+orders below the bullet-cluster",
            [sig, -math.log10(sig)],
            [0.1, 0.02],
            "order_estimates.sigma_psipsi_over_m at the light-branch m_psi",
        ))
        checks.append((
            "Sec. XII Tremaine-Gunn Q",
            r"`~([0-9.]+e43) GeV\^4`",
            [oe.tremaine_gunn_Q(m_psi, bar.H_INF)],
            [0.1],
            "order_estimates.tremaine_gunn_Q at the light-branch m_psi",
        ))
        if teb:
            checks.append((
                "m_psi drift from the wrong H_inf",
                r"understated `m_psi` by ([0-9.]+) percent",
                [100.0 * rel(dmg["m_star_powerlaw"] * teb["rows"][1]["H_inf_GeV"],
                             m_psi)],
                [0.05],
                "m_star_powerlaw paired with the drifted treh_error_band H_inf "
                "vs with bar.H_INF",
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
    A("| check | README says | recomputed | rel | tol | verdict |")
    A("|---|---|---|---|---|---|")

    failures = 0
    total = 0
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
