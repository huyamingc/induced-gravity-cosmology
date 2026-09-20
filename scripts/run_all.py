# -*- coding: utf-8 -*-
"""Run all verification and figure scripts."""
from __future__ import annotations

import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Console safety.  Several scripts print non-ASCII text (lambda_0, "<=", Unicode
# subscripts, minus signs) and the default Windows console code page is GBK,
# which raises UnicodeEncodeError *after* the .md artifact has been written --
# so the artifact is fine but the script exits 1 and the traceback looks like a
# failure.  runpy executes every script in THIS process, so one reconfigure here
# covers all of them.  Running a single script directly still needs
# PYTHONIOENCODING=utf-8 (see the project-root README.md, section "Numerical scripts").
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):  # pragma: no cover - non-TextIOWrapper
        pass

SCRIPTS = [
    "derive_from_action.py",           # independent derivation from the action alone
    "lock_n_convention.py",            # exact locked-N grid; the source of Table I
    # Audit the .tex once the source of its tables has run.  The figure scripts
    # below write only .pdf, so neither the .tex nor n_convention_results.json
    # (which this reads) can change after this point -- a second call used to
    # sit after them and was an exact duplicate.
    "audit_tex_numbers.py",
    # Figures: fig2_ns_r.py uses locked exact PS and writes both filenames
    "fig1_einstein_potential.py",
    "fig2_ns_r.py",
    "fig3_domain_wall.py",
    # --- Independent first-principles checks added in the revision (read-only w.r.t. the .tex;
    #     each writes only its own .md/.json next to this file).  Every file's docstring records
    #     the requirement it closes and which P0 it supports.
    "background_and_reheating.py",     # exact KG integration; N window; reheating channels
    "psi_production_bogoliubov.py",    # de Sitter Bogoliubov index; Omega_psi matching; dilution
    "residual_quintessence.py",        # two-fluid Boltzmann integration; Delta w budget
    # SLOW tail.  Measured end-to-end on the author's machine (Python 3.13,
    # numpy 2.4 / scipy 1.18 / matplotlib 3.11): psi_abundance_oscillating ~5.5 min,
    # dm_gap_closure_test 10-15 min, so a full run_all.py is 15-20 min (the spread is
    # machine load).  dm_gap_closure_test imports psi_abundance_oscillating, so it runs
    # after it.
    "psi_abundance_oscillating.py",
    "dm_gap_closure_test.py",          # light-branch abundance matching + free-streaming check
    # Reads dm_gap_closure_test.json, so it must run after it.  Cheap (analytic + one
    # small brentq scan), and it is the script behind the T_reh error-band table of Sec. V.
    "treh_error_band.py",              # T_reh uncertainty -> (N, n_s, r) and (g, m_psi) band
    # --- Auditors.  All three read the .json artefacts written above AND the .tex
    #     table bodies, so they must run LAST.  verify_numerics.py recomputes every
    #     printed claim from the current scripts; consistency_checks.py checks every
    #     table cell against the script that produces it and every independent route
    #     against the other.  Together they take well under a second.
    "verify_numerics.py",              # paper claims -> code
    "consistency_checks.py",           # paper tables -> code, and code -> code
    "audit_readme_numbers.py",         # README prose numbers -> code
]


def main() -> None:
    # runpy executes each script in THIS process, so any script ending in
    # sys.exit(code) aborts the whole run at that point.  verify_numerics.py,
    # consistency_checks.py and audit_readme_numbers.py all do exactly that, so
    # the loop used to stop at verify_numerics.py: consistency_checks.py never
    # ran and "ALL SCRIPTS DONE" was never reached.  Catch SystemExit per script,
    # count the non-zero ones, and fail the run if any audit reported a problem.
    failed = []
    for name in SCRIPTS:
        print("=" * 60)
        print("RUN", name)
        print("=" * 60)
        try:
            runpy.run_path(str(HERE / name), run_name="__main__")
        except SystemExit as exc:
            code = exc.code
            if code not in (0, None):
                failed.append((name, code))
                print("!! %s exited with code %r" % (name, code))
    print()
    if failed:
        for name, code in failed:
            print("FAILED: %s (exit %r)" % (name, code))
        print("ALL SCRIPTS DONE with %d failing script(s)" % len(failed))
        raise SystemExit(1)
    print("ALL SCRIPTS DONE")


if __name__ == "__main__":
    main()
