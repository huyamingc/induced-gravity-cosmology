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
# PYTHONIOENCODING=utf-8 (see the project-root README.md, section "数值脚本").
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):  # pragma: no cover - non-TextIOWrapper
        pass

SCRIPTS = [
    "verify_numerics.py",
    "consistency_checks.py",
    "extended_checks.py",
    "derive_from_action.py",
    "n_definition_check.py",
    "lock_n_convention.py",
    "evaluate_paper.py",
    "audit_tex_numbers.py",
    "quick_claims_check.py",
    # Figures: fig2_ns_r.py uses locked exact PS and writes both filenames
    "fig1_einstein_potential.py",
    "fig2_ns_r.py",
    "fig3_domain_wall.py",
    # Re-audit after figures regenerate
    "audit_tex_numbers.py",
    # --- Independent first-principles checks added in the revision (read-only w.r.t. the .tex;
    #     each writes only its own .md/.json next to this file).  Every file's docstring records
    #     the requirement it closes and which P0 it supports.
    "background_and_reheating.py",     # exact KG integration; N window; reheating channels
    "psi_production_bogoliubov.py",    # de Sitter Bogoliubov index; Omega_psi matching; dilution
    "residual_quintessence.py",        # two-fluid Boltzmann integration; Delta w budget
    # SLOW tail.  Measured end-to-end on the author's machine (Python 3.13,
    # numpy 2.4 / scipy 1.18 / matplotlib 3.11): psi_abundance_oscillating ~5.5 min,
    # dm_gap_closure_test ~19 min, so a full run_all.py is ~25 min.  dm_gap_closure_test
    # imports psi_abundance_oscillating, so it runs after it.
    "psi_abundance_oscillating.py",
    "dm_gap_closure_test.py",          # light-branch abundance matching + free-streaming check
]


def main() -> None:
    for name in SCRIPTS:
        print("=" * 60)
        print("RUN", name)
        print("=" * 60)
        runpy.run_path(str(HERE / name), run_name="__main__")
    print("\nALL SCRIPTS DONE")


if __name__ == "__main__":
    main()
