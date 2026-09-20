#!/usr/bin/env python3
"""
Cross-script consistency ledger -- the manuscript against the scripts, and the
scripts against each other
=============================================================================
Type:           PAPER
Paper Sec.:     III (Table I), V (Table 2), App. D0
Experiment:     consistency ledger
What it does:   Two directions that no other script covers.
                (A) paper -> code: every numeric cell of Table I and Table 2 is
                    parsed out of the .tex and recomputed from the script that
                    produces it.  A stale table cell is caught without anyone
                    having to remember to look.
                (B) code -> code: the same physical quantity computed by two
                    independent scripts is compared, with the residual reported
                    against the uncertainty the manuscript declares.
Inputs:         paper_prd_merged.tex, dm_gap_closure_test.json,
                treh_error_band.json
Outputs:        consistency_report.md, consistency_report.json
Dependencies:   must run AFTER every script that writes a .json artefact
                (run_all.py places it last for this reason).
=============================================================================
"""
import io
import json
import math
import os
import re
import sys
import time
from pathlib import Path

os.environ.setdefault("PYTHONUNBUFFERED", "1")
if hasattr(sys.stdout, "reconfigure"):          # Python 3.7+
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TEX = ROOT / "paper_prd_merged.tex"
OUT_MD = HERE / "consistency_report.md"
OUT_JSON = HERE / "consistency_report.json"

sys.path.insert(0, str(HERE))
import cosmo_model as cm                # noqa: E402
import derive_from_action as dfa        # noqa: E402
import lock_n_convention as lk          # noqa: E402
import background_and_reheating as bar  # noqa: E402
import dm_gap_closure_test as dg        # noqa: E402
import psi_abundance_oscillating as pa  # noqa: E402
import psi_production_bogoliubov as pb  # noqa: E402
import residual_quintessence as rq      # noqa: E402

XI = 11.1
N_FID = 50.0


# ---------------------------------------------------------------- helpers
def load_json(name):
    p = HERE / name
    if not p.exists():
        return None
    with io.open(p, encoding="utf-8") as fh:
        return json.load(fh)


def num(tok):
    """First number in a LaTeX table cell: '$7.24\\times10^{-8}$' -> 7.24e-08."""
    s = tok.replace("$", "").replace("{", "").replace("}", "")
    m = re.search(r"([\d.]+)\\times10\^?(-?\d+)", s)
    if m:
        return float(m.group(1)) * 10.0 ** int(m.group(2))
    m = re.search(r"-?\d+\.?\d*", s)
    return float(m.group(0)) if m else float("nan")


def rel(a, b):
    if b == 0.0:
        return abs(a - b)
    return abs(a - b) / abs(b)


def table_rows(text, label):
    """Numeric rows of the tabular that follows \\label{<label>}."""
    i = text.index("\\label{%s}" % label)
    j = text.index("\\end{tabular}", i)
    seg = text[i:j].replace("\\hline", "")
    rows = []
    for chunk in seg.split("\\\\"):
        if "&" not in chunk:
            continue
        cells = [c.strip() for c in chunk.split("&")]
        vals = [num(c) for c in cells]
        if not vals or math.isnan(vals[0]):
            continue
        rows.append(vals)
    return rows


# ---------------------------------------------------------------- ledger
LEDGER = []


def record(group, name, a_src, a_val, b_src, b_val, tol, kind="identity", note=""):
    """kind: 'identity' same formula+inputs (tol ~1e-9);
             'rounding' the .tex cell is quoted to 2-3 significant figures;
             'route'    two different implementations of the same physics;
             'precision' the two sources carry the same constant to different digits."""
    ok = rel(a_val, b_val) <= tol
    LEDGER.append({
        "group": group, "name": name, "kind": kind,
        "a_src": a_src, "a_val": a_val,
        "b_src": b_src, "b_val": b_val,
        "rel": rel(a_val, b_val), "tol": tol, "ok": ok, "note": note,
    })
    return ok


def main():
    t0 = time.time()
    print("[%s] START: cross-script consistency ledger" % time.strftime("%H:%M:%S"))
    text = TEX.read_text(encoding="utf-8")
    teb = load_json("treh_error_band.json")
    dmg = load_json("dm_gap_closure_test.json")

    p50 = lk.point(XI, N_FID)

    # ---- A. constants shared by several scripts -------------------------
    if teb:
        record("constants", "Phi_V [GeV]", "treh_error_band.json", teb["Phi_V_GeV"],
               "M_Pl/sqrt(xi)", dfa.M_PL / math.sqrt(XI), 1e-5, "precision",
               note="treh_error_band carries M_Pl = 2.435015e18, derive_from_action "
                    "2.435e18: a 6 ppm difference in a shared constant, not a physics gap")
        record("constants", "M_Pl [GeV]", "treh_error_band.json",
               teb["M_Pl_reduced_GeV"], "derive_from_action.M_PL", dfa.M_PL, 1e-5,
               "precision",
               note="same 6 ppm constant-precision difference")
    # lam0_of_N was rewired onto the exact A_s inversion (derive_from_action, via
    # cosmo_model.lambda0_for_As_locked), so the three entries below are now the SAME
    # formula evaluated in two places rather than two independent routes.  While
    # lam0_of_N was the fitted closed form 6.70e-8*(50/N)^2 they were "route" checks
    # with a 4.46e-4 spread.
    record("constants", "V0 [GeV^4]", "background.V0_of_N(50)", bar.V0_of_N(N_FID),
           "lock_n_convention.point V0", p50["V0"], 1e-9, "identity",
           note="same exact lambda0 and the same V0 = lambda0 M_Pl^4/(4 xi^2); the "
                "former 4.46e-4 route spread was entirely the fitted-vs-exact lambda0")
    record("constants", "V_end [GeV^4]", "background.V_end_of_N(50)",
           bar.V_end_of_N(N_FID), "lock_n_convention.point V_end",
           p50["V_end"], 1e-9, "identity",
           note="both sides now evaluate the same closed form V_end/V0 = (1-u_e)^2 "
                "with u_e = 1/(1+sqrt(2) beta_p), so this is an identity rather than a "
                "route check.  The former 1.4e-6 residual was the 6-digit literal "
                "0.285204; note that the dynamical integration itself sits 2.6e-6 away "
                "from that closed form, so freezing the integral would have been no "
                "better than the literal it replaced")
    record("constants", "H_inf [GeV]", "background.H_inf_of_N(50)",
           bar.H_inf_of_N(N_FID), "lock_n_convention.point H_inf",
           p50["H_inf"], 1e-9, "identity",
           note="algebraically identical (M_Pl sqrt(lambda0)/(2 sqrt3 xi)) once lambda0 agrees")

    # ---- end-of-inflation ratios: one source, guarded -------------------
    # These exist to stop a hand-copied literal from creeping back.  The numbers
    # 0.285204 / 0.19938 / 1.19938 used to sit in four scripts, and
    # psi_abundance_oscillating.RHO_END had silently kept an OLD-lambda0 value.
    sr = bar.slowroll_to_end()
    record("end of inflation", "V_end/V0", "cosmo_model.V_end_over_V0(11.1)",
           cm.V_end_over_V0(XI), "background.slowroll_to_end (RK4)",
           sr["V_end_frac_of_V0"], 5e-6, "route",
           note="closed form (1-u_e)^2 against the dynamical integration.  The "
                "integration converges onto the closed form only to 2.6e-6; that "
                "residual is the resolution of the integration, and it is the closed "
                "form that is exact")
    record("end of inflation", "rho_end/V_end", "cosmo_model.rho_end_over_V_end()",
           cm.rho_end_over_V_end(), "1 + eps_H/(3-eps_H), integrated",
           1.0 + sr["eps_H_end"] / (3.0 - sr["eps_H_end"]), 1e-9, "identity")
    pab = load_json("psi_abundance_oscillating.json")
    if pab:
        record("end of inflation", "rho_end [GeV^4]",
               "psi_abundance_oscillating.json", pab["constants"]["rho_end"],
               "(V_end/V0)(rho_end/V_end)V0(50), exact lambda0",
               cm.V_end_over_V0(XI) * cm.rho_end_over_V_end()
               * cm.V0(cm.lambda0_for_As_locked(N_FID, XI), XI), 1e-9, "identity",
               note="this script used to carry an isolated literal 1.63485e63 that "
                    "was computed under the OLD fitted lambda0 = 6.70e-8 and had "
                    "drifted 4.3e-4 from the current exact value; it now derives "
                    "from the exact lambda0")

    # ---- end-of-inflation scalars: one source, guarded ------------------
    # x_end is a closed form in cosmo_model, and the dilution factor now comes from
    # cosmo_model.dilution.  Both used to be literals in two or three scripts
    # (0.7636653 / 0.76367, 1.0204e-92, 1.0204e-101 twice).
    record("end of inflation", "x_end", "background.X_END", bar.X_END,
           "cosmo_model.x_end(11.1)", cm.x_end(XI), 1e-12, "identity")
    record("end of inflation", "x_end", "psi_production.X_END", pb.X_END,
           "cosmo_model.x_end(11.1)", cm.x_end(XI), 1e-12, "identity")
    record("end of inflation", "dilution at T_reh=1e9",
           "dm_gap.DIL_PER_GEV * 1e9", dg.DIL_PER_GEV * 1e9,
           "cosmo_model.dilution(1e9)", cm.dilution(1e9), 1e-12, "identity")
    record("end of inflation", "dilution at T_reh=1e9",
           "residual.A_END_OVER_A0_CUBED_1E9", rq.A_END_OVER_A0_CUBED_1E9,
           "cosmo_model.dilution(1e9)", cm.dilution(1e9), 1e-12, "identity")

    # ---- derived scales: paper-display references, guarded --------------
    # Several self-contained cross-check scripts quote a scale straight from the
    # manuscript and keep the manuscript's digits (H_inf = 1.6388e13 and so on).
    # That is deliberate, but it must not silently drift -- going stale unnoticed
    # is exactly how psi_abundance's RHO_END survived until the previous round.
    # Each tolerance is the rounding that script declares, not a wish for agreement.
    lam0_x = cm.lambda0_for_As_locked(N_FID, XI)
    scales = [
        ("lambda0", "psi_production.LAM0", pb.LAM0, lam0_x, 1e-9),
        ("V0 [GeV^4]", "psi_production.V0", pb.V0, cm.V0(lam0_x, XI), 1e-9),
        ("V0 [GeV^4]", "residual.V0", rq.V0, cm.V0(lam0_x, XI), 1e-5),
        ("H_inf [GeV]", "psi_production.H_INF", pb.H_INF,
         cm.H_inf(lam0_x, XI), 1e-9),
        ("H_inf [GeV]", "psi_abundance.H_INF", pa.H_INF,
         cm.H_inf(lam0_x, XI), 2e-7),
        ("H_inf [GeV]", "residual.H_INF", rq.H_INF, cm.H_inf(lam0_x, XI), 2e-7),
        ("m_chi [GeV]", "psi_production.M_CHI", pb.M_CHI,
         cm.m_chi(lam0_x, XI), 1e-9),
        ("m_chi [GeV]", "residual.M_CHI", rq.M_CHI, cm.m_chi(lam0_x, XI), 2e-5),
        ("Phi_V [GeV]", "psi_abundance.PHI_V", pa.PHI_V,
         cm.M_P / math.sqrt(XI), 1e-5),
    ]
    for name, label, val, exact, tol in scales:
        record("derived scales", "%s (%s)" % (name, label), label, val,
               "single source", exact, tol, "precision",
               note="this script quotes the manuscript's own digits; the tolerance "
                    "is that declared rounding, so a re-truncation or a value going "
                    "stale against a moved lambda0 fails here")

    # ---- A. Table I (tab:sens) vs lock_n_convention ---------------------
    # tolerances follow the rounding of the .tex cell, not a wish for agreement
    tol_i = [3e-2, 5e-3, 1e-4, 3e-3, 5e-3, 5e-2]   # N, lam0, ns, r, H_inf, T_reh*
    cols = ["N", "lambda0", "n_s", "r", "H_inf", "T_reh*"]
    n_tab1 = 0
    for vals in table_rows(text, "tab:sens"):
        if len(vals) < 6:
            continue
        N = vals[0]
        p = lk.point(XI, N)
        computed = [N, p["lambda0"], p["ns"], p["r"], p["H_inf"],
                    lk.T_reh_for_N_derived(p, N)]
        for k in range(1, 6):
            record("Table I", "%s at N=%d" % (cols[k], int(N)),
                   "paper_prd_merged.tex", vals[k],
                   "lock_n_convention", computed[k], tol_i[k], "rounding")
        n_tab1 += 1

    # ---- A. Table 2 (tab:trehband) vs treh_error_band -------------------
    # g, m_psi and lambda_fs are quoted to two significant figures in the table,
    # so a 5e-2 tolerance is the most that column can support
    tol_t = [1e-2, 1e-3, 1e-4, 3e-3, 5e-2, 5e-2, 5e-2]  # T,N,ns,r,g,m,lfs
    n_tab2 = 0
    if teb:
        rows2 = teb["rows"]
        for vals in table_rows(text, "tab:trehband"):
            if len(vals) < 7:
                continue
            T = vals[0]
            match = min(rows2, key=lambda r: abs(math.log(r["T_reh_GeV"] / T)))
            comp = [match["T_reh_GeV"], match["N"], match["ns"], match["r"],
                    match["g"], match["m_psi_GeV"], match["lambda_fs_Mpc"]]
            labels = ["T_reh", "N", "n_s", "r", "g", "m_psi", "lambda_fs"]
            for k in range(1, 7):
                record("Table 2", "%s at T_reh=%.1e" % (labels[k], T),
                       "paper_prd_merged.tex", vals[k],
                       "treh_error_band", comp[k], tol_t[k], "rounding")
            n_tab2 += 1

    # ---- B. the two independent N <-> T_reh routes ----------------------
    T_lock = lk.T_reh_for_N_derived(p50, N_FID)
    T_indep = bar.T_reh_star_first_principles(N_FID)
    record("routes", "T_reh*(50) [GeV]", "lock_n_convention.T_reh_for_N_derived",
           T_lock, "background.T_reh_star_first_principles", T_indep, 1e-2,
           "route", note="independent implementations of Eq. (matchN)")
    for N in (48.0, 52.0, 55.0):
        p = lk.point(XI, N)
        record("routes", "T_reh*(%d) [GeV]" % int(N),
               "lock_n_convention", lk.T_reh_for_N_derived(p, N),
               "background first principles", bar.T_reh_star_first_principles(N),
               1e-2, "route")

    # ---- B. N_max two ways ---------------------------------------------
    nw = bar.solve_N_max()
    fit = nw["fit"]
    from scipy.optimize import brentq  # noqa: E402
    direct = brentq(lambda N: bar.T_reh_star_from_table(N, fit) - bar.T_reh_inst(N),
                    50.0, 70.0, xtol=1e-12)
    record("N window", "N_max (proper rho_end ceiling)", "solve_N_max()",
           nw["N_max_proper"], "direct root find", direct, 1e-6)
    record("N window", "rho_end/V_end at N=50",
           "background.rho_end_of_N/V_end_of_N",
           bar.rho_end_of_N(N_FID) / bar.V_end_of_N(N_FID),
           "paper_prd_merged.tex (App. D0)", 1.1994, 1e-3, "rounding")
    record("N window", "K_end/V_end at N=50", "1.1994 - 1",
           bar.rho_end_of_N(N_FID) / bar.V_end_of_N(N_FID) - 1.0,
           "paper_prd_merged.tex (Sec. III)", 0.199, 1e-2, "rounding")

    # ---- B. dark-matter sector -----------------------------------------
    if dmg and teb:
        g_star = dmg["g_star_powerlaw"]
        m_star = dmg["m_star_powerlaw"]
        phi_V = teb["Phi_V_GeV"]
        H_inf = dg.H_INF
        T_model = dmg["T_reh_model"]
        record("dark matter", "m_psi [GeV] via Phi_V", "g_star * Phi_V",
               g_star * phi_V, "m_star * H_inf", m_star * H_inf, 1e-6)
        record("dark matter", "m_psi/H_inf (Eq. 33)",
               "(g/sqrt(xi)) * M_Pl/H_inf",
               (g_star / math.sqrt(XI)) * (dfa.M_PL / H_inf),
               "dm_gap_closure_test.json m_star_powerlaw", m_star, 1e-5,
               "precision", note="inherits the same 6 ppm M_Pl precision difference")
        record("dark matter", "Omega_psi target",
               "dm_gap_closure_test.OMEGA_TARGET", dg.OMEGA_TARGET,
               "manuscript Omega_DM", 0.265, 1e-3)
        # identity: recompute the anchor with the script's own function
        record("dark matter", "lambda_fs at the anchor [Mpc]",
               "dm_gap_closure_test.free_streaming_length",
               dg.free_streaming_length(m_star, T_model)["lambda_fs_Mpc"],
               "dm_gap_closure_test.json free_streaming",
               dmg["free_streaming"]["lambda_fs_Mpc"], 1e-9)
        # route: the nearest Tab. 2 grid row is not exactly the anchor point
        anchor_row = min(teb["rows"], key=lambda r: abs(
            math.log(r["T_reh_GeV"] / T_model)))
        record("dark matter", "anchor vs nearest Tab. 2 row [Mpc]",
               "dm_gap_closure_test.json", dmg["free_streaming"]["lambda_fs_Mpc"],
               "treh_error_band nearest grid row", anchor_row["lambda_fs_Mpc"],
               5e-3, "route",
               note="the Tab. 2 grid is coarser than the anchor, so the nearest row "
                    "differs at the few-10^-3 level; this is grid spacing, not a "
                    "disagreement")
        record("dark matter", "cold-DM bound margin [orders]",
               "log10(0.1/lambda_fs)",
               math.log10(dg.LYMAN_ALPHA_LFS_MPC
                          / dmg["free_streaming"]["lambda_fs_Mpc"]),
               "manuscript claim (17-18)", 17.5, 5e-2, "rounding")

        ms = [r["m_over_H_derived"] for r in teb["rows"]]
        ls = [r["lambda_fs_Mpc"] for r in teb["rows"]]
        slope = ((math.log(ls[-1]) - math.log(ls[0]))
                 / (math.log(ms[-1]) - math.log(ms[0])))
        record("dark matter", "dln(lambda_fs)/dln(m_psi)",
               "measured along the Table 2 trajectory", slope,
               "analytic -1/2", -0.5, 1.5e-1, "route",
               note="measured %.3f; the manuscript discloses -0.54 to -0.55 and "
                    "derives the sqrt(a_NR) origin of the near-half-integer slope"
                    % slope)

    # ---- B. reheating ceiling quoted in the manuscript ------------------
    record("reheating", "T_reh^inst(50) [GeV]", "manuscript Sec. IV (2.6e15)",
           2.6e15, "background.T_reh_inst(50)", bar.T_reh_inst(N_FID), 2e-2,
           "rounding")

    # ---- A. figures: the manuscript must include one that the run produced,
    # and the run must not leave an unreferenced one behind.  Nothing checked
    # this link before: a renamed or failed figure script surfaced only as a
    # LaTeX error at compile time, never as a run_all failure.  Both sides are
    # scored 1.0/0.0 so the same rel() machinery reports a mismatch either way
    # (included-but-absent, or present-but-never-included).
    figdir = ROOT / "figures"
    included = {os.path.basename(m.group(1)) for m in re.finditer(
        r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", text)}
    on_disk = {p.name for p in figdir.glob("*.pdf")}
    for name in sorted(included | on_disk):
        record("figures", name, "manuscript includes it",
               1.0 if name in included else 0.0, "written into figures/",
               1.0 if name in on_disk else 0.0, 1e-9, "identity",
               note="both sides must be 1: the manuscript must include the figure "
                    "AND the run must have produced it (and produce nothing else)")

    # ---- write ---------------------------------------------------------
    n_fail = sum(1 for r in LEDGER if not r["ok"])
    worst = max(r["rel"] for r in LEDGER)
    out = []
    A = out.append
    A("# Cross-script consistency ledger")
    A("")
    A("Two directions: **(A) paper -> code** (every numeric cell of Table I and")
    A("Table 2 recomputed from the script that produces it) and **(B) code -> code**")
    A("(the same quantity from two independent scripts).")
    A("")
    A("Sources: `paper_prd_merged.tex`, `dm_gap_closure_test.json`,")
    A("`treh_error_band.json`.")
    A("")
    A("## Summary")
    A("")
    A(f"- Checks run: **{len(LEDGER)}** ({n_tab1} Table I rows, {n_tab2} Table 2 rows)")
    A(f"- Failures: **{n_fail}**")
    A(f"- Worst relative deviation: **{worst:.3e}**")
    A("")
    A("Tolerance classes: `identity` = same formula and inputs (1e-9 expected);")
    A("`rounding` = the .tex cell is quoted to 2--3 significant figures;")
    A("`route` = two different implementations of the same physics;")
    A("`precision` = the two sources carry the same constant to different digits.")
    A("")
    kinds = sorted({r["kind"] for r in LEDGER})
    A("| tolerance class | checks | failures | worst deviation |")
    A("|---|---|---|---|")
    for k in kinds:
        sub = [r for r in LEDGER if r["kind"] == k]
        A("| %s | %d | %d | %.2e |"
          % (k, len(sub), sum(1 for r in sub if not r["ok"]),
             max(r["rel"] for r in sub)))
    A("")
    A("## Ledger")
    A("")
    A("| group | quantity | class | source A | source B | rel. dev | tol | verdict |")
    A("|---|---|---|---|---|---|---|---|")
    for r in LEDGER:
        A("| %s | %s | %s | %s | %s | %.2e | %.0e | %s |"
          % (r["group"], r["name"], r["kind"], r["a_src"], r["b_src"], r["rel"],
             r["tol"], "OK" if r["ok"] else "**FAIL**"))
    A("")
    noted = [r for r in LEDGER if r["note"]]
    if noted:
        A("## Recorded deviations (structure, not defects)")
        A("")
        for r in noted:
            A("- **%s** (%s): %.2e -- %s" % (r["name"], r["kind"], r["rel"], r["note"]))
        A("")
    if n_fail:
        A("## Failures needing attention")
        A("")
        for r in LEDGER:
            if not r["ok"]:
                A("- [%s] %s: %s = %.6e vs %s = %.6e (rel %.3e > tol %.0e)"
                  % (r["group"], r["name"], r["a_src"], r["a_val"],
                     r["b_src"], r["b_val"], r["rel"], r["tol"]))
        A("")
    else:
        A("No failure: every tabulated cell is reproduced by its script and every")
        A("independent route agrees within the tolerance recorded above.")
        A("")
    A(f"Runtime: {time.time() - t0:.2f} s")
    OUT_MD.write_text("\n".join(out), encoding="utf-8")
    with io.open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump({"checks": LEDGER, "failures": n_fail,
                   "runtime_s": time.time() - t0}, fh, indent=1)
    print("\n".join(out))
    print("[%s] DONE: %d checks, %d failure(s)"
          % (time.strftime("%H:%M:%S"), len(LEDGER), n_fail))
    return 0


if __name__ == "__main__":
    sys.exit(main())
