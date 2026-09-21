#!/usr/bin/env python3
"""
Numerical claims audit -- every quantitative claim in the manuscript recomputed
from the current scripts
=============================================================================
Type:           PAPER
Paper Sec.:     III, IV, V, VI, App. D
Experiment:     numerical claims audit
What it does:   The reverse direction of audit_tex_numbers.py.  That script asks
                "is a superseded value still present?";  this one asks "is the
                value that IS present correct?"  Each claim below carries the
                number as printed in main.tex, the script that
                computes it, and the tolerance justified by how the number is
                quoted.
Tolerance classes:
    identity  same formula and inputs                     (1e-9 .. 1e-6)
    rounding  the manuscript quotes 2-3 significant figs  (1e-2 .. 5e-2)
    route     two different implementations of one physics
                (declared; the residual is reported, not hidden)
    order     the manuscript prints the number with a leading "~", so the
                meaningful test is a decade match; ``tol`` is then the
                allowed distance in log10, not a relative deviation
Outputs:        verification_report.md, verification_report.json
Dependencies:   reads background_and_reheating.json, dm_gap_closure_test.json,
                residual_quintessence.json, psi_production_bogoliubov.json,
                treh_error_band.json, psi_abundance_oscillating.json and
                psi_mode_oscillating.json, so it must run after those scripts.
Limitation:     the ``paper`` values below are hardcoded mirrors of what
                main.tex prints; this script does NOT parse the
                tex.  If a number is edited in the manuscript without editing
                the matching claim here, this audit stays green while the
                mirror drifts -- audit_tex_numbers.py (superseded-value
                scan) and consistency_checks.py (table-cell recomputation)
                are the guards against that direction, and the pairing must
                be maintained by hand at every manuscript edit.
Note:           the ``order`` block carries the RG-running and closed-form
                claims that used to live in the standalone extended_checks.py
                and quick_claims_check.py.  Those two scripts were retired
                because the manuscript still prints these numbers: had they
                been removed without this block, Secs. III, IV, V, VIII and XII
                would have kept quoting values with no executing source left.
                Since the exact-background computation (psi_mode_oscillating.py)
                moved the abundance-matched coupling to g = 1.5e-7, the
                transition-only values (g = 1.0e-7, m_psi = 7.5e10) are kept as
                the labelled transition-only baseline and the exact-background
                numbers are checked against psi_mode_oscillating.json.
=============================================================================
"""
import io
import json
import math
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("PYTHONUNBUFFERED", "1")
if hasattr(sys.stdout, "reconfigure"):          # Python 3.7+
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT_MD = HERE / "verification_report.md"
OUT_JSON = HERE / "verification_report.json"

sys.path.insert(0, str(HERE))
import derive_from_action as dfa        # noqa: E402
import lock_n_convention as lk          # noqa: E402
import background_and_reheating as bar  # noqa: E402
import order_estimates as oe            # noqa: E402

XI = 11.1
CLAIMS = []


def load_json(name):
    p = HERE / name
    if not p.exists():
        return None
    with io.open(p, encoding="utf-8") as fh:
        return json.load(fh)


def rel(a, b):
    if b == 0.0:
        return abs(a - b)
    return abs(a - b) / abs(b)


def claim(section, quantity, where, paper, computed, tol, kind="rounding", note=""):
    if computed is None:
        CLAIMS.append({"section": section, "quantity": quantity, "where": where,
                       "paper": paper, "computed": None, "rel": float("nan"),
                       "tol": tol, "ok": None, "kind": kind,
                       "note": note or "source artefact missing"})
        return
    CLAIMS.append({"section": section, "quantity": quantity, "where": where,
                   "paper": paper, "computed": computed, "rel": rel(computed, paper),
                   "tol": tol, "ok": rel(computed, paper) <= tol, "kind": kind,
                   "note": note})


def claim_order(section, quantity, where, paper, computed, decades=0.5, note=""):
    """A claim the manuscript prints with a leading ``~``.

    ``decades`` is the allowed distance in log10, not a relative deviation:
    the meaningful statement is "same order of magnitude", so an answer good
    to a factor 2 (0.30 decades) passes and one off by a factor 10 does not.
    These are exactly the claims whose microscopic coefficients (b_s, the
    anomaly alpha_s, g_*) the model does not derive, so holding them to the
    ``rounding`` class would misstate what the paper claims.
    """
    if computed is None:
        CLAIMS.append({"section": section, "quantity": quantity, "where": where,
                       "paper": paper, "computed": None, "rel": float("nan"),
                       "tol": decades, "ok": None, "kind": "order",
                       "note": note or "source artefact missing"})
        return
    d = abs(math.log10(abs(computed)) - math.log10(abs(paper)))
    CLAIMS.append({"section": section, "quantity": quantity, "where": where,
                   "paper": paper, "computed": computed, "rel": d,
                   "tol": decades, "ok": d <= decades, "kind": "order",
                   "note": note})


def main():
    t0 = time.time()
    print("[%s] START: numerical claims audit" % time.strftime("%H:%M:%S"))
    br = load_json("background_and_reheating.json")
    dmg = load_json("dm_gap_closure_test.json")
    rq = load_json("residual_quintessence.json")
    bog = load_json("psi_production_bogoliubov.json")
    teb = load_json("treh_error_band.json")
    pm = load_json("psi_mode_oscillating.json")

    p50 = lk.point(XI, 50.0)
    p55 = lk.point(XI, 55.0)
    pao = load_json("psi_abundance_oscillating.json")
    H0 = rq["constants"]["H0_GeV"] if rq else dfa.H0_GeV_v2()
    # Light-branch dark-matter mass, m_psi = (m_psi/H_inf) H_inf.  The ratio is
    # dm_gap_closure_test.json's anchor and was measured at the fiducial H_inf, so
    # it must be multiplied by the fiducial H_inf (p50).  The H_inf column of
    # treh_error_band.json drifts with T_reh and N and belongs to a different
    # background, so it is the wrong partner here.  Sec. V needs this product for
    # the abundance-matched point and Sec. XII for its two consistency entries,
    # so it is formed once.
    H_INF = p50["H_inf"]
    M_PSI = (dmg["m_star_powerlaw"] * H_INF) if dmg else None
    # Exact-background anchors (psi_mode_oscillating.json): the abundance
    # matching through the oscillating-condensate era, which Sec. V now quotes
    # as the primary result.  The transition-only anchors above remain the
    # labelled baseline the manuscript still prints.
    G_EXACT = pm["matching"]["exact"]["g"] if pm else None
    M_PSI_EXACT = pm["matching"]["exact"]["m_psi_GeV"] if pm else None
    M_OVER_H_EXACT = pm["matching"]["exact"]["m_over_Hinf"] if pm else None

    # ---------------------------------------------------------- Sec. III
    S = "III inflation"
    claim(S, "lambda0 at N=50", "abstract / Table I", 6.70e-8, p50["lambda0"], 5e-3)
    claim(S, "n_s at N=50", "abstract / Table I", 0.9616, p50["ns"], 1e-4)
    claim(S, "r at N=50", "abstract / Table I", 0.00425, p50["r"], 3e-3)
    claim(S, "n_s at N=55", "Table I", 0.9650, p55["ns"], 1e-4)
    claim(S, "r at N=55", "Table I", 0.00355, p55["r"], 3e-3)
    claim(S, "H_inf at N=50", "text / Table I", 1.64e13, p50["H_inf"], 5e-3)
    claim(S, "m_chi at N=50", "text / Table I", 3.25e13, p50["m_chi"], 5e-3)
    claim(S, "V0", "Sec. III", 4.78e63, p50["V0"], 5e-3)
    claim(S, "beta_p", "Eq. (Nint)", 0.8104, p50["beta_p"], 5e-4)
    claim(S, "x_end", "App. D0", 0.76367, bar.X_END, 1e-4, "identity")
    claim(S, "eps_H(x_end)", "App. D0", 0.4988,
          br["slowroll"]["eps_H_end"] if br else None, 5e-4)
    claim(S, "V_end/V0", "App. D0", 0.28520,
          br["slowroll"]["V_end_frac_of_V0"] if br else None, 5e-4)
    claim(S, "K_end/V_end", "Sec. III / App. D0", 0.1994,
          br["slowroll"]["K_over_V_end"] if br else None, 5e-4)
    claim(S, "rho_end/V_end", "Sec. III", 1.1994,
          br["background"]["rho_end_over_V_end"] if br else None, 5e-4)
    rmax = max(lk.point(XI, float(n))["r"] for n in range(45, 56))
    claim(S, "max r over the 2-sigma band", "Sec. III", 0.0052, rmax, 2e-2)
    # Closed form, so this one is exact arithmetic, but the manuscript rounds
    # it to a single figure ("~-0.02") and the value is a slow-roll statement.
    fnl = oe.f_nl_local(p50["ns"])
    claim_order(S, "f_NL^local", "Sec. III (non-Gaussianity)", -0.02, fnl, 0.5,
                note="-5/12 (1-n_s) from the printed n_s; the manuscript quotes "
                     "~-0.02, which is the same decade but rounds up")

    nw = bar.solve_N_max()
    claim(S, "N_max (proper rho_end ceiling)", "Sec. III / App. D0", 55.6,
          nw["N_max_proper"], 2e-3, "route")
    claim(S, "N_max (naive V_end^{1/4} ceiling)", "App. D0", 55.9,
          nw["N_max_naive_Vend14"], 2e-3, "route")

    T50 = lk.T_reh_for_N_derived(p50, 50.0)
    T56 = lk.T_reh_for_N_derived(lk.point(XI, 56.0), 56.0)
    claim(S, "T_reh*(50)", "Sec. III / Table I", 1.1e8, T50, 5e-2)
    claim(S, "T_reh*(56)", "Sec. III", 7.8e15, T56, 5e-2)
    T51 = lk.T_reh_for_N_derived(lk.point(XI, 51.0), 51.0)
    claim(S, "T_reh*(51)", "App. A / Table I", 2.2e9, T51, 5e-2,
          note="guards the App. A verbatim listing, which once printed 1.4e9 with "
               "no source anywhere in the suite; Table I and the "
               "n_convention_results.json T_reh_selfcons entry give 2.166e9")
    claim(S, "T_reh~1e9 GeV maps to N", "abstract / Sec. III", 50.7,
          lk.N_derived_for_T(p50, 1.0e9), 2e-3, "route")
    claim(S, "T_reh~2.1e8 GeV maps to N", "Table 2 row 2", 50.23,
          lk.N_derived_for_T(p50, 2.1e8), 2e-3, "route")

    # ---------------------------------------------------------- Sec. IV
    S = "IV reheating"
    claim(S, "T_reh^inst at N=50", "Sec. IV", 2.6e15, bar.T_reh_inst(50.0), 2e-2)
    claim(S, "T_reh^inst at N=56", "Sec. III (ceiling that excludes N>=56)",
          2.5e15, bar.T_reh_inst(56.0), 2e-2,
          note="the ceiling falls as N^-1/2, so the N=56 value must lie BELOW the "
               "N=50 value quoted in Sec. IV")
    if br:
        rows = br["reheating"]["anomaly"]["rows"]
        phys = [r for r in rows if r["label"].startswith("physical")]
        if phys:
            claim(S, "alpha_s(m_chi)", "Sec. IV", 0.0262, phys[0]["alpha_s"], 5e-3)
            claim(S, "Gamma_chi->gg", "Sec. IV", 0.065, phys[0]["Gamma_GeV"], 5e-2)
            claim(S, "T_reh^(anomaly)", "Sec. IV", 2.1e8, phys[0]["T_reh_GeV"], 5e-2)
    # Mathieu parameter and thermalization rate.  The quark channel is the
    # largest and is still ~21 orders below the broad-resonance threshold
    # q >= 0.1 quoted next to it, so the conclusion does not rest on the
    # factor here; only the decade is claimed.
    claim_order(S, "parametric-resonance q (top)", "Sec. IV",
                2.8e-23, oe.mathieu_q(oe.M_TOP, p50["m_chi"]), 0.5,
                note="q ~ (m_f/m_chi)^2; the electrons and gauge bosons sit "
                     "further below still")
    claim_order(S, "Gamma_therm/H at T=1e9", "Sec. IV / Sec. XII",
                1.4e6, oe.gamma_therm_over_H(1.0e9), 0.5,
                note="printed as Gamma_therm ~ 1.4e6 GeV with the one-loop "
                     "running alpha_s(1e9)=0.0377 (same anchor chain as the "
                     "Sec. IV alpha_s(m_chi)=0.0262; the retired hand estimate "
                     "alpha_s~0.1 inflated this to ~1e7) against "
                     "H ~ sqrt(pi^2 g_*/90) T^2/M_Pl ~ 1.4 GeV")
    claim(S, "alpha_s at T=1e9 GeV", "Sec. XII #17", 0.038,
          oe.alpha_s_running(1.0e9), 5e-2,
          note="one-loop from alpha_s(m_Z)=0.1179, b_3=7; running the result "
               "back up to m_chi reproduces the printed 0.0262 to 1%")

    # ---------------------------------------------------------- Sec. V
    S = "V dark matter"
    if dmg and teb:
        claim(S, "Phi_V", "Sec. V", 7.31e17, teb["Phi_V_GeV"], 5e-3)
        claim(S, "g (transition-only baseline)", "Sec. V (labelled baseline)",
              1.0e-7, dmg["g_star_powerlaw"], 5e-2)
        claim(S, "m_psi (transition-only baseline)",
              "Sec. V (labelled baseline)", 7.5e10, M_PSI, 5e-2)
        claim(S, "m_psi/H_inf (transition-only baseline)",
              "Sec. V (labelled baseline)", 4.6e-3,
              dmg["m_star_powerlaw"], 5e-2)
        claim(S, "lambda_fs anchor", "abstract / Table 2", 9.5e-20,
              dmg["free_streaming"]["lambda_fs_Mpc"], 5e-2)
        claim(S, "lambda_fs bound margin [orders]", "Sec. V", 17.5,
              math.log10(0.1 / dmg["free_streaming"]["lambda_fs_Mpc"]), 5e-2)
        claim(S, "m_psi/H_inf from the mode equation", "Sec. V", 4.6e-3,
              dmg["m_star_powerlaw"], 5e-2, "identity",
              note="the anchor of dm_gap_closure_test.json, i.e. the same number "
                   "the light-branch scan converges to")
        s = teb["summary"]
        claim(S, "band: Delta N", "Sec. V", 1.5, s["dN"], 3e-2)
        claim(S, "band: Delta n_s", "Sec. V", 1.1e-3, s["dns"], 5e-2)
        claim(S, "band: Delta r / r", "Sec. V", 0.055, s["dr_over_r"], 5e-2)
        claim(S, "band: g ratio", "Sec. V", 10.0, s["g_ratio"], 1e-3)
        # Window endpoints of the T_reh band must agree with the printed
        # table rows; the prose window lower edge once drifted to 0.5e-7
        # while the tabulated T=1e9 row reads 4.7e-8.
        g_lo = next((r["g"] for r in teb["rows"]
                     if abs(r["T_reh_GeV"] - 1.0e9) / 1.0e9 < 1e-3), None)
        g_hi = next((r["g"] for r in teb["rows"]
                     if abs(r["T_reh_GeV"] - 1.0e8) / 1.0e8 < 1e-3), None)
        m_lo = next((r["m_psi_GeV"] for r in teb["rows"]
                     if abs(r["T_reh_GeV"] - 1.0e9) / 1.0e9 < 1e-3), None)
        claim(S, "g window lower edge (T=1e9 row)", "Sec. V / Table 2",
              4.7e-8, g_lo, 5e-2)
        claim(S, "g window upper edge (T=1e8 row)", "Sec. V / Table 2",
              1.5e-7, g_hi, 5e-2)
        claim(S, "m_psi window lower edge (T=1e9 row)", "Sec. V / Table 2",
              3.4e10, m_lo, 5e-2)
        m_hi = next((r["m_psi_GeV"] for r in teb["rows"]
                     if abs(r["T_reh_GeV"] - 1.0e8) / 1.0e8 < 1e-3), None)
        if pm and g_lo and g_hi and m_lo and m_hi:
            # Robustness-window paragraph: the 1.48-scaled window is quoted
            # as a range, not a single point (Sec. V revision).
            mult = (pm["g_shift"]["g_avg_over_paper_transition_only"]
                    * pm["g_shift"]["g_exact_over_g_avg"])
            claim(S, "exact-background g window, low edge (T=1e9 row x 1.48)",
                  "Sec. V", 0.7e-7, g_lo * mult, 5e-2)
            claim(S, "exact-background g window, high edge (T=1e8 row x 1.48)",
                  "Sec. V", 2.2e-7, g_hi * mult, 5e-2)
            claim(S, "exact-background m_psi window, low edge (T=1e9 row x 1.48)",
                  "Sec. V", 5.0e10, m_lo * mult, 5e-2)
            claim(S, "exact-background m_psi window, high edge (T=1e8 row x 1.48)",
                  "Sec. V", 1.6e11, m_hi * mult, 5e-2)
        slope = ((math.log(teb["rows"][-1]["lambda_fs_Mpc"])
                  - math.log(teb["rows"][0]["lambda_fs_Mpc"]))
                 / (math.log(teb["rows"][-1]["m_over_H_derived"])
                    - math.log(teb["rows"][0]["m_over_H_derived"])))
        claim(S, "dln(lambda_fs)/dln(m_psi)", "Sec. V", -0.55, slope, 3e-2,
              "route", note="the manuscript quotes -0.54 to -0.55")
    # Exact-background primary results: Sec. V and the abstract now quote the
    # matching through the oscillating-condensate era (psi_mode_oscillating).
    if pm:
        claim(S, "g (exact background, primary)", "abstract / Sec. V",
              1.5e-7, G_EXACT, 5e-2)
        claim(S, "m_psi (exact background, primary)", "abstract / Sec. V",
              1.1e11, M_PSI_EXACT, 5e-2)
        claim(S, "m_psi/H_inf (exact background)", "abstract / Sec. V",
              6.6e-3, M_OVER_H_EXACT, 5e-2)
        claim(S, "condensate enhancement R = n_exact/n_avg",
              "Sec. V / App. D0", 1.24, pm["scan"]["0.0046"]["R"], 3e-2)
        claim(S, "g_exact/g_avg (condensate effect)", "Sec. V",
              0.90, pm["g_shift"]["g_exact_over_g_avg"], 3e-2)
        claim(S, "transition baseline factor (n_p2/n_avg)",
              "Sec. V / App. D0", 2.7,
              1.4e-5 / pm["convergence"]["m=0.0046"]["avg"]["base_dphase0.04"],
              5e-2, "route",
              note="the published transition-only n at the matching point "
                   "(1.4e-5, Sec. V) over the exact-background avg baseline")
        claim(S, "transition-shape systematic in g (sqrt of the n-space factor)",
              "abstract / Sec. V", 1.6,
              pm["g_shift"]["g_avg_over_paper_transition_only"], 5e-2)
        claim(S, "Table-2 multiplier g_exact/g_transition_only",
              "Sec. V / Table 2", 1.48,
              pm["g_shift"]["g_avg_over_paper_transition_only"]
              * pm["g_shift"]["g_exact_over_g_avg"], 3e-2,
              note="g_avg/g_paper x g_exact/g_avg = 1.64 x 0.90; the exact "
                   "value multiplies the transition-only window of Table 2")
        _v1dev = max(abs(r["ratio"] - 1.0) for r in pm["reproduction_p2"])
        claim(S, "V1 anchor: p=2 reproduction at the 1e-9 level", "App. D0",
              1.0, 1.0 if _v1dev <= 1e-9 else 0.0, 1e-9, "identity",
              note="largest deviation of the three published-point "
                   "reproductions: %.2e" % _v1dev)
        claim(S, "massless limit (conformal invariance)", "App. D0",
              0.0, max(v["n_over_H3"] for v in pm["massless_limit"].values()),
              1e-12, "identity")
        _drift = max(c[b]["unitarity_drift_max"]
                     for c in pm["convergence"].values()
                     for b in ("avg", "exact"))
        claim(S, "unitarity drift within the printed 1e-12 bound", "App. D0",
              1.0, 1.0 if _drift <= 1e-12 else 0.0, 1e-9, "identity",
              note="the manuscript writes '<=1e-12'; recomputed max drift "
                   "%.2e" % _drift)
    # Direct detection and the kinematic-closure bound quoted in Sec. IV.
    claim_order(S, "sigma_psiN (graviton exchange)",
                "Sec. V (direct detection)", 1.0e-104, oe.sigma_psiN_cm2(), 0.5,
                note="G_N^2 m_N^2 (hbar c)^2 with G_N = 1/(8 pi M_Pl^2); the "
                     "manuscript prints a single significant figure and the "
                     "conclusion is that this is ~58 orders below XENONnT")
    g_kin = oe.g_kinematic(p50["m_chi"], XI)
    claim(S, "g at kinematic closure (m_psi = m_chi/2)", "Sec. IV / Sec. V",
          1.0, 1.0 if 1.0e-5 <= g_kin <= 1.0e-4 else 0.0, 1e-9, "identity",
          note="the manuscript states m_psi >~ m_chi/2, i.e. g >~ few x 1e-5; "
               "computed %.4e, so 'few x 1e-5' holds" % g_kin)
    if bog:
        claim(S, "Bogoliubov exponent audit present", "Sec. V", 1.0,
              1.0 if bog["exponent_audit"] else 0.0, 1e-9, "identity")

    # ---------------------------------------------------------- Sec. VI
    S = "VI dark energy"
    claim(S, "m_chi/H0 [orders]", "Sec. VI", 55.0,
          math.log10(p50["m_chi"] / H0), 2e-2)
    if rq:
        dw = rq["delta_w"]
        claim(S, "Delta w (quantum) exponent", "Sec. VI", 121.0,
              -math.log10(dw["quant_recomputed_H0^4/Vc"]), 2e-2)
        claim(S, "Delta w (Ricci) exponent", "Sec. VI", 111.0,
              -math.log10(dw["Ricci_recomputed"]), 2e-2)
        claim(S, "overclosure bound rho_cond/V_c", "Sec. VI", 0.39,
              dw["overclosure_bound_OmegaDM_over_OmegaL"], 5e-2)
        claim(S, "V_c [GeV^4]", "Sec. VI", 2.5e-47, rq["constants"]["V_C"], 2e-2,
              note="V_c = Omega_Lambda * rho_c with the manuscript's own "
                   "Omega_Lambda = 0.683; 2.7e-47 would require Omega_Lambda = 0.73")
        rho_c = rq["constants"]["RHO_C"]
        claim(S, "rho_DM^(0) [GeV^4]", "Sec. VI / App. D5c", 9.7e-48,
              0.265 * rho_c, 2e-2,
              note="= Omega_DM * rho_c with the manuscript's own Omega_DM = 0.265; "
                   "1.06e-47 would require Omega_DM = 0.288, and the two places that "
                   "quoted it disagreed with each other")

    # ---------------------------------------------------------- Sec. VIII
    # RG running.  Migrated from extended_checks.py, which was retired.  The
    # manuscript quotes these in Sec. VIII E; before this block NOTHING ran
    # them -- Delta lambda0 in particular had no executing source at all.
    S = "VIII embedding"
    lam0 = p50["lambda0"]
    claim_order(S, "Delta lambda0", "Sec. VIII E", 3.0e-14,
                oe.delta_lambda0(lam0), 0.5,
                note="9 lambda0^2/(8 pi^2) x 60; this quantity had no source "
                     "in the suite before the migration")
    claim_order(S, "Delta xi (quartic drive)", "Sec. VIII E (Eq. dxilam)",
                8.5e-7, oe.delta_xi_quartic(XI, lam0), 0.5)
    g_am = G_EXACT if G_EXACT else 1.5e-7     # abundance-matched coupling
    # Unrolled so that each printed value sits as a literal in args[3]:
    # the reverse-coverage auditor resolves claim mirrors positionally,
    # and a loop-unpacked printed name is invisible to it.
    claim_order(S, "Delta xi (Yukawa), g=1.5e-07", "Sec. VIII E (Eq. dxig)",
                9.1e-14, oe.delta_xi_yukawa(XI, g_am), 0.5)
    claim_order(S, "Delta xi (Yukawa), g=2.3e-05", "Sec. VIII E (Eq. dxig)",
                2.2e-9, oe.delta_xi_yukawa(XI, 2.3e-5), 0.5)
    claim_order(S, "Delta xi (Yukawa), g=1.0e-04", "Sec. VIII E (Eq. dxig)",
                4.2e-8, oe.delta_xi_yukawa(XI, 1.0e-4), 0.5)
    claim_order(S, "quartic / Yukawa dominance", "Sec. VIII E", 9.0e6,
                oe.dn_quartic_over_yukawa(XI, lam0, g_am), 0.5,
                note="at the abundance-matched g = 1.5e-7")
    claim_order(S, "prefactor (xi - 1/6) at xi=1", "Sec. VIII E",
                0.83, 1.0 - 1.0 / 6.0, 0.5)
    claim_order(S, "prefactor (xi - 1/6) at xi=100", "Sec. VIII E",
                99.8, 100.0 - 1.0 / 6.0, 0.5)
    # The printed range [6.5e-8, 7.8e-6] is reproduced only with the
    # anomaly-channel coupling g = 6.9e-5; the retired extended_checks.py
    # scanned 2.3e-5 and 1e-4 instead, so it never covered this sentence.
    dxi_tot = oe.delta_xi_total(XI, lam0, oe.G_ANOM)
    claim_order(S, "Delta xi(xi) at xi=1", "Sec. VIII E", 6.5e-8,
                oe.delta_xi_total(1.0, lam0, oe.G_ANOM), 0.5,
                note="the quoted range assumes g = 6.9e-5")
    claim_order(S, "Delta xi(xi) at xi=100", "Sec. VIII E", 7.8e-6,
                oe.delta_xi_total(100.0, lam0, oe.G_ANOM), 0.5)
    claim(S, "Delta xi(xi=11.1) <= 1e-6 (bound holds)", "Sec. VIII E",
          1.0, 1.0 if dxi_tot <= 1.0e-6 else 0.0, 1e-9, "identity",
          note="the manuscript states an upper bound; computed %.3e" % dxi_tot)

    # ---------------------------------------------------------- Sec. VII
    # Decoupling and the induced Higgs portal.  Migrated into order_estimates
    # in the anchoring round: before that, the threshold suppressions and
    # delta lambda_PhiH had NO executing source anywhere in the suite, and the
    # manuscript printed a stale <= 1e-10 threshold bound and a ~1e-15 portal.
    S = "VII embedding"
    claim_order(S, "threshold suppression (chi)", "Sec. VII", 2.2e-8,
                oe.threshold_suppression(p50["m_chi"]), 0.5,
                note="(m_chi/Lambda_J)^2 with Lambda_J = M_Pl/xi; an earlier "
                     "draft printed a single <= 1e-10 bound, which this "
                     "chi value exceeds by two orders")
    if M_PSI_EXACT:
        claim_order(S, "threshold suppression (psi)", "Sec. VII", 2.4e-13,
                    oe.threshold_suppression(M_PSI_EXACT), 0.5,
                    note="(m_psi/Lambda_J)^2 at the exact-background matched "
                         "m_psi = 1.1e11 GeV")
    claim(S, "delta lambda_PhiH (gravity-induced portal)", "Sec. VII",
          6.0e-17, oe.delta_lambda_PhiH(lam0), 2e-1,
          note="xi lambda0/(16 pi^2) (m_Phi/M_Pl)^2; an earlier draft "
               "printed ~1e-15, two orders too large")

    # ---------------------------------------------------------- Sec. XII
    # Closed-form claims from the Discussion list; same order class as above.
    S = "XII discussion"
    # The two entries of the "additional consistency checks" list that had NO
    # executing source anywhere in the suite -- current or retired.  Both were
    # found to carry a stale wimpzilla-scale mass (m ~ 1e13 GeV): 1e-67 cm^2/g
    # and 1e52 GeV^4 are what these closed forms give at m = 1e13, not at the
    # m_psi = 1.1e11 GeV that Sec. V actually matches (exact background).
    # The manuscript now quotes the values checked here.
    if M_PSI_EXACT:
        sig_ov_m = oe.sigma_psipsi_over_m(M_PSI_EXACT)
        claim_order(S, "sigma_psipsi/m_psi (graviton exchange)", "Sec. XII (i)",
                    1.0e-69, sig_ov_m, 0.5,
                    note="G_N^2 m_psi (hbar c)^2 / G_PER_GEV with the same G_N and "
                         "the same (hbar c)^2 convention as sigma_psiN; this entry "
                         "had no source before this migration")
        claim(S, "bullet-cluster margin [orders]", "Sec. XII (i)", 69.0,
              -math.log10(sig_ov_m), 2e-2,
              note="the manuscript writes '~69 orders below the bullet-cluster "
                   "bound'; the printed cross section must move with it")
        claim_order(S, "Tremaine-Gunn Q", "Sec. XII (ii)",
                    1.4e44, oe.tremaine_gunn_Q(M_PSI_EXACT, H_INF), 0.5,
                    note="Q = rho/sigma_v^3 = m_psi^4 with rho = m_psi H_inf^3 and "
                         "sigma_v = H_inf/m_psi; also had no source before")
    claim_order(S, "spectral running alpha_s", "Sec. XII (xiv)", -8.0e-4,
                oe.alpha_s_attractor(), 0.5,
                note="attractor closed form -2/N^2 at the locked N=50")
    claim_order(S, "Delta w quartic chain (Discussion #14)", "Sec. XII (#14)",
                2.0e-231, oe.dw_quartic_ricci(lam0, XI, p50["m_chi"]), 0.5,
                note="(lambda0/xi^2)(H0/m_chi)^4 with the Sec. VI Ricci "
                     "displacement; the draft chain 1e-6 * (1e-46)^2 = 1e-98 "
                     "was internally inconsistent")
    claim(S, "N at which r = 0.01", "Sec. III (LiteBIRD)", 32.0,
          oe.n_at_r(0.01, XI), 5e-2,
          note="exact slow-roll inversion; the manuscript writes N ~ 32, "
               "outside the Planck band")

    # ---------------------------------------------------------- reverse sweep
    # Registrations from the bidirectional provenance audit (see
    # review_workspace/provenance_audit_report.md): these magnitudes were
    # printed in the manuscript but sat outside every audit's field of view.
    # Each value was re-derived by hand before registration; the tolerance
    # states how the manuscript quotes it (rounding vs decade).  Two genuine
    # numerical errors found by the same sweep (App. D m_chi/H floor and the
    # D5c suppression factor) were fixed in the manuscript first; the claims
    # below check the corrected print.
    PhiV = teb["Phi_V_GeV"] if teb else None
    if rq:
        V_END = rq["constants"]["V_end"]
        RHO_END = rq["constants"]["rho_end"]
        H0C = rq["constants"]["H0_GeV"]
    M_CHI = p50["m_chi"]
    BETA = p50["beta_p"]

    S = "II model scales"
    if teb:
        claim(S, "m_Phi = sqrt(2 lam0) Phi_V", "Sec. II / V / VIII", 2.7e14,
              math.sqrt(2.0 * lam0) * PhiV, 5e-2)
    if rq:
        claim(S, "U^(1/4) = V0^(1/4)", "Sec. III / App. A", 8.3e15,
              p50["V0"] ** 0.25, 5e-2)
        claim_order(S, "V_end^(1/4)", "Sec. IV", 6.1e15, V_END ** 0.25, 0.5)
        claim_order(S, "Gamma(chi -> tt) / m_chi", "Sec. II E", 5.0e-35,
                    (173.0 / (oe.M_P * math.sqrt(6.0 + 1.0 / XI))) ** 2
                    / (16.0 * math.pi), 0.5)
        claim_order(S, "Gamma(chi -> tt) absolute", "Sec. II E", 1.5e-21,
                    (173.0 / (oe.M_P * math.sqrt(6.0 + 1.0 / XI))) ** 2
                    * M_CHI / (16.0 * math.pi), 0.5)
        claim(S, "V_end (absolute)", "Sec. III", 1.36e63, V_END, 5e-3)
        claim(S, "rho_end (absolute)", "Sec. IV", 1.64e63, RHO_END, 5e-3)
        claim(S, "K_end (absolute)", "Sec. III", 2.7e62,
              br["slowroll"]["K_over_V_end"] * V_END, 5e-2)
        claim(S, "Lambda_J = M_Pl / xi", "Sec. II E / Sec. VII", 2.2e17,
              oe.lambda_J(XI), 5e-2)
        claim(S, "m_chi / H_0 today", "Sec. VI / App. D", 2.3e55,
              M_CHI / H0C, 5e-2)
        claim_order(S, "m_chi / H at the model T_reh", "App. D (quasi-static)",
                    5.0e14,
                    M_CHI / math.sqrt((math.pi ** 2 / 30.0) * 106.75
                                      * (2.1e8) ** 4 / (3.0 * oe.M_P ** 2)),
                    0.5, note="replaces the erroneous 'never below 1e50'")
        claim_order(S, "wall tension sigma_wall", "Sec. II F", 1.0e50,
                    (4.0 / 3.0) * math.sqrt(lam0 / 2.0) * PhiV ** 3
                    if PhiV else None, 0.5)

    S = "IV reheating"
    if rq:
        claim(S, "initial gravitational pulse rho_rad(a_end)", "Sec. IV",
              4.0e50, rq["constants"]["rho_pulse"], 1e-1,
              note="the seed bath residual_quintessence.py actually integrates; "
                   "the manuscript quotes the order")

    S = "V dark matter"
    if teb and rq:
        claim(S, "heavy-branch closure mass m_psi = 3.1 H_inf", "Sec. V",
              5.1e13, 3.1 * H_INF, 5e-2)
        claim(S, "heavy-branch closure coupling g = m_psi / Phi_V", "Sec. V",
              6.9e-5, 3.1 * H_INF / PhiV, 5e-2)
        claim(S, "de Sitter closure g = H_inf / (2 pi Phi_V)", "Sec. V",
              3.6e-6, H_INF / (2.0 * math.pi * PhiV), 5e-2)
        g_row = next((r["g"] for r in teb["rows"]
                      if abs(r["T_reh_GeV"] - 1.0e9) / 1.0e9 < 1e-3), None)
        m_row = next((r["m_psi_GeV"] for r in teb["rows"]
                      if abs(r["T_reh_GeV"] - 1.0e9) / 1.0e9 < 1e-3), None)
        if g_row and dmg:
            # light branch: n_psi prop m_psi prop g, so the n of the T=1e9
            # Table-2 row follows from the T_reh_model anchor of dm_gap
            rho_psi = (m_row * dmg["n_star_powerlaw"]
                       * (g_row / dmg["g_star_powerlaw"]) * H_INF ** 3)
            claim(S, "rho_psi at production (T_reh = 1e9 row)", "Sec. V",
                  9.5e44, rho_psi, 5e-2)
            claim(S, "rho_psi / rho_end", "Sec. V", 5.8e-19,
                  rho_psi / RHO_END, 5e-2)
        if pm:
            claim(S, "n_psi / H_inf^3 (exact background)", "Sec. V", 9.5e-6,
                  pm["matching"]["exact"]["n_over_H3"], 5e-2)
        if dmg:
            claim(S, "n_psi / H_inf^3 (transition-only baseline)", "Sec. V",
                  1.4e-5, dmg["n_star_powerlaw"], 5e-2)
            claim_order(S, "adiabaticity numerator m_Phi Phi_osc / Phi_V^2",
                        "Sec. V", 4.4e-4, 0.8 * math.sqrt(2.0 * lam0), 0.5,
                        note="Phi_osc = 0.8 Phi_V (Eq. adiab); the printed "
                             "4.4e-4 corresponds to a slightly larger amplitude")
            if G_EXACT:
                claim_order(S, "adiabaticity parameter at g_exact", "Sec. V",
                            3.0e3, 0.8 * math.sqrt(2.0 * lam0) / G_EXACT, 0.5)
                claim(S, "adiabaticity ratio light / heavy branch", "Sec. V",
                      470.0, (3.1 * H_INF / PhiV) / G_EXACT, 5e-2,
                      note="the numerator cancels: this is g_heavy / g_exact")
    if pao and dmg:
        rows_n = {r["m_over_Hinf"]: r["n_num"] for r in pao["compare"]}
        m_match = dmg["m_star_powerlaw"]
        if 3.0 in rows_n:
            omega3 = 0.265 * (3.0 * rows_n[3.0]) / (m_match * dmg["n_star_powerlaw"])
            claim(S, "heavy-branch Omega_psi at m/H = 3 (low end)", "Sec. V",
                  1.6e3, omega3, 1e-1)
            claim_order(S, "heavy-branch overproduction ratio", "Sec. V",
                        6.0e3, omega3 / 0.265, 0.5)
        treh_hi = pao["T_reh_required"].get("3")
        if treh_hi:
            claim_order(S, "heavy-branch required T_reh", "Sec. V", 3.0e4,
                        treh_hi, 0.5)
    if teb:
        g10 = next((r["g"] for r in teb["rows"]
                    if abs(r["T_reh_GeV"] - 1.0e10) / 1.0e10 < 1e-3), None)
        m10 = next((r["m_psi_GeV"] for r in teb["rows"]
                    if abs(r["T_reh_GeV"] - 1.0e10) / 1.0e10 < 1e-3), None)
        if g10:
            claim(S, "g at the T_reh = 1e10 row", "Sec. V", 1.5e-8, g10, 5e-2)
        if m10:
            claim(S, "m_psi at the T_reh = 1e10 row", "Sec. V", 1.1e10, m10, 5e-2)
    if dmg:
        claim_order(S, "lambda_fs order (abstract rounding)", "abstract / Sec. V",
                    1.0e-19, dmg["free_streaming"]["lambda_fs_Mpc"], 0.5)
    if M_PSI_EXACT:
        claim_order(S, "m_chi order (decade rounding)", "abstract / Secs. II-VIII",
                    1.0e13, M_CHI if M_CHI else None, 0.6,
                    note="decade rounding of the locked 3.25e13 GeV")

    S = "VI dark energy"
    if rq:
        V_C = rq["constants"]["V_C"]
        RHO_DM = 0.265 * rq["constants"]["RHO_C"]
        claim_order(S, "V_0 / V_c hierarchy", "abstract / VI / XII", 1.0e110,
                    p50["V0"] / V_C, 0.5)
        claim_order(S, "radiative correction delta V_c (psi)", "Sec. VI", 1.0e42,
                    M_PSI_EXACT ** 4 / (16.0 * math.pi ** 2)
                    if M_PSI_EXACT else None, 0.5,
                    note="m^4/(16 pi^2); the sweep moved the printed decade "
                         "from 1e41, which the exact m_psi exceeds by 0.94 dec")
        claim_order(S, "radiative correction delta V_c (chi)", "Sec. VI", 1.0e52,
                    M_CHI ** 4 / (16.0 * math.pi ** 2), 0.5)
        claim_order(S, "tuning vs V_c (psi end)", "Sec. VI", 1.0e89,
                    M_PSI_EXACT ** 4 / (16.0 * math.pi ** 2) / V_C
                    if M_PSI_EXACT else None, 0.5)
        claim_order(S, "tuning vs V_c (chi end)", "Sec. VI", 1.0e98,
                    M_CHI ** 4 / (16.0 * math.pi ** 2) / V_C, 0.5)
        lo = RHO_DM * (1.0e-2 / 2.36e-13) ** 3
        hi = RHO_DM * (1.0e15 / 2.36e-13) ** 3
        claim_order(S, "rho_cond(a_reh) bound, low edge", "App. D5c", 5.2e-16,
                    lo, 0.5)
        claim_order(S, "rho_cond(a_reh) bound, high edge", "App. D5c", 5.2e35,
                    hi, 0.5)
        claim(S, "rho_cond(a_end) = m_chi^2 M_Pl^2 / 2", "App. D5c", 3.0e63,
              0.5 * M_CHI ** 2 * oe.M_P ** 2, 1e-1)
        claim_order(S, "required suppression factor, low", "App. D5c", 4.0e27,
                    0.5 * M_CHI ** 2 * oe.M_P ** 2 / hi, 0.5,
                    note="the sweep replaced the printed 2e26, which disagreed "
                         "with this sentence's own rho_cond(a_end) by 29x")
        claim_order(S, "required suppression factor, high", "App. D5c", 4.0e78,
                    0.5 * M_CHI ** 2 * oe.M_P ** 2 / lo, 0.5)
        sc = rq.get("stable_condensate", {}).get("no_pulse", {})
        if sc:
            claim_order(S, "stable-condensate overclosure vs rho_DM", "App. D5c",
                        1.0e18, sc["rho_cond_a0"] / RHO_DM, 0.5)
            claim_order(S, "stable-condensate rho_cond(a0)/V_c", "App. D5c",
                        7.0e17, sc["rho_cond_a0_over_Vc"], 0.5)
        claim_order(S, "V_c order (decade rounding)", "Sec. VI / VIII", 1.0e-47,
                    V_C, 0.5)

    S = "VIII embedding"
    claim(S, "Coleman-Weinberg VEV shift", "Sec. VIII E", 4.0e-10,
          lam0 / (16.0 * math.pi ** 2), 1e-1)
    claim(S, "quartic VEV-shift coefficient lam0 / xi^2", "Sec. XII (field space)",
          5.4e-10, lam0 / XI ** 2, 5e-2)
    claim_order(S, "2-loop gamma_Phi", "Sec. VIII E", 1.0e-19,
                (lam0 / (16.0 * math.pi ** 2)) ** 2, 0.5)
    claim(S, "c_3 cubic term of V_E", "App. A", 1.8e8,
          BETA * M_CHI ** 2 / (2.0 * oe.M_P), 5e-2)

    S = "App D residual budget"
    if rq and teb and PhiV:
        dA = XI * 9.2 * H0C ** 2 / (lam0 * PhiV ** 2)
        claim(S, "Delta A / A today (Ricci)", "Sec. VI / App. D5b", 6.3e-111,
              dA, 1e-1)
        claim_order(S, "delta Phi / Phi_V (Ricci, today)", "App. D5", 1.0e-111,
                    (H0C / M_CHI) ** 2, 0.5)
        claim_order(S, "Delta w_Ricci mantissa", "App. D5 / Sec. VI", 2.0e-111,
                    (H0C / M_CHI) ** 2, 0.5)
        claim_order(S, "Delta alpha/alpha scenario (ii)", "Sec. VI", 1.0e-110,
                    dA, 0.5)
        claim_order(S, "Delta w quantum fluctuations", "App. D5b", 2.0e-121,
                    H0C ** 4 / V_C, 0.5,
                    note="the manuscript's own chain: rho_quant ~ 3H0^4/(16 pi^2)"
                         " << V_c, so Delta w <= H0^4/V_c (the prefactor is "
                         "dropped inside the inequality)")
        claim_order(S, "domain-wall suppression e^{-3N}", "Sec. II / XII",
                    1.0e-66, math.exp(-3.0 * 50.7), 0.5,
                    note="the sweep moved the printed decade from 1e-65")
        claim_order(S, "orders below DESI 3-sigma", "App. D5b", 1.0e108,
                    10.0 ** math.log10(1.0e-3 / ((H0C / M_CHI) ** 2)), 0.5)

    # ---------------------------------------------------------- summary
    graded = [c for c in CLAIMS if c["ok"] is not None]
    n_fail = sum(1 for c in graded if not c["ok"])
    n_skip = len(CLAIMS) - len(graded)

    out = []
    A = out.append
    A("# Numerical claims audit")
    A("")
    A("Every quantitative claim in `main.tex` recomputed from the current")
    A("scripts.  This is the reverse direction of `audit_tex_numbers.py`: that script")
    A("looks for *superseded* values, this one checks that the values that *are*")
    A("printed are correct.")
    A("")
    A("## Summary")
    A("")
    A(f"- Claims checked: **{len(graded)}** ({n_skip} skipped: source artefact absent)")
    A(f"- Failures: **{n_fail}**")
    A("")
    A("| section | claims | failures |")
    A("|---|---|---|")
    for sec in dict.fromkeys(c["section"] for c in CLAIMS):
        sub = [c for c in CLAIMS if c["section"] == sec and c["ok"] is not None]
        if sub:
            A("| %s | %d | %d |" % (sec, len(sub),
                                    sum(1 for c in sub if not c["ok"])))
    A("")
    A("## Claims")
    A("")
    A("| section | quantity | where | printed | recomputed | dev | tol | verdict |")
    A("|---|---|---|---|---|---|---|---|")
    for c in CLAIMS:
        if c["computed"] is None:
            A("| %s | %s | %s | %.6g | - | - | %.0e | SKIP |"
              % (c["section"], c["quantity"], c["where"], c["paper"], c["tol"]))
        else:
            A("| %s | %s | %s | %.6g | %.6g | %.2e%s | %.0e | %s |"
              % (c["section"], c["quantity"], c["where"], c["paper"], c["computed"],
                 c["rel"], " dec" if c["kind"] == "order" else "", c["tol"],
                 "OK" if c["ok"] else "**FAIL**"))
    A("")
    if n_fail:
        A("## Failing claims")
        A("")
        for c in graded:
            if not c["ok"]:
                A("- **%s** (%s): printed %.6g, recomputed %.6g, rel %.3e > tol %.0e"
                  % (c["quantity"], c["where"], c["paper"], c["computed"],
                     c["rel"], c["tol"]))
                if c["note"]:
                    A("  - %s" % c["note"])
        A("")
    else:
        A("No failing claim.")
        A("")
    noted = [c for c in CLAIMS if c["note"] and c["ok"]]
    if noted:
        A("## Notes on claims that pass")
        A("")
        for c in noted:
            A("- **%s**: %s" % (c["quantity"], c["note"]))
        A("")
    A(f"Runtime: {time.time() - t0:.2f} s")

    OUT_MD.write_text("\n".join(out), encoding="utf-8")
    with io.open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump({"claims": CLAIMS, "failures": n_fail,
                   "runtime_s": time.time() - t0}, fh, indent=1)
    print("\n".join(out))
    print("[%s] DONE: %d claims, %d failure(s)"
          % (time.strftime("%H:%M:%S"), len(graded), n_fail))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
