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
                number as printed in paper_prd_merged.tex, the script that
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
                residual_quintessence.json, psi_production_bogoliubov.json and
                treh_error_band.json, so it must run after those scripts.
Note:           the ``order`` block carries the RG-running and closed-form
                claims that used to live in the standalone extended_checks.py
                and quick_claims_check.py.  Those two scripts were retired
                because the manuscript still prints these numbers: had they
                been removed without this block, Secs. III, IV, V, VIII and XII
                would have kept quoting values with no executing source left.
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

    p50 = lk.point(XI, 50.0)
    p55 = lk.point(XI, 55.0)
    H0 = rq["constants"]["H0_GeV"] if rq else dfa.H0_GeV_v2()

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
                1.0e7, oe.gamma_therm_over_H(1.0e9), 0.5,
                note="Gamma_therm ~ alpha_s^2 T ~ 1e7 GeV against "
                     "H ~ sqrt(pi^2 g_*/90) T^2/M_Pl ~ 1.4 GeV")

    # ---------------------------------------------------------- Sec. V
    S = "V dark matter"
    if dmg and teb:
        claim(S, "Phi_V", "Sec. V", 7.31e17, teb["Phi_V_GeV"], 5e-3)
        claim(S, "g (light branch)", "abstract / Sec. V", 1.0e-7,
              dmg["g_star_powerlaw"], 5e-2)
        claim(S, "m_psi", "abstract / Sec. V", 7.5e10,
              dmg["m_star_powerlaw"] * teb["rows"][1]["H_inf_GeV"], 5e-2)
        claim(S, "m_psi/H_inf", "abstract / Sec. V", 4.6e-3,
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
        slope = ((math.log(teb["rows"][-1]["lambda_fs_Mpc"])
                  - math.log(teb["rows"][0]["lambda_fs_Mpc"]))
                 / (math.log(teb["rows"][-1]["m_over_H_derived"])
                    - math.log(teb["rows"][0]["m_over_H_derived"])))
        claim(S, "dln(lambda_fs)/dln(m_psi)", "Sec. V", -0.55, slope, 3e-2,
              "route", note="the manuscript quotes -0.54 to -0.55")
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
    for g_y, printed in ((1.0e-7, 4.2e-14), (2.3e-5, 2.2e-9), (1.0e-4, 4.2e-8)):
        claim_order(S, "Delta xi (Yukawa), g=%.1e" % g_y,
                    "Sec. VIII E (Eq. dxig)", printed,
                    oe.delta_xi_yukawa(XI, g_y), 0.5)
    claim_order(S, "quartic / Yukawa dominance", "Sec. VIII E", 2.0e7,
                oe.dn_quartic_over_yukawa(XI, lam0, 1.0e-7), 0.5,
                note="at the abundance-matched g = 1.0e-7")
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

    # ---------------------------------------------------------- Sec. XII
    # Closed-form claims from the Discussion list; same order class as above.
    S = "XII discussion"
    claim_order(S, "spectral running alpha_s", "Sec. XII (xiv)", -8.0e-4,
                oe.alpha_s_attractor(), 0.5,
                note="attractor closed form -2/N^2 at the locked N=50")
    claim(S, "N at which r = 0.01", "Sec. III (LiteBIRD)", 32.0,
          oe.n_at_r(0.01, XI), 5e-2,
          note="exact slow-roll inversion; the manuscript writes N ~ 32, "
               "outside the Planck band")

    # ---------------------------------------------------------- summary
    graded = [c for c in CLAIMS if c["ok"] is not None]
    n_fail = sum(1 for c in graded if not c["ok"])
    n_skip = len(CLAIMS) - len(graded)

    out = []
    A = out.append
    A("# Numerical claims audit")
    A("")
    A("Every quantitative claim in `paper_prd_merged.tex` recomputed from the current")
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
    return 0


if __name__ == "__main__":
    sys.exit(main())
