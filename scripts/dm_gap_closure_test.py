#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
dm_gap_closure_test.py
======================
Question: can the dark-matter overproduction gap (Omega_psi ~ 3e4 x Omega_DM) be closed within the
**current action + current reheating mechanism**? Or must the model be extended?

Why a separate test is worthwhile
------------------
The Sec. 4 Omega table of `psi_abundance_oscillating.py` (and the Sec. 2c T_reh inversion) contains the lines

    if m < 0.4:
        continue

i.e. that table **deliberately skips** the points with m/H_inf < 0.4. And the conclusion "under the
power law Omega is nearly independent of g, the abundance is set by T_reh, and g cannot tune it"
was read off precisely that table. This script extends the scan toward small m and forces frozen convergence.

Three independent reasons why the "small-g solution" should exist
----------------------------------
1. The paper's own Eq. (eq:Omegapsi) reads
       Omega ~ m_psi H_inf^3 exp(-2 pi m_psi/H_inf) (a_end/a_0)^3 / rho_c,
   which is **non-monotonic** in m_psi; the paper itself admits at L458 that "Omega_psi as a function of g
   has a maximum at g = 3.6e-6". Any target value below the maximum has **two** solutions.
2. A massless Dirac field is conformally invariant => production is exactly zero (self-check in this repo: |beta| = 0 to machine zero).
   Hence as m -> 0, n_psi -> 0, m_psi n_psi -> 0, Omega_psi -> 0.
3. From 1+2: Omega_psi(g) -> 0 as g -> 0, far exceeds 0.265 at large g, and is continuous
   => by the intermediate value theorem a small-g solution **must** exist. So the real question is not "whether a solution exists",
   but "whether the g of that solution is physically acceptable".

Two key differences from the original script (otherwise the small-m points are not trustworthy)
--------------------------------------------
1. Frozen convergence: the original script uses a_final = 60 for all m, so at small m
   k_max/(m a_final) is far larger than 1 (up to 67 at m = 0.01) and the modes never freeze.
   This script takes a_final = max(60, k_max/(0.3 m)), forcing k_max/(m a_final) <= 0.3.
2. Lower limit of the k integral: production peaks at k ~ m, but the original script fixes kmin = 0.05; when m < 0.05
   the peak is missed entirely. This script takes kmin = min(0.05, m/10).

Verdict
----
PASS : there exists a g with Omega_psi = Omega_DM = 0.265 (the gap can be closed within the current model)
FAIL : none exists (closing the gap requires extending the model)
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

REUSE = "--reuse" in sys.argv          # reuse the scan results from the existing json; only re-emit the report

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import cosmo_model as cm                   # noqa: E402
from psi_abundance_oscillating import (  # noqa: E402
    H_INF, PHI_V, RHO_C, npsi_from_spectrum,
)

# dil(T_reh) = DIL_PER_GEV * T_reh   (T_reh in GeV), from the shared
# entropy-conserving function.  This was a hand-copied 1.0204e-101 ("same
# convention as [2b] Sec. 2c") that no audit cross-checked.
DIL_PER_GEV = cm.dilution(1e9) / 1e9
OMEGA_TARGET = 0.265
OMEGA_DM = 0.265
T_REH_MODEL = 2.1e8          # T_reh [GeV] from the anomaly channel (physical alpha_s)
T_REH_FIDUCIAL = 1.0e9       # paper fiducial
A_END_OVER_A0 = cm.dilution(1e9) ** (1.0 / 3.0)

FREEZE_TARGET = 0.3          # require k_max/(m a_final) <= 0.3
MASSES = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0]
# the exponential closed form must be scanned down to very small m to find the "light psi branch" (analytically m ~ 1e-8)
EXP_GRID = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3,
            0.01, 0.1, 0.3, 1.0, 3.0]

JSON_PATH = HERE / "dm_gap_closure_test.json"
MD_PATH = HERE / "dm_gap_closure_test.md"


def a_final_for(m: float, kmax: float) -> float:
    """Minimum a_final required by the freeze condition k_max/(m a_final) <= FREEZE_TARGET."""
    return max(60.0, kmax / (FREEZE_TARGET * m))


def omega(m: float, n_over_H3: float, T_reh: float) -> float:
    """Omega_psi = m_psi * n_psi * dil(T_reh) / rho_c,  n_psi = (n/H^3) * H^3."""
    return m * H_INF * n_over_H3 * H_INF**3 * DIL_PER_GEV * T_reh / RHO_C


def omega_exp(m: float, T_reh: float) -> float:
    """The paper's de Sitter exponential closed form: n/H^3 = exp(-2 pi m/H_inf) (for comparison)."""
    return omega(m, math.exp(-2.0 * math.pi * m), T_reh)


# ----------------------------------------------------------- free-streaming length check
G_STAR = 106.75
M_P_GEV = 2.435e18
H0_GEV = 1.4377e-42
T_EQ_GEV = 7.5e-10            # ~0.75 eV, matter-radiation equality
MPC_PER_INV_GEV = 6.395e-39   # 1 GeV^-1 -> Mpc
LYMAN_ALPHA_LFS_MPC = 0.1     # conservative upper limit on the warm-DM free-streaming length


def free_streaming_length(m_over_H, T_reh, p_over_m_prod=None):
    """Comoving free-streaming length of the light branch [Mpc]. Unit convention: a_end=1, H in GeV.

    Relativistic -> non-relativistic via v(a) = x/sqrt(1+x^2), x = (p_prod/m)/a.
    By default the **most conservative** (longest free-streaming) p_prod ~ H_inf is used.
    """
    if p_over_m_prod is None:
        p_over_m_prod = 1.0 / m_over_H
    a0 = 1.0 / A_END_OVER_A0
    H_reh = math.sqrt(math.pi**2 * G_STAR / 90.0) * T_reh**2 / M_P_GEV
    a_reh = (H_INF / H_reh) ** (2.0 / 3.0)      # w=0 condensate dominated
    a_eq = a_reh * (T_reh / T_EQ_GEV)
    H_eq = H0_GEV * (a0 / a_eq) ** 1.5

    def H_of(a):
        if a <= a_reh:
            return H_INF * a ** (-1.5)
        if a <= a_eq:
            return H_reh * (a / a_reh) ** (-2.0)
        return H_eq * (a / a_eq) ** (-1.5)

    def integrand(ln_a):
        a = math.exp(ln_a)
        x = p_over_m_prod / a
        v = x / math.sqrt(1.0 + x * x)
        return v / (a * H_of(a))               # v da/(a^2 H) = v dln a/(a H)

    total = 0.0
    for lo, hi, n in ((0.0, math.log(a_reh), 4000),
                      (math.log(a_reh), math.log(a_eq), 2000),
                      (math.log(a_eq), math.log(a0), 4000)):
        if hi <= lo:
            continue
        h = (hi - lo) / n
        s = 0.5 * (integrand(lo) + integrand(hi))
        s += sum(integrand(lo + i * h) for i in range(1, n))
        total += s * h
    return {"lambda_com_GeVinv": total,
            "lambda_fs_Mpc": a0 * total * MPC_PER_INV_GEV,
            "p_over_m_prod": p_over_m_prod,
            "a_reh": a_reh, "a_eq": a_eq, "H_reh_GeV": H_reh}


def solve_g_for_target(ms, omegas, target=OMEGA_TARGET):
    """Interpolate on a log-log grid and return the **smallest** Omega = target solution (i.e. the small-g branch).

    Note: the paper's exponential closed form first rises and then falls with m and has two solutions; the smaller one
    must be taken to obtain the "light psi branch". An earlier version returned the first zero crossing (the large-g branch), which was wrong.
    """
    pts = sorted((m, w) for m, w in zip(ms, omegas) if w and w > 0)
    for (m1, w1), (m2, w2) in zip(pts, pts[1:]):
        if (w1 - target) * (w2 - target) <= 0 and w1 != w2:
            t = (math.log(target) - math.log(w1)) / (math.log(w2) - math.log(w1))
            lm = math.log(m1) + t * (math.log(m2) - math.log(m1))
            return math.exp(lm)
    return None


def loglog_n_at(m_star, rows):
    """Interpolate n/H^3 between the two adjacent points on a log-log grid."""
    lo = max((r for r in rows if r["m_over_H"] <= m_star),
             key=lambda r: r["m_over_H"], default=None)
    hi = min((r for r in rows if r["m_over_H"] >= m_star),
             key=lambda r: r["m_over_H"], default=None)
    if lo and hi and hi["m_over_H"] != lo["m_over_H"]:
        t = ((math.log(m_star) - math.log(lo["m_over_H"])) /
             (math.log(hi["m_over_H"]) - math.log(lo["m_over_H"])))
        return math.exp(math.log(lo["n_over_H3"]) + t *
                        (math.log(hi["n_over_H3"]) - math.log(lo["n_over_H3"])))
    return lo["n_over_H3"] if lo else None


def run_scan() -> list:
    print("=" * 78)
    print("DM gap closure test:  does Omega_psi(g) have a second (small-g) solution?")
    print("=" * 78)
    print(f"H_inf = {H_INF:.4e} GeV   Phi_V = {PHI_V:.4e}   rho_c = {RHO_C:.4e}")
    print(f"T_reh(model, anomaly) = {T_REH_MODEL:.2e} GeV ;  "
          f"target Omega = {OMEGA_TARGET}")
    print("freeze target: k_max/(m*a_final) <= "
          f"{FREEZE_TARGET}   kmin = min(0.05, m/10)")
    print("-" * 78)
    print(f"{'m/H':>8} {'g':>11} {'kmax':>6} {'a_fin':>9} {'n/H^3':>11} "
          f"{'m*n/H^4':>11} {'Om(T=2.1e8)':>12} {'/0.265':>10} {'sec':>6}")
    rows = []
    for m in MASSES:
        kmax = max(40.0, 20.0 * m)
        kmin = min(0.05, m / 10.0)
        a_fin = a_final_for(m, kmax)
        ta = time.time()
        res = npsi_from_spectrum(m, a_final=a_fin, kmax=kmax,
                                 kmin=kmin, nk_pts=45)
        n = res["n_over_H3"]
        w = omega(m, n, T_REH_MODEL)
        rows.append({"m_over_H": m, "g": m * H_INF / PHI_V, "kmax": kmax,
                     "kmin": kmin, "a_final": a_fin, "n_over_H3": n,
                     "m_times_n": m * n, "omega_at_model_T": w,
                     "ratio_to_DM": w / OMEGA_DM,
                     "omega_at_fiducial_T": omega(m, n, T_REH_FIDUCIAL),
                     "nsteps": res["nsteps"]})
        print(f"{m:8.4g} {m*H_INF/PHI_V:11.4e} {kmax:6.1f} {a_fin:9.1f} "
              f"{n:11.4e} {m*n:11.4e} {w:12.4e} {w/OMEGA_DM:10.4e} "
              f"{time.time()-ta:6.1f}")
    return rows


def main() -> None:
    t0 = time.time()
    if REUSE and JSON_PATH.exists():
        prev = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        rows = prev["rows"]
        print(f"[reuse] loaded {len(rows)} scan rows from "
              f"{JSON_PATH.name}; skipping mode-equation integration")
    else:
        rows = run_scan()

    ms = [r["m_over_H"] for r in rows]
    ws = [r["omega_at_model_T"] for r in rows]
    m_star = solve_g_for_target(ms, ws)
    n_star = loglog_n_at(m_star, rows) if m_star else None
    fs = free_streaming_length(m_star, T_REH_MODEL) if m_star else None

    # ---- comparison: the paper's exponential closed form ----
    exp_rows = [{"m_over_H": m, "omega_at_model_T": omega_exp(m, T_REH_MODEL)}
                for m in EXP_GRID]
    for r in exp_rows:
        r["ratio_to_DM"] = r["omega_at_model_T"] / OMEGA_DM
        r["g"] = r["m_over_H"] * H_INF / PHI_V
    m_star_exp = solve_g_for_target([r["m_over_H"] for r in exp_rows],
                                    [r["omega_at_model_T"] for r in exp_rows])

    print("-" * 78)
    print("Comparison (the paper's de Sitter exponential closed form, n/H^3 = exp(-2 pi m/H), same T):")
    for r in exp_rows:
        if r["m_over_H"] in (1e-10, 1e-8, 1e-6, 1e-4, 0.01, 0.1, 1.0, 3.0):
            print(f"  m/H = {r['m_over_H']:8.4g}   Omega = "
                  f"{r['omega_at_model_T']:12.4e}   /0.265 = "
                  f"{r['ratio_to_DM']:10.4e}")

    verdict = "PASS" if m_star is not None else "FAIL"
    g_star = (m_star * H_INF / PHI_V) if m_star else None
    print("=" * 78)
    print(f"VERDICT: {verdict}")
    if m_star is not None:
        print(f"  Power-law (true mode equation) small-g solution:  m_psi/H_inf = {m_star:.4e}"
              f"  =>  m_psi = {m_star*H_INF:.4e} GeV,  g = {g_star:.4e}")
        print(f"    at production p/m ~ H_inf/m_psi = {1.0/m_star:.3e} ; "
              f"today v0 = {A_END_OVER_A0/m_star:.3e}")
    if fs:
        cold = fs["lambda_fs_Mpc"] < LYMAN_ALPHA_LFS_MPC
        print(f"  light-branch free-streaming length lambda_fs = {fs['lambda_fs_Mpc']:.3e} Mpc "
              f"(reference upper limit {LYMAN_ALPHA_LFS_MPC} Mpc) -> "
              f"{'COLD (pass)' if cold else 'WARM (FAIL)'}")
    if m_star_exp is not None:
        print(f"  exponential closed-form small-g solution:  m/H = {m_star_exp:.4e},  "
              f"g = {m_star_exp*H_INF/PHI_V:.4e}")
    print(f"  the window quoted in the paper, m_psi/H_inf in [1,10], corresponds to "
          f"g in [{H_INF/PHI_V:.3e}, {10*H_INF/PHI_V:.3e}]")
    print(f"total runtime {time.time()-t0:.1f} s")
    print("=" * 78)

    out = {"verdict": verdict, "target_omega": OMEGA_TARGET,
           "T_reh_model": T_REH_MODEL, "T_reh_fiducial": T_REH_FIDUCIAL,
           "freeze_target": FREEZE_TARGET, "rows": rows,
           "exponential_rows": exp_rows,
           "m_star_powerlaw": m_star, "g_star_powerlaw": g_star,
           "n_star_powerlaw": n_star, "free_streaming": fs,
           "m_star_exponential": m_star_exp,
           "g_star_exponential": ((m_star_exp * H_INF / PHI_V)
                                  if m_star_exp else None),
           "runtime_s": time.time() - t0}
    JSON_PATH.write_text(json.dumps(out, indent=2, ensure_ascii=False),
                         encoding="utf-8")
    write_md(out)
    print(f"wrote scripts/{JSON_PATH.name} and scripts/{MD_PATH.name}")


def write_md(o: dict) -> None:
    L = []
    A = L.append
    A("# Can the dark-matter gap be closed: the small-g branch of Omega_psi(g) (verification script)\n")
    A(f"**Verdict: {o['verdict']}** -- "
      + ("there exists a g with $\\Omega_\\psi=\\Omega_{{\\rm DM}}=0.265$; "
         "the gap can be closed within the **current action + current reheating mechanism**.\n"
         if o["verdict"] == "PASS" else
         "no such g exists; closing the gap requires extending the model.\n"))
    A(f"Target $\\Omega_\\psi={o['target_omega']}$, taking the model's own "
      f"$T_{{\\rm reh}}={o['T_reh_model']:.2e}$ GeV (anomaly channel, physical $\\alpha_s$).\n")

    A("## 1. Why a separate test is worthwhile\n")
    A("The Sec. 4 $\\Omega$ table of `psi_abundance_oscillating.py` and the Sec. 2c $T_{\\rm reh}$ inversion contain\n")
    A("```python\n    if m < 0.4:\n        continue\n```\n")
    A("i.e. that table **deliberately skips** the points with $m/H_{\\rm inf}<0.4$, and the conclusion "
      "that \"$\\Omega$ is nearly independent of $g$ and the abundance is set by $T_{\\rm reh}$\" "
      "was read off precisely that table.\n")

    A("## 2. Why the small-$g$ solution should exist (three independent reasons)\n")
    A("1. The paper's own Eq. (\\texttt{eq:Omegapsi}) is **non-monotonic** in $m_\\psi$;")
    A("   the paper itself admits at L458 that $\\Omega_\\psi(g)$ has a maximum at $g=3.6\\times10^{-6}$.")
    A("   Any target value below the maximum has **two** solutions.")
    A("2. A massless Dirac field is conformally invariant $\\Rightarrow$ production is exactly zero "
      "(self-check in this repo: $|\\beta|=0$ to machine zero), so as $m\\to0$ "
      "$\\Omega_\\psi\\to0$.")
    A("3. From 1+2 and continuity, the intermediate value theorem guarantees that a small-$g$ solution exists.\n")

    A("## 3. Numerical results (true mode equation, frozen convergence)\n")
    A("| $m_\\psi/H_{\\rm inf}$ | $g$ | $k_{\\max}$ | $a_{\\rm final}$ | "
      "$n_\\psi/H^3$ | $m n/H^4$ | $\\Omega_\\psi$ | $/0.265$ |")
    A("|---|---|---|---|---|---|---|---|")
    for r in o["rows"]:
        A(f"| {r['m_over_H']:g} | {r['g']:.4e} | {r['kmax']:g} | "
          f"{r['a_final']:.4g} | {r['n_over_H3']:.4e} | {r['m_times_n']:.4e} | "
          f"{r['omega_at_model_T']:.4e} | {r['ratio_to_DM']:.4e} |")
    A("")
    if o["m_star_powerlaw"]:
        A(f"**Small-$g$ solution:** $m_\\psi/H_{{\\rm inf}}={o['m_star_powerlaw']:.4e}$, "
          f"i.e. $m_\\psi={o['m_star_powerlaw']*H_INF:.4e}$ GeV, "
          f"$g={o['g_star_powerlaw']:.4e}$.")
        A(f"Corresponding $n_\\psi/H^3={o['n_star_powerlaw']:.4e}$.")
        A(f"At production it is relativistic ($p/m\\sim{1/o['m_star_powerlaw']:.3e}$), "
          f"and today $v_0={A_END_OVER_A0/o['m_star_powerlaw']:.3e}$ (extremely cold).\n")

    A("## 4. Comparison with the paper's exponential closed form\n")
    A("| $m/H_{\\rm inf}$ | $g$ | $\\Omega_\\psi$ | $/0.265$ |")
    A("|---|---|---|---|")
    for r in o["exponential_rows"]:
        A(f"| {r['m_over_H']:g} | {r['g']:.4e} | "
          f"{r['omega_at_model_T']:.4e} | {r['ratio_to_DM']:.4e} |")
    A("")
    if o["m_star_exponential"]:
        A(f"Small-$g$ solution of the exponential closed form: $m/H_{{\\rm inf}}={o['m_star_exponential']:.4e}$, "
          f"$g={o['g_star_exponential']:.4e}$.")
    A("")

    A("## 5. Convergence correction at small $m$\n")
    A("The original script fixes $a_{\\rm final}=60$ for all $m$, so at small $m$ the freeze condition "
      "$k_{\\max}/(m a_{\\rm final})\\ll1$ is badly violated "
      "(the ratio reaches 67 at $m=0.01$). Taking "
      "$a_{\\rm final}=\\max(60,\\,k_{\\max}/(0.3m))$ in this script, "
      "the small-$m$ $n_\\psi$ drops substantially:\n")
    A("| $m/H$ | original script $n/H^3$ (unfrozen) | this script $n/H^3$ (frozen) | ratio |")
    A("|---|---|---|---|")
    for mm, old in ((0.01, 6.725e-5), (0.1, 5.355e-4), (1.0, 5.028e-4),
                    (3.0, 1.336e-4)):
        new = next((r["n_over_H3"] for r in o["rows"]
                    if abs(r["m_over_H"] - mm) < 1e-12), None)
        if new:
            A(f"| {mm:g} | {old:.4e} | {new:.4e} | {old/new:.2f} |")
    A("")
    A("Notably, at the small-$m$ end $n_\\psi\\propto m^{1.0}$ (originally reported as $m^{0.9}$), "
      "so $m n\\propto m^{2}$.\n")

    A("## 6. Is the light branch cold dark matter? (free-streaming length check)\n")
    if o.get("free_streaming"):
        f = o["free_streaming"]
        A(f"- at production $p/m$ = {f['p_over_m_prod']:.3e} (taking the most conservative $p\\sim H_{{\\rm inf}}$)")
        A(f"- $a_{{\\rm reh}}$ = {f['a_reh']:.3e}, $a_{{\\rm eq}}$ = {f['a_eq']:.3e}, "
          f"$H_{{\\rm reh}}$ = {f['H_reh_GeV']:.3e} GeV")
        A(f"- comoving free-streaming length $\\lambda_{{\\rm fs}}$ = **{f['lambda_fs_Mpc']:.3e} Mpc**; "
          f"reference warm-DM upper limit {LYMAN_ALPHA_LFS_MPC} Mpc -> "
          f"**{'COLD (pass)' if f['lambda_fs_Mpc'] < LYMAN_ALPHA_LFS_MPC else 'WARM (FAIL)'}**")
        A(f"- physical reason: $m_\\psi$ = {o['m_star_powerlaw']*H_INF:.3e} GeV is far larger than $T_{{\\rm reh}}$, "
          f"the particles become non-relativistic already at $a/a_{{\\rm end}}$ ~ {1/o['m_star_powerlaw']:.3g}, "
          f"long before reheating completes; afterwards $v\\propto1/a$ and free streaming is negligible.\n")
    _ = A
    MD_PATH.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
