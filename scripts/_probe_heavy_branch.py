"""Diagnostic probe: does the Omega = 0.265 crossing exist on the HEAVY-BRANCH
side under the true mode equation, and where?

Background: the power-law scan in dm_gap_closure_test.py stops at m/H = 3,
where Omega = 1.6e3 (still 6e3 times above 0.265), so the heavy-side matching
point under the true mode equation had never been located. This probe only
fills in a few larger m/H points, reusing the same mode equation and the same
normalization. Diagnostic retained for reference; it is not part of
run_all.py, writes no files, and modifies no existing artifact.
"""
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from dm_gap_closure_test import (  # noqa: E402
    H_INF, PHI_V, OMEGA_DM, T_REH_MODEL,
    a_final_for, omega, npsi_from_spectrum,
)

MASSES = [4.0, 5.0, 6.0, 8.0, 12.0]

print(f"H_inf = {H_INF:.4e}   Phi_V = {PHI_V:.4e}   "
      f"T_reh = {T_REH_MODEL:.2e}   target = {OMEGA_DM}")
print(f"{'m/H':>7} {'g':>11} {'kmax':>7} {'a_fin':>8} {'n/H^3':>11} "
      f"{'m*n/H^4':>11} {'Om(T)':>11} {'/0.265':>10} {'s':>6}")
for m in MASSES:
    kmax = max(40.0, 20.0 * m)
    kmin = min(0.05, m / 10.0)
    a_fin = a_final_for(m, kmax)
    t0 = time.time()
    r = npsi_from_spectrum(m, a_final=a_fin, kmax=kmax, kmin=kmin, nk_pts=45)
    n = r["n_over_H3"]
    w = omega(m, n, T_REH_MODEL)
    print(f"{m:7.3g} {m * H_INF / PHI_V:11.4e} {kmax:7.1f} {a_fin:8.1f} "
          f"{n:11.4e} {m * n:11.4e} {w:11.4e} {w / OMEGA_DM:10.4e} "
          f"{time.time() - t0:6.1f}", flush=True)
