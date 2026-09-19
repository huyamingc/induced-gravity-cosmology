"""临时探针：真实模方程下 Omega=0.265 在**重支侧**的交点是否存在、在哪？

背景：dm_gap_closure_test.py 的幂律扫描停在 m/H = 3，那里 Omega = 1.6e3
（仍比 0.265 高 6e3 倍），所以"重支"在真实模方程下的匹配点从未被定位。
本探针只补几个 m/H 较大的点，复用同一套模方程与归一化。
用完即删，不改动任何既有产物。
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
