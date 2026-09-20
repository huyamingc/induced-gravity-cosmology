# Shared constants and model functions for induced-gravity cosmology checks.
# Fiducial parameters follow paper_prd_merged.tex (analytic A5 normalization).
from __future__ import annotations

import functools
import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

# --- fundamental constants (reduced Planck mass) ---
M_P = 2.435e18  # GeV
A_S = 2.100e-9  # Planck 2018
N_S_PLANCK = 0.9649
N_S_SIGMA = 0.0042
XI_FID = 11.1
N_FID = 50
H0_KMS_MPC = 67.4
MPC_KM = 3.086e19
HBAR_GEV_S = 6.582e-25

# Paper (merged draft) convention: beta_p = 2/sqrt(6+1/xi)
# Our convention: beta_o = sqrt(6+1/xi), r = 2 beta_o^2 / N^2 = 8/(beta_p^2 N^2)


def beta_o(xi: float) -> float:
    return math.sqrt(6.0 + 1.0 / xi)


def beta_p(xi: float) -> float:
    return 2.0 / beta_o(xi)


def lambda0_analytic(xi: float, N: float, As: float = A_S) -> float:
    """Large-field analytic: A_s = lambda0 N^2 / (12 pi^2 xi^2 (6+1/xi))."""
    b2 = 6.0 + 1.0 / xi
    return 12.0 * math.pi**2 * xi**2 * b2 * As / N**2


def lambda0_for_As_locked(N: float, xi: float, As: float = A_S) -> float:
    """Exact potential slow-roll A_s inversion (paper locked-N primary).

    Imported lazily to avoid circular import with derive_from_action.
    """
    from derive_from_action import lambda0_for_As as _lam

    return _lam(N, xi, As)[0]


def lambda0_paper_A5_wrong(xi: float, N: float, As: float = A_S) -> float:
    """Incorrect intermediate formula printed in some drafts: 48 pi^2 xi^2 As / N^2."""
    return 48.0 * math.pi**2 * xi**2 * As / N**2


def ns_LO(N: float) -> float:
    return 1.0 - 2.0 / N


def ns_NLO_starobinsky(N: float) -> float:
    """NLO with Starobinsky-limit coefficient -3/(2N^2)."""
    return 1.0 - 2.0 / N - 1.5 / N**2


def r_of(xi: float, N: float) -> float:
    return 2.0 * (6.0 + 1.0 / xi) / N**2


def H_inf(lam0: float, xi: float) -> float:
    return math.sqrt(lam0) * M_P / (2.0 * math.sqrt(3.0) * xi)


def V0(lam0: float, xi: float) -> float:
    return lam0 * M_P**4 / (4.0 * xi**2)


def U_quarter(lam0: float, xi: float) -> float:
    return V0(lam0, xi) ** 0.25


def m_chi(lam0: float, xi: float) -> float:
    return math.sqrt(2.0 * lam0) * M_P / (xi * beta_o(xi))


def m_Phi(lam0: float, xi: float) -> float:
    """Jordan-frame mass at vacuum: m = sqrt(2 lambda0) * Phi0 = sqrt(2 lambda0) M_P / sqrt(xi)."""
    return math.sqrt(2.0 * lam0) * M_P / math.sqrt(xi)


# Alias kept for older call sites
m_Phi_jordan = m_Phi


def H0_GeV() -> float:
    return H0_KMS_MPC / MPC_KM * HBAR_GEV_S


def Vc(Omega_Lambda: float = 0.683) -> float:
    H0 = H0_GeV()
    return Omega_Lambda * 3.0 * H0**2 * M_P**2


def epsilon_exact(x: float, xi: float) -> float:
    """eps = 2 beta_p^2 e^{-2x} / (1-e^{-x})^2, x = beta_p chi / M_P."""
    b = beta_p(xi)
    u = math.exp(-x)
    return 2.0 * b**2 * u**2 / (1.0 - u) ** 2


def x_end(xi: float) -> float:
    """Solve eps=1: u = 1/(1+sqrt(2)*beta_p), x=-ln u."""
    b = beta_p(xi)
    u = 1.0 / (1.0 + math.sqrt(2.0) * b)
    return -math.log(u)


def u_end(xi: float) -> float:
    """u_e = e^{-x_end} = 1/(1+sqrt(2) beta_p), fixed by eps_V = 1."""
    return math.exp(-x_end(xi))


def V_end_over_V0(xi: float) -> float:
    """Exact end-of-inflation potential ratio, V_end/V0 = (1-u_e)^2.

    A pure number: it contains no lambda0, no M_Pl and no N.  This function is
    the single source for that ratio across the repository.  The 6-digit literal
    0.285204 that used to be copied into four scripts was a truncation of this
    value (relative error 1.4e-6) and has been removed.
    """
    return (1.0 - u_end(xi)) ** 2


@functools.lru_cache(maxsize=None)
def K_over_V_end() -> float:
    """Kinetic-to-potential ratio at the end of inflation, K/V = eps_H/(3-eps_H).

    This is NOT a closed form.  Slow roll has already failed at eps_V = 1, so
    the slow-roll relation K/V = 1/3 is wrong by a factor ~1.7 there and the
    ratio has to come from the dynamical integration of the e-fold equations.
    We take it from background_and_reheating.slowroll_to_end (RK4, 2e4 steps;
    imported lazily to avoid a circular import), which returns
    eps_H,end = 0.4987655588641701.

    Like V_end/V0 this is a pure number: the e-fold ODE in M_Pl = 1 units is
    invariant under V -> c V, so the ratio carries no lambda0 and no N.
    """
    from background_and_reheating import slowroll_to_end

    eps_h = slowroll_to_end()["eps_H_end"]
    return eps_h / (3.0 - eps_h)


def rho_end_over_V_end() -> float:
    """rho_end/V_end = 1 + K/V at the end of inflation."""
    return 1.0 + K_over_V_end()


# Cosmological constants entering the reheating dilution factor.  They live here
# so that the scripts which quote a dilution no longer each carry their own copy.
T0_GEV = 2.3491e-13       # CMB temperature today
G_STAR_REH = 106.75       # g_* at reheating
G_STAR_S0 = 3.91          # g_{*s} today


def rho_end(N: float = N_FID, xi: float = XI_FID) -> float:
    """Energy density at the end of inflation: V0 (V_end/V0)(rho_end/V_end)."""
    lam0 = lambda0_for_As_locked(N, xi)
    return V0(lam0, xi) * V_end_over_V0(xi) * rho_end_over_V_end()


def dilution(T_reh: float, N: float = N_FID, xi: float = XI_FID) -> float:
    """(a_end/a_0)^3 at reheating temperature T_reh, entropy-conserving.

    Matter-like condensate domination between a_end and a_reh, radiation after:
        rho_end (a_end/a_reh)^3 = rho_rad(T_reh)
        a_reh/a_0 = (g_s,0/g_s,reh)^(1/3) (T_0/T_reh)

    psi_abundance_oscillating, dm_gap_closure_test and residual_quintessence all
    quote this one quantity.  It used to exist as three separate literals
    (1.0204e-92, and 1.0204e-101 in two places), which is why it is a function here.
    """
    rho_e = rho_end(N, xi)
    return (((math.pi**2 / 30.0) * G_STAR_REH * T_reh**4 / rho_e)
            * (G_STAR_S0 * T0_GEV**3) / (G_STAR_REH * T_reh**3))


def VE_of_varphi(varphi, lam0: float, xi: float, Vc_val: float) -> np.ndarray:
    varphi = np.asarray(varphi, dtype=float)
    return (lam0 * M_P**4 / (4.0 * xi**2)) * (1.0 - np.exp(-2.0 * varphi)) ** 2 + Vc_val * np.exp(-4.0 * varphi)


def F_of_varphi(varphi, xi: float) -> np.ndarray:
    """F = xi Phi^2 = xi Phi0^2 e^{2 varphi} = M_P^2 e^{2 varphi}."""
    varphi = np.asarray(varphi, dtype=float)
    return M_P**2 * np.exp(2.0 * varphi)


def N_match(T_reh: float, V_end: float) -> float:
    """N = 61 + (1/4) ln(V_end/M_P^4) - (1/3) ln(V_end^{1/4}/T_reh)."""
    return 61.0 + 0.25 * math.log(V_end / M_P**4) - (1.0 / 3.0) * math.log(V_end**0.25 / T_reh)


def T_reh_from_N(N: float, V_end: float) -> float:
    """Invert N_match for T_reh."""
    # N - 61 - (1/4)ln(V/M^4) = -(1/3) ln(V^{1/4}/T)
    rhs = -(N - 61.0 - 0.25 * math.log(V_end / M_P**4))
    # rhs = (1/3) ln(V^{1/4}/T) => 3 rhs = ln(V^{1/4}/T) => V^{1/4}/T = exp(3 rhs)
    return V_end**0.25 / math.exp(3.0 * rhs)


@dataclass
class FiducialPoint:
    xi: float
    N: float
    lam0: float
    beta_o: float
    beta_p: float
    r: float
    ns_lo: float
    ns_nlo: float
    H_inf: float
    V0: float
    U_quarter: float
    m_chi: float
    m_Phi: float
    Phi0: float
    Vc: float
    H0: float
    u_end: float
    V_end: float
    V_end_frac: float
    m_chi_over_H0: float


def fiducial(xi: float = XI_FID, N: float = N_FID, use_numeric_lam0: bool = False, use_locked_lam0: bool = False) -> FiducialPoint:
    """Fiducial scales.

    use_numeric_lam0: legacy Table-I draft lambda0=6.78e-8 (attractor-matched draft).
    use_locked_lam0: exact slow-roll A_s inversion at this (N, xi) -- paper primary.
    Default False keeps large-field analytic lambda0 (App. A5) for comparison only.
    """
    if use_locked_lam0:
        lam0 = lambda0_for_As_locked(N, xi)
    elif use_numeric_lam0:
        lam0 = 6.78e-8
    else:
        lam0 = lambda0_analytic(xi, N)
    V_end_frac = (1.0 - 1.0 / (1.0 + math.sqrt(2.0) * beta_p(xi))) ** 2
    V0_ = V0(lam0, xi)
    H0 = H0_GeV()
    return FiducialPoint(
        xi=xi,
        N=N,
        lam0=lam0,
        beta_o=beta_o(xi),
        beta_p=beta_p(xi),
        r=r_of(xi, N),
        ns_lo=ns_LO(N),
        ns_nlo=ns_NLO_starobinsky(N),
        H_inf=H_inf(lam0, xi),
        V0=V0_,
        U_quarter=U_quarter(lam0, xi),
        m_chi=m_chi(lam0, xi),
        m_Phi=m_Phi(lam0, xi),
        Phi0=M_P / math.sqrt(xi),
        Vc=Vc(),
        H0=H0,
        u_end=1.0 / (1.0 + math.sqrt(2.0) * beta_p(xi)),
        V_end=V0_ * V_end_frac,
        V_end_frac=V_end_frac,
        m_chi_over_H0=m_chi(lam0, xi) / H0,
    )
