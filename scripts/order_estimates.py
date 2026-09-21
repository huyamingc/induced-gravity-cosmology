#!/usr/bin/env python3
"""
Order-of-magnitude estimate library -- the single source for every claim the
manuscript prints with a leading "~"
=============================================================================
Type:           CORE
Paper Sec.:     III, IV, V, VIII E, XII
Experiment:     order-of-magnitude claim support
What it does:   Holds the closed forms behind the order-of-magnitude numbers
                printed in main.tex: primordial non-Gaussianity,
                the spectral running, the Mathieu parameter of parametric
                resonance, the thermalization rate (with the one-loop running
                coupling anchored to the same alpha_s chain Sec. IV prints),
                the direct-detection cross section, the kinematic-closure
                coupling, the Lambda_J threshold suppressions and the
                gravity-induced Higgs portal (Sec. VII), the one-loop RG
                shifts of xi and lambda0, and the quartic Delta-w chain of
                Discussion #14.

                This module exists so that each of those formulas has exactly
                ONE implementation.  They used to live in two standalone
                scripts, extended_checks.py and quick_claims_check.py, which
                were retired.  The manuscript still prints every one of these
                numbers, so the formulas were moved here rather than deleted,
                and verify_numerics.py checks the paper against this module.
                Any script that needs one of these quantities must import it
                from here instead of re-deriving it.

                The INPUTS (lambda0, m_chi, n_s, ...) are the exact locked-N
                values.  Only the verdicts are order-of-magnitude, because the
                microscopic coefficients b_s, the anomaly alpha_s and g_* are
                not derived from the model.  Nothing here should be quoted as
                a precise prediction.

Constants:      M_P, G_STAR_REH, XI_FID and N_FID come from cosmo_model, which
                remains the single source for the exact quantities.  This file
                does not redefine any of them.
Conventions:    DLNMU = 60 is the inflationary-to-electroweak running range
                used in Eqs. (dxilam)/(dxig) of the manuscript.
Outputs:        none (library).  Running this file directly prints a
                human-readable self-check; the machine check on the manuscript
                numbers themselves lives in verify_numerics.py.
Dependencies:   cosmo_model (constants), derive_from_action + scipy (n_at_r)
=============================================================================
"""
from __future__ import annotations

import math
import os
import sys

os.environ.setdefault("PYTHONUNBUFFERED", "1")
if hasattr(sys.stdout, "reconfigure"):          # Python 3.7+
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)

from cosmo_model import G_STAR_REH, H0_GeV, M_P, N_FID, XI_FID  # noqa: E402
from cosmo_model import m_Phi, m_chi  # noqa: E402

HBAR_C = 1.973269804e-14       # GeV cm
G_PER_GEV = 1.782661921e-24    # grams per GeV (c^2 included)
DLNMU = 60.0                   # running range of Eqs. (dxilam)/(dxig)
M_Z = 91.19                    # GeV, Z pole: the alpha_s anchor of Sec. IV
ALPHA_S_MZ = 0.1179            # alpha_s(m_Z), the value the manuscript prints
B3_THERM = 7.0                 # SM b_3 (n_f = 5), the same coefficient Sec. IV uses
G_ANOM = 6.9e-5                # anomaly-framework matching point of Sec. III
M_TOP = 173.0                  # GeV


# --------------------------------------------------------------------------
# Sec. III -- inflation
# --------------------------------------------------------------------------
def f_nl_local(ns: float) -> float:
    """Local primordial non-Gaussianity, f_NL = -5/12 (1 - n_s).

    The manuscript quotes ~-0.02; the exact arithmetic on its own printed
    n_s = 0.9616 gives -0.0160, i.e. the same decade.
    """
    return -5.0 / 12.0 * (1.0 - ns)


# --------------------------------------------------------------------------
# Sec. IV -- reheating
# --------------------------------------------------------------------------
def alpha_s_running(Q: float, Q_ref: float = M_Z,
                    alpha_s_ref: float = ALPHA_S_MZ,
                    b3: float = B3_THERM) -> float:
    """One-loop running alpha_s(Q) = alpha_s(Q_ref)/(1 + b3 alpha_s/(2pi) ln(Q/Q_ref)).

    Anchored to the SAME chain the manuscript prints in Sec. IV
    (alpha_s(m_Z) = 0.1179, b_3 = 7).  Cross-check: running this anchor up to
    m_chi = 3.25e13 GeV reproduces the printed alpha_s(m_chi) = 0.0262 to
    better than 1%, so the thermalization coupling below is not a second,
    inconsistent estimate.

    The manuscript once used a hand-estimated alpha_s ~ 0.1 here, which
    inflated Gamma_therm/H by about one decade (1e7 instead of ~1e6).
    """
    return alpha_s_ref / (1.0 + (b3 * alpha_s_ref / (2.0 * math.pi))
                          * math.log(Q / Q_ref))


# Thermalization coupling at T = 1e9 GeV: one-loop running of the Sec. IV
# anchor.  Computed once here so that gamma_therm_over_H's default argument
# is a single, locked constant (0.0377).
ALPHA_S_THERM = alpha_s_running(1.0e9)   # = 0.0377


def mathieu_q(m_f: float, m_chi: float) -> float:
    """Mathieu parameter of parametric resonance, q ~ (m_f/m_chi)^2.

    The manuscript quotes q ~ 2.8e-23 for the top quark, over twenty orders
    below the broad-resonance threshold q >= 0.1.
    """
    return (m_f / m_chi) ** 2


def gamma_therm_over_H(T: float, alpha_s: float = ALPHA_S_THERM,
                       g_star: float = G_STAR_REH) -> float:
    """Thermalization rate over Hubble, Gamma_therm/H, at temperature T.

    Gamma_therm ~ alpha_s^2 T and H = sqrt(pi^2 g_*/90) T^2/M_Pl, with the
    default alpha_s = ALPHA_S_THERM the one-loop running value 0.0377 at
    T = 1e9 GeV.  The manuscript quotes ~1e6 at T = 1e9 GeV (Sec. IV and
    Sec. XII); the retired hand estimate alpha_s ~ 0.1 gave ~1e7.
    """
    gamma = alpha_s**2 * T
    H = math.sqrt(math.pi**2 * g_star / 90.0) * T**2 / M_P
    return gamma / H


# --------------------------------------------------------------------------
# Sec. V -- dark matter
# --------------------------------------------------------------------------
def sigma_psiN_cm2(m_N: float = 1.0) -> float:
    """psi-nucleon cross section from graviton exchange, in cm^2.

    sigma_psiN ~ G_N^2 m_N^2 (hbar c)^2 with G_N = 1/(8 pi M_Pl^2).  The
    manuscript quotes ~1e-104 cm^2, far below the XENONnT/LZ floor.
    """
    G_N = 1.0 / (8.0 * math.pi * M_P**2)
    return G_N**2 * m_N**2 * HBAR_C**2


def g_kinematic(m_chi: float, xi: float = XI_FID) -> float:
    """Coupling at which m_psi = m_chi/2, i.e. g = m_chi/(2 Phi_V).

    Phi_V = M_Pl/sqrt(xi).  The manuscript states that kinematic closure of
    chi -> psi psibar requires g >~ few x 1e-5 at xi = 11.1.
    """
    return m_chi / (2.0 * M_P / math.sqrt(xi))


def r_of_locked_at_N(N: float, xi: float = XI_FID) -> float:
    """Exact potential slow-roll r at a locked e-fold number N."""
    import derive_from_action as dfa

    bp = dfa.beta_p(xi)
    return dfa.observables_at_x(dfa.x_star_for_N(N, bp), bp, 1.0)["r_ps"]


def n_at_r(r_target: float, xi: float = XI_FID,
           n_lo: float = 20.0, n_hi: float = 48.0) -> float:
    """Locked-N e-fold number at which r equals r_target (exact slow roll).

    The manuscript quotes r = 0.01 at N ~ 32 for xi = 11.1, which is what
    makes a LiteBIRD r > 0.01 detection exclude the model class.
    """
    from scipy.optimize import brentq

    return brentq(lambda n: r_of_locked_at_N(n, xi) - r_target, n_lo, n_hi)


# --------------------------------------------------------------------------
# Sec. VII -- SM embedding: decoupling and the induced Higgs portal
# --------------------------------------------------------------------------
def lambda_J(xi: float = XI_FID) -> float:
    """Induced-gravity UV cutoff Lambda_J = M_Pl/xi.

    The manuscript quotes Lambda_J ~ 2.2e17 GeV at xi = 11.1.
    """
    return M_P / xi


def threshold_suppression(m: float, xi: float = XI_FID) -> float:
    """Heavy-field threshold suppression (m/Lambda_J)^2.

    Sec. VII: (m_chi/Lambda_J)^2 ~ 2.2e-8 and (m_psi/Lambda_J)^2 ~ 1.2e-13.
    An earlier draft printed a single bound <= 1e-10, which the chi value
    exceeds by two orders; the two decoupled values are what the manuscript
    now prints.
    """
    return (m / lambda_J(xi)) ** 2


def delta_lambda_PhiH(lam0: float, xi: float = XI_FID) -> float:
    """Gravity-induced Higgs portal, delta_lambda_PhiH ~ xi lambda0/(16 pi^2) (m_Phi/M_Pl)^2.

    Sec. VII: with m_Phi = sqrt(2 lambda0) M_Pl/sqrt(xi) ~ 2.7e14 GeV this
    gives ~6e-17.  An earlier draft printed ~1e-15, which overstated the
    portal by about two orders (the conclusion of negligibility is unchanged).
    """
    return xi * lam0 / (16.0 * math.pi**2) * (m_Phi(lam0, xi) / M_P) ** 2


# --------------------------------------------------------------------------
# Sec. VIII E -- one-loop RG running
# --------------------------------------------------------------------------
def delta_lambda0(lam0: float, dlnmu: float = DLNMU) -> float:
    """One-loop shift of lambda0, Delta lambda0 ~ 9 lam0^2/(8 pi^2) dlnmu.

    The manuscript quotes ~3e-14 over dlnmu ~ 60.  Beta_lambda0 ~ lam0^2 > 0,
    so the self-coupling decreases monotonically from the UV to the IR.
    """
    return 9.0 * lam0**2 / (8.0 * math.pi**2) * dlnmu


def delta_xi_quartic(xi: float, lam0: float, dlnmu: float = DLNMU) -> float:
    """Quartic-curvature drive: Delta xi ~ (xi - 1/6) 3 lam0/(16 pi^2) dlnmu.

    This is the direct curvature-operator renormalization from the
    quartic-curvature cross term.  The manuscript quotes ~8.5e-7; it is
    g-independent and therefore sets the floor of the running.
    """
    return (xi - 1.0 / 6.0) * 3.0 * lam0 / (16.0 * math.pi**2) * dlnmu


def delta_xi_yukawa(xi: float, g: float, dlnmu: float = DLNMU) -> float:
    """Yukawa drive: Delta xi ~ (xi - 1/6) g^2/(16 pi^2) dlnmu.

    The manuscript quotes 4.2e-14 (g = 1.0e-7), 2.2e-9 (g = 2.3e-5) and
    4.2e-8 (g = 1.0e-4).
    """
    return (xi - 1.0 / 6.0) * g**2 / (16.0 * math.pi**2) * dlnmu


def delta_xi_total(xi: float, lam0: float, g: float,
                   dlnmu: float = DLNMU) -> float:
    """Total Delta xi; the two drives are additive and carry the same sign.

    The manuscript states Delta xi <= 1e-6 at xi = 11.1 and gives the range
    [6.5e-8, 7.8e-6] over xi in [1, 100] at the anomaly coupling g = 6.9e-5.
    """
    return (delta_xi_quartic(xi, lam0, dlnmu)
            + delta_xi_yukawa(xi, g, dlnmu))


def dn_quartic_over_yukawa(xi: float, lam0: float, g: float,
                           dlnmu: float = DLNMU) -> float:
    """Ratio of the quartic drive to the Yukawa drive at coupling g.

    The manuscript quotes ~2e7 at the abundance-matched g = 1.0e-7.
    """
    return (delta_xi_quartic(xi, lam0, dlnmu)
            / delta_xi_yukawa(xi, g, dlnmu))


# --------------------------------------------------------------------------
# Sec. XII -- Discussion closed forms
# --------------------------------------------------------------------------
def sigma_psipsi_over_m(m_psi: float, m_N: float = 1.0) -> float:
    """psi--psi self-interaction cross section per unit mass, in cm^2/g.

    sigma_psipsi ~ G_N^2 m_psi^2 (hbar c)^2 (graviton t-channel exchange) and
    the mass is converted with G_PER_GEV, so this is sigma/m in cm^2/g:

        sigma/m = G_N^2 m_psi (hbar c)^2 / G_PER_GEV .

    ``m_N`` is an unused-by-default reference mass kept in the signature only to
    mirror sigma_psiN_cm2; the scaling is linear in m_psi, so the result is the
    SAME G_N and the same (hbar c)^2 convention the manuscript uses for
    sigma_psiN in Sec. V.

    The manuscript prints ~7e-70 cm^2/g (about 69 orders below the bullet-cluster
    bound of 1 cm^2/g) at the light-branch m_psi.  NOTE: an earlier draft
    printed ~1e-67 cm^2/g, which is what this formula gives for m_psi = 1e13 GeV,
    not for the m_psi = 7.5e10 GeV that Sec. V actually selects.  Sec. XII had
    kept the old wimpzilla-scale mass; it is now tied to the matching value.
    """
    G_N = 1.0 / (8.0 * math.pi * M_P**2)
    sigma_cm2 = G_N**2 * m_psi**2 * HBAR_C**2
    m_grams = m_psi * G_PER_GEV
    return sigma_cm2 / m_grams


def tremaine_gunn_Q(m_psi: float, H_inf: float) -> float:
    """Non-thermal phase-space density Q = rho/sigma_v^3 at production, GeV^4.

    With the parametric number density of Eq. (npsi), n_psi ~ H_inf^3, and the
    velocity dispersion the manuscript quotes, sigma_v ~ H_inf/m_psi:

        rho     = m_psi n_psi ~ m_psi H_inf^3
        Q       = rho/sigma_v^3 ~ m_psi H_inf^3 (m_psi/H_inf)^3 = m_psi^4 .

    The closed form is written out rather than collapsed to ``m_psi**4`` so the
    two cancellations that produce it stay visible.

    The manuscript prints ~3e43 GeV^4 at the light-branch m_psi.  NOTE: an
    earlier draft printed ~1e52 GeV^4, which is m_psi^4 for m_psi = 1e13 GeV --
    the same stale mass as in sigma_psipsi_over_m.  Either way Q is far above the
    dwarf-galaxy lower limit, which is the only statement the entry makes.
    """
    n_psi = H_inf**3                 # parametric, Eq. (npsi), up to O(1)
    rho = m_psi * n_psi
    sigma_v = H_inf / m_psi
    return rho / sigma_v**3


def alpha_s_attractor(N: float = float(N_FID)) -> float:
    """Spectral running from the attractor closed form, alpha_s ~ -2/N^2.

    The manuscript quotes ~-8e-4 at the locked N = 50.
    """
    return -2.0 / N**2


def dw_quartic_ricci(lam0: float, xi: float = XI_FID,
                     m_chi_val: float = 0.0) -> float:
    """Quartic correction to w at the Ricci-frozen field value (Discussion #14).

    The chain is delta w/w ~ (lambda0/xi^2) varphi^2 with the frozen
    displacement varphi ~ (H0/m_chi)^2 ~ 2e-111 (the Sec. VI Ricci driving),
    so the correction is (lambda0/xi^2) (H0/m_chi)^4 ~ 2e-231.  An earlier
    draft printed the inconsistent chain 1e-6 * (1e-46)^2 = 1e-98: the
    coefficient was lambda0/xi^2 = 5.4e-10, not 1e-6, and the frozen value is
    the squared ratio, not H0/m_chi itself.
    """
    if m_chi_val <= 0.0:
        m_chi_val = m_chi(lam0, xi)
    varphi = (H0_GeV() / m_chi_val) ** 2
    return (lam0 / xi**2) * varphi**2


def _light_branch():
    """(m_psi, H_inf) at the abundance-matched light branch, or (None, None).

    m_psi is an OUTPUT of dm_gap_closure_test.py, so it is read from that
    script's artefact instead of being duplicated here as a constant.  H_inf
    comes from cosmo_model, the single source for the exact quantities.
    """
    import json
    from pathlib import Path

    from cosmo_model import H_inf as _H_inf
    from cosmo_model import lambda0_for_As_locked as _lam0

    p = Path(__file__).resolve().parent / "dm_gap_closure_test.json"
    if not p.exists():
        return None, None
    with p.open(encoding="utf-8") as fh:
        d = json.load(fh)
    H = _H_inf(_lam0(N_FID, XI_FID), XI_FID)
    return d["m_star_powerlaw"] * H, H


# --------------------------------------------------------------------------
# self-check
# --------------------------------------------------------------------------
def _selfcheck() -> None:
    """Print every quantity this module provides at the manuscript fiducial."""
    from cosmo_model import lambda0_for_As_locked

    lam0 = lambda0_for_As_locked(N_FID, XI_FID)
    mchi = math.sqrt(2.0 * lam0) * M_P / (XI_FID * math.sqrt(6.0 + 1.0 / XI_FID))
    ns = 0.9616
    print("[order_estimates] fiducial: xi=%s N=%s lam0=%.6e m_chi=%.6e"
          % (XI_FID, N_FID, lam0, mchi))
    print("  Sec. III  f_NL^local          = %+.4e   (paper ~-0.02)" % f_nl_local(ns))
    print("  Sec. IV   q (top)             = %.4e   (paper 2.8e-23)" % mathieu_q(M_TOP, mchi))
    print("  Sec. IV   Gamma_therm/H       = %.4e   (paper ~1e6)" % gamma_therm_over_H(1.0e9))
    print("  Sec. IV   alpha_s(1e9 GeV)    = %.4e   (paper ~0.038; -> alpha_s(m_chi)=%.4e)"
          % (ALPHA_S_THERM, alpha_s_running(3.25e13, 1.0e9, ALPHA_S_THERM)))
    print("  Sec. V    sigma_psiN [cm^2]   = %.4e   (paper ~1e-104)" % sigma_psiN_cm2())
    print("  Sec. V    g kinematic         = %.4e   (paper few x 1e-5)" % g_kinematic(mchi))
    print("  Sec. VII  (m_chi/Lambda_J)^2  = %.4e   (paper 2.2e-8)"
          % threshold_suppression(mchi))
    print("  Sec. VII  delta lambda_PhiH   = %.4e   (paper ~6e-17)"
          % delta_lambda_PhiH(lam0))
    print("  Sec. VIII Delta lambda0       = %.4e   (paper ~3e-14)" % delta_lambda0(lam0))
    print("  Sec. VIII Delta xi quartic    = %.4e   (paper ~8.5e-7)" % delta_xi_quartic(XI_FID, lam0))
    for g in (1.0e-7, 2.3e-5, 1.0e-4):
        print("  Sec. VIII Delta xi Yukawa g=%.1e = %.4e" % (g, delta_xi_yukawa(XI_FID, g)))
    print("  Sec. VIII quartic/Yukawa      = %.4e   (paper ~2e7)"
          % dn_quartic_over_yukawa(XI_FID, lam0, 1.0e-7))
    lo = delta_xi_total(1.0, lam0, G_ANOM)
    hi = delta_xi_total(100.0, lam0, G_ANOM)
    print("  Sec. VIII Delta xi(xi) range  = [%.3e, %.3e]  (paper [6.5e-8, 7.8e-6], g=6.9e-5)"
          % (lo, hi))
    print("  Sec. VIII Delta xi(xi=11.1)   = %.4e   (paper <= 1e-6)"
          % delta_xi_total(XI_FID, lam0, G_ANOM))
    print("  Sec. XII  alpha_s             = %+.4e  (paper ~-8e-4)" % alpha_s_attractor())
    print("  Sec. XII  Delta w quartic     = %.4e   (paper ~2e-231)"
          % dw_quartic_ricci(lam0, XI_FID, mchi))
    m_psi, H_i = _light_branch()
    if m_psi is None:
        print("  Sec. XII  (dm_gap_closure_test.json absent -- "
              "sigma_psipsi/m_psi and Tremaine-Gunn Q skipped)")
    else:
        print("  Sec. VII  (m_psi/Lambda_J)^2  = %.4e   (paper 1.2e-13)"
              % threshold_suppression(m_psi))
        print("  Sec. XII  sigma_psipsi/m_psi  = %.4e cm^2/g (paper ~7e-70)"
              % sigma_psipsi_over_m(m_psi))
        print("  Sec. XII  Tremaine-Gunn Q     = %.4e GeV^4  (paper ~3e43)"
              % tremaine_gunn_Q(m_psi, H_i))
    print("  Sec. III  N at r = 0.01       = %.3f    (paper ~32)" % n_at_r(0.01))
    print("[order_estimates self-check complete]")


if __name__ == "__main__":
    _selfcheck()
