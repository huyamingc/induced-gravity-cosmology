# Induced-Gravity Cosmology (EPJC submission repository)

Paper: **Induced-Gravity Cosmology: Inflation, Dark Energy, and Dark Matter from a Common Scalar Origin**
Target journal: **Eur. Phys. J. C** (Springer Nature `sn-jnl`; SCOAP3 covers the whole journal, APC = 0).

## Repository layout

```text
paper_prd_merged.tex     # the manuscript (EPJC / sn-jnl)
paper_prd_merged.pdf     # compiled output (pdflatex x2)
sn-jnl.cls               # Springer Nature journal class
sn-mathphys-num.bst      # numbered physics bibliography style
figures/*.pdf            # figures used by the manuscript (fig1-fig3)
scripts/*.py             # numerical verification and figure scripts
README.md
LICENSE                  # MIT for the code; manuscript text (c) the author
.gitignore
.gitattributes
```

## Compilation

```powershell
cd <repository root>
pdflatex -interaction=nonstopmode paper_prd_merged.tex
pdflatex -interaction=nonstopmode paper_prd_merged.tex
```

Dependencies: `sn-jnl.cls`, `sn-mathphys-num.bst`, `figures/*.pdf`.
The bibliography is inlined as `thebibliography`; **no** bibtex/biber pass is needed.

## Document class and format

- `\documentclass[pdflatex,sn-mathphys-num]{sn-jnl}`
- EPJC: abstract 150-250 words; `Declarations` (Funding / Competing interests / Data / Code / Author contributions)
- Independent researcher; the affiliation records the city/country of residence (Guiyang, Guizhou, China). ORCID: 0009-0003-1406-0485.
- Funding statement: no funding was received for this study.

## Numerical scripts

Environment (local; do not commit the venv):

```powershell
python -m venv scripts\.venv
scripts\.venv\Scripts\python.exe -m pip install numpy scipy matplotlib mpmath
$env:PYTHONUNBUFFERED=1
& scripts\.venv\Scripts\python.exe scripts\run_all.py
```

`scripts/run_all.py` is the single entry point. It runs every verification and
figure script in order and regenerates the reports (`scripts/*.md`, `*.json`)
and the three figures. Measured end-to-end runtime on the author's machine
(Python 3.13, numpy 2.4 / scipy 1.18 / matplotlib 3.11) is **15--20 minutes**
in total, of which `psi_abundance_oscillating.py` takes ~5.5 min and
`dm_gap_closure_test.py` 10--15 min (the spread is machine load); `treh_error_band.py` finishes in under a
second. On a Windows console, either run through
`run_all.py` (it reconfigures the output streams to UTF-8) or set
`$env:PYTHONIOENCODING='utf-8'` before running a single script.

Scripts and their role. Every script below solves the stated equations
numerically and its numbers may be quoted. The paper also prints a block of
order-of-magnitude claims (`f_NL`, `alpha_s`, `q`, `Gamma_therm/H`,
`sigma_psiN`, the RG shifts `Delta xi` and `Delta lambda0`); those are checked
inside `verify_numerics.py` against a decade tolerance, described further down.
There is no longer a separate script whose role is to support a scaling rather
than a value.

| Script | Role | Content |
|---|---|---|
| `lock_n_convention.py` | **exact / source of Table I** | the authoritative locked-N grid: lambda0 inverted from A_s at fixed N, exact potential slow roll, and the T_reh*(N) matching. Every cell of Table I (`tab:sens`) is produced here |
| `derive_from_action.py` | exact / independent derivation | rebuilds the model from the action alone, presupposing none of the paper's formulas; holds the exact A_s inversion and the N(T_reh) matching used everywhere else |
| `background_and_reheating.py` | exact / cross-check | exact KG integration of the e-fold equations, the N window, and the two reheating channels; recomputes T_reh*(N) from first principles as an independent check on the Table I column. Shares the same exact lambda0 as `lock_n_convention.py` |
| `cosmo_model.py` | exact / shared constants | common constants, the memoised `lambda0_for_As_locked` wrapper used by the figure scripts, and the **single source** of `x_end`, the end-of-inflation ratios `V_end_over_V0`, `K_over_V_end`, `rho_end_over_V_end`, and the reheating dilution `rho_end` / `dilution` |
| `order_estimates.py` | exact inputs / order-of-magnitude closes | the **single source** for the closed forms behind every claim the paper prints with a leading `~`: `f_NL`, `alpha_s`, the Mathieu `q`, `Gamma_therm/H`, `sigma_psiN`, the kinematic-closure `g`, the one-loop RG shifts of `xi` and `lambda0`, and the two Sec. XII consistency entries `sigma_psipsi/m_psi` and the Tremaine-Gunn `Q`. `verify_numerics.py` checks the manuscript against this module, and any new script that needs one of these quantities must import it rather than re-derive it |
| `psi_production_bogoliubov.py` | exact / cross-check | de Sitter exponent 2 pi audit, g matching |
| `psi_abundance_oscillating.py` | exact / cross-check | cross-transition mode equation (power-law spectrum) |
| `dm_gap_closure_test.py` | exact / cross-check | small-g light branch and free-streaming length |
| `treh_error_band.py` | exact / error propagation | propagates the T_reh uncertainty to (N, n_s, r) and to the (g, m_psi) window |
| `residual_quintessence.py` | exact / cross-check | two-fluid integration, Delta w budget |
| `verify_numerics.py` | audit | recomputes every quantitative claim printed in the paper from the current scripts, including the order-of-magnitude block; the reverse direction of `audit_tex_numbers.py` |
| `consistency_checks.py` | audit | parses the Table I and Table 2 bodies out of the `.tex` and checks each cell against the script that produces it, then compares independent routes against each other |
| `audit_tex_numbers.py` | audit | scans the `.tex` for superseded values and checks that any survivor sits in a comparison or historical context |
| `audit_readme_numbers.py` | audit | recomputes the prose numbers in this README from the scripts. Its patterns must still match, so rewording a sentence without updating the audit is itself a failure |
| `fig1_einstein_potential.py`, `fig2_ns_r.py`, `fig3_domain_wall.py` | figure | generate `figures/*.pdf`. Each figure has exactly one writer: `fig2_ns_r.py` alone produces `figures/fig2_ns_r.pdf` |

The three auditors run last in `run_all.py` (they read the `.json` artefacts, the
manuscript tables and this README) and together take under a second.
`verify_numerics.py` and `consistency_checks.py` are the automated form of the
manual cross-checks that closed the N <-> T_reh matching review; run them after any
edit to the manuscript or to the matching code.

`run_all.py` catches the `SystemExit` that each auditor raises, so a failing audit
makes the whole run exit non-zero instead of silently truncating it. Until that was
fixed the loop stopped at `verify_numerics.py`: `consistency_checks.py` never ran as
part of a full regression, and an audit appended after it would have been dead code.

Two different things were previously conflated under one label, and the
distinction matters when a number is quoted:

- **A literal copied from an exact value.** `0.285204`, `0.19938` and `1.19938`
  are the 6- and 5-digit truncations of `(1-u_e)^2`, `K/V` and `1 + K/V` at the
  end of inflation. They had been copied into four scripts, and
  `psi_abundance_oscillating.py` separately carried `rho_end = 1.63485e63`,
  computed under the *old* fitted lambda0 = 6.70e-8 and drifted 4.3e-4 away from
  the current exact value. Those literals are gone: `cosmo_model.py` evaluates
  the closed form `V_end/V0 = (1-u_e)^2` and the memoised `K/V` once, and every
  script imports them. Neither is an approximation level any more -- the closed
  form is exact, and the dynamical integration sits 2.6e-6 from it -- so any
  future drift would be a bug, and `consistency_checks.py` fails if a literal
  reappears.
- **A claim that is order-of-magnitude in the paper.** `f_NL ~ -0.02`,
  `q ~ (m_t/m_chi)^2`, `Gamma_therm/H ~ 1e7`, `sigma_psiN ~ 1e-104 cm^2`, the
  anomaly width and the RG shifts `Delta xi`, `Delta lambda0` are
  order-of-magnitude *as physics*: they are closed forms, or they depend on
  microscopic coefficients (`b_s`, `alpha_s`, `g_*`) that the model does not
  derive. `verify_numerics.py` checks them in an `order` tolerance class, which
  measures the distance in decades rather than a relative deviation, so an
  answer good to a factor 2 passes and one off by a factor 10 does not. Their
  inputs are exact; only the verdicts are order-of-magnitude, so none of them
  should be quoted for a precise value.

Five scripts were retired in this pass, and the reason differs by script.
`extended_checks.py` and `quick_claims_check.py` held the order-of-magnitude
claims above. The paper still prints every one of those numbers, so deleting
them outright would have left Secs. III, IV, V, VIII and XII quoting values with
no executing source. Their formulas were therefore consolidated into one shared
module, `order_estimates.py`, which `verify_numerics.py` now uses to check the
manuscript against; what was dropped was the superseded content -- an
`Omega_DM = 0.12` abundance scaling and a pre-lattice `Gamma_anom` table that the
mode-equation analysis has since replaced, plus a table of review verdicts the
revision had already absorbed into the manuscript. A shared module rather than a
copy inside the audit is the point: the next script that needs one of these
quantities imports it instead of writing a seventh copy.
`evaluate_paper.py`, `n_definition_check.py` and `_probe_heavy_branch.py` had no
dependants at all: nothing imports them, the manuscript never named them, and the
first had recommended creating `audit_tex_numbers.py`, which is what made its own
report obsolete. All five are archived outside the submission tree rather than
deleted.

Two further notes on the auditors. `audit_tex_numbers.py` used to be called
twice, the second time "after the figures regenerate"; the figure scripts write
only `.pdf`, so neither the `.tex` nor `n_convention_results.json` can change
between the two calls, and the duplicate is gone. And `Delta lambda0 ~ 3e-14`
in Sec. VIII E had no executing source anywhere in the suite until
`order_estimates.py` supplied one.

A second pass closed the two remaining gaps. Sec. XII's "additional consistency
checks" list printed `sigma_psipsi/m_psi ~ 1e-67 cm^2/g` and `Q ~ m_psi^4 ~
1e52 GeV^4`, and no script -- current or retired -- had ever computed either
one; every other entry in that list states a parametric bound for which the
manuscript gives no closed form, so none of them carries an executable claim.
Both turned
out to carry a stale wimpzilla-scale mass: those two numbers are exactly what
the printed closed forms give at `m = 1e13 GeV`, not at the `m_psi = 7.5e10 GeV`
that Sec. V actually matches. `order_estimates.py` now computes them
(`sigma_psipsi_over_m`, `tremaine_gunn_Q`) from the same `m_star_powerlaw`
anchor Sec. V uses, `verify_numerics.py` checks them, and the manuscript quotes
the matching values (`~7e-70 cm^2/g`, `~69` orders below the bullet-cluster
bound; `~3e43 GeV^4`).

The light-branch mass is now formed once, as `m_star_powerlaw * H_inf`. It had
taken its `H_inf` from the `T_reh` sweep table of `treh_error_band.json`, whose
`H_inf` column drifts with `T_reh` and `N` and belongs to a different
background; the anchor ratio is measured at the fiducial `H_inf`, so pairing it
with a drifted value understated `m_psi` by 0.43 percent.

Fig. 2 has a single writer. `fig2_ns_r.py` wrote both `fig2_ns_r.pdf` and a
byte-identical `fig2_ns_r_lockedN.pdf` alias, and `lock_n_convention.py` redrew
the same figure again under the alias name, while the manuscript includes only
`figures/fig2_ns_r.pdf`. The alias is gone, `lock_n_convention.py` no longer
plots (it is a numerical script, and its `n_convention_results.json` is
unchanged), and `figures/` no longer accumulates an untracked duplicate on
every run.


## Dark-matter convention (paper Sec. V)

With the true mode-equation spectrum, the abundance matching lands on the
**light branch**: g ~ 1.0e-7, m_psi ~ 7.5e10 GeV (cold dark matter). The heavy
branch of the exponential closed form overproduces at the model's own T_reh.
The absolute normalization awaits a lattice/Floquet computation.

`treh_error_band.py` shows that this matching point is a window, not a number:
on the light branch g and m_psi scale as T_reh^(-1/2), so a factor-100 band in
T_reh moves them by a factor 10, while it moves n_s by only ~1.1e-3 and r by
~5.5 percent. The branch stays cold dark matter across the whole band. The same
script verifies the exact relation m_psi/H_inf = (g/sqrt(xi)) (M_Pl/H_inf),
which ties the dark-matter mass to the inflationary scale through the shared VEV.

The N <-> T_reh matching behind every tabulated T_reh* value uses comoving-entropy
conservation across the radiation era (`derive_from_action.N_match_derived` with
`entropy_matching=True`) together with rho_end = K_end + V_end at the end of
inflation, rather than the constant-g_* scaling a ~ rho^(-1/4). The independent
first-principles recomputation in `background_and_reheating.py` now agrees with
that table to 0.092 percent at N=50 and 0.087 percent at N=55, i.e. a single
N-independent constant offset; the two routes previously differed by a factor 3.4
because of that scaling plus a wrong x_* relation in the cross-check.

Both routes now share the same lambda0: `background_and_reheating.lam0_of_N`
delegates to the exact slow-roll A_s inversion in `derive_from_action`, through
`cosmo_model.lambda0_for_As_locked`, memoised. What remains is an implementation
difference, not a difference in approximation level, and the fitted closed form
6.70e-8 (50/N)^2 that this script used before is gone.

## Submission packaging

Suggested zip contents for submission: `paper_prd_merged.tex`, `sn-jnl.cls`,
`sn-mathphys-num.bst`, `figures/*.pdf`. The cover letter may note SCOAP3
funding (EPJC APC = 0).

**Funding**: the Declarations state that the author is an independent
researcher and that no funding was received for this study.

## Git

Do not commit: `scripts/.venv/`, `__pycache__/`, LaTeX intermediates (see
`.gitignore`). `.gitattributes` pins `* -text` so that every platform checks
out byte-exact files.
