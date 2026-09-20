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

`scripts/_probe_heavy_branch.py` is a standalone diagnostic that is NOT part of
`run_all.py`: it probes a few large m/H points on the heavy branch of the
mode-equation abundance matching, reusing the machinery of
`dm_gap_closure_test.py`. It writes no files.

Scripts and their role. `exact` means the script solves the stated equations
numerically and its numbers may be quoted; `order-of-magnitude` means it only
supports a claimed scaling, never a precise value.

| Script | Role | Content |
|---|---|---|
| `lock_n_convention.py` | **exact / source of Table I** | the authoritative locked-N grid: lambda0 inverted from A_s at fixed N, exact potential slow roll, and the T_reh*(N) matching. Every cell of Table I (`tab:sens`) is produced here |
| `derive_from_action.py` | exact / independent derivation | rebuilds the model from the action alone, presupposing none of the paper's formulas; holds the exact A_s inversion and the N(T_reh) matching used everywhere else |
| `background_and_reheating.py` | exact / cross-check | exact KG integration of the e-fold equations, the N window, and the two reheating channels; recomputes T_reh*(N) from first principles as an independent check on the Table I column. Shares the same exact lambda0 as `lock_n_convention.py` |
| `cosmo_model.py` | exact / shared constants | common constants and the memoised `lambda0_for_As_locked` wrapper used by the figure scripts |
| `psi_production_bogoliubov.py` | exact / cross-check | de Sitter exponent 2 pi audit, g matching |
| `psi_abundance_oscillating.py` | exact / cross-check | cross-transition mode equation (power-law spectrum) |
| `dm_gap_closure_test.py` | exact / cross-check | small-g light branch and free-streaming length |
| `treh_error_band.py` | exact / error propagation | propagates the T_reh uncertainty to (N, n_s, r) and to the (g, m_psi) window |
| `residual_quintessence.py` | exact / cross-check | two-fluid integration, Delta w budget |
| `n_definition_check.py` | exact / cross-check | how the definition of N changes r and n_s |
| `extended_checks.py` | **order-of-magnitude** | RG running, Gamma_anom and Omega_psi scaling: supports the claimed orders of magnitude only, not precise values |
| `quick_claims_check.py` | **order-of-magnitude** | small closed-form claims (f_NL, alpha_s, q, sigma_psiN, Gamma_th/H) |
| `verify_numerics.py` | audit | recomputes every quantitative claim printed in the paper from the current scripts; the reverse direction of `audit_tex_numbers.py` |
| `consistency_checks.py` | audit | parses the Table I and Table 2 bodies out of the `.tex` and checks each cell against the script that produces it, then compares independent routes against each other |
| `audit_tex_numbers.py` | audit | scans the `.tex` for superseded values and checks that any survivor sits in a comparison or historical context |
| `evaluate_paper.py` | audit | structural and cross-reference evaluation of the manuscript |
| `fig1_einstein_potential.py`, `fig2_ns_r.py`, `fig3_domain_wall.py` | figure | generate `figures/*.pdf` |

The two consistency auditors, `verify_numerics.py` and `consistency_checks.py`,
run last in `run_all.py` (they read the `.json` artefacts and the manuscript
tables) and together take under a second. They are the automated form of the
manual cross-checks that closed the N <-> T_reh matching review: run them after
any edit to the manuscript or to the matching code.

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
that table to 0.094 percent at N=50 and 0.090 percent at N=55, i.e. a single
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
