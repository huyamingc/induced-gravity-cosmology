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
(Python 3.13, numpy 2.4 / scipy 1.18 / matplotlib 3.11) is about **25 minutes**
in total, of which `psi_abundance_oscillating.py` takes ~5.5 min and
`dm_gap_closure_test.py` ~19 min; `treh_error_band.py` finishes in under a
second. On a Windows console, either run through
`run_all.py` (it reconfigures the output streams to UTF-8) or set
`$env:PYTHONIOENCODING='utf-8'` before running a single script.

`scripts/_probe_heavy_branch.py` is a standalone diagnostic that is NOT part of
`run_all.py`: it probes a few large m/H points on the heavy branch of the
mode-equation abundance matching, reusing the machinery of
`dm_gap_closure_test.py`. It writes no files.

Key result scripts:

| Script | Content |
|---|---|
| `background_and_reheating.py` | exact KG integration, N window, reheating channels |
| `psi_production_bogoliubov.py` | de Sitter exponent 2 pi audit, g matching |
| `psi_abundance_oscillating.py` | cross-transition mode equation (power-law spectrum) |
| `dm_gap_closure_test.py` | small-g light branch and free-streaming length |
| `treh_error_band.py` | propagates the T_reh uncertainty to (N, n_s, r) and to the (g, m_psi) window |
| `residual_quintessence.py` | two-fluid integration, Delta w budget |

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
that table to 1.4 percent; the two routes previously differed by a factor 3.4
because of that scaling plus a wrong x_* relation in the cross-check.

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
