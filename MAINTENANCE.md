# MAINTENANCE.md — maintainer notes

This file is public in the repository as part of the maintenance history but
stays out of the journal submission zip. `README.md` is the front door for
readers; this file keeps the historical narrative that is useful when editing
the code but noisy for readers.

---

## run_all.py and SystemExit

`run_all.py` executes every script with `runpy` in one process.
`verify_numerics.py`, `consistency_checks.py`, `audit_tex_numbers.py` and
`audit_readme_numbers.py`
end in `sys.exit()`. The driver now catches `SystemExit` per script so a
failing audit makes the whole run exit non-zero instead of silently truncating
the loop. Before that fix the loop stopped at `verify_numerics.py`:
`consistency_checks.py` never ran inside a full regression. `verify_numerics.py`
and `audit_tex_numbers.py` each returned/ignored a constant success code until
this round -- `verify_numerics` returned 0 even with failing claims, and
`audit_tex_numbers` had no exit code at all, so its FLAG lived only in
`tex_number_audit.md` and `run_all.py` could not see it.

## Exact vs order-of-magnitude (why the split exists)

Two different things were previously conflated under one label:

1. **A literal copied from an exact value.** `0.285204`, `0.19938` and
   `1.19938` were 6- and 5-digit truncations of `(1-u_e)^2`, `K/V` and
   `1 + K/V` at the end of inflation. They had been copied into four scripts,
   and `psi_abundance_oscillating.py` separately carried
   `rho_end = 1.63485e63`, computed under the *old* fitted
   `lambda0 = 6.70e-8` and drifted 4.3e-4 away from the current exact value
   (rel `0.00043` vs `cosmo_model.rho_end()`). Those literals are gone:
   `cosmo_model.py` evaluates the closed form once. The dynamical integration
   sits `2.6e-06` from the closed form `V_end/V0 = (1-u_e)^2`; that residual is
   integration resolution, and the closed form is what is exact.
   `consistency_checks.py` fails if a literal reappears.

2. **A claim that is order-of-magnitude in the paper.** `f_NL ~ -0.02`,
   `q ~ (m_t/m_chi)^2`, `Gamma_therm/H ~ 1e7`, `sigma_psiN ~ 1e-104 cm^2`,
   the anomaly width and the RG shifts `Delta xi`, `Delta lambda0` are
   order-of-magnitude *as physics*. `verify_numerics.py` checks them in an
   `order` tolerance class measured in decades.

## Retired scripts

Five scripts were retired; formulas that the manuscript still quotes were
consolidated into `order_estimates.py`.

- `extended_checks.py` and `quick_claims_check.py` held the order-of-magnitude
  claims. Deleting them outright would have left Secs. III, IV, V, VIII and XII
  quoting values with no executing source. What was dropped: an
  `Omega_DM = 0.12` abundance scaling, a pre-lattice `Gamma_anom` table, and a
  table of review verdicts already absorbed into the manuscript.
- `evaluate_paper.py`, `n_definition_check.py` and `_probe_heavy_branch.py` had
  no dependants; the first recommended creating `audit_tex_numbers.py`, which
  made its own report obsolete.

All five are archived outside the submission tree
(`review_workspace/retired_scripts/` locally; that directory is gitignored).

## audit_tex_numbers called twice

`audit_tex_numbers.py` used to be called twice in `run_all.py`, the second time
"after the figures regenerate". Figure scripts write only `.pdf`, so neither
the `.tex` nor `n_convention_results.json` can change between the calls. The
duplicate is gone.

`Delta lambda0 ~ 3e-14` in Sec. VIII E had no executing source until
`order_estimates.py` supplied one.

## Sec. XII stale mass pair

Sec. XII's "additional consistency checks" printed
`~7e-70 cm^2/g` (sigma_psipsi/m_psi; ~69 orders below the bullet-cluster bound)
and `Q ~ m_psi^4 ~ 3e43 GeV^4` after correction. Earlier drafts printed
`1e-67 cm^2/g` and `1e52 GeV^4`. A provenance sweep showed both land on those
older printed values at `m = 1e13 GeV` (a stale wimpzilla-scale mass), not at
the `m_psi = 7.5e10 GeV` that Sec. V selected at the time.

Both entries scale with the dark-matter mass, so the exact-background round
moved them again: at the current `m_psi = 1.1e11 GeV` the manuscript now prints
`~1e-69 cm^2/g` and `Q ~ 1.4e44 GeV^4` (the bullet-cluster margin stays at
~69 orders). `verify_numerics.py` forms them from `M_PSI_EXACT`, the matched
mass of `psi_mode_oscillating.json`.

`order_estimates.py` computes `sigma_psipsi_over_m` and `tremaine_gunn_Q`
from the same mass anchor Sec. V uses; the manuscript
quotes the matching values. Every other entry in that Sec. XII list states a
parametric bound for which the manuscript gives no closed form, so none of the
others carries an executable claim.

## m_psi paired with the wrong H_inf

`verify_numerics.py` used to form `m_psi` as `m_star_powerlaw` times the
drifting `H_inf` column of `treh_error_band.json`, which belongs to a different
background. The anchor ratio is measured at the fiducial `H_inf`, so the
correct partner is `p50["H_inf"]` / `bar.H_INF`. The old pairing understated
`m_psi` by ~0.43 percent.

The light-branch mass is now formed once as `m_star_powerlaw * H_inf`.

## Fig. 2 single writer

`fig2_ns_r.py` used to write both `fig2_ns_r.pdf` and a byte-identical
`fig2_ns_r_lockedN.pdf` alias, and `lock_n_convention.py` redrew the same
figure under the alias name, while the manuscript includes only
`figures/fig2_ns_r.pdf`. The alias is gone; `lock_n_convention.py` no longer
plots (its `n_convention_results.json` is unchanged).

## N ↔ T_reh route history

Earlier README prose said the first-principles $T_{\rm reh}^*$ agreed with the
tabulated one "to 1.4 percent"; that figure predated the `Omega_r` correction.
The measured spread is **0.092% at N=50** and **0.087% at N=55**. Both routes
now share the same exact `lambda0`. A second pass found that
`background_and_reheating` once hardcoded the stale percentage; it is now
computed from the two values.

## README audit patterns

`audit_readme_numbers.py` only checks prose that remains in the public
`README.md`:

- T_reh* agreement at N=50/55 (`0.092` / `0.087` percent)
- dark-matter abundance-matched coupling (`g ~ 1.5e-7`, `m_psi ~ 1.1e11 GeV`),
  now recomputed from `psi_mode_oscillating.json` (exact background) rather
  than from the transition-only `dm_gap_closure_test.json` anchor
- T_reh band effect (factor `10`; `~1.1e-3` in n_s; `~5.5` percent in r)

The README quotes the transition-only baseline with `=` rather than `~`
(`g = 1.0e-7, m_psi = 7.5e10 GeV`), so the `g ~ ...` pattern stays unique and
the audit cannot silently match the wrong pair. A missing upstream JSON is now
a reported failure instead of a silently skipped check.

Historical numbers that used to live in README (RHO_END drift, Sec. XII
correction pair, `m_psi` 0.43% understatement, closed-form `2.6e-06` residual)
are intentionally **not** audited against README anymore. If you put them back
into the public README, restore the matching patterns in
`scripts/audit_readme_numbers.py`.

## First audit-gap round: five tex fixes + window co-occurrence

A review pass found three structural blind spots that let semantic
contradictions survive a green audit suite: (1) a blacklist-only scan cannot
catch values that are new or unsourced, (2) verbatim code listings inside the
tex are scanned by nothing, and (3) single-point existence checks do not
enforce cross-section semantic consistency. The fixes:

- `audit_tex_numbers.py` gained a "Window co-occurrence" section: every quote
  of the `45--56` band must sit near the `N_max = 55.6` cap; the abstract and
  the `tab:sens` caption are required to carry it.
- `verify_numerics.py` gained the `T_reh*(51)` claim (2.2e9, 5% tolerance),
  which pins the App. A verbatim listing to the JSON ground truth.
- Five tex fixes: the Majorana mass term is even under Z2_psi (not excluded
  by it; the Dirac nature rests on the UV charge assignment), the N=51
  bracket now cites Table I (2.2e9) and maps the anomaly 1e9 GeV to the
  interpolated N=50.7, the verbatim listing 1.4e9 became 2.2e9, the abstract
  gained the N_max cap plus the DESI DR2 exposure, and the Table I caption
  states the viable interval [45, 55.6].

## Second round: stale section pointer, D0 script list, the 2.04 factor

- `Secs.~IV and VIII` pointed at the domain-wall chapter after the section
  reorganization; the chi--gauge discussion lives in Secs. IV and VII. A full
  sweep of all 60 `Sec.~` cross-references found this one error; the literal
  is now blacklisted in `audit_tex_numbers.py`.
- Appendix D0 listed six compute scripts but omitted `lock_n_convention.py`,
  the source of every Table I cell (README calls it "exact / source of
  Table I"). It is now listed with that role.
- The matching-route factor was quoted as 2.04 (0.23 e-folds), but the
  tabulated ratio across the 0fe113c route change is 4.4e7 -> 1.1e8 GeV:
  2.48 / 2.50 / 2.44 at N = 49/50/51. `ln(2.04)/3.05 = 0.234` shows 2.04 and
  0.23 were one bound pair under the e^{3.05 N} scaling, so both were
  replaced together: about 2.5 and about 0.30 (ln(2.5)/3.018 from the JSON
  grid). `factor $2.04$` and `0.23$ $e$-folds` are blacklisted.
- The conclusions' first 45--56 quote now carries the statistical-band vs
  N_max distinction, and the abstract was compressed 249 -> 243 words.

## Anchoring round: every hand-written magnitude tied to an executing source

- Root cause of the reviewer's number disputes: `order_estimates.py` (the
  declared "single source for every `~` claim") never absorbed the `Lambda_J`
  threshold suppressions, the Higgs portal, or the Discussion #14 quartic
  chain — those numbers had no executing source anywhere, and
  `verify_numerics.py` had no claim to catch them. The g-window and `6.6e-8`
  mismatches were prose rounding drift off correct script values.
- Manuscript fixes (18 spots): threshold bound `<= 1e-10` -> per-field
  `2.2e-8` / `1.2e-13`; Discussion #14 chain `1e-6 * (1e-46)^2 = 1e-98` ->
  `5.4e-10 * (2e-111)^2 ~ 2e-231` (the Sec. VI Ricci displacement); abstract
  "both consequences of `F = xi Phi^2`" -> the `Z2^psi` symmetry stated as an
  assumed anomaly-free UV remnant; g window `0.5--1.5e-7` -> `(0.47--1.5)e-7`
  matching the `T = 1e9` table row; portal `~1e-15` -> `~6e-17`; `lambda0
  6.6e-8` -> `6.70e-8` in the lattice paragraph; `Gamma_therm/H 1e7` -> `1e6`
  with the one-loop running `alpha_s(1e9) = 0.038` (cross-checks back to
  `alpha_s(m_chi) = 0.0262`); "only non-degenerate observational handle"
  softened at three spots; `45--56` capped at `N_max = 55.6` in Sec. X,
  Sec. XII and App. C; "selects" -> "allows on n_s alone"; the 2.5 route
  factor attributed to the combined entropy-matching / rho_end / Omega_r
  effect (g_*s dilution alone gives 3.0); the upper-bound g sentence carries
  its `T_reh = 2.1e8 GeV`; the Fig. 2 caption now reads "uses the locked-N
  grid from lock_n_convention.py".
- Script changes: `ALPHA_S_THERM` `0.1` -> `0.0377` computed by the new
  one-loop `alpha_s_running()`; new `lambda_J()`, `threshold_suppression()`,
  `delta_lambda_PhiH()`, `dw_quartic_ricci()`; `verify_numerics.py`
  Gamma_therm claim `1e7` -> `1e6` plus 7 new claims (67 -> 75, all green);
  `audit_tex_numbers.py` OLD_NOTE extended with 18 anchoring-round regression
  patterns (regex self-tested against sample old text).
- Suite after the round: run_all exit 0; verify 75/0; consistency 93/0;
  audit_tex FLAG=0; audit_readme 7/0; PDF 55 pages, no undefined refs,
  0 overfull boxes; abstract 241 words.

## Exact-background round: the dark-matter normalization is no longer pending

`psi_mode_oscillating.py` (Type: PAPER) integrates the Dirac mode equation on
the exact homogeneous background through the oscillating-condensate era, with
three controlled backgrounds: `p2` (the published transition-only convention,
kept as the validation anchor), `avg` (same slow roll, analytic
matter-domination extension) and `exact` (numerically resolved condensate
oscillations). It reports two systematics separately:

- the condensate oscillations raise production by `R = n_exact/n_avg ~ 1.24` at
  the matching point, which alone would lower the coupling by `R^-1/2 ~ 0.90`;
- replacing the instantaneous splice by the real slow-roll descent lowers the
  transition-only production baseline by a factor ~2.7, and dominates.

The net multiplier on the transition-only baseline is 1.48, so the matched
coupling moves from `g = 1.0e-7` to `g = 1.5e-7` (`m_psi = 1.1e11 GeV` at the
model's own `T_reh = 2.1e8 GeV`). **This reverses the manuscript's earlier
upper-bound reading** of the transition-only value (old L469 argued the
oscillating branch "can only add to the production", hence must lower `g`); the
sentence now states the measured direction, and the transition-only numbers are
retained throughout as the *labelled baseline*, never as the primary result.

Manuscript changes: ~45 spots, all moving the "absolute normalization awaits a
lattice" narrative to "closed at the homogeneous level, with a 3D lattice
simulation as the non-perturbative cross-check". The `tab:trehband` caption now
says the table is the transition-only baseline and that the exact-background
computation multiplies its `g` and `m_psi` columns by 1.48 (the
`T_reh^-1/2` scaling is unchanged), so the table did **not** need regenerating.
Appendix D0 gained the computation's description and the script's entry in its
list; its unitarity bound was corrected from `1e-13` to `1e-12` (the measured
worst drift is 7.7e-13, which the old written bound did not contain).

`verify_numerics.py` grew from 75 to 85 claims: the exact-background primary
values (`g`, `m_psi`, `m_psi/H_inf`), `R`, `g_exact/g_avg`, the transition
baseline factor, the Table 2 multiplier, and the three validation results (V1
reproduction at the `1e-9` level, exactly zero massless production, unitarity
within the printed bound) are now recomputed from `psi_mode_oscillating.json`.
Its docstring also records the structural limitation that the manuscript values
are hardcoded mirrors -- this script does not parse the tex, so
`audit_tex_numbers.py` and `consistency_checks.py` remain the guards against
drift in that direction.

Script review: all 20 scripts in `scripts/` were re-read and checked (logic,
physics formulas, output consistency, CLAUDE.md compliance). No bug was found
in any of them; the only defects fixed were the three audit/driver gaps above
(constant exit codes, unchecked missing JSON, the new script absent from
`run_all.py`) plus the missing `Type:` fields in the three figure scripts.

## Provenance audit round: what actually has a source, and what does not

The question asked was "can every number and every argument in the paper be
traced to a root?"  The honest answer needed a reverse check, because every
existing guard runs the other way.  `verify_numerics.py` holds 83 hand-written
`claim(...)` registrations and says so itself: *"the `paper` values below are
hardcoded mirrors of what main.tex prints; this script does NOT
parse the tex."*  So there was **no guarantee that a number printed in the
manuscript is inside any audit's field of view** -- if it was never registered,
nothing looks at it.

Method (tools in `review_workspace/`, see `provenance_audit_report.md`):

- `provenance_coverage_probe.py` extracts every numeric token from the tex
  (masking `\ref`/`\cite`/`\label` so cross-references are not counted) and
  classifies each scientific magnitude by which mechanism covers it:
  registered claim / table cell / external input / historical wording /
  untraced.  Scope is `a\times10^{b}` plus bare `10^{b}` with `|b| >= 3`, i.e.
  the manuscript's magnitude assertions; auxiliary coefficients (0.01, 0.02)
  and bookkeeping integers are out of scope by construction.
- `verify_untraced.py` recomputes each untraced magnitude from
  `lock_n_convention.point(11.1, 50)` and the artefact JSONs.

Result: 17 claim-covered + 34 table cells + 6 external + 6 historical against
**21 untraced distinct magnitudes (28 occurrences)**; five of the 21 are probe
false negatives (the claim's `paper` value is built from a format string or
carries a sign), leaving **16 genuinely sourceless magnitudes**.

Sixteen is not the alarming part -- `verify_untraced.py` shows most of them are
arithmetically right (the `10^{110}` hierarchy, `c_3 = -1.8e8 GeV`,
`Delta A/A = 6.3e-111`, `m_chi/H_0 = 2.3e55`, `rho_DM^(0) = 9.7e-48` all
reproduce to better than a factor 2).  The alarming part is that **they have no
executing source, so they cannot drift loudly**: change `lambda_0`, `m_chi` or
`T_reh` and they simply stay where they are.  Two of them have already gone
wrong that way:

- **L956** claims `m_chi/H` "never drops below `10^{50}`".  At reheating --
  where the ratio is smallest -- it is `2.3e13` at `T_reh = 1e9` GeV and still
  only `2.3e5` at the `1e15` GeV edge D5c allows.  `10^{55}` is *today's*
  value (L982/L1006); the sentence overstates by ~37 decades.
- **L1004** quotes a required suppression factor `~2e26--2e77` while the same
  sentence gives `rho_cond(a_end) ~ 3e63` and the line above gives
  `rho_cond(a_reh) < 5.2e-16--5.2e35`.  Those three numbers give
  `4.2e27--4.2e78`; the printed pair is exactly `29x` smaller, which is what
  one gets from `rho_cond(a_end) ~ 1.1e62`.  One of the two has to move.

Searching the ops record settled the fallback question: the sixteen magnitudes
appear **nowhere** in `MAINTENANCE.md`, `modification_log_scripts.md` (96 KB) or
`modification_log_tex.md` except `rho_paper = 9.5e44` in
`review_workspace/diag_sec12.py`.  The logs only ever handled disputed or
superseded values; a reverse coverage sweep had never been run, so numbers
nobody disputed were never logged either.

Fix path, in order of value: (1) register the closed-form ones in
`order_estimates.py` + `verify_numerics.py` so the README table covers them;
(2) correct L956 and L1004; (3) promote the probe into a failing audit with an
explicit exemption table (literature inputs, historical wording, auxiliary
coefficients), each exemption carrying its reason.  Nothing in `scripts/` or the
tex was modified in this round -- the audit is read-only.

## Registration round: the reverse sweep is now a failing audit

All three fix-path items executed and gated (run_all chain exit 0):

- `scripts/audit_provenance.py` (new, Type AUDIT): the probe promoted.  447
  magnitudes classified, `EXEMPT` holds 32 reasons, `UNTRACED > 0` exits 1;
  wired into `run_all.py` after `verify_numerics.py` and into the README fast
  path.  It writes `provenance_coverage.json/.md`, and `audit_readme_numbers`
  now checks the README's coverage sentence against that JSON (zero tolerance).
- `verify_numerics.py` 85 -> 139 claims: the closed-form and artefact-backed
  magnitudes registered in six new sections; the Delta-xi(Yukawa) loop unrolled
  so each printed value is a literal mirror; `Gamma_therm` mirror aligned to the
  printed 1.4e6.  Nothing pre-existing was retuned.
- Manuscript: the two numerical errors fixed (the `m_chi/H >= 1e50` floor and
  the D5c suppression factor), plus four decade roundings (Sec. VI radiative
  and tuning decades, `e^{-3N}`).  PDF 57 pages, clean.

Do not add an `EXEMPT` entry without a reason, and do not silence the audit by
loosening the 2% modulus match -- register the value instead.  The probe era
files in `review_workspace/` (`evidence_check.py`, `verify_untraced.py`,
`provenance_coverage_probe.py`, `provenance_audit_report.md`) are the audit's
working papers; `audit_provenance.py` supersedes them as the standing check.

## README dedup round: B5 single-source-of-truth cleanup

Audited README.md against the CLAUDE.md B5.1/B5.3 rules and removed six
hand-written duplications (Declarations summary in Packaging, the
audit_provenance mechanism described both in its roles row and in the Reverse
provenance section, the "EXEMPT is the only door" rule in three places, a
hand-written "139" in Forward provenance, the order-of-magnitude claim list in
two extra places, and the Maintenance rule duplicating the New-number
workflow) plus compressed the Dark-matter convention section per the B5.2
shrink table. All three RNUM anchor sentences were preserved verbatim and the
generated provenance block was untouched (verified: zero diff lines inside
BEGIN/END PROVENANCE after re-running provenance_map.py). Fast-path audits
all exit 0; RNUM 10 checks / 0 failures. See
`review_workspace/modification_log_scripts.md` (README dedup round) for the
per-item record.

## External-review round: pending-wording closure and window-as-range

An external review flagged three hard items; verification against the tex
confirmed all three and found two more instances of the same stale wording.
All fixed:

1. Introduction dark-matter bullet no longer says the relic abundance
   "requires a preheating simulation" (L58; also L802 "pending preheating
   simulation" and the Modifications (15) "pending dark-matter normalization"
   -- a full-text sweep found exactly these three residues).
2. Modifications list: item (19) moved after (18); numbering now runs 1..19
   in order, content unchanged.
3. README audit_tex_numbers row now points at the script's pattern list
   instead of enumerating a per-round subset of it.
4. Robustness-window paragraph: the x1.48-scaled window is now quoted as a
   range, g = (0.7--2.2)e-7 and m_psi = (5.0--16)e10 across the window, with
   the matched point g = 1.5e-7 tied to the model's own T_reh = 2.1e8 GeV.
   The four endpoints are registered as rounding claims in verify_numerics.py
   (teb row values times the g_shift multiplier); claims 139 -> 143.
5. README "Exact values vs order-of-magnitude claims" promoted to its own
   section directly above the Data-provenance block.

Gates: all six audits exit 0 (143/0, 93/0, FLAG=0, 450 magnitudes / 0
untraced, 143 call sites / 58 symbols, RNUM 10/0); PDF recompiled, 57 pages,
0 errors / 0 undefined / 0 overfull. Reverse-sweep delta accounting: +3
tokens (2.1e8 -> CLAIM via the T_reh^(anomaly) claim, bare 1e-7 -> CLAIM via
the transition-only baseline, bare 1e10 -> existing EXEMPT entry).
Not adopted, with reasons: extra x1.48 columns in Table 2 (three-way change
for a caption that already states the multiplier), merging Project map into
the roles table (B5.1 assigns them different facts).

## 1.6-vs-2.7 round: the two transition-shape factors are now cross-referenced

Reader-facing inconsistency: the abstract and Discussion quote the
transition-shape systematic as x1.6 (g space) while Sec. V / App. D0 quote a
factor ~2.7 (production-n space). Values are exactly consistent
(1.6447^2 = 2.705; 0.9016 x 1.6447 = 1.4829 = net multiplier), but no sentence
stated the relation. Fix: the App. D0 definition sentence now adds
"(a factor ~1.6 in the matched g, since g prop n_psi^{-1/2} on the light
branch)"; the abstract is unchanged. The x1.6 itself is now a registered
claim (verify_numerics 144 claims / 0 failures; README counters synced to
144 call sites). All six audits exit 0; PDF 57 pages, clean.

## Submission-cleanup round: revision history removed, abstract within the 150-250 limit

External pre-submission review verified against the tex (all six items real,
two worse than reported: the abstract was 271 words against the official
150-250 limit, and there were five "earlier drafts/versions" residues, not
three). Executed after confirmation:

- Removed the 19-item "Modifications relative to earlier drafts" paragraph
  (5487 chars), the Sec. IV revision-record parenthesis, and the "quoted in
  earlier versions" tail in Sec. III; rewrote the Liddle-Leach and the
  App. A shortcut-formula sentences as direct statements (the audit_tex
  blacklist context keywords "formula"/"Note." preserved and verified).
- Abstract: 271 -> 247 words (official limit 250), \textbf removed, one
  secondary sentence (Omega=Phi/PhiV identity, SM fifth-force suppression)
  dropped.
- Deduplication: the Acknowledgements "Supplemental Material:" sentence
  (near-verbatim duplicate of Code availability, and a misuse of the SI
  term) removed; Data/Code availability kept (officially required).
- Gates: six audits exit 0 (144/0, 93/0, FLAG=0, 428 magnitudes / 0
  untraced, 144 sites / 58 symbols, RNUM 10/0); PDF recompiled, 55 pages,
  0 errors / 0 undefined / 0 overfull; no "earlier draft/version" wording
  remains.
- Submission-file conclusions from the EPJC guidelines: Highlights.docx not
  required; DECLARATION_OF_INTERESTS.docx not required (system form +
  in-paper Statements and Declarations); cover letter is a system text box
  (draft at review_workspace/cover_letter_draft.txt with the three open
  points: conditional Z2^psi UV completion, pending lattice cross-check,
  DESI live exposure).

## Cover-letter + abstract margin round (2026-09-21 evening)

Readiness re-check on this machine (Python 3.13 system install; no `scripts/.venv`):
all six audits exit 0 (verify 144/0, consistency 93/0, tex FLAG=0, provenance 0
untraced, pmap OK, readme 10/0); pdflatex x2 clean, 55 pages.

Fixes applied after the check:

- `review_workspace/cover_letter_draft.txt`: DESI open-point citation
  `Discussion #18` was wrong (#18 is PPN beta) -> `Sec. X and Conclusions`;
  `our manuscript` -> `my manuscript`; `We wish` -> `I wish`; `GUIyang` ->
  `Guiyang`.
- Abstract compressed 250 -> ~227 words (presubmit split count) to leave margin
  under the EPJC 250-word cap. Locked tokens kept: `6.70e-8`, `0.00425`,
  `0.962`, `N~50.7`, `45--56` with `N_max~55.6`, `g~1.5e-7`,
  `m_psi~1.1e11`, `lambda_fs~1e-19`, `x1.6`, DESI DR2, `Delta xi<=1e-6`.
  audit_tex abstract PASS; window co-occurrence mandatory-missing-cap=0.

Not done this round (still open for submission): suggested-reviewers list.

## Full run_all.py gate (2026-09-21 19:54–20:05, this machine)

System Python 3.13.14 (numpy 2.4.6 / scipy 1.18 / matplotlib 3.11 / mpmath
1.3 / numba 0.66); no `scripts/.venv`. Command:
`scripts/run_all.py` via `PYTHONUNBUFFERED=1` + `PYTHONIOENCODING=utf-8`.
Log: `review_workspace/run_all_full_20260921_195456.log`.

**RUN_ALL_EXIT=0.** End-to-end about 11 min (faster than the 20–25 min
estimate because `psi_mode_oscillating` resumed a V1/V3/V4 snapshot and skipped
the validation battery; the spectrum scan + matching still recomputed:
exact-background `g=1.4828e-07`, `m_psi=1.0837e+11` GeV, `R=1.2420`).

Gates after the full run:
- verify_numerics **144/0**
- consistency_checks **93/0**
- audit_tex_numbers **FLAG=0**, mandatory-missing-cap=0
- audit_provenance **428 magnitudes / 0 untraced**
- provenance_map: README block rewritten (144 sites / 58 symbols)
- audit_readme_numbers **10/0**
- figures fig1/fig2/fig3 regenerated; dm_gap_closure_test **PASS**
  (light branch `g=1.0273e-07`, `m_psi=7.5085e+10`, `lambda_fs=9.52e-20` Mpc)

Noise only: `audit_provenance.py:19` SyntaxWarning (`\c` in a docstring);
does not affect exit codes.

## Figure-label polish (cosmetic only)

Fixed three in-figure annotation issues found during PDF visual QA; no
numerical content changed:

- `fig1_einstein_potential.py`: legend moved to `lower left`; inflation
  marker text moved to mid-right (`varphi=1: inflation`) so it no longer
  collides with the legend.
- `fig2_ns_r.py`: locked-N marker labels staggered in a diagonal fan
  (N=48..52) so the dense cluster is readable; N=55 stays to the right.
- `fig3_domain_wall.py`: `conformal artifact` note moved into the empty
  middle of the Einstein panel; legend to `lower right`.

Figures regenerated; pdflatex x2 clean (55 pages); consistency_checks
figure identity checks and verify_numerics still exit 0.

## Manuscript rename round: paper_prd_merged -> main (2026-09-21 late)

The PRD-era development filename `paper_prd_merged` was renamed to the
conventional `main` before EPJC submission. `git mv` renamed the tracked
`paper_prd_merged.tex` -> `main.tex` and `paper_prd_merged.pdf` -> `main.pdf`
(history preserved as renames); the untracked LaTeX intermediates
(`.aux/.log/.out`) were deleted and regenerated under the new name.

- 33 textual references updated by a single idempotent replace
  (`paper_prd_merged` -> `main`): 12 scripts under `scripts/` (28 spots,
  including the four functional `TEX = ROOT / "..."` readers
  `audit_tex_numbers.py`, `consistency_checks.py`, `audit_provenance.py`
  and the docstring mentions), `README.md` (9 spots: Project map, Layout,
  Run commands, provenance prose), `CLAUDE.md` (1 spot, App.1 Project map).
- No physics content, claim registration, JSON schema, or output filename
  of any compute script changed; the rename is purely nominal. The
  `review_workspace/` archives (apply_round*.py, modification logs,
  retired scripts) intentionally keep the old name -- they are historical
  records of rounds that ran under it.
- Gates after the rename: fast-path six audits re-run, all exit 0
  (144/0, 93/0, tex FLAG=0, provenance 0 untraced, PMAP block
  regenerated idempotently, RNUM 10/0); pdflatex x2 on `main.tex`, clean,
  same page count as before the rename.
- README/provenance: Project-map paths and Run commands synced; no numeric-zone
  changes; generated block rewritten by PMAP.

## Backups

All ten `.bak` files that used to sit in `scripts/` (3) and
`review_workspace/` (7) were deleted in the exact-background round, at the
user's request, after the full script review came back clean. They were
snapshots of files whose history is already in git (`HEAD` is clean, and every
round that created one was committed), so git is the rollback path now; if a
future round wants an on-disk snapshot, take it with `git stash` or a branch
rather than a `.bak` copy that no check knows about.

- `review_workspace/modification_log_scripts.md` -- CLAUDE.md 6.6 style log (scripts)
- `review_workspace/modification_log_tex.md` -- CLAUDE.md 6.6 style log (manuscript)

## Red lines still in force (CLAUDE.md §6)

Do not change: physics constants, CORE formulas, CSV/JSON output column names,
RNG seed rules, existing claim registrations without re-running the suite, or
output filenames the manuscript already cites. `MAINTENANCE.md` is public in
the repository but stays out of the journal submission zip; do not commit
`CLAUDE.md` or `review_workspace/`.
