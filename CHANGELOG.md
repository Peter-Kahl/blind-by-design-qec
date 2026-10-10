# Changelog

All notable changes to this repository are recorded here.

## 0.6.1 — 2026-10-08

Removes the dependence of the history-used lower bounds on the grouping of syndromes into classes. Results of the other four scripts are unchanged.

### `prop52_numerical_check.py`

- `optimal_repair` takes an optional list of individual branches and returns a further lower bound: from the same dual point, the dual constraint is made feasible separately for every one of the 256 syndromes, so the bound holds for policies with one repair per syndrome and assumes nothing about proportionality within classes. Upper bounds were already valid for such policies, since one repair per class is one of them.
- Section 4 prints these bounds for both history-used repairs. In the reference run the sign-blind bound agrees with the grouped one to within 0.1 per cent at every angle; the known-sign bound at θ = 0.05 is negative, hence uninformative, and only the upper end of that bracket is used.
- New summary check: the sign-blind lower bound with the history used holds over all syndromes without grouping, within 1 per cent. The check that sign-blind repair is more than twice the known-sign optimum now uses the smaller of the two lower bounds.

### All scripts

- Version raised to 0.6.1. No other changes to the four other scripts.

### `README.md`

- Setup uses `python3 -m venv`, with a note for systems (including macOS) where `python` is not on the path until the environment is activated. Run time of `prop52_numerical_check.py` corrected to about 20 to 30 minutes. Records a cross-platform rerun of `prop52_numerical_check.py` on macOS (Python 3.14.7, NumPy 2.5.3): all checks pass, with brackets consistent with the reference. The reference outputs of `blind_by_design_continuous_misspecification.py` now keep the elapsed-time line the script prints, as the run produced it.

## 0.6.0 — 2026-10-08

Adds repairs that use the syndrome history, and lower bounds on every optimal repair. Results of the other four scripts are unchanged.

### `prop52_numerical_check.py`

- Section 4 rewritten. Four optimal repairs over all channels are bracketed at each of four angles: known sign and both signs, each with the syndrome history discarded (one repair after the standard correction) and with the history used (one repair per class of syndromes whose branch maps agree up to weight; 20 classes for the L = 3 instrument).
- Upper bounds: the defect achieved by the solver's repairs, made exact channels. Lower bounds: the dual objective at a dual point built from the solver's dual variables and made exactly feasible, by weak duality. Both are floating-point computations.
- The optimisation is vectorised (the composed Choi matrix is a fixed linear map of the repair's Choi matrix), which makes the history-conditioned problems tractable.
- Reference-run findings: with the history discarded, the fitted unitary repair is optimal within the bracket and no repair is within 3 per cent of the best sign-blind repair; using the history improves sign-blind repair by 1 to 29 per cent, and known-sign repair far more, so the known-sign advantage grows.
- Summary checks revised accordingly.

### All scripts

- Version raised to 0.6.0. No other changes to the four other scripts.

### Repository

- `README.md`, `CITATION.cff` and `.zenodo.json` updated; all reference outputs regenerated in the same environment.

## 0.5.0 — 2026-10-08

Replaces the unitary-only bound on uniform repair with an optimisation over all repair channels. Results of the other four scripts are unchanged.

### `prop52_numerical_check.py`

- Section 4 now computes, by semidefinite programming over all repair channels, the best repair for a known sign (pointwise recoverability) and the best single repair for both signs (uniform recoverability). Each is reported as the solver's objective value and as the defect achieved by the solver's repair after it is made an exact channel (clipped to positive and normalised to preserve trace), evaluated with the bounded diamond-norm routine.
- Reference-run findings: the fitted unitary repair equals the best repair over all channels to solver accuracy; the best single repair for both signs equals no repair to solver accuracy, at every angle checked.
- New summary checks: the uniform optimum exceeds twice the pointwise optimum; no repair is within 2 per cent of the uniform optimum; the fitted unitary repair is within 2 per cent of the pointwise optimum.
- The unitary-only lower bound ½‖Λ₊ − Λ₋‖⋄ is still printed for comparison.

### All scripts

- Version raised to 0.5.0. No other changes to the four other scripts.

### Repository

- `README.md`, `CITATION.cff` and `.zenodo.json` updated; all reference outputs regenerated in the same environment.

## 0.4.0 — 2026-10-08

Adds a comparison of pointwise and uniform repair, and states the status of the numerical bounds more cautiously. Results of the other four scripts are unchanged.

### `prop52_numerical_check.py`

- New section 4: for four angles, the defect of the repair fitted to +θ applied at +θ and at −θ, the defect with no repair, and a lower bound, ½‖Λ₊ − Λ₋‖⋄, on what any single unitary repair can achieve at both signs at once. The recoverability defect of Proposition 5.2 is pointwise; since the syndrome record cannot reveal the sign, a repair chosen from the record must be uniform. New summary check: the uniform unitary bound exceeds the pointwise defect at every angle.
- The bounds are now described as numerically constructed estimates in floating point, not certified by interval arithmetic or validated eigenvalue bounds ('rigorous up to floating-point rounding' withdrawn).
- Sections renumbered (repetition code is now 5, sentinel 6).

### All scripts

- Version raised to 0.4.0. No other changes to the four other scripts.

### Repository

- `README.md`, `CITATION.cff` and `.zenodo.json` updated; all reference outputs regenerated in the same environment.

## 0.3.0 — 2026-10-08

Adds a comparison of what the environment receives with what the record keeps, and corrects the description of the numerical bounds. Results of the other four scripts are unchanged.

### `prop52_numerical_check.py`

- New section 3: for four angles, the largest trace distance found between the complementary outputs of two pure encoded states (the syndrome register with coherences kept), compared with the record's leak (coherences discarded) and with the bound 2√ε_dec. In the reference run the complementary output distinguishes encoded states roughly 14,000 times better than the record at θ = 0.05, growing as θ² against the record's θ⁶, and comes within a factor of about 1.4 of the bound. Local search from fixed seeds; the values are lower bounds on the maxima. New summary check: record leak ≤ complement distance ≤ 2√ε_dec.
- Description of the bounds corrected: upper bounds are repaired dual points and lower bounds explicit entangled-input witnesses, both rigorous up to floating-point rounding (not interval arithmetic). 'Certified' removed.
- The docstring now states that the checks use one fitted unitary repair, so they test a necessary consequence of Proposition 5.2, not the bound at the infimum over repairs.
- Sections renumbered (repetition code is now 4, sentinel 5).

### All scripts

- Version raised to 0.3.0. No other changes to the four other scripts.

### Repository

- `README.md`, `requirements.txt` (SciPy noted), `CITATION.cff` and `.zenodo.json` updated; all reference outputs regenerated in the same environment as 0.2.0.

## 0.2.0 — 2026-10-07

Adds the numerical check of Proposition 5.2 (approximate classical-output blindness, §5.3.1 of the paper). Results of the four existing scripts are unchanged.

### New: `prop52_numerical_check.py`

- Checks the inequality of Proposition 5.2, leak ≤ 2√ε, on the L = 3 toric-code instrument of Özgüler (2026) over seven rotation angles (θ = 0.05 to 0.5), with ε computed both against the identity and after a fitted unitary decoder.
- Checks part (ii) of the proposition over N = 1, 2 and 3 rounds at θ = 0.1 and 0.3: the exact N-round leak against its bound, and ‖Λ^N − id‖⋄ ≤ N ε.
- Shows that the converse fails: in the three-qubit repetition code every three-round history law is independent of the encoded state, although the logical channel is not the identity.
- Reports, for contrast, the sign-discrimination profile of a sentinel qubit.
- Diamond norms are computed by semidefinite programme (Watrous 2013) with CVXPY and Clarabel. Reported values are certified upper bounds; every verdict is tested at a certified lower bound, so neither depends on solver tolerance. The SDP is validated against a closed form.
- The leak is computed exactly from the parity form of the history effects; the script reports the deviation from that form.
- Reuses the instrument built by `ozguler_prop1_check.py`; no code is duplicated.
- Ends with a PASS/FAIL summary, as `blind_by_design_toy_repetition_code.py` does.

### All scripts

- Version number raised to 0.2.0 in each header and in the printed run configuration. No other changes to the four existing scripts.

### Reference outputs

- All six reference outputs produced with version 0.2.0 in one environment: Python 3.12.3, NumPy 2.4.4, Stim 1.16.0, PyMatching 2.4.0, CVXPY 1.9.3, Clarabel 0.11.1.
- The five outputs carried over from 0.1.0 are identical to their 0.1.0 versions except for the version line.
- Added `output/prop52_numerical_check_output.txt`.

### Repository

- `README.md`: new script described; Propositions 5.1 and 5.2 summarised; correspondence table, requirements, running instructions, reference-output table, limits, structure listing and citations updated; version correspondence now points Version 1.0 of the paper to 0.2.0, superseding the correspondence stated for 0.1.0 below (0.1.0 accompanied a draft without Proposition 5.2).
- `requirements.txt`: CVXPY added; reference environment updated.
- `CITATION.cff` and `.zenodo.json`: version, date, description and keywords updated; references to Kretschmann, Schlingemann and Werner (2008) and Watrous (2013) added.

## 0.1.0 — 2026-10-04

First versioned release, corresponding to Version 1.0 of the paper. Numerical results are unchanged from the unversioned scripts except where noted.

### All scripts

- Version number, repository URL and paper section added to each header; header format aligned across scripts.
- Each script now prints a run configuration (script version, package versions, parameters) before its results.
- Minimum Python version stated as 3.8 for the NumPy scripts; `audit_checks.py` requires a Python version supported by Stim.

### `blind_by_design_toy_repetition_code.py`

- Output labelled by check, with a PASS/FAIL summary.
- New checks: effects are multiples of the identity on the code space; effects sum to the identity.
- History probabilities now computed from products of Kraus operators. The largest difference at ±θ is unchanged in substance (floating-point rounding) but now reads about 10⁻¹⁸ instead of about 2 × 10⁻¹⁶.

### `ozguler_prop1_check.py`

- New sanity checks: code states orthonormal and stabilised by every plaquette and star check; logical operators commute with the plaquette checks.
- Second logical operator renamed from `logical_y` to `logical_x2` (it is the second logical X operator, not a Y operator).
- Output now states that L = 2 lies outside the proposition's scope and that the parity form is not expected there.
- Numerical results unchanged.

### `blind_by_design_continuous_misspecification.py`

- Command-line options for every parameter (`--kappa`, `--gamma`, `--dt`, `--eta-true`, `--eta-assumed`, `--runs`, `--seed`, `--durations`); defaults reproduce §6.4.2.
- `--kappa 1` reproduces the parameter-sensitivity figures quoted in the paper.
- Input validation: the true and assumed efficiencies must lie on the candidate grid.
- Documentation states that the two default durations share a seed, so their records are not independent, and that re-estimation reruns the filter retrospectively on the same records.
- Numerical results unchanged.

### `audit_checks.py`

- Labels explained in the header (`B***` is the paper's baseline, `A5` its correlated error, `C_B` a discarded comparison family, e* the D6 boundary edge).
- Check 3 now lists D6's matching-graph edges and the boundary edge of its partner D9.
- New check 4 (`--lemma-checks`): all 520 single-channel deletions and 131 seeded random reweightings of the baseline, confirming that e* is never created and no edge is added.
- Fixed: comparison of matching edges no longer fails when no edges are shared; sorting of edge keys containing boundary edges no longer raises an error in demonstration mode.

### Repository

- `README.md` rewritten: correspondence with the paper's sections, reference outputs, reproducibility, limits.
- Added `CHANGELOG.md`, `CITATION.cff`, `.zenodo.json`, `docs/reliance_note_ozguler2026.md` and the reference outputs in `output/`.
- `requirements.txt` corrected (it had been copied from another repository) to list NumPy, Stim and PyMatching.
- Repository-structure listing corrected to the image actually present.
