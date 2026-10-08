# Changelog

All notable changes to this repository are recorded here.

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
