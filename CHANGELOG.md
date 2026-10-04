# Changelog

All notable changes to this repository are recorded here.

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
